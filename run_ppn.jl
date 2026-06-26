#!/usr/bin/env julia

include(joinpath(@__DIR__, "NovaRunTools.jl"))
using .NovaRunTools
using Dates
using Printf

const ROOT = dirname(@__DIR__)
const PROJECT_ROOT = dirname(ROOT)
const DEFAULT_CASES_CONFIG = joinpath(ROOT, "config", "iliadis2002_model_cases.json")
const NETWORK_MODES = Set(["generated", "fixed"])

Base.@kwdef mutable struct Options
    cases_config::String = DEFAULT_CASES_CONFIG
    case::Union{Nothing,String} = nothing
    all_cases::Bool = false
    jobs::Int = get(ENV, "MAX_JOBS", nothing) === nothing ? 1 : parse(Int, ENV["MAX_JOBS"])
    baseline_only::Bool = false
    build_only::Bool = false
    dry_run::Bool = false
    flux::Bool = false
    network_mode::String = "generated"
    starlib_option::Int = 0
end

function usage()
    println("""
Usage:
  ./run_ppn [options]

Options:
  --case NAME           Build/run one nova case from config/iliadis2002_model_cases.json
  --all                 Build/run all configured nova cases
  --config PATH         Cases config file (default: config/iliadis2002_model_cases.json)
  --jobs N              Maximum concurrent ppn jobs (default: MAX_JOBS or 1)
  --baseline-only       Build/run only the baseline for each selected case
  --build-only          Build run directories but do not execute ppn.exe
  --dry-run             Resolve and print planned runs without writing directories
  --flux                Enable flux_*.DAT output in ppn_frame.input
  --network-mode MODE   Network handling: generated writes a fresh networksetup.txt (ININET=0);
                        fixed prepares ININET=3 inputs for build-only inspection. Default: generated
  --starlib-option N    STARLIB_OPTION value written to ppn_physics.input. Default: 0
  -h, --help            Show this help
""")
end

function parse_args(args)
    opts = Options()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--case"
            i += 1
            i <= length(args) || error("--case requires a value")
            opts.case = args[i]
        elseif arg == "--all"
            opts.all_cases = true
        elseif arg == "--config"
            i += 1
            i <= length(args) || error("--config requires a value")
            opts.cases_config = abspath(args[i])
        elseif arg == "--jobs"
            i += 1
            i <= length(args) || error("--jobs requires a value")
            opts.jobs = parse(Int, args[i])
            opts.jobs >= 1 || error("--jobs must be >= 1")
        elseif arg == "--baseline-only"
            opts.baseline_only = true
        elseif arg == "--build-only"
            opts.build_only = true
        elseif arg == "--dry-run"
            opts.dry_run = true
        elseif arg == "--flux"
            opts.flux = true
        elseif arg == "--network-mode"
            i += 1
            i <= length(args) || error("--network-mode requires a value")
            opts.network_mode = lowercase(args[i])
            opts.network_mode in NETWORK_MODES || error("--network-mode must be one of: $(join(sort(collect(NETWORK_MODES)), ", "))")
        elseif arg == "--starlib-option"
            i += 1
            i <= length(args) || error("--starlib-option requires a value")
            opts.starlib_option = parse(Int, args[i])
            opts.starlib_option >= 0 || error("--starlib-option must be >= 0")
        elseif arg in ("-h", "--help")
            usage()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end
    opts.case === nothing && !opts.all_cases && (opts.case = default_case(opts.cases_config))
    return opts
end

network_mode_ininet(mode) = mode == "generated" ? 0 : 3

function default_case(cases_config)
    config = NovaRunTools.parse_json_file(cases_config)
    cases = config["cases"]
    isempty(cases) && error("No cases configured in $cases_config")
    return first(cases)["nova_case"]
end

function project_path(path)
    isabspath(path) ? path : joinpath(ROOT, path)
end

function selected_cases(opts)
    config = NovaRunTools.parse_json_file(opts.cases_config)
    cases = config["cases"]
    if opts.all_cases
        return cases
    end
    selected = [case for case in cases if case["nova_case"] == opts.case]
    isempty(selected) && error("Case $(opts.case) is not configured in $(opts.cases_config)")
    return selected
end

function load_reactions(case_config)
    reaction_file = get(case_config, "reaction_file", nothing)
    reaction_file !== nothing || error("Case $(case_config["nova_case"]) is missing reaction_file")
    config = NovaRunTools.parse_json_file(project_path(reaction_file))
    return config["reactions"]
end

function case_input_paths(case_name)
    case_dir = joinpath(ROOT, "nova_cases", case_name)
    trajectory = joinpath(case_dir, "trajectory_$(case_name).input")
    initial = joinpath(case_dir, "iniab_$(case_name).dat")
    isfile(trajectory) || error("Missing trajectory file: $trajectory")
    isfile(initial) || error("Missing initial abundance file: $initial")
    return case_dir, trajectory, initial
end

function force_namelist_int!(path, key, value)
    text = read(path, String)
    pattern = Regex("(?im)^(\\s*$(key)\\s*=\\s*)\\d+")
    new_text = replace(text, pattern => SubstitutionString("\\g<1>$value"); count=1)
    if new_text == text
        lines = readlines(path, keep=true)
        inserted = false
        out = String[]
        for line in lines
            if !inserted && strip(line) == "/"
                push!(out, "        $key = $value\n")
                inserted = true
            end
            push!(out, line)
        end
        inserted || error("Could not find ppn_physics namelist terminator in $path")
        new_text = join(out)
    end
    write(path, new_text)
end

force_ininet!(path, value) = force_namelist_int!(path, "ININET", value)
force_starlib_option!(path, value) = force_namelist_int!(path, "STARLIB_OPTION", value)

function set_flux_option!(path, enabled)
    text = read(path, String)
    value = enabled ? "1" : "0"
    new_text = replace(text, r"(?im)^(\s*iplot_flux_option\s*=\s*)\d+" => SubstitutionString("\\g<1>$value"); count=1)
    write(path, new_text)
end

function replace_rate_factor_column(line, value)
    eol = endswith(line, "\n") ? "\n" : ""
    body = chomp(line)
    # Reaction rows end with rfac followed by bind_energy_diff; preserve the final Q-value field.
    pattern = r"(\s+)[-+]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?(\s+[-+]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?\s*)$"
    replace(body, pattern => SubstitutionString("\\g<1>$value\\g<2>"); count=1) * eol
end

function reset_network_factors!(path)
    lines = readlines(path, keep=true)
    network = parse_networksetup(path)
    for row in network
        lines[row.line_no] = replace_rate_factor_column(lines[row.line_no], "1.000E+00")
    end
    write(path, join(lines))
end

function set_network_factors!(path, factors_by_index)
    lines = readlines(path, keep=true)
    network = parse_networksetup(path)
    rows_by_index = Dict(row.index => row for row in network)
    missing = Int[]
    for (idx, factor) in factors_by_index
        row = get(rows_by_index, idx, nothing)
        if row === nothing
            push!(missing, idx)
            continue
        end
        lines[row.line_no] = replace_rate_factor_column(lines[row.line_no], @sprintf("%.3E", factor))
    end
    isempty(missing) || error("Missing networksetup indices: $(join(missing, ", "))")
    write(path, join(lines))
end

function write_physics_rate_factors!(path, factors_by_index)
    lines = readlines(path, keep=true)
    new_lines = String[]
    inserted = false
    i = 1
    for line in lines
        if !inserted && strip(line) == "/"
            for idx in sort(collect(keys(factors_by_index)))
                i <= 10 || error("PPN physics_knobs supports at most 10 rate factors per run")
                push!(new_lines, "        rate_index($i) = $idx\n")
                push!(new_lines, "        rate_factor($i) = $(@sprintf("%.12E", factors_by_index[idx]))\n")
                i += 1
            end
            inserted = true
        end
        push!(new_lines, line)
    end
    inserted || error("Could not find ppn_physics namelist terminator in $path")
    write(path, join(new_lines))
end

function row_by_name(network, reaction)
    if haskey(reaction, "index")
        row = NovaRunTools.row_by_index(network, Int(reaction["index"]))
        row === nothing && return nothing
        return row
    end
    haskey(reaction, "name") || return nothing
    candidates, _ = NovaRunTools.matching_rows_for_reaction(reaction, network)
    active = [row for row in candidates if row.active]
    !isempty(active) && return sort(active, by=row -> row.index)[1]
    isempty(candidates) ? nothing : sort(candidates, by=row -> row.index)[1]
end

function npdata_source_path()
    project_npdata = joinpath(PROJECT_ROOT, "NPDATA")
    if ispath(project_npdata)
        return realpath(project_npdata)
    end

    return realpath(joinpath(ROOT, "..", "..", "physics", "NPDATA"))
end

function run_label(reaction, factor)
    label = if haskey(reaction, "name")
        reaction["name"]
    elseif haskey(reaction, "article_reaction")
        reaction["article_reaction"]
    else
        "index_$(reaction["index"])"
    end
    clean_label = replace(label, r"[^A-Za-z0-9_+.-]" => "_")
    factor_label = replace(@sprintf("%.3E", factor), "+" => "p", "-" => "m")
    return clean_label, "factor_$factor_label"
end

function copy_common_inputs!(run_dir, trajectory, initial; flux=false, network_mode="generated", starlib_option=0)
    mkpath(run_dir)
    npdata_src = npdata_source_path()
    parent_npdata = joinpath(dirname(run_dir), "NPDATA")
    if !ispath(parent_npdata)
        symlink(npdata_src, parent_npdata)
    end
    local_npdata = joinpath(run_dir, "NPDATA")
    if !ispath(local_npdata)
        symlink(npdata_src, local_npdata)
    end
    for file in ("ppn_solver.input", "ppn_physics.input", "ppn_frame.input",
                 "networksetup.txt", "isotopedatabase.txt", "isotopedatabase_cf.txt",
                 "isotopedatabase_all.txt")
        src = joinpath(ROOT, file)
        isfile(src) && cp(src, joinpath(run_dir, file); force=true)
    end
    cp(trajectory, joinpath(run_dir, "trajectory.input"); force=true)
    cp(initial, joinpath(run_dir, "initial_abundance.dat"); force=true)
    exe_link = joinpath(run_dir, "ppn.exe")
    ispath(exe_link) && rm(exe_link; force=true)
    symlink(joinpath(ROOT, "ppn.exe"), exe_link)

    physics_path = joinpath(run_dir, "ppn_physics.input")
    force_ininet!(physics_path, network_mode_ininet(network_mode))
    force_starlib_option!(physics_path, starlib_option)
    set_flux_option!(joinpath(run_dir, "ppn_frame.input"), flux)
    reset_network_factors!(joinpath(run_dir, "networksetup.txt"))
end

function write_run_manifest!(run_dir, case_config, reaction, factor, indices; starlib_option=0)
    open(joinpath(run_dir, "run_manifest.json"), "w") do io
        println(io, "{")
        println(io, "  \"nova_case\": \"$(case_config["nova_case"])\",")
        println(io, "  \"iliadis_model\": \"$(get(case_config, "iliadis_model", "unknown"))\",")
        println(io, "  \"reaction\": \"$(get(reaction, "article_reaction", get(reaction, "name", "baseline")))\",")
        println(io, "  \"factor\": $factor,")
        println(io, "  \"applied_factor\": $(get(reaction, "applied_factor", factor)),")
        println(io, "  \"baseline_factor\": $(get(reaction, "baseline_factor", 1.0)),")
        println(io, "  \"network_mode\": \"$(get(case_config, "network_mode", "generated"))\",")
        println(io, "  \"starlib_option\": $starlib_option,")
        println(io, "  \"network_indices\": [$(join(indices, ", "))]")
        println(io, "}")
    end
end

function baseline_factor(reaction)
    return Float64(get(reaction, "baseline_factor", 1.0))
end

function configured_baseline_factors(network, reactions)
    factors_by_index = Dict{Int,Float64}()
    for reaction in reactions
        base = baseline_factor(reaction)
        base == 1.0 && continue
        row = row_by_name(network, reaction)
        row === nothing && continue
        linked = [Int(i) for i in get(reaction, "linked_indices", Any[])]
        for idx in vcat([row.index], linked)
            factors_by_index[idx] = base
        end
    end
    return factors_by_index
end

function build_case(case_config; baseline_only=false, dry_run=false, flux=false, network_mode="generated", starlib_option=0)
    case_name = case_config["nova_case"]
    case_config["network_mode"] = network_mode
    _, trajectory, initial = case_input_paths(case_name)
    output_dir = joinpath(ROOT, get(case_config, "output_dir", joinpath("nova_runs", case_name)))
    baseline_dir = joinpath(output_dir, "baseline")
    run_dirs = String[]

    network = parse_networksetup(joinpath(ROOT, "networksetup.txt"))
    reactions = load_reactions(case_config)
    baseline_factors_by_index = configured_baseline_factors(network, reactions)

    if dry_run
        println("Would build baseline: $baseline_dir using network mode $network_mode starlib_option $starlib_option")
        if !isempty(baseline_factors_by_index)
            details = ["$idx => $factor" for (idx, factor) in sort(collect(baseline_factors_by_index))]
            println("Would apply baseline rate factors: $(join(details, ", "))")
        end
    else
        ispath(baseline_dir) && rm(baseline_dir; recursive=true)
        copy_common_inputs!(baseline_dir, trajectory, initial; flux=flux, network_mode=network_mode, starlib_option=starlib_option)
        if !isempty(baseline_factors_by_index)
            set_network_factors!(joinpath(baseline_dir, "networksetup.txt"), baseline_factors_by_index)
            write_physics_rate_factors!(joinpath(baseline_dir, "ppn_physics.input"), baseline_factors_by_index)
        end
        write_run_manifest!(baseline_dir, case_config, Dict{String,Any}("article_reaction" => "baseline"), 1.0, collect(keys(baseline_factors_by_index)); starlib_option=starlib_option)
        println("Built baseline: $baseline_dir")
    end
    push!(run_dirs, baseline_dir)
    baseline_only && return run_dirs

    skipped = Tuple{String,String}[]

    for reaction in reactions
        row = row_by_name(network, reaction)
        label = get(reaction, "article_reaction", get(reaction, "name", "unknown"))
        if row === nothing
            push!(skipped, (label, "no matching networksetup row"))
            continue
        end
        linked = [Int(i) for i in get(reaction, "linked_indices", Any[])]
        indices = vcat([row.index], linked)
        factors = get(reaction, "factors", Any[])
        isempty(factors) && (factors = [get(reaction, "factor", 1.0)])
        for factor_any in factors
            factor = Float64(factor_any)
            base = baseline_factor(reaction)
            applied_factor = factor * base
            reaction_dir, factor_dir = run_label(reaction, factor)
            run_dir = joinpath(output_dir, reaction_dir, factor_dir)
            if dry_run
                println("Would build $run_dir using indices $(join(indices, ", ")) factor $factor applied_factor $applied_factor network mode $network_mode starlib_option $starlib_option")
            else
                ispath(run_dir) && rm(run_dir; recursive=true)
                copy_common_inputs!(run_dir, trajectory, initial; flux=flux, network_mode=network_mode, starlib_option=starlib_option)
                factors_by_index = copy(baseline_factors_by_index)
                for idx in indices
                    factors_by_index[idx] = applied_factor
                end
                set_network_factors!(joinpath(run_dir, "networksetup.txt"), factors_by_index)
                write_physics_rate_factors!(joinpath(run_dir, "ppn_physics.input"), factors_by_index)
                manifest_reaction = copy(reaction)
                manifest_reaction["applied_factor"] = applied_factor
                write_run_manifest!(run_dir, case_config, manifest_reaction, factor, sort(collect(keys(factors_by_index))); starlib_option=starlib_option)
            end
            push!(run_dirs, run_dir)
        end
    end

    if !isempty(skipped) && !dry_run
        mkpath(output_dir)
        open(joinpath(output_dir, "skipped_reactions.txt"), "w") do io
            for (label, reason) in skipped
                println(io, "$label\t$reason")
            end
        end
        println("Skipped $(length(skipped)) reaction(s) for $case_name; see $(joinpath(output_dir, "skipped_reactions.txt"))")
    elseif !isempty(skipped)
        println("Would skip $(length(skipped)) reaction(s) for $case_name:")
        for (label, reason) in skipped
            println("  $label: $reason")
        end
    end
    return run_dirs
end

function run_and_log(cmd::Cmd, logfile)
    open(logfile, "a") do io
        process = run(pipeline(cmd, stdout=io, stderr=io), wait=false)
        wait(process)
        success(process) || error("Command failed: $cmd")
    end
end

function run_one(run_dir)
    case_root = basename(run_dir) == "baseline" ? dirname(run_dir) : dirname(dirname(dirname(run_dir)))
    logs_dir = joinpath(case_root, "logs")
    mkpath(logs_dir)
    name = replace(relpath(run_dir, case_root), '/' => '_')
    logfile = joinpath(logs_dir, "$name.log")
    open(logfile, "w") do io
        println(io, "Running $run_dir")
        println(io, "Started $(Dates.now())")
    end
    println("Running $name")
    start = time()
    run_and_log(Cmd(`./ppn.exe`, dir=run_dir), logfile)
    elapsed = round(Int, time() - start)
    open(logfile, "a") do io
        println(io, "Finished $(Dates.now()) in $elapsed seconds")
    end
    println("Finished $name in $(elapsed)s")
end

function run_parallel(run_dirs, jobs)
    queue = Channel{String}(length(run_dirs))
    foreach(dir -> put!(queue, dir), run_dirs)
    close(queue)
    failures = Channel{Tuple{String,Any}}(length(run_dirs))
    workers = Task[]
    for _ in 1:min(jobs, length(run_dirs))
        push!(workers, @async begin
            for run_dir in queue
                try
                    run_one(run_dir)
                catch err
                    put!(failures, (run_dir, err))
                end
            end
        end)
    end
    foreach(wait, workers)
    close(failures)
    failed = collect(failures)
    if !isempty(failed)
        for (run_dir, err) in failed
            println("FAILED $run_dir")
            println("  $err")
        end
        error("$(length(failed)) run(s) failed")
    end
end

function main()
    opts = parse_args(ARGS)
    if opts.network_mode == "fixed" && !opts.dry_run && !opts.build_only
        error("--network-mode fixed is build-only right now because PPN's rnetw2008 exits for ININET=3; use --network-mode generated for executable runs")
    end
    ENV["OMP_NUM_THREADS"] = get(ENV, "OMP_NUM_THREADS", "1")
    ENV["OPENBLAS_NUM_THREADS"] = get(ENV, "OPENBLAS_NUM_THREADS", "1")
    ENV["MKL_NUM_THREADS"] = get(ENV, "MKL_NUM_THREADS", "1")
    ENV["BLIS_NUM_THREADS"] = get(ENV, "BLIS_NUM_THREADS", "1")
    ENV["VECLIB_MAXIMUM_THREADS"] = get(ENV, "VECLIB_MAXIMUM_THREADS", "1")

    all_run_dirs = String[]
    for case_config in selected_cases(opts)
        append!(all_run_dirs, build_case(case_config;
            baseline_only=opts.baseline_only,
            dry_run=opts.dry_run,
            flux=opts.flux,
            network_mode=opts.network_mode,
            starlib_option=opts.starlib_option))
    end
    opts.dry_run && return
    opts.build_only && return
    println("Starting $(length(all_run_dirs)) PPN run(s) with up to $(opts.jobs) job(s).")
    run_parallel(all_run_dirs, opts.jobs)
end

main()
