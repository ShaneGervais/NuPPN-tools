#!/usr/bin/env julia

include(joinpath(@__DIR__, "NovaRunTools.jl"))
using .NovaRunTools
using Dates
using Printf
using Random

Base.@kwdef mutable struct MCOptions
    nova::String = "nova_test"
    samples::Int = 100
    jobs::Int = get(ENV, "MAX_JOBS", nothing) === nothing ? 4 : parse(Int, ENV["MAX_JOBS"])
    seed::Int = 12345
    sigma_ln::Float64 = log(2.0)
    starlib_option::Int = 1
    runs_name::String = "runs_star_mc"
    build_only::Bool = false
    dry_run::Bool = false
    include_baseline::Bool = true
end

function usage()
    println("""
Usage:
  julia tools/run_ppn_star_mc.jl [options]

Options:
  --nova NAME             Nova directory under nova_cases/ (default: nova_test)
  --samples N             Number of Monte Carlo samples to build (default: 100)
  --jobs N                Maximum concurrent ppn jobs (default: MAX_JOBS or 4)
  --seed N                Random seed (default: 12345)
  --sigma-ln X            Lognormal sigma for sampled rate factors (default: log(2))
  --starlib-option N      STARLIB option: 1 = MC10/MC13 STL01, 2 = TALYS/ATOMKI STL02 (default: 1)
  --runs-name NAME        Output runs directory under the nova case (default: runs_star_mc)
  --no-baseline           Do not create/run the baseline directory
  --build-only            Build run directories but do not execute ppn.exe
  --dry-run               Resolve and print planned runs without writing directories
  -h, --help              Show this help
""")
end

function parse_args(args)
    opts = MCOptions()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            opts.nova = args[i]
        elseif arg == "--samples"
            i += 1
            i <= length(args) || error("--samples requires a value")
            opts.samples = parse(Int, args[i])
            opts.samples >= 1 || error("--samples must be >= 1")
        elseif arg == "--jobs"
            i += 1
            i <= length(args) || error("--jobs requires a value")
            opts.jobs = parse(Int, args[i])
            opts.jobs >= 1 || error("--jobs must be >= 1")
        elseif arg == "--seed"
            i += 1
            i <= length(args) || error("--seed requires a value")
            opts.seed = parse(Int, args[i])
        elseif arg == "--sigma-ln"
            i += 1
            i <= length(args) || error("--sigma-ln requires a value")
            opts.sigma_ln = parse(Float64, args[i])
            opts.sigma_ln >= 0 || error("--sigma-ln must be >= 0")
        elseif arg == "--starlib-option"
            i += 1
            i <= length(args) || error("--starlib-option requires a value")
            opts.starlib_option = parse(Int, args[i])
            opts.starlib_option in (1, 2) || error("--starlib-option must be 1 or 2")
        elseif arg == "--runs-name"
            i += 1
            i <= length(args) || error("--runs-name requires a value")
            opts.runs_name = args[i]
            isempty(strip(opts.runs_name)) && error("--runs-name cannot be empty")
        elseif arg == "--no-baseline"
            opts.include_baseline = false
        elseif arg == "--build-only"
            opts.build_only = true
        elseif arg == "--dry-run"
            opts.dry_run = true
        elseif arg in ("-h", "--help")
            usage()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end
    return opts
end

json_escape(text) = replace(string(text), "\\" => "\\\\", "\"" => "\\\"")

function set_starlib_option!(path, option)
    lines = readlines(path, keep=true)
    changed = false
    for i in eachindex(lines)
        if occursin(r"(?i)^\s*starlib_option\s*=", lines[i])
            lines[i] = replace(lines[i], r"(?i)^(\s*starlib_option\s*=\s*)\d+" => SubstitutionString("\\g<1>$option"); count=1)
            changed = true
            break
        end
    end
    if !changed
        insert_at = findfirst(line -> strip(line) == "/", lines)
        insert_at === nothing && error("Could not find ppn_physics namelist terminator in $path")
        insert!(lines, insert_at, "        starlib_option = $option\n")
    end
    write(path, join(lines))
end

function write_physics_mc_input(ppn_dir, run_dir, starlib_option, indices, factors)
    NovaRunTools.write_physics_input(ppn_dir, run_dir, indices, factors)
    set_starlib_option!(joinpath(run_dir, "ppn_physics.input"), starlib_option)
end

function resolve_mc_reactions(config, network)
    by_index = Dict{Int,Vector{Dict{String,Any}}}()
    index_order = Int[]

    for reaction in config["reactions"]
        index = NovaRunTools.resolve_reaction_index(reaction, network)
        row = NovaRunTools.row_by_index(network, index)
        row === nothing && error("$(reaction["name"]): index $index was not found after resolution")
        if !haskey(by_index, index)
            by_index[index] = Dict{String,Any}[]
            push!(index_order, index)
        end
        push!(by_index[index], reaction)

        reverse_index = NovaRunTools.validate_reverse_index(reaction, network, index)
        if reverse_index !== nothing
            if !haskey(by_index, reverse_index)
                by_index[reverse_index] = Dict{String,Any}[]
                push!(index_order, reverse_index)
            end
            push!(by_index[reverse_index], reaction)
        end
    end

    return index_order, by_index
end

function sample_factors(rng, indices, sigma_ln)
    Dict(index => exp(sigma_ln * randn(rng)) for index in indices)
end

function write_mc_manifest(path, opts, sample_id, factors_by_index, reactions_by_index)
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"nova\": \"$(json_escape(opts.nova))\",")
        println(io, "  \"sample\": $sample_id,")
        println(io, "  \"seed\": $(opts.seed),")
        println(io, "  \"sigma_ln\": $(opts.sigma_ln),")
        println(io, "  \"starlib_option\": $(opts.starlib_option),")
        println(io, "  \"sampling_model\": \"lognormal_factor = exp(sigma_ln * randn())\",")
        println(io, "  \"rates\": [")
        sorted_indices = sort(collect(keys(factors_by_index)))
        for (i, index) in enumerate(sorted_indices)
            reactions = reactions_by_index[index]
            names = [get(reaction, "name", "unknown") for reaction in reactions]
            articles = [get(reaction, "article_reaction", get(reaction, "name", "unknown")) for reaction in reactions]
            comma = i == length(sorted_indices) ? "" : ","
            println(io, "    {")
            println(io, "      \"index\": $index,")
            println(io, "      \"factor\": $(@sprintf("%.16E", factors_by_index[index])),")
            println(io, "      \"reaction_names\": [$(join(["\"" * json_escape(name) * "\"" for name in names], ", "))],")
            println(io, "      \"article_reactions\": [$(join(["\"" * json_escape(name) * "\"" for name in articles], ", "))]")
            println(io, "    }$comma")
        end
        println(io, "  ]")
        println(io, "}")
    end
end

function relative_run_name(runs_dir, run_dir)
    replace(relpath(run_dir, runs_dir), '/' => '_')
end

function list_ppn_executables(runs_dir; include_baseline=true)
    exes = String[]
    for (root, _, files) in walkdir(runs_dir)
        if "ppn.exe" in files
            if include_baseline || basename(root) != "baseline"
                push!(exes, joinpath(root, "ppn.exe"))
            end
        end
    end
    sort!(exes)
    isempty(exes) && error("No ppn.exe files found under $runs_dir")
    return exes
end

function run_and_log(cmd::Cmd, logfile)
    open(logfile, "a") do io
        process = run(pipeline(cmd, stdout=io, stderr=io), wait=false)
        wait(process)
        success(process) || error("Command failed: $cmd")
    end
end

function run_one_ppn(exe, runs_dir, logs_dir)
    run_dir = dirname(exe)
    name = relative_run_name(runs_dir, run_dir)
    logfile = joinpath(logs_dir, "$name.log")
    start_time = time()

    open(logfile, "w") do io
        println(io, "===========================================")
        println(io, "Running $run_dir")
        println(io, "Started $(Dates.now())")
        println(io, "===========================================")
    end

    println("Running $name")
    run_and_log(Cmd(`./ppn.exe`, dir=run_dir), logfile)

    elapsed = round(Int, time() - start_time)
    open(logfile, "a") do io
        println(io, "Finished $run_dir in $elapsed seconds")
    end
    println("Finished $name in $(elapsed)s")
end

function run_parallel(exes, runs_dir, logs_dir, jobs)
    queue = Channel{String}(length(exes))
    for exe in exes
        put!(queue, exe)
    end
    close(queue)

    failures = Channel{Tuple{String,Any}}(length(exes))
    workers = Task[]
    for _ in 1:min(jobs, length(exes))
        push!(workers, @async begin
            for exe in queue
                try
                    run_one_ppn(exe, runs_dir, logs_dir)
                catch err
                    put!(failures, (exe, err))
                end
            end
        end)
    end
    foreach(wait, workers)
    close(failures)

    failed = collect(failures)
    if !isempty(failed)
        println("Failures:")
        for (exe, err) in failed
            println("  $exe")
            println("    $err")
        end
        error("$(length(failed)) ppn jobs failed")
    end
end

function build_mc_runs(opts)
    base = NovaRunTools.nova_dir(opts.nova)
    ppn_dir = joinpath(base, "ppn")
    runs_dir = joinpath(base, opts.runs_name)
    plan_path = joinpath(base, "config", "reaction_plan.json")
    network_path = joinpath(ppn_dir, "networksetup.txt")

    isdir(ppn_dir) || error("Missing ppn directory: $ppn_dir")
    isfile(plan_path) || error("Missing reaction plan: $plan_path")
    isfile(network_path) || error("Missing networksetup: $network_path")

    config = NovaRunTools.parse_json_file(plan_path)
    network = NovaRunTools.parse_networksetup(network_path)
    indices, reactions_by_index = resolve_mc_reactions(config, network)

    println("Resolved $(length(config["reactions"])) reaction-plan entries to $(length(indices)) unique network index/indices.")
    duplicates = [(index, reactions) for (index, reactions) in reactions_by_index if length(reactions) > 1]
    if !isempty(duplicates)
        println("Shared network indices will receive one sampled factor each:")
        for (index, reactions) in sort(duplicates, by=x -> x[1])
            names = join([get(reaction, "name", "unknown") for reaction in reactions], ", ")
            println("  index $index: $names")
        end
    end

    if opts.dry_run
        println("Would write $(opts.samples) sample run(s) under $runs_dir with starlib_option=$(opts.starlib_option).")
        return runs_dir
    end

    if ispath(runs_dir)
        println("Deleting existing MC runs directory: $runs_dir")
        rm(runs_dir; recursive=true)
    end
    mkpath(runs_dir)

    if opts.include_baseline
        baseline_dir = joinpath(runs_dir, "baseline")
        NovaRunTools.copy_ppn(ppn_dir, baseline_dir)
        set_starlib_option!(joinpath(baseline_dir, "ppn_physics.input"), opts.starlib_option)
        write_mc_manifest(joinpath(baseline_dir, "mc_manifest.json"), opts, 0, Dict{Int,Float64}(), Dict{Int,Vector{Dict{String,Any}}}())
        println("Built STARLIB MC baseline in $baseline_dir")
    end

    rng = MersenneTwister(opts.seed)
    for sample_id in 1:opts.samples
        factors_by_index = sample_factors(rng, indices, opts.sigma_ln)
        run_dir = joinpath(runs_dir, @sprintf("sample_%06d", sample_id))
        NovaRunTools.copy_ppn(ppn_dir, run_dir)
        sorted_indices = sort(collect(keys(factors_by_index)))
        factors = [factors_by_index[index] for index in sorted_indices]
        write_physics_mc_input(ppn_dir, run_dir, opts.starlib_option, sorted_indices, factors)
        write_mc_manifest(joinpath(run_dir, "mc_manifest.json"), opts, sample_id, factors_by_index, reactions_by_index)
    end
    println("Built $(opts.samples) STARLIB MC sample run(s) in $runs_dir")
    return runs_dir
end

function main()
    opts = parse_args(ARGS)

    ENV["OMP_NUM_THREADS"] = get(ENV, "OMP_NUM_THREADS", "1")
    ENV["OPENBLAS_NUM_THREADS"] = get(ENV, "OPENBLAS_NUM_THREADS", "1")
    ENV["MKL_NUM_THREADS"] = get(ENV, "MKL_NUM_THREADS", "1")
    ENV["BLIS_NUM_THREADS"] = get(ENV, "BLIS_NUM_THREADS", "1")
    ENV["VECLIB_MAXIMUM_THREADS"] = get(ENV, "VECLIB_MAXIMUM_THREADS", "1")

    global_start = time()
    runs_dir = build_mc_runs(opts)
    opts.dry_run && return
    if opts.build_only
        println("Build-only requested; not running ppn.exe.")
        return
    end

    base = NovaRunTools.nova_dir(opts.nova)
    logs_dir = joinpath(base, "logs_$(opts.runs_name)")
    mkpath(logs_dir)
    exes = list_ppn_executables(runs_dir; include_baseline=opts.include_baseline)
    println("Starting $(length(exes)) STARLIB MC ppn run(s) with up to $(opts.jobs) concurrent job(s).")
    run_parallel(exes, runs_dir, logs_dir, opts.jobs)

    elapsed = round(Int, time() - global_start)
    println("==================================================")
    println("ALL REQUESTED STARLIB MC RUNS COMPLETE")
    println("Total time: $elapsed seconds")
    println("==================================================")
end

main()
