#!/usr/bin/env julia

# Decay every PPN run in a nova case and mirror runs/ into decay_runs/.
# Julia port of decay_all_runs.py.

include(joinpath(@__DIR__, "NovaRunTools.jl"))
using .NovaRunTools
using Dates
using Printf

const DEFAULT_DECAY_TIME_SECONDS = 2.0 * 60.0 * 60.0
const TEXT_INPUTS = sort([
    "Makefile",
    "isotopedatabase.txt",
    "isotopedatabase_all.txt",
    "isotopedatabase_cf.txt",
    "networksetup.txt",
    "ppn_frame.input",
    "ppn_physics.input",
    "ppn_solver.input",
])

struct DecayJob
    source_run::String
    decay_run::String
    source_iso_massf::String
    final_cycle::Int
end

mutable struct DecayResult
    job::DecayJob
    returncode::Int
    output::String
    log::String
    status::String
end

Base.@kwdef mutable struct SweepOptions
    nova::String = "nova_test"
    runs_name::String = "runs"
    decay_runs_name::String = "decay_runs"
    decay_time::Float64 = DEFAULT_DECAY_TIME_SECONDS
    jobs::Int = get(ENV, "MAX_JOBS", nothing) === nothing ? 4 : parse(Int, ENV["MAX_JOBS"])
    dry_run::Bool = false
    no_run::Bool = false
end

function usage()
    println("""
Usage:
  julia tools/decay_ppn_sweep.jl [options]

Mirror nova_cases/<nova>/runs into nova_cases/<nova>/decay_runs and run PPN
decay for each run's final iso_massf file.

Options:
  --nova NAME             Nova directory under nova_cases/ (default: nova_test)
  --runs-name NAME        Source runs directory under the nova case (default: runs)
  --decay-runs-name NAME  Output decay runs directory under the nova case (default: decay_runs)
  --decay-time SECONDS    Decay time in seconds (default: 7200, i.e. 2 hours)
  -j, --jobs N            Number of ppn.exe decay runs to execute in parallel (default: MAX_JOBS or 4)
  --dry-run               Print the jobs that would be built and executed without writing files
  --no-run                Build decay run directories and manifest but do not execute ppn.exe
  -h, --help              Show this help
""")
end

function parse_args(args)
    opts = SweepOptions()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            opts.nova = args[i]
        elseif arg == "--runs-name"
            i += 1
            i <= length(args) || error("--runs-name requires a value")
            opts.runs_name = args[i]
            isempty(strip(opts.runs_name)) && error("--runs-name cannot be empty")
        elseif arg == "--decay-runs-name"
            i += 1
            i <= length(args) || error("--decay-runs-name requires a value")
            opts.decay_runs_name = args[i]
            isempty(strip(opts.decay_runs_name)) && error("--decay-runs-name cannot be empty")
        elseif arg == "--decay-time"
            i += 1
            i <= length(args) || error("--decay-time requires a value")
            opts.decay_time = parse(Float64, args[i])
        elseif arg in ("-j", "--jobs")
            i += 1
            i <= length(args) || error("--jobs requires a value")
            opts.jobs = parse(Int, args[i])
            opts.jobs >= 1 || error("--jobs must be >= 1")
        elseif arg == "--dry-run"
            opts.dry_run = true
        elseif arg == "--no-run"
            opts.no_run = true
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

cycle_name(cycle::Int) = @sprintf("%05d", cycle)

function final_cycle_from_xtime(run_dir)
    xtime = joinpath(run_dir, "x-time.dat")
    isfile(xtime) || error("missing x-time.dat: $xtime")

    final_cycle = nothing
    for line in eachline(xtime)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        final_cycle = parse(Int, split(stripped)[1])
    end
    final_cycle === nothing && error("no cycle rows found in $xtime")
    return final_cycle
end

function discover_source_runs(runs_dir)
    runs = String[]
    for (root, _, files) in walkdir(runs_dir)
        "ppn.exe" in files && isfile(joinpath(root, "ppn.exe")) && push!(runs, root)
    end
    sort!(runs)
    return runs
end

function parse_iso_massf_rows(path)
    isfile(path) || error("missing iso_massf input: $path")

    rows = NamedTuple[]
    for line in eachline(path)
        parts = split(line)
        isempty(parts) && continue
        all(isdigit, parts[1]) || continue
        push!(
            rows,
            (
                z=round(Int, parse(Float64, parts[2])),
                a=round(Int, parse(Float64, parts[3])),
                abundance=parse(Float64, parts[5]),
                label_parts=parts[6:end],
            ),
        )
    end
    isempty(rows) && error("no abundance rows found in $path")
    return rows
end

function format_iso_label(a, label_parts, source_path)
    if label_parts == ["NEUT"]
        return "NEUT"
    elseif label_parts == ["PROT"]
        return "PROT"
    end

    element = if length(label_parts) >= 2
        lowercase(label_parts[1])
    else
        m = match(r"^([A-Za-z]+)(\d+)$", label_parts[1])
        m === nothing && error("cannot convert isotope label from $source_path: $label_parts")
        lowercase(m.captures[1])
    end
    @sprintf("%-2s%3d", element, a)
end

function write_post_abundance(source_iso_massf, out_path)
    rows = parse_iso_massf_rows(source_iso_massf)
    open(out_path, "w") do io
        for row in rows
            iso = format_iso_label(row.a, row.label_parts, source_iso_massf)
            @printf(io, "%3d %-5s         %16.10E\n", row.z, iso, row.abundance)
        end
    end
end

function copy_or_link(src, dst)
    if islink(src)
        target = realpath(src)
        symlink(relpath(target, dirname(dst)), dst)
    else
        cp(src, dst)
    end
end

function ensure_symlink(link_path, target)
    (ispath(link_path) || islink(link_path)) && return
    symlink(relpath(realpath(target), dirname(link_path)), link_path)
end

function copy_decay_inputs(source_run, decay_run)
    for name in TEXT_INPUTS
        src = joinpath(source_run, name)
        ispath(src) && copy_or_link(src, joinpath(decay_run, name))
    end

    input_files = sort(filter(name -> endswith(name, ".input"), readdir(source_run)))
    for name in input_files
        dst = joinpath(decay_run, name)
        ispath(dst) || copy_or_link(joinpath(source_run, name), dst)
    end

    ppn_exe = joinpath(source_run, "ppn.exe")
    isfile(ppn_exe) || error("missing ppn.exe in source run: $ppn_exe")
    symlink(relpath(realpath(ppn_exe), decay_run), joinpath(decay_run, "ppn.exe"))

    npdata = joinpath(dirname(source_run), "NPDATA")
    ispath(npdata) && ensure_symlink(joinpath(dirname(decay_run), "NPDATA"), npdata)
end

function update_namelist(text, replacements::Vector{Pair{String,String}})
    lines = split(text, '\n')
    !isempty(lines) && lines[end] == "" && pop!(lines)

    replacement_keys = Dict(lowercase(k) => v for (k, v) in replacements)
    found = Set{String}()
    output = String[]

    for line in lines
        stripped = strip(line)
        if stripped == "/"
            for (key, value) in replacements
                lowercase(key) in found || push!(output, "        $key = $value")
            end
            push!(output, line)
            continue
        end

        m = match(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
        if m !== nothing && haskey(replacement_keys, lowercase(m.captures[2]))
            key = m.captures[2]
            push!(found, lowercase(key))
            push!(output, "$(m.captures[1])$key = $(replacement_keys[lowercase(key)])")
        else
            push!(output, line)
        end
    end

    return join(output, '\n') * "\n"
end

fortran_float(value) = replace(@sprintf("%.10E", value), "E" => "d")

function patch_decay_inputs(decay_run, decay_time)
    frame = joinpath(decay_run, "ppn_frame.input")
    physics = joinpath(decay_run, "ppn_physics.input")
    isfile(frame) || error("missing ppn_frame.input in decay run: $frame")
    isfile(physics) || error("missing ppn_physics.input in decay run: $physics")

    write(
        frame,
        update_namelist(
            read(frame, String),
            [
                "nsource" => "0",
                "iabuini" => "11",
                "ini_filename" => "'post_abundance.DAT'",
                "iplot_flux_option" => "0",
                "i_flux_integrated" => "0",
            ],
        ),
    )
    write(
        physics,
        update_namelist(
            read(physics, String),
            [
                "decay" => ".true.",
                "decay_time" => fortran_float(decay_time),
                "detailed_balance" => ".false.",
            ],
        ),
    )
end

function build_jobs(runs_dir, decay_runs_dir, source_runs)
    jobs = DecayJob[]
    for source_run in source_runs
        final_cycle = final_cycle_from_xtime(source_run)
        source_iso_massf = joinpath(source_run, "iso_massf$(cycle_name(final_cycle)).DAT")
        isfile(source_iso_massf) || error("missing final iso_massf file: $source_iso_massf")
        rel = relpath(source_run, runs_dir)
        push!(jobs, DecayJob(source_run, joinpath(decay_runs_dir, rel), source_iso_massf, final_cycle))
    end
    return jobs
end

function build_decay_run(nova_case, job::DecayJob, decay_time)
    mkpath(job.decay_run)
    copy_decay_inputs(job.source_run, job.decay_run)
    write_post_abundance(job.source_iso_massf, joinpath(job.decay_run, "post_abundance.DAT"))
    patch_decay_inputs(job.decay_run, decay_time)
    open(joinpath(job.decay_run, "decay_source.txt"), "w") do io
        println(io, "nova_case = $nova_case")
        println(io, "source_run = $(job.source_run)")
        println(io, "source_iso_massf = $(job.source_iso_massf)")
        println(io, "source_final_cycle = $(job.final_cycle)")
        println(io, "decay_time_seconds = $(decay_time)")
        println(io, "generated = $(Dates.now())")
    end
end

function set_single_threaded_env!()
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    )
        ENV[key] = "1"
    end
end

function run_ppn_decay(job::DecayJob)
    log_path = joinpath(job.decay_run, "decay.log")
    process = open(log_path, "w") do io
        p = run(pipeline(Cmd(`./ppn.exe`, dir=job.decay_run), stdout=io, stderr=io), wait=false)
        wait(p)
        p
    end

    output = joinpath(job.decay_run, "iso_massfdecay.DAT")
    returncode = process.exitcode
    status = if returncode == 0 && isfile(output)
        "ok"
    elseif returncode == 0
        "missing_output"
    else
        "failed"
    end
    return DecayResult(job, returncode, output, log_path, status)
end

function csv_field(value)
    s = string(value)
    occursin(r"[,\"\n]", s) || return s
    "\"" * replace(s, "\"" => "\"\"") * "\""
end

function write_manifest(path, jobs::Vector{DecayJob}, results::Union{Nothing,Vector{DecayResult}}=nothing)
    result_by_run = Dict{String,DecayResult}()
    if results !== nothing
        for result in results
            result_by_run[result.job.decay_run] = result
        end
    end

    open(path, "w") do io
        println(io, "source_run,decay_run,source_iso_massf,final_cycle,output,log,returncode,status")
        for job in jobs
            result = get(result_by_run, job.decay_run, nothing)
            fields = if result === nothing
                (job.source_run, job.decay_run, job.source_iso_massf, job.final_cycle, "", "", "", "built")
            else
                (
                    job.source_run,
                    job.decay_run,
                    job.source_iso_massf,
                    job.final_cycle,
                    result.output,
                    result.log,
                    result.returncode,
                    result.status,
                )
            end
            println(io, join(csv_field.(fields), ","))
        end
    end
end

function run_jobs_parallel(jobs::Vector{DecayJob}, jobs_count::Int)
    queue = Channel{DecayJob}(length(jobs))
    for job in jobs
        put!(queue, job)
    end
    close(queue)

    results = Channel{DecayResult}(length(jobs))
    workers = Task[]
    for _ in 1:min(jobs_count, length(jobs))
        push!(workers, @async begin
            for job in queue
                put!(results, run_ppn_decay(job))
            end
        end)
    end
    foreach(wait, workers)
    close(results)
    return collect(results)
end

function main()
    opts = parse_args(ARGS)

    nova_case = NovaRunTools.nova_dir(opts.nova)
    isdir(nova_case) || error("nova_case does not exist: $nova_case")

    runs_dir = joinpath(nova_case, opts.runs_name)
    decay_runs_dir = joinpath(nova_case, opts.decay_runs_name)
    isdir(runs_dir) || error("runs directory does not exist: $runs_dir")

    source_runs = discover_source_runs(runs_dir)
    isempty(source_runs) && error("no ppn.exe files found below $runs_dir")
    jobs = build_jobs(runs_dir, decay_runs_dir, source_runs)

    println("found $(length(jobs)) source runs")
    println("decay time: $(opts.decay_time) seconds")
    println("runs dir: $runs_dir")
    println("decay runs dir: $decay_runs_dir")

    if opts.dry_run
        for job in jobs
            println("$(job.source_run) -> $(job.decay_run)")
        end
        return
    end

    if ispath(decay_runs_dir)
        println("WARNING: removing existing decay_runs directory: $decay_runs_dir")
        rm(decay_runs_dir; recursive=true)
    end
    mkpath(decay_runs_dir)

    root_npdata = joinpath(runs_dir, "NPDATA")
    ispath(root_npdata) && ensure_symlink(joinpath(decay_runs_dir, "NPDATA"), root_npdata)

    for job in jobs
        build_decay_run(nova_case, job, opts.decay_time)
    end

    manifest = joinpath(decay_runs_dir, "decay_manifest.csv")
    write_manifest(manifest, jobs)
    println("built $(length(jobs)) decay run directories")

    if opts.no_run
        println("wrote manifest: $manifest")
        println("--no-run requested; not executing ppn.exe")
        return
    end

    set_single_threaded_env!()
    results = run_jobs_parallel(jobs, opts.jobs)
    for (idx, result) in enumerate(results)
        rel = relpath(result.job.decay_run, decay_runs_dir)
        println("[$idx/$(length(jobs))] $(result.status): $rel")
    end

    sort!(results; by=result -> result.job.decay_run)
    write_manifest(manifest, jobs, results)
    failures = filter(result -> result.status != "ok", results)
    println("wrote manifest: $manifest")
    println("successful decays: $(length(results) - length(failures))")
    println("failed decays: $(length(failures))")
    if !isempty(failures)
        for result in first(failures, 20)
            println("  $(result.status): $(result.job.decay_run) log=$(result.log)")
        end
        exit(1)
    end
end

main()
