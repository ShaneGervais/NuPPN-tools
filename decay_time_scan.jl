#!/usr/bin/env julia

# Run real PPN decay for one source run across a grid of decay times, so a
# downstream script can pick whichever time best matches a reference baseline
# (e.g. Iliadis 2002 Table 4 JCH1). Shares its per-job building blocks with
# decay_ppn_sweep.jl, but sweeps decay time for one run instead of mirroring
# every run at one decay time.

include(joinpath(@__DIR__, "NovaRunTools.jl"))
using .NovaRunTools
using Dates
using Printf

const DEFAULT_DECAY_TIME_GRID_SECONDS = [
    0.0, 60.0, 300.0, 900.0, 1800.0,
    3600.0, 7200.0, 14400.0, 28800.0, 43200.0,
    86400.0, 604800.0, 2592000.0, 15552000.0, 31557600.0, 94672800.0,
]
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

struct TimeJob
    decay_time::Float64
    decay_run::String
end

mutable struct TimeResult
    job::TimeJob
    returncode::Int
    output::String
    log::String
    status::String
end

Base.@kwdef mutable struct ScanOptions
    nova::String = "nova_test"
    run_name::String = "baseline"
    out_name::String = "decay_time_scan"
    times::Vector{Float64} = DEFAULT_DECAY_TIME_GRID_SECONDS
    jobs::Int = get(ENV, "MAX_JOBS", nothing) === nothing ? 4 : parse(Int, ENV["MAX_JOBS"])
    dry_run::Bool = false
end

function usage()
    println("""
Usage:
  julia tools/decay_time_scan.jl [options]

Run real PPN decay for one source run (default: runs/baseline) across a grid
of decay times, so a downstream comparison script (e.g. baseline_decay_checker.py)
can pick whichever decay time best matches a reference baseline.

Options:
  --nova NAME          Nova directory under nova_cases/ (default: nova_test)
  --run-name NAME      Source run under runs/ to decay (default: baseline)
  --out-name NAME      Output directory under the nova case (default: decay_time_scan)
  --times s1,s2,...    Comma-separated decay times in seconds (default: a 16-point
                       grid from 0 s to 3 years)
  -j, --jobs N         Number of ppn.exe decay runs to execute in parallel (default: MAX_JOBS or 4)
  --dry-run            Print the (time, decay_run) pairs that would be built without writing files
  -h, --help           Show this help
""")
end

function parse_times(text)
    seconds = Float64[]
    for chunk in split(text, ',')
        chunk = strip(chunk)
        isempty(chunk) && continue
        push!(seconds, parse(Float64, chunk))
    end
    isempty(seconds) && error("--times must contain at least one value")
    return seconds
end

function parse_args(args)
    opts = ScanOptions()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            opts.nova = args[i]
        elseif arg == "--run-name"
            i += 1
            i <= length(args) || error("--run-name requires a value")
            opts.run_name = args[i]
            isempty(strip(opts.run_name)) && error("--run-name cannot be empty")
        elseif arg == "--out-name"
            i += 1
            i <= length(args) || error("--out-name requires a value")
            opts.out_name = args[i]
            isempty(strip(opts.out_name)) && error("--out-name cannot be empty")
        elseif arg == "--times"
            i += 1
            i <= length(args) || error("--times requires a value")
            opts.times = parse_times(args[i])
        elseif arg in ("-j", "--jobs")
            i += 1
            i <= length(args) || error("--jobs requires a value")
            opts.jobs = parse(Int, args[i])
            opts.jobs >= 1 || error("--jobs must be >= 1")
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

function format_decay_time(seconds)
    seconds == 0.0 && return "0 s"
    minute, hour, day, year = 60.0, 3600.0, 86400.0, 31557600.0
    if seconds < hour
        return @sprintf("%.3g min", seconds / minute)
    elseif seconds < day
        return @sprintf("%.3g h", seconds / hour)
    elseif seconds < year
        return @sprintf("%.3g d", seconds / day)
    end
    return @sprintf("%.3g yr", seconds / year)
end

time_dir_name(seconds) = @sprintf("t%010ds", round(Int, seconds))

function build_jobs(out_dir, times)
    jobs = TimeJob[]
    for seconds in times
        push!(jobs, TimeJob(seconds, joinpath(out_dir, time_dir_name(seconds))))
    end
    return jobs
end

function build_time_run(source_run, source_iso_massf, job::TimeJob)
    mkpath(job.decay_run)
    copy_decay_inputs(source_run, job.decay_run)
    write_post_abundance(source_iso_massf, joinpath(job.decay_run, "post_abundance.DAT"))
    patch_decay_inputs(job.decay_run, job.decay_time)
    open(joinpath(job.decay_run, "decay_source.txt"), "w") do io
        println(io, "source_run = $source_run")
        println(io, "source_iso_massf = $source_iso_massf")
        println(io, "decay_time_seconds = $(job.decay_time)")
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

function run_ppn_decay(job::TimeJob)
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
    return TimeResult(job, returncode, output, log_path, status)
end

function csv_field(value)
    s = string(value)
    occursin(r"[,\"\n]", s) || return s
    "\"" * replace(s, "\"" => "\"\"") * "\""
end

function write_manifest(path, jobs::Vector{TimeJob}, results::Union{Nothing,Vector{TimeResult}}=nothing)
    result_by_run = Dict{String,TimeResult}()
    if results !== nothing
        for result in results
            result_by_run[result.job.decay_run] = result
        end
    end

    open(path, "w") do io
        println(io, "decay_time_seconds,decay_time_label,decay_run,output,log,returncode,status")
        for job in jobs
            result = get(result_by_run, job.decay_run, nothing)
            fields = if result === nothing
                (job.decay_time, format_decay_time(job.decay_time), job.decay_run, "", "", "", "built")
            else
                (
                    job.decay_time,
                    format_decay_time(job.decay_time),
                    job.decay_run,
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

function run_jobs_parallel(jobs::Vector{TimeJob}, jobs_count::Int)
    queue = Channel{TimeJob}(length(jobs))
    for job in jobs
        put!(queue, job)
    end
    close(queue)

    results = Channel{TimeResult}(length(jobs))
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

    source_run = joinpath(nova_case, "runs", opts.run_name)
    isdir(source_run) || error("source run does not exist: $source_run")

    final_cycle = final_cycle_from_xtime(source_run)
    source_iso_massf = joinpath(source_run, "iso_massf$(cycle_name(final_cycle)).DAT")
    isfile(source_iso_massf) || error("missing final iso_massf file: $source_iso_massf")

    out_dir = joinpath(nova_case, opts.out_name)
    jobs = build_jobs(out_dir, opts.times)

    println("source run: $source_run")
    println("final cycle: $final_cycle")
    println("decay times: $(length(jobs))")
    println("output dir: $out_dir")

    if opts.dry_run
        for job in jobs
            println("$(format_decay_time(job.decay_time)) ($(job.decay_time) s) -> $(job.decay_run)")
        end
        return
    end

    if ispath(out_dir)
        println("WARNING: removing existing scan directory: $out_dir")
        rm(out_dir; recursive=true)
    end
    mkpath(out_dir)

    for job in jobs
        build_time_run(source_run, source_iso_massf, job)
    end

    manifest = joinpath(out_dir, "decay_time_scan_manifest.csv")
    write_manifest(manifest, jobs)
    println("built $(length(jobs)) decay time directories")

    set_single_threaded_env!()
    results = run_jobs_parallel(jobs, opts.jobs)
    for (idx, result) in enumerate(results)
        println("[$idx/$(length(jobs))] $(result.status): $(format_decay_time(result.job.decay_time))")
    end

    sort!(results; by=result -> result.job.decay_time)
    write_manifest(manifest, jobs, results)
    successes = filter(result -> result.status == "ok", results)
    failures = filter(result -> result.status != "ok", results)
    println("wrote manifest: $manifest")
    println("successful decay times: $(length(successes))")
    println("failed decay times: $(length(failures))")
    if !isempty(failures)
        for result in failures
            println("  $(result.status): $(format_decay_time(result.job.decay_time)) log=$(result.log)")
        end
    end
    isempty(successes) && error("all decay times failed; see logs above")
end

main()
