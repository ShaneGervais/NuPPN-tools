#!/usr/bin/env julia

using Dates
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_NOVA = "nova_test"
const SECONDS_PER_MINUTE = 60.0
const SECONDS_PER_HOUR = 60.0 * SECONDS_PER_MINUTE
const SECONDS_PER_DAY = 24.0 * SECONDS_PER_HOUR
const SECONDS_PER_WEEK = 7.0 * SECONDS_PER_DAY
const SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY
const DEFAULT_DECAY_SECONDS = SECONDS_PER_WEEK

Base.@kwdef mutable struct Options
    nova::String = DEFAULT_NOVA
    decay_seconds::Float64 = DEFAULT_DECAY_SECONDS
    jobs::Int = get(ENV, "MAX_JOBS", nothing) === nothing ? 4 : parse(Int, ENV["MAX_JOBS"])
    run::Bool = false
end

function usage()
    println("""
Usage:
  julia --project=NovaProject_0.2/NovaJL NovaProject_0.2/tools/prepare_n_run_decay.jl [options]

Options:
  --nova NAME              Nova directory under novae/ (default: nova_test)
  --decay-time DURATION    Decay time, e.g. 1y3w9d18h7m2s
  --decay-seconds VALUE    Decay time in seconds (default: 604800, one week)
  --jobs N                 Maximum concurrent decay jobs when --run is used
  --run                    Run all prepared decay ppn.exe jobs
  -h, --help               Show this help

Assumptions:
  novae/<nova>/decay_ppn already exists and has been manually prepared from ppn/.
  The script only updates decay_ppn/ppn_physics.input for decay mode and decay_time.
  If novae/<nova>/decay_runs exists, it is deleted before preparing the new decay sweep.
""")
end

function parse_duration_seconds(text)
    s = lowercase(strip(text))
    isempty(s) && error("--decay-time cannot be empty")
    total = 0.0
    consumed = falses(length(s))
    found = false
    for m in eachmatch(r"([0-9]+(?:\.[0-9]+)?)(y|w|d|h|m|s)", s)
        found = true
        value = parse(Float64, m.captures[1])
        unit = m.captures[2]
        multiplier = if unit == "y"
            SECONDS_PER_YEAR
        elseif unit == "w"
            SECONDS_PER_WEEK
        elseif unit == "d"
            SECONDS_PER_DAY
        elseif unit == "h"
            SECONDS_PER_HOUR
        elseif unit == "m"
            SECONDS_PER_MINUTE
        else
            1.0
        end
        total += value * multiplier
        for i in m.offset:(m.offset + length(m.match) - 1)
            consumed[i] = true
        end
    end
    found || error("could not parse --decay-time=$(repr(text)); expected e.g. 1w or 1y3w9d18h7m2s")
    all(consumed) || error("invalid --decay-time=$(repr(text)); expected only number+unit chunks with units y,w,d,h,m,s")
    total > 0 || error("--decay-time must be > 0")
    return total
end

function parse_args(args)
    opts = Options()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            opts.nova = args[i]
        elseif arg == "--decay-seconds"
            i += 1
            i <= length(args) || error("--decay-seconds requires a value")
            opts.decay_seconds = parse(Float64, args[i])
            opts.decay_seconds > 0 || error("--decay-seconds must be > 0")
        elseif arg == "--decay-time"
            i += 1
            i <= length(args) || error("--decay-time requires a value")
            opts.decay_seconds = parse_duration_seconds(args[i])
        elseif arg == "--jobs"
            i += 1
            i <= length(args) || error("--jobs requires a value")
            opts.jobs = parse(Int, args[i])
            opts.jobs >= 1 || error("--jobs must be >= 1")
        elseif arg == "--run"
            opts.run = true
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

function nova_dir(name)
    joinpath(PROJECT_ROOT, "novae", name)
end

function replace_or_insert_namelist_value(text, key, value)
    lines = split(text, '\n'; keepempty=true)
    pattern = Regex("^(\\s*$(key)\\s*=\\s*)[^!\\n]*(.*)\$")
    for i in eachindex(lines)
        m = match(pattern, lines[i])
        m === nothing && continue
        lines[i] = string(m.captures[1], value, m.captures[2])
        return join(lines, '\n')
    end
    for i in eachindex(lines)
        if occursin(r"^\s*/\s*$", lines[i])
            insert!(lines, i, "        $key = $value")
            return join(lines, '\n')
        end
    end
    error("could not find namelist terminator while inserting $key")
end

function fortran_float(value)
    replace(@sprintf("%.10e", value), "e" => "d")
end

function update_decay_template!(template_dir, decay_seconds)
    physics_path = joinpath(template_dir, "ppn_physics.input")
    isfile(physics_path) || error("missing $physics_path")

    physics = read(physics_path, String)
    physics = replace_or_insert_namelist_value(physics, "decay", ".true.")
    physics = replace_or_insert_namelist_value(physics, "decay_time", fortran_float(decay_seconds))
    physics = replace_or_insert_namelist_value(physics, "detailed_balance", ".false.")
    write(physics_path, physics)
end

function prepare_decay_template(base, opts)
    template_dir = joinpath(base, "decay_ppn")
    isdir(template_dir) || error("missing manually prepared decay_ppn template: $template_dir")
    isfile(joinpath(template_dir, "ppn.exe")) || error("missing ppn.exe in decay template: $template_dir")
    isfile(joinpath(template_dir, "ppn_frame.input")) || error("missing ppn_frame.input in decay template: $template_dir")
    update_decay_template!(template_dir, opts.decay_seconds)
    return template_dir
end

function parse_iso_massf(path)
    rows = NamedTuple[]
    for line in eachline(path)
        startswith(strip(line), "#") && continue
        parts = split(line)
        length(parts) < 6 && continue
        try
            z = round(Int, parse(Float64, parts[2]))
            a = round(Int, parse(Float64, parts[3]))
            x = parse(Float64, parts[5])
            iso_tokens = parts[6:end]
            iso = join(iso_tokens, " ")
            push!(rows, (z=z, a=a, x=x, iso=iso))
        catch
            continue
        end
    end
    isempty(rows) && error("no isotope rows parsed from $path")
    return rows
end

function iso_field_for_initial(row)
    if row.iso == "PROT" || row.iso == "NEUT"
        return row.iso
    end
    parts = split(row.iso)
    if length(parts) == 1
        m = match(r"^([A-Za-z]+)(\d+)$", parts[1])
        m === nothing && return uppercase(parts[1])
        return @sprintf("%-2s%3d", lowercase(m.captures[1]), parse(Int, m.captures[2]))
    end
    return @sprintf("%-2s%3d", lowercase(parts[1]), parse(Int, parts[2]))
end

function write_post_abundance(path, rows)
    open(path, "w") do io
        for row in rows
            iso = iso_field_for_initial(row)
            @printf(io, "%3d %-5s         %.10E\n", row.z, iso, row.x)
        end
    end
end

function final_iso_massf(run_dir)
    files = String[]
    for name in readdir(run_dir)
        occursin(r"^iso_massf\d+\.DAT$", name) && push!(files, joinpath(run_dir, name))
    end
    isempty(files) && return nothing
    sort!(files, by=path -> parse(Int, match(r"iso_massf(\d+)\.DAT$", basename(path)).captures[1]))
    return last(files)
end

function list_completed_runs(runs_dir)
    run_dirs = String[]
    for (root, _, files) in walkdir(runs_dir)
        "ppn.exe" in files || continue
        iso = final_iso_massf(root)
        iso === nothing && continue
        push!(run_dirs, root)
    end
    sort!(run_dirs)
    return run_dirs
end

function sync_template_to_run(template_dir, decay_run_dir)
    if !isdir(decay_run_dir)
        mkpath(dirname(decay_run_dir))
        cp(template_dir, decay_run_dir)
        return
    end
    for name in readdir(template_dir)
        src = joinpath(template_dir, name)
        dst = joinpath(decay_run_dir, name)
        if isdir(src)
            cp(src, dst; force=true, follow_symlinks=true)
        else
            cp(src, dst; force=true)
        end
    end
end

function ensure_npdata_link(parent_dir, npdata_source)
    link_path = joinpath(parent_dir, "NPDATA")
    if ispath(link_path)
        return
    end
    symlink(npdata_source, link_path)
end

function prepare_decay_runs(base, template_dir)
    runs_dir = joinpath(base, "runs")
    decay_runs_dir = joinpath(base, "decay_runs")
    npdata_source = joinpath(base, "NPDATA")
    isdir(runs_dir) || error("missing runs directory: $runs_dir")
    isdir(npdata_source) || error("missing NPDATA directory: $npdata_source")
    if ispath(decay_runs_dir)
        println("Deleting existing decay run sweep: $decay_runs_dir")
        rm(decay_runs_dir; recursive=true)
    end
    mkpath(decay_runs_dir)
    ensure_npdata_link(decay_runs_dir, npdata_source)

    prepared = 0
    for run_dir in list_completed_runs(runs_dir)
        rel = relpath(run_dir, runs_dir)
        decay_run_dir = joinpath(decay_runs_dir, rel)
        mkpath(decay_run_dir)
        sync_template_to_run(template_dir, decay_run_dir)
        ensure_npdata_link(dirname(decay_run_dir), npdata_source)
        ensure_npdata_link(decay_run_dir, npdata_source)

        iso_path = final_iso_massf(run_dir)
        rows = parse_iso_massf(iso_path)
        post_path = joinpath(decay_run_dir, "post_abundance.DAT")
        write_post_abundance(post_path, rows)

        open(joinpath(decay_run_dir, "decay_source.txt"), "w") do io
            println(io, "source_run = $run_dir")
            println(io, "source_iso_massf = $iso_path")
            println(io, "generated = $(Dates.now())")
        end

        prepared += 1
    end
    return prepared
end

function relative_run_name(decay_runs_dir, run_dir)
    replace(relpath(run_dir, decay_runs_dir), '/' => '_')
end

function run_and_log(cmd::Cmd, logfile)
    open(logfile, "a") do io
        process = run(pipeline(cmd, stdout=io, stderr=io), wait=false)
        wait(process)
        success(process) || error("Command failed: $cmd")
    end
end

function run_one_decay(run_dir, decay_runs_dir, logs_dir)
    name = relative_run_name(decay_runs_dir, run_dir)
    logfile = joinpath(logs_dir, "$name.log")
    open(logfile, "w") do io
        println(io, "===========================================")
        println(io, "Running decay PPN in $run_dir")
        println(io, "Started $(Dates.now())")
        println(io, "===========================================")
    end
    println("Running decay $name")
    run_and_log(Cmd(`./ppn.exe`, dir=run_dir), logfile)
    isfile(joinpath(run_dir, "iso_massfdecay.DAT")) || error("decay output missing for $run_dir")
    open(logfile, "a") do io
        println(io, "Finished $(Dates.now())")
    end
end

function run_decay_jobs(base, jobs)
    decay_runs_dir = joinpath(base, "decay_runs")
    logs_dir = joinpath(base, "logs", "decay_ppn")
    mkpath(logs_dir)
    run_dirs = String[]
    for (root, _, files) in walkdir(decay_runs_dir)
        "ppn.exe" in files && "post_abundance.DAT" in files && push!(run_dirs, root)
    end
    # If no runs found, check if it's because the template was copied but not the executable
    if isempty(run_dirs)
        println("Warning: No prepared decay run directories found with ppn.exe and post_abundance.DAT.")
        println("Checking contents of $decay_runs_dir:")
        for (root, dirs, files) in walkdir(decay_runs_dir)
            println("  $root: $files")
        end
    end
    sort!(run_dirs)
    isempty(run_dirs) && error("no prepared decay run directories under $decay_runs_dir")

    queue = Channel{String}(length(run_dirs))
    for run_dir in run_dirs
        put!(queue, run_dir)
    end
    close(queue)

    failures = Channel{Tuple{String,Any}}(length(run_dirs))
    workers = Task[]
    for _ in 1:min(jobs, length(run_dirs))
        push!(workers, @async begin
            for run_dir in queue
                try
                    run_one_decay(run_dir, decay_runs_dir, logs_dir)
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
        println("Failures:")
        for (run_dir, err) in failed
            println("  $run_dir")
            println("    $err")
        end
        error("$(length(failed)) decay jobs failed")
    end
    return length(run_dirs)
end

function main()
    opts = parse_args(ARGS)
    base = nova_dir(opts.nova)

    template_dir = prepare_decay_template(base, opts)
    prepared = prepare_decay_runs(base, template_dir)

    println("Updated decay template: $template_dir")
    println("Decay time: $(opts.decay_seconds) seconds")
    println("Prepared $prepared decay run(s).")

    if opts.run
        count = run_decay_jobs(base, opts.jobs)
        println("Completed $count decay run(s).")
    else
        println("Run decay jobs with: julia --project=NovaProject_0.2/NovaJL NovaProject_0.2/tools/prepare_n_run_decay.jl --nova $(opts.nova) --run")
    end
end

main()
