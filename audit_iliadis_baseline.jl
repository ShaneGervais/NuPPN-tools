#!/usr/bin/env julia

# Compare the local nova_test baseline model against the JCH1 reference
# quantities in Iliadis et al. (2002), Tables 1, 2, and 4.

using CSV
using DataFrames
using Dates
using Statistics

include(joinpath(@__DIR__, "NovaRunTools.jl"))
const NRT = NovaRunTools

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_NOVA = "nova_test"
const ARTICLE_PDF = "/home/sgervais/Documents/PhD/articles/Iliadis-ApJS-2002.pdf"
const SECONDS_PER_YEAR = 365.25 * 24 * 3600

const JCH1_PROPERTIES = Dict(
    "WD mass (Msun)" => 1.15,
    "mixing_percent" => 50.0,
    "Tpeak_T9" => 0.231,
    "Lpeak_1e4_Lsun" => 26.0,
    "Macc_1e-5_Msun" => 3.2,
    "Mdot_1e-10_Msun_per_yr" => 2.0,
    "Mej_1e-5_Msun" => 2.6,
)

const JCH1_INITIAL = Dict(
    "H-1" => 3.5e-1, "H-2" => 2.4e-5, "HE-3" => 1.5e-5, "HE-4" => 1.4e-1,
    "LI-6" => 3.2e-10, "LI-7" => 4.7e-9, "BE-9" => 8.3e-11, "B-10" => 5.3e-10,
    "B-11" => 2.4e-9, "C-12" => 6.1e-3, "C-13" => 1.8e-5, "N-14" => 5.5e-4,
    "N-15" => 2.2e-6, "O-16" => 2.6e-1, "O-17" => 1.9e-6, "O-18" => 1.1e-5,
    "F-19" => 2.0e-7, "NE-20" => 1.6e-1, "NE-21" => 3.0e-3, "NE-22" => 2.2e-3,
    "NA-23" => 3.2e-2, "MG-24" => 2.8e-2, "MG-25" => 7.9e-3, "MG-26" => 5.0e-3,
    "AL-27" => 5.4e-3, "SI-28" => 3.3e-4, "SI-29" => 1.7e-5, "SI-30" => 1.2e-5,
    "P-31" => 1.1e-6, "S-32" => 2.0e-4, "S-33" => 4.5e-7, "S-34" => 2.6e-6,
    "CL-35" => 4.9e-7, "CL-37" => 1.7e-7, "AR-36" => 3.9e-5, "AR-38" => 1.9e-6,
    "K-39" => 4.8e-7, "CA-40" => 3.0e-5,
)

const JCH1_FINAL = Dict(
    "H-1" => 1.6e-1, "HE-3" => 0.0, "HE-4" => 3.6e-1, "BE-7" => 1.9e-10,
    "C-12" => 3.5e-2, "C-13" => 6.6e-2, "N-14" => 1.2e-1, "N-15" => 2.5e-5,
    "O-16" => 3.4e-3, "O-17" => 1.7e-3, "O-18" => 1.7e-7, "F-18" => 2.1e-6,
    "F-19" => 7.5e-9, "NE-20" => 1.6e-1, "NE-21" => 4.2e-6, "NE-22" => 2.8e-4,
    "NA-22" => 1.1e-4, "NA-23" => 1.6e-4, "MG-24" => 8.1e-6, "MG-25" => 8.0e-4,
    "MG-26" => 3.1e-5, "AL-26" => 9.9e-5, "AL-27" => 6.5e-4, "SI-28" => 6.9e-2,
    "SI-29" => 8.7e-4, "SI-30" => 1.1e-2, "P-31" => 5.2e-3, "S-32" => 4.2e-3,
    "S-33" => 9.3e-7, "S-34" => 1.3e-6, "CL-35" => 2.2e-6, "CL-37" => 1.5e-7,
    "AR-36" => 1.2e-5, "AR-37" => 2.6e-5, "AR-38" => 3.4e-6, "K-39" => 6.6e-7,
    "CA-40" => 3.0e-5, "CA-41" => 4.2e-9, "CA-42" => 0.0,
)

const SHORT_LIVED_DECAYS = Dict(
    "N-13" => "C-13",
    "O-14" => "N-14",
    "O-15" => "N-15",
    "F-17" => "O-17",
    "NE-19" => "F-19",
    "NA-21" => "NE-21",
    "MG-23" => "NA-23",
    "AL-25" => "MG-25",
    "SI-27" => "AL-27",
    "P-29" => "SI-29",
    "P-30" => "SI-30",
    "S-31" => "P-31",
    "CL-34" => "S-34",
    "AR-35" => "CL-35",
    "K-37" => "AR-37",
    "K-38" => "AR-38",
    "CA-39" => "K-39",
    "SC-40" => "CA-40",
    "SC-41" => "CA-41",
)

const EVOLUTION_ISOTOPES = [
    "H-1", "HE-4", "BE-7", "C-12", "C-13", "N-14", "N-15", "O-16", "O-17",
    "F-18", "NE-20", "NE-21", "NE-22", "NA-22", "NA-23", "MG-25", "AL-26",
    "SI-28", "SI-30", "P-31",
]

function parse_args(args)
    nova = DEFAULT_NOVA
    i = 1
    while i <= length(args)
        if args[i] == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            nova = args[i]
        elseif args[i] in ("-h", "--help")
            println("Usage: julia --project=NovaProject_0.2/NovaJL NovaProject_0.2/tools/audit_iliadis_baseline.jl [--nova nova_test]")
            exit(0)
        else
            error("unknown argument $(args[i])")
        end
        i += 1
    end
    return nova
end

function normalize_isotope(el, a)
    el = uppercase(String(el))
    el == "PROT" && return "H-1"
    el == "NEUT" && return "NEUT"
    return "$(el)-$(parse(Int, string(a)))"
end

function parse_initial_abundance(path)
    values = Dict{String,Float64}()
    for line in eachline(path)
        parts = split(strip(line))
        length(parts) < 3 && continue
        try
            z = parse(Int, parts[1])
            isotope = z == 1 && uppercase(parts[2]) == "PROT" ? "H-1" : normalize_isotope(parts[2], parts[3])
            x = parse(Float64, parts[end])
            values[isotope] = x
        catch
            continue
        end
    end
    return values
end

function parse_final_abundances(path)
    df = CSV.read(path, DataFrame; missingstring=["", "missing"])
    values = Dict{String,Float64}()
    for row in eachrow(df)
        iso = string(row.isotope)
        iso == "PROT" && (iso = "H-1")
        iso == "NEUT" && (iso = "NEUT")
        values[iso] = parse(Float64, strip(string(row.X)))
    end
    return values
end

function normalize_xtime_isotope(text)
    s = strip(replace(text, r"^\d+-" => ""))
    s == "PROT" && return "H-1"
    s == "NEUT" && return "NEUT"
    parts = split(s)
    if length(parts) == 2
        return "$(uppercase(parts[1]))-$(parse(Int, parts[2]))"
    end
    m = match(r"^([A-Z]+)\s*(\d+)$", replace(s, " " => ""))
    m === nothing && return s
    return "$(m.captures[1])-$(parse(Int, m.captures[2]))"
end

function parse_xtime(path; isotopes=EVOLUTION_ISOTOPES)
    isfile(path) || error("missing x-time output: $path")
    wanted = Set(isotopes)
    header = nothing
    for line in eachline(path)
        if startswith(line, "#|")
            header = line
            break
        end
    end
    header === nothing && error("could not find x-time header in $path")

    chunks = strip.(split(header, "|")[2:end])
    species_cols = Dict{String,Int}()
    for (i, chunk) in enumerate(chunks)
        i <= 6 && continue
        iso = normalize_xtime_isotope(chunk)
        iso in wanted && (species_cols[iso] = i)
    end

    cols = Dict{Symbol,Vector{Float64}}(
        :cycle => Float64[],
        :time_yr => Float64[],
        :time_s => Float64[],
        :temperature_T9 => Float64[],
        :density_cgs => Float64[],
        :ye => Float64[],
    )
    for iso in sort(collect(wanted))
        cols[Symbol(replace(iso, "-" => "_"))] = Float64[]
    end

    for line in eachline(path)
        startswith(line, "#") && continue
        parts = split(strip(line))
        length(parts) < length(chunks) && continue
        try
            time_yr = parse(Float64, parts[2])
            push!(cols[:cycle], parse(Float64, parts[1]))
            push!(cols[:time_yr], time_yr)
            push!(cols[:time_s], time_yr * SECONDS_PER_YEAR)
            push!(cols[:temperature_T9], parse(Float64, parts[3]))
            push!(cols[:density_cgs], parse(Float64, parts[4]))
            push!(cols[:ye], parse(Float64, parts[6]))
            for iso in sort(collect(wanted))
                col = Symbol(replace(iso, "-" => "_"))
                idx = get(species_cols, iso, 0)
                push!(cols[col], idx == 0 ? 0.0 : parse(Float64, parts[idx]))
            end
        catch
            continue
        end
    end
    return DataFrame(cols)
end

function parse_trajectory(path)
    time = Float64[]
    temp = Float64[]
    rho = Float64[]
    ageunit = "SEC"
    tunit = "T9K"
    rhounit = "CGS"
    for line in eachline(path)
        line = strip(line)
        isempty(line) && continue
        startswith(line, "#") && continue
        if occursin("=", line)
            key, value = strip.(split(line, "=", limit=2))
            key = uppercase(key)
            value = uppercase(value)
            key == "AGEUNIT" && (ageunit = value)
            key == "TUNIT" && (tunit = value)
            key == "RHOUNIT" && (rhounit = value)
            continue
        end
        parts = split(line)
        length(parts) < 3 && continue
        try
            t = parse(Float64, parts[1])
            T = parse(Float64, parts[2])
            r = parse(Float64, parts[3])
            ageunit in ("YRS", "YR", "YEAR", "YEARS") && (t *= 365.25 * 24 * 3600)
            tunit == "T8K" && (T /= 10)
            rhounit == "LOG" && (r = 10.0^r)
            push!(time, t); push!(temp, T); push!(rho, r)
        catch
            continue
        end
    end
    return DataFrame(time_s=time, temperature_T9=temp, density_cgs=rho)
end

function isotope_value(row, isotope)
    col = Symbol(replace(isotope, "-" => "_"))
    hasproperty(row, col) ? getproperty(row, col) : 0.0
end

function first_below_xtime(xtime, isotope, threshold)
    col = Symbol(replace(isotope, "-" => "_"))
    hasproperty(xtime, col) || return missing
    idx = findfirst(<(threshold), xtime[!, col])
    idx === nothing && return missing
    return xtime.time_s[idx] - first(xtime.time_s)
end

function proton_consumption_table(xtime, final_decay)
    ipeak = argmax(xtime.temperature_T9)
    stages = [
        ("initial", 1),
        ("Tpeak", ipeak),
        ("final", nrow(xtime)),
    ]
    rows = NamedTuple[]
    for (stage, idx) in stages
        h = xtime.H_1[idx]
        push!(rows, (
            metric="H-1 at $stage",
            local_value=h,
            iliadis_jch1=stage == "final" ? JCH1_FINAL["H-1"] : missing,
            local_over_iliadis=stage == "final" ? h / JCH1_FINAL["H-1"] : missing,
            time_since_start_s=xtime.time_s[idx] - first(xtime.time_s),
            temperature_T9=xtime.temperature_T9[idx],
            notes=stage == "Tpeak" ? "H-1 remaining at peak trajectory temperature" : "",
        ))
    end
    for threshold in (0.16, 0.1, 0.01, 1e-6)
        t = first_below_xtime(xtime, "H-1", threshold)
        push!(rows, (
            metric="first H-1 below $threshold",
            local_value=threshold,
            iliadis_jch1=missing,
            local_over_iliadis=missing,
            time_since_start_s=t,
            temperature_T9=missing,
            notes=ismissing(t) ? "threshold not crossed" : "threshold crossing time from x-time.dat",
        ))
    end
    push!(rows, (
        metric="final decay-collapsed H-1",
        local_value=get(final_decay, "H-1", 0.0),
        iliadis_jch1=JCH1_FINAL["H-1"],
        local_over_iliadis=get(final_decay, "H-1", 0.0) / JCH1_FINAL["H-1"],
        time_since_start_s=last(xtime.time_s) - first(xtime.time_s),
        temperature_T9=last(xtime.temperature_T9),
        notes="decay collapse does not replenish H-1",
    ))
    return DataFrame(rows)
end

function evolution_table(xtime, final_decay)
    ipeak = argmax(xtime.temperature_T9)
    points = [("initial", 1), ("Tpeak", ipeak), ("final", nrow(xtime))]
    rows = NamedTuple[]
    for isotope in EVOLUTION_ISOTOPES
        ref = get(JCH1_FINAL, isotope, missing)
        final_x = get(final_decay, isotope, 0.0)
        for (stage, idx) in points
            x = isotope_value(xtime[idx, :], isotope)
            push!(rows, (
                isotope=isotope,
                stage=stage,
                local_value=x,
                final_decay_collapsed=final_x,
                iliadis_jch1_final=ref,
                final_local_over_iliadis=ismissing(ref) || ref == 0 || final_x <= 0 ? missing : final_x / ref,
                time_since_start_s=xtime.time_s[idx] - first(xtime.time_s),
                temperature_T9=xtime.temperature_T9[idx],
            ))
        end
    end
    return DataFrame(rows)
end

function decay_convention_table(raw_values, decay_values)
    rows = NamedTuple[]
    for isotope in sort(collect(keys(JCH1_FINAL)))
        ref = JCH1_FINAL[isotope]
        raw = get(raw_values, isotope, 0.0)
        collapsed = get(decay_values, isotope, 0.0)
        raw_ratio = ref > 0 && raw > 0 ? raw / ref : missing
        collapsed_ratio = ref > 0 && collapsed > 0 ? collapsed / ref : missing
        raw_err = ismissing(raw_ratio) ? missing : abs(log10(raw_ratio))
        collapsed_err = ismissing(collapsed_ratio) ? missing : abs(log10(collapsed_ratio))
        improvement = ismissing(raw_err) || ismissing(collapsed_err) ? missing : raw_err - collapsed_err
        push!(rows, (
            isotope=isotope,
            raw_local=raw,
            decay_collapsed_local=collapsed,
            iliadis_jch1=ref,
            raw_over_iliadis=raw_ratio,
            decay_collapsed_over_iliadis=collapsed_ratio,
            abs_log10_error_raw=raw_err,
            abs_log10_error_decay_collapsed=collapsed_err,
            decay_improvement_log10=improvement,
        ))
    end
    df = DataFrame(rows)
    sort!(df, [:decay_improvement_log10], rev=true)
    return df
end

function network_rate_table(network_path)
    isfile(network_path) || error("missing networksetup: $network_path")
    rows = NRT.parse_networksetup(network_path)
    records = NamedTuple[]
    for row in rows
        target_a = row.reactant === nothing ? missing : row.reactant[1]
        max_a = maximum(skipmissing([
            row.reactant === nothing ? missing : row.reactant[1],
            row.projectile === nothing ? missing : row.projectile[1],
            row.product_1 === nothing ? missing : row.product_1[1],
            row.product_2 === nothing ? missing : row.product_2[1],
        ]))
        bucket = ismissing(target_a) ? "unknown" : target_a < 20 ? "A<20" : target_a <= 40 ? "20<=A<=40" : "A>40"
        controlled = row.rtype in NRT.CONTROLLED_RTYPES
        active_crosses_a40 = controlled && row.active && (target_a > 40 || max_a > 40)
        push!(records, (
            index=row.index,
            active=row.active,
            source=row.source,
            rtype=row.rtype,
            target_a=target_a,
            max_a=max_a,
            target_bucket=bucket,
            controlled_charged_particle=controlled,
            active_within_a40=controlled && row.active ? (target_a <= 40 && max_a <= 40) : missing,
            active_crosses_a40=active_crosses_a40,
            reaction=strip(row.line),
        ))
    end
    return DataFrame(records)
end

function network_rate_summary(network_rates, recreation_audit_path)
    rows = combine(groupby(network_rates, [:active, :source, :rtype, :target_bucket]), nrow => :count)
    sort!(rows, [:active, :source, :rtype, :target_bucket])
    if isfile(recreation_audit_path)
        audit = CSV.read(recreation_audit_path, DataFrame)
        bad = audit[in.(audit.status, Ref(["fallback", "isomer_split_unresolved", "local_surrogate", "missing_network_row", "source_mismatch"])), :]
        rows.recreation_audit_open_items .= nrow(bad)
    else
        rows.recreation_audit_open_items .= missing
    end
    return rows
end

function decay_collapsed(values)
    out = copy(values)
    for (parent, daughter) in SHORT_LIVED_DECAYS
        x = get(out, parent, 0.0)
        x <= 0 && continue
        out[daughter] = get(out, daughter, 0.0) + x
        out[parent] = 0.0
    end
    return out
end

function comparison_table(local_values, reference_values; local_label="local")
    rows = NamedTuple[]
    for isotope in sort(collect(keys(reference_values)))
        ref = reference_values[isotope]
        local_x = get(local_values, isotope, 0.0)
        ratio = ref > 0 && local_x > 0 ? local_x / ref : missing
        abslog = ratio === missing ? missing : abs(log10(ratio))
        push!(rows, (
            isotope=isotope,
            local_value=local_x,
            iliadis_jch1=ref,
            local_over_iliadis=ratio,
            abs_log10_error=abslog,
            local_label=local_label,
        ))
    end
    df = DataFrame(rows)
    sort!(df, [:abs_log10_error], rev=true)
    return df
end

function trajectory_table(traj)
    ipeak = argmax(traj.temperature_T9)
    t0 = first(traj.time_s)
    tf = last(traj.time_s)
    return DataFrame([
        (property="Tpeak_T9", local_value=traj.temperature_T9[ipeak], iliadis_jch1=JCH1_PROPERTIES["Tpeak_T9"], local_over_iliadis=traj.temperature_T9[ipeak] / JCH1_PROPERTIES["Tpeak_T9"], notes="Table 1 JCH1 Tpeak = 231 MK"),
        (property="rho_at_Tpeak_cgs", local_value=traj.density_cgs[ipeak], iliadis_jch1=missing, local_over_iliadis=missing, notes="Iliadis Table 1 does not list density"),
        (property="time_at_Tpeak_s", local_value=traj.time_s[ipeak] - t0, iliadis_jch1=missing, local_over_iliadis=missing, notes="seconds since first trajectory point"),
        (property="duration_s", local_value=tf - t0, iliadis_jch1=missing, local_over_iliadis=missing, notes="network trajectory duration"),
        (property="n_trajectory_points", local_value=Float64(nrow(traj)), iliadis_jch1=missing, local_over_iliadis=missing, notes=""),
    ])
end

function write_summary(path, trajectory, initial, final_raw, final_decay, proton, network_rates, decay_conv)
    generated_at = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
    function top_rows(df; n=8)
        valid = df[.!ismissing.(df.abs_log10_error), :]
        first(valid, min(n, nrow(valid)))
    end
    open(path, "w") do io
        println(io, "# Iliadis 2002 JCH1 Baseline Audit")
        println(io)
        println(io, "Generated: $generated_at")
        println(io)
        println(io, "Reference: `$ARTICLE_PDF`")
        println(io)
        println(io, "## Trajectory")
        println(io)
        for row in eachrow(trajectory)
            println(io, "- $(row.property): local $(row.local_value), Iliadis JCH1 $(row.iliadis_jch1)")
        end
        println(io)
        println(io, "## Largest Initial Composition Differences")
        println(io)
        println(io, "| isotope | local | Iliadis JCH1 | local/Iliadis |")
        println(io, "| --- | ---: | ---: | ---: |")
        for row in eachrow(top_rows(initial))
            println(io, "| $(row.isotope) | $(row.local_value) | $(row.iliadis_jch1) | $(row.local_over_iliadis) |")
        end
        println(io)
        println(io, "## Largest Final Baseline Differences After Simple Decay Collapse")
        println(io)
        println(io, "| isotope | local | Iliadis JCH1 | local/Iliadis |")
        println(io, "| --- | ---: | ---: | ---: |")
        for row in eachrow(top_rows(final_decay))
            println(io, "| $(row.isotope) | $(row.local_value) | $(row.iliadis_jch1) | $(row.local_over_iliadis) |")
        end
        println(io)
        println(io, "The decay-collapsed comparison moves short-lived parents listed in the script into their stable daughters before comparing to Table 4.")
        println(io, "Large remaining baseline differences mean sensitivity ratios can disagree even when the varied rate source is correct.")
        println(io)
        println(io, "## Proton Consumption")
        println(io)
        hfinal = proton[proton.metric .== "final decay-collapsed H-1", :][1, :]
        println(io, "- Final local H-1 is $(hfinal.local_value), versus Iliadis JCH1 $(hfinal.iliadis_jch1).")
        for row in eachrow(proton[occursin.("first H-1 below", proton.metric), :])
            println(io, "- $(row.metric): $(row.time_since_start_s) s after start.")
        end
        println(io)
        println(io, "## Decay Convention")
        println(io)
        improved = decay_conv[.!ismissing.(decay_conv.decay_improvement_log10) .& (decay_conv.decay_improvement_log10 .> 0.25), :]
        println(io, "- Isotopes materially improved by simple decay collapse: $(nrow(improved)).")
        println(io, "- Decay collapse does not address H-1 exhaustion, N-15 overproduction, or most Ne-Na/Mg-Al offsets.")
        println(io)
        println(io, "## Network Rule Check")
        println(io)
        crossings = network_rates[network_rates.active_crosses_a40 .== true, :]
        active_heavy = network_rates[(network_rates.active .== true) .& (network_rates.controlled_charged_particle .== true) .& (network_rates.target_a .> 40), :]
        println(io, "- Active controlled charged-particle rows that cross above A=40: $(nrow(crossings)).")
        println(io, "- Active controlled charged-particle rows with target A>40: $(nrow(active_heavy)).")
        println(io)
        println(io, "Raw final comparison is also written for checking whether decay convention is the dominant issue.")
    end
end

function audit_baseline(; nova=DEFAULT_NOVA)
    base = joinpath(PROJECT_ROOT, "novae", nova)
    ppn = joinpath(base, "ppn")
    runs = joinpath(base, "runs")
    results = joinpath(base, "analysis", "results")
    mkpath(results)

    initial_path = joinpath(ppn, "initial_abundance.dat")
    trajectory_path = joinpath(ppn, "trajectory.input")
    network_path = joinpath(ppn, "networksetup.txt")
    final_path = joinpath(runs, "baseline", "final_abundances.csv")
    xtime_path = joinpath(runs, "baseline", "x-time.dat")
    isfile(final_path) || error("missing baseline final abundances: $final_path")

    local_initial = parse_initial_abundance(initial_path)
    local_final_raw = parse_final_abundances(final_path)
    local_final_decay = decay_collapsed(local_final_raw)
    traj = parse_trajectory(trajectory_path)
    xtime = parse_xtime(xtime_path)

    trajectory = trajectory_table(traj)
    initial = comparison_table(local_initial, JCH1_INITIAL; local_label="initial_abundance.dat")
    final_raw = comparison_table(local_final_raw, JCH1_FINAL; local_label="baseline final raw")
    final_decay = comparison_table(local_final_decay, JCH1_FINAL; local_label="baseline final decay-collapsed")
    proton = proton_consumption_table(xtime, local_final_decay)
    evolution = evolution_table(xtime, local_final_decay)
    decay_conv = decay_convention_table(local_final_raw, local_final_decay)
    network_rates = network_rate_table(network_path)
    network_summary = network_rate_summary(network_rates, joinpath(results, "iliadis_recreation_rate_audit.csv"))

    trajectory_path_out = joinpath(results, "iliadis_baseline_trajectory_audit.csv")
    initial_path_out = joinpath(results, "iliadis_baseline_initial_abundance_audit.csv")
    final_raw_path_out = joinpath(results, "iliadis_baseline_final_abundance_raw_audit.csv")
    final_decay_path_out = joinpath(results, "iliadis_baseline_final_abundance_decay_audit.csv")
    proton_path_out = joinpath(results, "iliadis_baseline_proton_consumption_audit.csv")
    evolution_path_out = joinpath(results, "iliadis_baseline_evolution_audit.csv")
    decay_conv_path_out = joinpath(results, "iliadis_baseline_decay_convention_audit.csv")
    network_rates_path_out = joinpath(results, "iliadis_baseline_network_rate_audit.csv")
    network_summary_path_out = joinpath(results, "iliadis_baseline_network_rate_summary.csv")
    summary_path = joinpath(results, "iliadis_baseline_audit.md")

    CSV.write(trajectory_path_out, trajectory)
    CSV.write(initial_path_out, initial)
    CSV.write(final_raw_path_out, final_raw)
    CSV.write(final_decay_path_out, final_decay)
    CSV.write(proton_path_out, proton)
    CSV.write(evolution_path_out, evolution)
    CSV.write(decay_conv_path_out, decay_conv)
    CSV.write(network_rates_path_out, network_rates)
    CSV.write(network_summary_path_out, network_summary)
    write_summary(summary_path, trajectory, initial, final_raw, final_decay, proton, network_rates, decay_conv)

    println("Wrote $trajectory_path_out")
    println("Wrote $initial_path_out")
    println("Wrote $final_raw_path_out")
    println("Wrote $final_decay_path_out")
    println("Wrote $proton_path_out")
    println("Wrote $evolution_path_out")
    println("Wrote $decay_conv_path_out")
    println("Wrote $network_rates_path_out")
    println("Wrote $network_summary_path_out")
    println("Wrote $summary_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    audit_baseline(nova=parse_args(ARGS))
end
