#!/usr/bin/env julia

# Audit the PPN rate-source choices against the Iliadis et al. (2002)
# Table 8 sensitivity reactions and the local NuGrid/NETGEN data bundle.

using CSV
using DataFrames
using Dates
using Printf

include(joinpath(@__DIR__, "NovaRunTools.jl"))
const NRT = NovaRunTools

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_NOVA = "nova_test"

const ARTICLE_PDF = "/home/sgervais/Documents/PhD/articles/Iliadis-ApJS-2002.pdf"

const EXCEPTION_REFERENCES = Dict(
    "8B_pg_9C" => ("Wiescher et al. 1989", ["JINAC", "VITAL"]),
    "17F_pg_18Ne" => ("Bardayan et al. 2000", ["JINAC", "VITAL"]),
    "17O_pg_18F" => ("J. Blackmon et al. 2002, in preparation", ["JINAC", "NACRR", "VITAL"]),
    "17O_pa_14N" => ("J. Blackmon et al. 2002, in preparation", ["JINAC", "NACRR", "VITAL"]),
    "18F_pg_19Ne" => ("Coc et al. 2000", ["JINAC", "VITAL"]),
    "18F_pa_15O" => ("Coc et al. 2000", ["JINAC", "NACRR", "VITAL"]),
    "13C_pg_14N" => ("Caughlan & Fowler 1988", ["JINAC", "NACRR", "VITAL"]),
    "14N_pg_15O" => ("Caughlan & Fowler 1988", ["JINAC", "NACRR", "VITAL"]),
    "16O_pg_17F" => ("Caughlan & Fowler 1988", ["JINAC", "NACRR", "VITAL"]),
    "18O_pa_15N" => ("Caughlan & Fowler 1988", ["JINAC", "NACRR", "VITAL"]),
    "19F_pg_20Ne" => ("Caughlan & Fowler 1988", ["JINAC", "NACRR", "VITAL"]),
)

const FACTOR_COLUMNS = ["100", "10", "2", "0.5", "0.1", "0.01"]

function usage()
    println("""
    Usage:
      julia --project=NovaProject_0.2/NovaJL NovaProject_0.2/tools/audit_iliadis_recreation.jl [--nova nova_test]

    Outputs:
      novae/<nova>/analysis/results/iliadis_recreation_rate_audit.csv
      novae/<nova>/analysis/results/iliadis_recreation_source_summary.csv
      novae/<nova>/analysis/results/iliadis_recreation_rules.md
    """)
end

function parse_args(args)
    nova = DEFAULT_NOVA
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            nova = args[i]
        elseif arg == "-h" || arg == "--help"
            usage()
            exit(0)
        else
            error("unknown argument $arg")
        end
        i += 1
    end
    return nova
end

function species_label(species)
    species === nothing && return ""
    a, el = species
    a == 0 && return "gamma"
    el == "H" && a == 1 && return "p"
    el == "N" && a == 1 && return "n"
    el == "HE" && a == 4 && return "alpha"
    return "$(a)$(titlecase(lowercase(el)))"
end

function species_from_netgen_line(line)
    text = strip(replace(strip(line), r"^#\s*" => ""))
    isempty(text) && return nothing
    occursin(r"^[+=]", text) && return nothing
    startswith(text, "Qrad") && return :qrad
    startswith(text, "Qnu") && return nothing
    startswith(text, "Type/Ne") && return nothing
    startswith(text, "T8") && return nothing
    tokens = split(text)
    length(tokens) < 2 && return nothing
    try
        parse(Int, tokens[1])
    catch
        return nothing
    end
    name = String(tokens[2])
    name == "PROT" && return (1, "H")
    name == "NEUT" && return (1, "N")
    name == "OOOOO" && return (0, "G")
    length(tokens) >= 3 || return nothing
    mass = match(r"^\d+", tokens[3])
    mass === nothing && return nothing
    return (parse(Int, mass.match), name)
end

function reaction_name_from_species(species)
    length(species) >= 4 || return nothing
    target, projectile, product1, product2 = species[1], species[2], species[3], species[4]
    projectile == (1, "H") || return nothing
    if product1 == (0, "G") || product2 == (0, "G")
        product = product1 == (0, "G") ? product2 : product1
        return "$(species_label(target))_pg_$(species_label(product))"
    elseif product1 == (4, "HE") || product2 == (4, "HE")
        product = product1 == (4, "HE") ? product2 : product1
        return "$(species_label(target))_pa_$(species_label(product))"
    end
    return nothing
end

function parse_netgen_reaction_names(path)
    names = Set{String}()
    species = Tuple{Int,String}[]
    for line in eachline(path)
        startswith(strip(line), "#") || continue
        parsed = species_from_netgen_line(line)
        if parsed === :qrad
            name = reaction_name_from_species(species)
            name !== nothing && push!(names, name)
            empty!(species)
        elseif parsed isa Tuple{Int,String}
            push!(species, parsed)
        end
    end
    return names
end

function reaction_policy(reaction_name)
    target, _, rtype, _ = NRT.reaction_from_name(reaction_name)
    if haskey(EXCEPTION_REFERENCES, reaction_name)
        ref, local_sources = EXCEPTION_REFERENCES[reaction_name]
        return ref, local_sources
    elseif 1 <= target[1] <= 19
        return "Angulo et al. 1999 / NACRE recommended library", ["NACRR", "JINAC", "VITAL"]
    elseif 20 <= target[1] <= 40 && rtype in Set(["(p,g)", "(p,a)"])
        return "Iliadis et al. 2001 NETGEN library", ["ILI01", "JINAC", "NACRR", "VITAL"]
    else
        return "outside Iliadis 2002 Table 8 local policy", ["JINAC", "VITAL"]
    end
end

function active_rows(candidates)
    [row for row in candidates if row.active]
end

function format_rows(rows)
    isempty(rows) && return ""
    join(["$(row.index):$(row.source):$(row.active ? "T" : "F")" for row in sort(rows, by=row -> row.index)], "; ")
end

function factors_present(df)
    out = Dict{String,String}()
    for reaction in unique(df.reaction)
        sub = df[df.reaction .== reaction, :]
        factors = String[]
        for factor in FACTOR_COLUMNS
            col = Symbol(factor)
            col in propertynames(sub) || continue
            any(!ismissing(x) for x in sub[!, col]) && push!(factors, factor)
        end
        out[reaction] = join(factors, ";")
    end
    return out
end

function load_plan_by_reaction(path)
    isfile(path) || return Dict{String,Any}()
    config = NRT.parse_json_file(path)
    plan = Dict{String,Any}()
    for reaction in config["reactions"]
        plan[string(reaction["name"])] = reaction
    end
    return plan
end

function audit_reactions(; nova=DEFAULT_NOVA)
    base = joinpath(PROJECT_ROOT, "novae", nova)
    network_path = joinpath(base, "ppn", "networksetup.txt")
    table8_path = joinpath(base, "analysis", "iliadis_table8_jch1_full.csv")
    plan_path = joinpath(base, "config", "reaction_plan.json")
    netgen_iliadis_path = joinpath(base, "NPDATA", "netgen", "netgen_iliadis2001_log_100.txt")
    results_dir = joinpath(base, "analysis", "results")
    mkpath(results_dir)

    isfile(network_path) || error("missing networksetup: $network_path")
    isfile(table8_path) || error("missing Iliadis Table 8 CSV: $table8_path")
    isfile(netgen_iliadis_path) || error("missing Iliadis NETGEN table: $netgen_iliadis_path")

    table8 = CSV.read(table8_path, DataFrame; missingstring=["", "missing"])
    network = NRT.parse_networksetup(network_path)
    plan = load_plan_by_reaction(plan_path)
    netgen_iliadis_names = parse_netgen_reaction_names(netgen_iliadis_path)
    factor_map = factors_present(table8)

    rows = NamedTuple[]
    for reaction_name in sort(unique(String.(table8.reaction)))
        target, _, rtype, _ = NRT.reaction_from_name(reaction_name)
        candidates, product_was_remapped = NRT.matching_rows_for_reaction(Dict("name" => reaction_name), network)
        active = active_rows(candidates)
        selected = isempty(active) ? nothing : sort(active, by=row -> row.index)[1]
        policy_ref, desired_sources = reaction_policy(reaction_name)
        table_rows = nrow(table8[table8.reaction .== reaction_name, :])
        plan_entry = get(plan, reaction_name, nothing)
        plan_index = plan_entry === nothing || !haskey(plan_entry, "index") ? missing : Int(plan_entry["index"])
        plan_source = plan_entry === nothing || !haskey(plan_entry, "source") ? missing : string(plan_entry["source"])
        reverse_index = plan_entry === nothing || !haskey(plan_entry, "reverse_index") ? missing : Int(plan_entry["reverse_index"])
        reverse_row = ismissing(reverse_index) ? nothing : NRT.row_by_index(network, reverse_index)

        active_source = selected === nothing ? missing : selected.source
        status = "ok"
        notes = String[]

        if isempty(candidates)
            status = "missing_network_row"
            push!(notes, "no matching networksetup row")
        elseif selected === nothing
            status = "inactive"
            push!(notes, "matching rows exist but none are active")
        elseif !(selected.source in desired_sources)
            status = "source_mismatch"
            desired = join(desired_sources, ">")
            push!(notes, "active source $(selected.source) is not in desired local source order $desired")
        elseif target[1] >= 20 && target[1] <= 40 && rtype in Set(["(p,g)", "(p,a)"]) && selected.source != "ILI01"
            status = "fallback"
            push!(notes, "Iliadis 2001 is preferred for this mass range but active source is $(selected.source)")
        elseif target[1] < 20 && haskey(EXCEPTION_REFERENCES, reaction_name)
            status = "local_surrogate"
            push!(notes, "article cites $(policy_ref); local PPN source labels do not identify that library explicitly")
        end

        if product_was_remapped
            push!(notes, "product was remapped while matching network row")
        end
        if product_was_remapped && occursin(r"Al[gm]$", reaction_name)
            status = "isomer_split_unresolved"
            push!(notes, "Table 8 separates 26Al ground/metastable products, but the local Iliadis NETGEN row is generic 26Al")
        end
        if selected !== nothing && selected.source == "ILI01" && !(reaction_name in netgen_iliadis_names)
            status = "iliadis_label_without_table_match"
            push!(notes, "network source is ILI01, but reaction name was not found in parsed Iliadis NETGEN headers")
        end
        if !ismissing(reverse_index)
            if reverse_row === nothing
                status = "bad_reverse_index"
                push!(notes, "configured reverse_index $reverse_index was not found")
            elseif !reverse_row.active
                push!(notes, "configured reverse_index $reverse_index is inactive")
            elseif reverse_row.source != active_source
                push!(notes, "reverse source $(reverse_row.source) differs from active forward source $(active_source)")
            end
        end

        push!(rows, (
            reaction = reaction_name,
            target_a = target[1],
            rtype = rtype,
            table8_isotope_rows = table_rows,
            table8_factors = get(factor_map, reaction_name, ""),
            article_rate_reference = policy_ref,
            desired_local_source_order = join(desired_sources, ">"),
            in_iliadis2001_netgen = reaction_name in netgen_iliadis_names,
            plan_index = plan_index,
            plan_source = plan_source,
            plan_reverse_index = reverse_index,
            active_index = selected === nothing ? missing : selected.index,
            active_source = active_source,
            active_reverse_index = reverse_row === nothing ? missing : reverse_row.index,
            active_reverse_source = reverse_row === nothing ? missing : reverse_row.source,
            candidates = format_rows(candidates),
            status = status,
            notes = join(notes, " | "),
        ))
    end

    audit = DataFrame(rows)
    sort!(audit, [:target_a, :reaction])

    summary = combine(groupby(audit, [:status, :active_source]), nrow => :count)
    sort!(summary, [:status, :active_source])

    audit_path = joinpath(results_dir, "iliadis_recreation_rate_audit.csv")
    summary_path = joinpath(results_dir, "iliadis_recreation_source_summary.csv")
    abundance_path = joinpath(results_dir, "iliadis_recreation_abundance_deviation_audit.csv")
    rules_path = joinpath(results_dir, "iliadis_recreation_rules.md")

    CSV.write(audit_path, audit)
    CSV.write(summary_path, summary)
    write_abundance_deviation_audit(abundance_path, audit, joinpath(base, "analysis", "ppnData+IliTable8Data.csv"))
    write_rules_markdown(rules_path, audit, network_path, table8_path, netgen_iliadis_path)

    println("Wrote $audit_path")
    println("Wrote $summary_path")
    println("Wrote $abundance_path")
    println("Wrote $rules_path")
    println()
    println(summary)
    return audit
end

function parse_optional_float(value)
    if value === missing
        return missing
    elseif value isa Number
        return Float64(value)
    end
    text = strip(string(value))
    isempty(text) && return missing
    return try
        parse(Float64, text)
    catch
        missing
    end
end

function write_abundance_deviation_audit(path, rate_audit, comparison_path)
    if !isfile(comparison_path)
        @warn "comparison CSV not found; skipping abundance deviation audit" comparison_path
        return
    end

    comparison = CSV.read(comparison_path, DataFrame; missingstring=["", "missing"])
    rate_lookup = Dict(string(row.reaction) => row for row in eachrow(rate_audit))
    rows = NamedTuple[]

    for row in eachrow(comparison)
        reaction = string(row.reaction)
        audit_row = get(rate_lookup, reaction, nothing)
        for factor in FACTOR_COLUMNS
            ppn_col = Symbol(factor * "_ppn")
            ili_col = Symbol(factor * "_iliadis")
            ppn_col in propertynames(comparison) || continue
            ili_col in propertynames(comparison) || continue
            ppn = parse_optional_float(row[ppn_col])
            ili = parse_optional_float(row[ili_col])
            (ismissing(ppn) || ismissing(ili) || ppn <= 0 || ili <= 0) && continue
            ratio = ppn / ili
            push!(rows, (
                reaction = reaction,
                isotope = string(row.isotope),
                factor = factor,
                ppn = ppn,
                iliadis = ili,
                ppn_over_iliadis = ratio,
                abs_log10_error = abs(log10(ratio)),
                signed_log10_error = log10(ratio),
                source_status = audit_row === nothing ? missing : audit_row.status,
                active_source = audit_row === nothing ? missing : audit_row.active_source,
                active_index = audit_row === nothing ? missing : audit_row.active_index,
                source_notes = audit_row === nothing ? "" : audit_row.notes,
            ))
        end
    end

    deviations = DataFrame(rows)
    if nrow(deviations) > 0
        sort!(deviations, :abs_log10_error, rev=true)
    end
    CSV.write(path, deviations)
end

function write_rules_markdown(path, audit, network_path, table8_path, netgen_iliadis_path)
    status_counts = combine(groupby(audit, :status), nrow => :count)
    sort!(status_counts, :status)
    fallback = audit[in.(audit.status, Ref(["fallback", "isomer_split_unresolved", "local_surrogate", "source_mismatch", "missing_network_row", "inactive", "iliadis_label_without_table_match"])), :]
    generated_at = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")

    open(path, "w") do io
        println(io, "# Iliadis 2002 PPN Recreation Rate Rules")
        println(io)
        println(io, "Generated: $generated_at")
        println(io)
        println(io, "Inputs:")
        println(io, "- Article: `$ARTICLE_PDF`")
        println(io, "- Table 8 CSV: `$table8_path`")
        println(io, "- PPN network: `$network_path`")
        println(io, "- Iliadis NETGEN table: `$netgen_iliadis_path`")
        println(io)
        println(io, "## Policy")
        println(io)
        println(io, "1. Use the local NACRE recommended source (`NACRR`) for A=1..19 charged-particle rates when no paper-specific exception is available as a named PPN source.")
        println(io, "2. Use `ILI01` for A=20..40 `(p,g)` and `(p,a)` reactions whenever the local `netgen_iliadis2001_log_100.txt` table provides the reaction.")
        println(io, "3. Use `JINAC`, then `NACRR`, then `VITAL` only as explicit fallbacks when the preferred local source is unavailable.")
        println(io, "4. Keep the sensitivity recreation scoped to reactions with target mass A<=40, matching the published nova sensitivity range.")
        println(io, "5. Treat the paper-specific A<20 sources as local surrogates unless their exact library has been imported: `17O+p` from Blackmon et al., `18F+p` from Coc et al., `17F(p,g)` from Bardayan et al., selected rates from CF88, and `8B(p,g)` from Wiescher et al.")
        println(io, "6. For reactions with configured inverse branches, audit the reverse index and source. Do not interpret a forward-only rate factor as an exact paper reproduction if the inverse branch is active and uses a different source.")
        println(io)
        println(io, "These rules are encoded for PPN setup in `novae/nova_test/config/network_edits.json`; this report audits the resulting current `networksetup.txt`.")
        println(io)
        println(io, "## Status Counts")
        println(io)
        println(io, "| status | count |")
        println(io, "| --- | ---: |")
        for row in eachrow(status_counts)
            println(io, "| $(row.status) | $(row.count) |")
        end
        println(io)
        println(io, "## Rows Requiring Attention")
        println(io)
        if nrow(fallback) == 0
            println(io, "No rows require attention under the current policy.")
        else
            println(io, "| reaction | status | active | desired | note |")
            println(io, "| --- | --- | --- | --- | --- |")
            for row in eachrow(fallback)
                active = ismissing(row.active_index) ? "" : "$(row.active_index):$(row.active_source)"
                note = replace(row.notes, "|" => "\\|")
                println(io, "| $(row.reaction) | $(row.status) | $active | $(row.desired_local_source_order) | $note |")
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    nova = parse_args(ARGS)
    audit_reactions(nova=nova)
end
