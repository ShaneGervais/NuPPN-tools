module NovaRunTools

using Printf

export PROJECT_ROOT,
    setup_network,
    create_factored_runs,
    parse_networksetup,
    apply_network_edits,
    normalize_networksetup,
    validate_plan

Base.@kwdef mutable struct JsonParser
    text::String
    pos::Int = 1
end

function skip_ws!(p::JsonParser)
    while p.pos <= lastindex(p.text) && p.text[p.pos] in (' ', '\n', '\r', '\t')
        p.pos = nextind(p.text, p.pos)
    end
end

function parse_json_string!(p::JsonParser)
    p.text[p.pos] == '"' || error("Expected JSON string at byte $(p.pos)")
    p.pos = nextind(p.text, p.pos)
    out = IOBuffer()
    while p.pos <= lastindex(p.text)
        ch = p.text[p.pos]
        if ch == '"'
            p.pos = nextind(p.text, p.pos)
            return String(take!(out))
        elseif ch == '\\'
            p.pos = nextind(p.text, p.pos)
            esc = p.text[p.pos]
            if esc == '"' || esc == '\\' || esc == '/'
                print(out, esc)
            elseif esc == 'b'
                print(out, '\b')
            elseif esc == 'f'
                print(out, '\f')
            elseif esc == 'n'
                print(out, '\n')
            elseif esc == 'r'
                print(out, '\r')
            elseif esc == 't'
                print(out, '\t')
            else
                error("Unsupported JSON escape \\$esc")
            end
        else
            print(out, ch)
        end
        p.pos = nextind(p.text, p.pos)
    end
    error("Unterminated JSON string")
end

function parse_json_number!(p::JsonParser)
    start = p.pos
    while p.pos <= lastindex(p.text) && p.text[p.pos] in Set(['-', '+', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.', 'e', 'E'])
        p.pos = nextind(p.text, p.pos)
    end
    raw = p.text[start:prevind(p.text, p.pos)]
    occursin(r"[.eE]", raw) ? parse(Float64, raw) : parse(Int, raw)
end

function parse_json_array!(p::JsonParser)
    p.pos = nextind(p.text, p.pos)
    values = Any[]
    while true
        skip_ws!(p)
        if p.text[p.pos] == ']'
            p.pos = nextind(p.text, p.pos)
            return values
        end
        push!(values, parse_json_value!(p))
        skip_ws!(p)
        if p.text[p.pos] == ','
            p.pos = nextind(p.text, p.pos)
        elseif p.text[p.pos] == ']'
            p.pos = nextind(p.text, p.pos)
            return values
        else
            error("Expected ',' or ']' at byte $(p.pos)")
        end
    end
end

function parse_json_object!(p::JsonParser)
    p.pos = nextind(p.text, p.pos)
    obj = Dict{String,Any}()
    while true
        skip_ws!(p)
        if p.text[p.pos] == '}'
            p.pos = nextind(p.text, p.pos)
            return obj
        end
        key = parse_json_string!(p)
        skip_ws!(p)
        p.text[p.pos] == ':' || error("Expected ':' at byte $(p.pos)")
        p.pos = nextind(p.text, p.pos)
        obj[key] = parse_json_value!(p)
        skip_ws!(p)
        if p.text[p.pos] == ','
            p.pos = nextind(p.text, p.pos)
        elseif p.text[p.pos] == '}'
            p.pos = nextind(p.text, p.pos)
            return obj
        else
            error("Expected ',' or '}' at byte $(p.pos)")
        end
    end
end

function parse_json_value!(p::JsonParser)
    skip_ws!(p)
    ch = p.text[p.pos]
    ch == '"' && return parse_json_string!(p)
    ch == '{' && return parse_json_object!(p)
    ch == '[' && return parse_json_array!(p)
    if startswith(p.text[p.pos:end], "true")
        p.pos += 4
        return true
    elseif startswith(p.text[p.pos:end], "false")
        p.pos += 5
        return false
    elseif startswith(p.text[p.pos:end], "null")
        p.pos += 4
        return nothing
    elseif ch in Set(['-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
        return parse_json_number!(p)
    end
    error("Unexpected JSON value at byte $(p.pos)")
end

function parse_json_file(path)
    parser = JsonParser(text=read(path, String))
    value = parse_json_value!(parser)
    skip_ws!(parser)
    parser.pos <= lastindex(parser.text) && error("Unexpected trailing JSON content")
    return value
end

const PROJECT_ROOT = dirname(@__DIR__)
const VALID_SOURCES = Set([
    "BASEL", "JINAR", "JINAC", "JINAV", "ILI01", "NACRL", "NACRR", "NACRU",
    "VITAL", "RVRSE", "KADON", "NETB1", "ODA94", "LMP00", "FFW85", "JBJ16",
    "NKK04", "MPG06", "NTRNO", "OOOOO",
])
const SOURCE_ALLOWED_RTYPES = Dict(
    "NACRL" => Set(["(n,a)", "(p,g)", "(p,n)", "(p,a)", "(a,g)", "(g,a)", "(a,n)", "(a,p)"]),
    "NACRR" => Set(["(n,a)", "(p,g)", "(p,n)", "(p,a)", "(a,g)", "(g,a)", "(a,n)", "(a,p)"]),
    "NACRU" => Set(["(n,a)", "(p,g)", "(p,n)", "(p,a)", "(a,g)", "(g,a)", "(a,n)", "(a,p)"]),
    "ILI01" => Set(["(p,g)", "(p,a)"]),
    "KADON" => Set(["(n,g)"]),
    "NETB1" => Set(["(-,g)", "(+,g)"]),
    "ODA94" => Set(["(-,g)", "(+,g)"]),
    "LMP00" => Set(["(-,g)", "(+,g)"]),
    "FFW85" => Set(["(-,g)", "(+,g)"]),
    "JBJ16" => Set(["(-,g)", "(+,g)"]),
    "NKK04" => Set(["(-,g)", "(+,g)"]),
)
const ACTIVE_NACRE_SOURCE = "NACRR"
const NOVA_MAX_A = 40
const CONTROLLED_RTYPES = Set(["(p,g)", "(p,a)", "(a,g)", "(a,n)", "(a,p)"])
const WEAK_RTYPES = Set(["(-,g)", "(+,g)"])
const ALPHA_REVERSE_RTYPES = Set(["(p,a)", "(a,p)", "(n,a)", "(a,n)"])
const REVERSE_RTYPE = Dict("(p,a)" => "(a,p)", "(a,p)" => "(p,a)", "(n,a)" => "(a,n)", "(a,n)" => "(n,a)")
const ELEMENT_SYMBOLS = Set([
    "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE",
    "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
    "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR",
    "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
    "SB", "TE", "I", "XE", "CS", "BA", "LA", "CE", "PR", "ND",
    "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
    "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",
    "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
    "PA", "U",
])

const NETWORK_RE = r"^\s*(\d+)\s+([TF])\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+->\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+\S+\s+(\S+)\s+(\S+)\s+(\d+)"

Base.@kwdef struct Row
    index::Int
    active::Bool
    reactant::Union{Nothing,Tuple{Int,String}}
    projectile::Union{Nothing,Tuple{Int,String}}
    product_1::Union{Nothing,Tuple{Int,String}}
    product_2::Union{Nothing,Tuple{Int,String}}
    source::String
    rtype::String
    line_no::Int
    line::String
end

nova_dir(name) = joinpath(PROJECT_ROOT, "nova_cases", name)

function getkey(d::Dict, key::AbstractString, default=nothing)
    haskey(d, key) ? d[key] : default
end

function species_id(text)
    text = strip(text)
    text == "PROT" && return (1, "H")
    text == "NEUT" && return (1, "N")
    text == "OOOOO" && return (0, "G")
    m = match(r"([A-Z]+)\s*(\d+)$", text)
    m === nothing && return nothing
    return (parse(Int, m.captures[2]), m.captures[1])
end

function species_from_name(text)
    m = match(r"(\d+)([A-Za-z]+)$", text)
    m === nothing && error("Cannot parse isotope name $(repr(text))")
    symbol = uppercase(m.captures[2])
    if !(symbol in ELEMENT_SYMBOLS) && length(symbol) > 1 && last(symbol) in ('G', 'M')
        base = symbol[1:end-1]
        base in ELEMENT_SYMBOLS && (symbol = base)
    end
    return (parse(Int, m.captures[1]), symbol)
end

function reaction_from_name(name)
    parts = split(name, "_")
    length(parts) == 3 || error("Reaction name must look like '20Ne_pg_21Na': $name")
    target, channel, product = parts
    channels = Dict(
        "pg" => ("(p,g)", (1, "H")),
        "pa" => ("(p,a)", (1, "H")),
        "ag" => ("(a,g)", (4, "HE")),
        "an" => ("(a,n)", (4, "HE")),
        "ap" => ("(a,p)", (4, "HE")),
    )
    haskey(channels, channel) || error("Unsupported reaction channel $(repr(channel)) in $name")
    rtype, projectile = channels[channel]
    return species_from_name(target), projectile, rtype, species_from_name(product)
end

function parse_networksetup(path)
    rows = Row[]
    for (line_no0, line) in enumerate(eachline(path))
        m = match(NETWORK_RE, line)
        m === nothing && continue
        c = m.captures
        push!(rows, Row(
            index=parse(Int, c[1]),
            active=c[2] == "T",
            reactant=species_id(c[4]),
            projectile=species_id(c[6]),
            product_1=species_id(c[8]),
            product_2=species_id(c[10]),
            source=c[11],
            rtype=c[12],
            line_no=line_no0,
            line=strip(line),
        ))
    end
    return rows
end

products(row::Row) = sort([row.product_1, row.product_2], by=x -> string(x))
reaction_key(row::Row) = (row.reactant, row.projectile, row.rtype, products(row))
function row_by_index(network, index::Integer)
    found = findfirst(row -> row.index == index, network)
    found === nothing ? nothing : network[found]
end
row_max_a(row::Row) = maximum([sp[1] for sp in (row.reactant, row.projectile, row.product_1, row.product_2) if sp !== nothing]; init=0)
row_has_species_above(row::Row, max_a) = any(sp !== nothing && sp[1] > max_a for sp in (row.reactant, row.projectile, row.product_1, row.product_2))

function preferred_sources(target_a)
    1 <= target_a < 20 && return ["NACRR", "NACRL", "NACRU", "JINAC", "VITAL"]
    20 <= target_a <= 40 && return ["ILI01", "JINAC", "VITAL"]
    return String[]
end

function source_priority(row::Row)
    target_a = row.reactant === nothing ? row_max_a(row) : row.reactant[1]
    if row.rtype in WEAK_RTYPES
        priorities = ["ODA94", "NETB1", "FFW85", "LMP00", "JINAC", "BASEL", "JINAR", "JINAV"]
    elseif row.rtype == "(n,g)"
        priorities = target_a < 20 ?
            ["NACRR", "NACRL", "NACRU", "KADON", "JINAC", "BASEL", "JINAR", "JINAV"] :
            ["ILI01", "NACRR", "NACRL", "NACRU", "KADON", "JINAC", "BASEL", "JINAR", "JINAV"]
    elseif 1 <= target_a < 20
        priorities = ["NACRR", "NACRL", "NACRU", "JINAC", "BASEL", "JINAR", "JINAV", "VITAL", "RVRSE"]
    elseif 20 <= target_a <= NOVA_MAX_A
        priorities = ["ILI01", "JINAC", "BASEL", "JINAR", "JINAV", "NACRR", "NACRL", "NACRU", "KADON", "VITAL", "RVRSE"]
    else
        priorities = ["JINAC", "BASEL", "JINAR", "JINAV", "KADON", "RVRSE", "VITAL"]
    end
    found = findfirst(==(row.source), priorities)
    return found === nothing ? length(priorities) + 1 : found
end

function matching_rows_for_reaction(reaction, network)
    target, projectile, rtype, product = reaction_from_name(getkey(reaction, "network_name", reaction["name"]))
    candidates = [row for row in network if row.reactant == target && row.projectile == projectile && row.rtype == rtype && (row.product_1 == product || row.product_2 == product)]
    product_was_remapped = false
    if isempty(candidates)
        candidates = [row for row in network if row.reactant == target && row.projectile == projectile && row.rtype == rtype]
        product_was_remapped = !isempty(candidates)
    end
    return candidates, product_was_remapped
end

function select_configured_row(reaction, candidates)
    forced_index = getkey(reaction, "index")
    forced_index === nothing && return nothing
    forced = [row for row in candidates if row.index == Int(forced_index)]
    isempty(forced) ? nothing : first(forced)
end

function set_networksetup_active(line, active)
    flag = active ? "T" : "F"
    replace(line, r"^(\s*\d+\s+)[TF](\s+)" => SubstitutionString("\\1$(flag)\\2"); count=1)
end

function set_networksetup_source(line, source)
    source in VALID_SOURCES || error("Unsupported network source label $(repr(source))")
    replace(line, r"(\s+\S+\s+)(\S+)(\s+\([^)]+\)\s+\d+)" => SubstitutionString("\\1$(source)\\3"); count=1)
end

as_vector(value) = value isa AbstractArray ? collect(value) : [value]

function row_name(row::Row)
    channel = Dict("(p,g)" => "pg", "(p,a)" => "pa", "(a,g)" => "ag", "(a,n)" => "an", "(a,p)" => "ap")
    haskey(channel, row.rtype) || return nothing
    projectile = row.projectile
    if row.rtype in Set(["(p,g)", "(p,a)"]) && projectile != (1, "H")
        return nothing
    elseif row.rtype in Set(["(a,g)", "(a,n)", "(a,p)"]) && projectile != (4, "HE")
        return nothing
    end
    product = row.product_1 == (0, "G") ? row.product_2 : row.product_1
    product === nothing && return nothing
    product[1] == 0 && return nothing
    return "$(row.reactant[1])$(titlecase(lowercase(row.reactant[2])))_$(channel[row.rtype])_$(product[1])$(titlecase(lowercase(product[2])))"
end

function match_value(actual, matcher)
    matcher === nothing && return true
    values = as_vector(matcher)
    return actual in values
end

function row_matches(row::Row, match_config)
    isempty(match_config) && return true

    if haskey(match_config, "index") && !match_value(row.index, match_config["index"])
        return false
    end
    if haskey(match_config, "indices") && !(row.index in match_config["indices"])
        return false
    end
    if haskey(match_config, "name") && !match_value(row_name(row), match_config["name"])
        return false
    end
    if haskey(match_config, "names") && !(row_name(row) in match_config["names"])
        return false
    end
    if haskey(match_config, "source") && !match_value(row.source, match_config["source"])
        return false
    end
    if haskey(match_config, "rtype") && !match_value(row.rtype, match_config["rtype"])
        return false
    end
    if haskey(match_config, "active") && row.active != Bool(match_config["active"])
        return false
    end
    if haskey(match_config, "target_a_min") && (row.reactant === nothing || row.reactant[1] < Int(match_config["target_a_min"]))
        return false
    end
    if haskey(match_config, "target_a_max") && (row.reactant === nothing || row.reactant[1] > Int(match_config["target_a_max"]))
        return false
    end
    if haskey(match_config, "max_a_min") && row_max_a(row) < Int(match_config["max_a_min"])
        return false
    end
    if haskey(match_config, "max_a_max") && row_max_a(row) > Int(match_config["max_a_max"])
        return false
    end
    return true
end

function write_network_lines(path, lines, changed)
    changed > 0 && write(path, join(lines))
    return changed
end

function set_row_active!(lines, row::Row, active)
    lines[row.line_no] = set_networksetup_active(lines[row.line_no], active)
end

function set_row_source!(lines, row::Row, source)
    lines[row.line_no] = set_networksetup_source(lines[row.line_no], source)
end

function apply_simple_set!(lines, rows, set_config)
    changed = 0
    for row in rows
        if haskey(set_config, "active") && row.active != Bool(set_config["active"])
            set_row_active!(lines, row, Bool(set_config["active"]))
            changed += 1
        end
        if haskey(set_config, "source") && row.source != set_config["source"]
            set_config["source"] in VALID_SOURCES || error("Unsupported source $(set_config["source"])")
            set_row_source!(lines, row, set_config["source"])
            changed += 1
        end
    end
    return changed
end

function apply_choose_source!(lines, rows, sources; active=true)
    changed = 0
    source_list = String[string(s) for s in as_vector(sources)]
    for source in source_list
        source in VALID_SOURCES || error("Unsupported source $source")
    end

    groups = Dict{Any,Vector{Row}}()
    for row in rows
        push!(get!(groups, reaction_key(row), Row[]), row)
    end

    for (_, group_rows) in groups
        selected = nothing
        for source in source_list
            candidates = [row for row in group_rows if row.source == source]
            if !isempty(candidates)
                selected = sort(candidates, by=row -> (row.active ? 0 : 1, row.index))[1]
                break
            end
        end
        selected === nothing && continue

        for row in group_rows
            should_be_active = active && row.index == selected.index
            if row.active != should_be_active
                set_row_active!(lines, row, should_be_active)
                changed += 1
            end
        end
    end
    return changed
end

function reaction_edit_match_config(edit)
    match_config = Dict{String,Any}()
    if haskey(edit, "name")
        match_config["name"] = edit["name"]
    elseif haskey(edit, "names")
        match_config["names"] = edit["names"]
    else
        error("reaction_edits entries require name or names")
    end

    for key in ("index", "indices", "source", "rtype", "active", "target_a_min", "target_a_max", "max_a_min", "max_a_max")
        if haskey(edit, key)
            match_config[key] = edit[key]
        end
    end
    return match_config
end

function reaction_edit_set_config(edit)
    set_config = Dict{String,Any}()
    if haskey(edit, "set")
        merge!(set_config, edit["set"])
    end
    if haskey(edit, "new_source")
        set_config["source"] = edit["new_source"]
    elseif haskey(edit, "set_source")
        set_config["source"] = edit["set_source"]
    end
    if haskey(edit, "set_active")
        set_config["active"] = edit["set_active"]
    elseif haskey(edit, "new_active")
        set_config["active"] = edit["new_active"]
    end
    if isempty(set_config)
        error("reaction_edits entries require set, set_source/new_source, or set_active/new_active")
    end
    return set_config
end

function apply_reaction_edits(path, reaction_edits)
    lines = readlines(path, keep=true)
    total_changed = 0

    for (i, edit) in enumerate(reaction_edits)
        network = parse_networksetup(path)
        match_config = reaction_edit_match_config(edit)
        matched = [row for row in network if row_matches(row, match_config)]
        if isempty(matched)
            description = get(edit, "description", "reaction edit $i")
            error("network_edits: $description matched no rows")
        end

        set_config = reaction_edit_set_config(edit)
        changed = apply_simple_set!(lines, matched, set_config)
        if changed > 0
            write(path, join(lines))
            total_changed += changed
        end

        description = get(edit, "description", "reaction edit $i")
        indices = join([string(row.index) for row in matched], ", ")
        println("network_edits: $description matched $(length(matched)) rows [$indices], changed $changed fields")
    end

    return total_changed
end

function apply_network_edits(path, edits_config)
    total_changed = 0

    if haskey(edits_config, "reaction_edits")
        total_changed += apply_reaction_edits(path, edits_config["reaction_edits"])
    end

    lines = readlines(path, keep=true)
    rules = get(edits_config, "rules", Any[])

    for (i, rule) in enumerate(rules)
        network = parse_networksetup(path)
        match_config = get(rule, "match", Dict{String,Any}())
        matched = [row for row in network if row_matches(row, match_config)]

        changed = 0
        if haskey(rule, "choose_source")
            set_config = get(rule, "set", Dict{String,Any}())
            active = get(set_config, "active", true)
            changed += apply_choose_source!(lines, matched, rule["choose_source"]; active=Bool(active))
        end
        if haskey(rule, "set")
            set_config = copy(rule["set"])
            if haskey(rule, "choose_source") && haskey(set_config, "active")
                delete!(set_config, "active")
            end
            changed += apply_simple_set!(lines, matched, set_config)
        end

        if changed > 0
            write(path, join(lines))
            total_changed += changed
        end
        description = get(rule, "description", "rule $i")
        println("network_edits: $description matched $(length(matched)) rows, changed $changed fields")
    end

    return total_changed
end

function validate_source_for_row(reaction, source, row::Row)
    allowed = get(SOURCE_ALLOWED_RTYPES, source, nothing)
    if allowed !== nothing && !(row.rtype in allowed)
        error("$(reaction["name"]): source $source cannot provide reaction type $(row.rtype)")
    end
    if source in Set(["NACRL", "NACRR", "NACRU"]) && source != ACTIVE_NACRE_SOURCE
        error("$(reaction["name"]): requested $source, but this PPN build currently loads $ACTIVE_NACRE_SOURCE for all NACRE labels")
    end
    if source == "ILI01" && row.reactant[1] < 20
        error("$(reaction["name"]): ILI01 is not available for target A=$(row.reactant[1]) in the local Iliadis NETGEN table")
    end
end

function apply_configured_sources(path, config, network)
    lines = readlines(path, keep=true)
    changed = 0
    for reaction in config["reactions"]
        source = getkey(reaction, "source")
        source === nothing && continue
        source in VALID_SOURCES || error("$(reaction["name"]): unsupported source $(repr(source))")
        candidates, _ = matching_rows_for_reaction(reaction, network)
        isempty(candidates) && error("$(reaction["name"]): could not find reaction in $path")
        selected = select_configured_row(reaction, candidates)
        selected === nothing && error("$(reaction["name"]): configured index $(getkey(reaction, "index")) was not found among matching rows in $path")
        validate_source_for_row(reaction, source, selected)
        if selected.source != source
            lines[selected.line_no] = set_networksetup_source(lines[selected.line_no], source)
            changed += 1
        end
        for row in candidates
            should_be_active = row.line_no == selected.line_no
            if row.active != should_be_active
                lines[row.line_no] = set_networksetup_active(lines[row.line_no], should_be_active)
                changed += 1
            end
        end
    end
    changed > 0 && write(path, join(lines))
    return changed
end

function configured_active_indices_by_key(network, config)
    active_by_key = Dict{Any,Int}()
    for reaction in config["reactions"]
        candidates, _ = matching_rows_for_reaction(reaction, network)
        selected = isempty(candidates) ? nothing : select_configured_row(reaction, candidates)
        selected !== nothing && (active_by_key[reaction_key(selected)] = selected.index)
        reverse_index = getkey(reaction, "reverse_index")
        reverse_index === nothing && continue
        reverse = row_by_index(network, Int(reverse_index))
        reverse === nothing && error("$(reaction["name"]): reverse_index $reverse_index was not found")
        active_by_key[reaction_key(reverse)] = reverse.index
    end
    return active_by_key
end

function enforce_activation_policy(path, config)
    lines = readlines(path, keep=true)
    network = parse_networksetup(path)
    active_by_key = configured_active_indices_by_key(network, config)
    groups = Dict{Any,Vector{Row}}()
    for row in network
        push!(get!(groups, reaction_key(row), Row[]), row)
    end
    changed = 0
    for (key, rows) in groups
        selected_index = get(active_by_key, key, nothing)
        if selected_index !== nothing
            selected = row_by_index(rows, selected_index)
            selected === nothing && error("Configured index $selected_index was not found in its duplicate group")
        elseif any(row_has_species_above(row, NOVA_MAX_A) for row in rows)
            selected = nothing
        else
            selected = sort(rows, by=row -> (source_priority(row), row.active ? 0 : 1, row.index))[1]
        end
        for row in rows
            should_be_active = selected !== nothing && row.index == selected.index
            if row.active != should_be_active
                lines[row.line_no] = set_networksetup_active(lines[row.line_no], should_be_active)
                changed += 1
            end
        end
    end
    changed > 0 && write(path, join(lines))
    return changed
end

function normalize_networksetup(path, config)
    changed = 0
    while true
        network = parse_networksetup(path)
        source_changes = apply_configured_sources(path, config, network)
        changed += source_changes
        source_changes == 0 && break
    end
    changed += enforce_activation_policy(path, config)
    return changed
end

function validate_reverse_index(reaction, network, forward_index)
    !haskey(reaction, "reverse_index") && return nothing
    reverse_index = Int(reaction["reverse_index"])
    forward = row_by_index(network, forward_index)
    reverse = row_by_index(network, reverse_index)
    forward === nothing && error("$(reaction["name"]): forward index $forward_index was not found")
    reverse === nothing && error("$(reaction["name"]): reverse_index $reverse_index was not found")
    forward.rtype in ALPHA_REVERSE_RTYPES || error("$(reaction["name"]): reverse_index is only supported for alpha-transfer reactions, not $(forward.rtype)")
    reverse.rtype == REVERSE_RTYPE[forward.rtype] || error("$(reaction["name"]): reverse_index $reverse_index has type $(reverse.rtype), expected $(REVERSE_RTYPE[forward.rtype])")
    Set([reverse.reactant, reverse.projectile]) == Set([forward.product_1, forward.product_2]) || error("$(reaction["name"]): reverse_index $reverse_index is not the reverse of forward index $forward_index")
    Set([reverse.product_1, reverse.product_2]) == Set([forward.reactant, forward.projectile]) || error("$(reaction["name"]): reverse_index $reverse_index is not the reverse of forward index $forward_index")
    return reverse_index
end

function validate_plan(network, config)
    selected = []
    for reaction in config["reactions"]
        candidates, product_was_remapped = matching_rows_for_reaction(reaction, network)
        product_was_remapped = product_was_remapped || get(reaction, "product_was_remapped", false)
        isempty(candidates) && error("$(reaction["name"]): no matching row in networksetup.txt")
        row = select_configured_row(reaction, candidates)
        row === nothing && error("$(reaction["name"]): configured index $(getkey(reaction, "index")) was not found among matching rows")
        row.active || error("$(reaction["name"]): configured index $(row.index) is inactive")
        reverse_index = validate_reverse_index(reaction, network, row.index)
        push!(selected, (reaction["name"], row, reverse_index, product_was_remapped))
    end
    groups = Dict{Any,Vector{Row}}()
    for row in network
        push!(get!(groups, reaction_key(row), Row[]), row)
    end
    duplicate_active = [active for active in ([row for row in rows if row.active] for rows in values(groups)) if length(active) > 1]
    return selected, duplicate_active
end

function setup_network(nova; check=false, validate_factor_plan=false)
    base = nova_dir(nova)
    plan_path = joinpath(base, "config", "reaction_plan.json")
    edits_path = joinpath(base, "config", "network_edits.json")
    network_path = joinpath(base, "ppn", "networksetup.txt")
    isfile(network_path) || error("Missing $network_path")
    if isfile(edits_path)
        edits_config = parse_json_file(edits_path)
        !check && apply_network_edits(network_path, edits_config)
    elseif !check
        println("No network edit plan found at $edits_path; leaving networksetup.txt unchanged.")
    end
    network = parse_networksetup(network_path)
    println("Network: $network_path")
    isfile(edits_path) && println("Edits:   $edits_path")
    println("Parsed $(length(network)) reaction rows.")

    if validate_factor_plan
        isfile(plan_path) || error("Missing $plan_path")
        plan_config = parse_json_file(plan_path)
        selected, duplicate_active = validate_plan(network, plan_config)
        println("Plan:    $plan_path")
        println("Checked $(length(selected)) configured factor reactions.")
        for (name, row, reverse_index, product_was_remapped) in selected
            reverse = ""
            if reverse_index !== nothing
                reverse_row = row_by_index(network, reverse_index)
                reverse = ", reverse $reverse_index ($(reverse_row.source))"
            end
            note = product_was_remapped ? " [product remapped]" : ""
            println("  $name: index $(row.index) ($(row.source)) active$reverse$note")
        end
        isempty(duplicate_active) ? println("No duplicate reaction groups have multiple active rows.") :
            println("WARNING: $(length(duplicate_active)) duplicate reaction groups have multiple active rows.")
    end
end

function resolve_reaction_index(reaction, network)
    candidates, product_was_remapped = matching_rows_for_reaction(reaction, network)
    product_was_remapped = product_was_remapped || get(reaction, "product_was_remapped", false)
    isempty(candidates) && error("$(reaction["name"]): could not resolve reaction in networksetup.txt")
    selected = select_configured_row(reaction, candidates)
    selected === nothing && error("$(reaction["name"]): configured index $(getkey(reaction, "index")) was not found")
    println("$(reaction["name"]): using index $(selected.index) ($(selected.source))")
    product_was_remapped && println("  note: product was remapped by network boundaries: $(selected.line)")
    return selected.index
end

function copy_ppn(ppn_dir, dest)
    ispath(dest) && rm(dest; recursive=true)
    mkpath(dirname(dest))
    cp(ppn_dir, dest; force=true)
    npdata_src = realpath(joinpath(ppn_dir, "..", "NPDATA"))
    parent_npdata = joinpath(dirname(dest), "NPDATA")
    ispath(parent_npdata) || symlink(npdata_src, parent_npdata)
    local_npdata = joinpath(dest, "NPDATA")
    ispath(local_npdata) || symlink(npdata_src, local_npdata)
end

function write_physics_input(ppn_dir, run_dir, indices, factors)
    template = joinpath(ppn_dir, "ppn_physics.input")
    output = joinpath(run_dir, "ppn_physics.input")
    lines = readlines(template, keep=true)
    new_lines = String[]
    for line in lines
        if strip(line) == "/"
            for (i, (idx, factor)) in enumerate(zip(indices, factors))
                push!(new_lines, "        rate_index($i) = $idx\n")
                push!(new_lines, "        rate_factor($i) = $factor\n")
            end
        end
        push!(new_lines, line)
    end
    write(output, join(new_lines))
end

function create_factored_runs(nova; baseline_only=false, dry_run=false, runs_name="runs")
    base = nova_dir(nova)
    ppn_dir = joinpath(base, "ppn")
    runs_dir = joinpath(base, runs_name)
    plan_path = joinpath(base, "config", "reaction_plan.json")
    network_path = joinpath(ppn_dir, "networksetup.txt")
    config = parse_json_file(plan_path)
    network = parse_networksetup(network_path)
    default_factors = get(config, "default_factors", Any[])
    baseline_dir = joinpath(runs_dir, "baseline")
    if dry_run
        println("Would build baseline run in $baseline_dir")
    else
        mkpath(runs_dir)
        copy_ppn(ppn_dir, baseline_dir)
        println("Built baseline run in $baseline_dir")
    end
    baseline_only && return
    created = 0
    for reaction in config["reactions"]
        index = resolve_reaction_index(reaction, network)
        reverse_index = validate_reverse_index(reaction, network, index)
        factors = get(reaction, "factors", default_factors)
        for factor in factors
            run_dir = joinpath(runs_dir, reaction["name"], "fact_$factor")
            indices = [index]
            factor_values = [factor]
            if reverse_index !== nothing
                push!(indices, reverse_index)
                push!(factor_values, factor)
            end
            if !dry_run
                copy_ppn(ppn_dir, run_dir)
                write_physics_input(ppn_dir, run_dir, indices, factor_values)
            end
            created += 1
        end
    end
    println("$(dry_run ? "Would build" : "Built") $created factored runs in $runs_dir")
end

end
