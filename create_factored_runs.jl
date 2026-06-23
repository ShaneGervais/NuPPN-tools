#!/usr/bin/env julia

include(joinpath(@__DIR__, "NovaRunTools.jl"))
using .NovaRunTools

function parse_args(args)
    nova = "nova_test"
    baseline_only = false
    dry_run = false
    i = 1
    while i <= length(args)
        if args[i] == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            nova = args[i]
        elseif args[i] == "--baseline-only"
            baseline_only = true
        elseif args[i] == "--dry-run"
            dry_run = true
        elseif args[i] in ("-h", "--help")
            println("Usage: julia tools/create_factored_runs.jl [--nova nova_test] [--baseline-only] [--dry-run]")
            exit(0)
        else
            error("Unknown argument: $(args[i])")
        end
        i += 1
    end
    return nova, baseline_only, dry_run
end

nova, baseline_only, dry_run = parse_args(ARGS)
create_factored_runs(nova; baseline_only=baseline_only, dry_run=dry_run)
