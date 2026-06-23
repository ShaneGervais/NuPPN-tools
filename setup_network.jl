#!/usr/bin/env julia

include(joinpath(@__DIR__, "NovaRunTools.jl"))
using .NovaRunTools

function parse_args(args)
    nova = "nova_test"
    check = false
    validate_factor_plan = false
    i = 1
    while i <= length(args)
        if args[i] == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            nova = args[i]
        elseif args[i] == "--check"
            check = true
        elseif args[i] == "--validate-factor-plan"
            validate_factor_plan = true
        elseif args[i] in ("-h", "--help")
            println("Usage: julia tools/setup_network.jl [--nova nova_test] [--check] [--validate-factor-plan]")
            exit(0)
        else
            error("Unknown argument: $(args[i])")
        end
        i += 1
    end
    return nova, check, validate_factor_plan
end

nova, check, validate_factor_plan = parse_args(ARGS)
setup_network(nova; check=check, validate_factor_plan=validate_factor_plan)
