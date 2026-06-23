# NuGrid-tools
Tools for running NuGrid PPN for multiple reaction rate factorization and general output analysis

# Dependencies
- Julia
- NuGrid PPN framework

Assuming novae are found in a directory e.g. nova_cases where in a particular nova i.e. nova_cases/nova has a directory containing ppn inputs ppn/, a reaction sweep according to a JSON file reaction can be run as 

`julia tools/run_ppn_sweep.jl` (add `--jobs` for parallel processing)

Outputs will be a runs/ directory containing all factored reactions and their ppn outputs. 

Typical sweeps of 180 reactions factored by 5 different amounts at once takes about 200 seconds. 

Modify source code as needed. 

# Install Julia

`curl -fsSL https://install.julialang.org | sh`
