1. julia tools/setup_network.jl --nova <case>          # only if rates/network changed
2. julia tools/run_ppn_sweep.jl --nova <case>           # nucleosynthesis, all reactions×factors
3. julia tools/decay_time_scan.jl --nova <case>         # find the right decay time
4. python3 baseline_decay_checker.py                    # confirms it against Table 4
5. julia tools/decay_ppn_sweep.jl --nova <case> --decay-time <chosen>   # decay everything
6. python3 make_ili02_table_with_ppn.py
7. python3 compare_PPN_v_ILI.py
8. python3 analyse_PPN_ILI.py
