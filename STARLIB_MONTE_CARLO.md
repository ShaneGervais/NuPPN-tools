# STARLIB Monte Carlo Sweep Plan

This project currently has deterministic factor sweeps, where one reaction rate is multiplied by fixed factors such as `2.0` or `0.5`. STARLIB makes a better statistical workflow possible because STARLIB rates include uncertainty information.

## Goal

Build many PPN runs where STARLIB is enabled and reaction rates from the existing `reaction_plan.json` are perturbed simultaneously with random multiplicative factors. The outputs can then be analyzed as an ensemble instead of one-at-a-time sensitivity tests.

## Practical First Implementation

1. Build a baseline run with STARLIB enabled.
2. Build sample directories:

   ```text
   runs_star_mc/
     baseline/
     sample_000001/
     sample_000002/
     ...
   ```

3. For each sample, resolve the reaction-plan entries to PPN `networksetup.txt` indices.
4. Draw one random factor per unique network index.
5. Write the sampled factors into `ppn_physics.input` using PPN's existing `rate_index(i)` / `rate_factor(i)` mechanism.
6. Write `starlib_option = 1` or `2` into `ppn_physics.input`.
7. Save a `mc_manifest.json` containing the seed, STARLIB option, reaction names, network indices, and sampled factors.
8. Run all sample `ppn.exe` jobs in parallel.

## STARLIB Options

From the current NuPPN Fortran source:

- `starlib_option = 1`: uses `starlib_mc10_mc13_082022.txt` and `sunet_mc10_mc13_082022.dat`; internally labelled `STL01`.
- `starlib_option = 2`: uses `atomki_rates.dat` and `sunet_taly_012025_atom.dat`; internally labelled `STL02`.

## Sampling Model

The initial script samples lognormal factors:

```text
factor = exp(sigma_ln * randn())
```

This is a pragmatic placeholder. The better long-term version should read STARLIB's native uncertainty fields and sample each rate from its own STARLIB uncertainty distribution.

## Caveats

- The current PPN network has generic `AL 26`, not separate `26Alg` and `26Alm`.
- If multiple plan entries map to the same network index, the MC script samples one factor for that index and records all mapped reaction names.
- PPN has a finite number of `rate_index/rate_factor` slots. If a sample tries to vary too many unique indices for the compiled settings, reduce the reaction set or extend the compiled PPN knob dimensions.
- The first implementation perturbs selected PPN rates multiplicatively while STARLIB is enabled. It does not yet parse native STARLIB distribution parameters.

## Example

```bash
julia tools/run_ppn_star_mc.jl \
  --nova ne_nova_1.15_12_X_weiss_mixed \
  --samples 1000 \
  --jobs 12 \
  --starlib-option 1 \
  --runs-name runs_star_mc
```

For a build-only smoke test:

```bash
julia tools/run_ppn_star_mc.jl \
  --nova ne_nova_1.15_12_X_weiss_mixed \
  --samples 5 \
  --build-only
```
