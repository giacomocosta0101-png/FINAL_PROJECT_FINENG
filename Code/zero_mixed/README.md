# Zero-Mixed

Scope: Zero-Mixed model — one calibration per active set, so `2^d`
calibrations in total for `d` risk classes.

For each subset `S` of the `d` risk classes, the model holds a separate
copy of marginal parameters and a separate Gaussian copula correlation
matrix on the active block. Calibration, simulation and bootstrap are
organised around this `2^d` decomposition.

## Files in this folder

Data preparation:

- `zero_mixed_case_data.m` — precomputes shared active-set metadata
  for an input `X`: the list of active sets `1..2^d`, the row indices
  belonging to each set, the sample size per set and the empirical
  probability of each set. Reused by both calibration and simulation
  to avoid re-enumerating subsets.
- `get_active_sets.m` — utility to enumerate the `2^d` active subsets
  of `{1, ..., d}`.

Calibration:

- `zero_mixed_calibration.m` — main calibration entry point. For each
  active set `S`, fits marginal parameters `(mu_S, sigma_S)` on the
  observations belonging to `S` and, when `|S| >= 2`, calibrates the
  Gaussian copula correlation matrix `R_S` on the active block via
  `calibration_rho_zero_mixed`.
- `calibration_rho_zero_mixed.m` — copula calibration on a single
  active block. Maximises the Gaussian copula log-likelihood using the
  spherical Cholesky parametrisation. Equivalent to the IFM second step
  in `Comb_and_Semi/calibrate_model.m`, but restricted to the rows of
  the active set and without the atom contribution (the conditioning on
  the active set has already happened at the data-splitting stage).
- `zero_mixed_unpacking.m` — flattens the `1 x 2^d` cell of calibrated
  structs into convenient vectors / matrices for downstream use
  (bootstrap CI tables, simulation, etc.).

Simulation:

- `copula_sim.m` — Gaussian copula sample on the active block of a
  given active set. Internal helper.
- `zero_mixed_sim.m` — main simulator. For each row of the output,
  draws an active set according to the calibrated empirical
  probabilities, then conditional on the active set simulates the
  active marginals via the Gaussian copula.
- `zero_mixed_sim_fixed.m` — variant where the active-set counts are
  fixed at their calibrated values (no sampling variability on
  occurrence). Used inside `zero_mixed_bootstrap_fixed`.

Bootstrap:

- `zero_mixed_bootstrap.m` — `B` parametric replicas of the full
  calibration with re-sampled active-set counts. Returns componentwise
  Bonferroni-corrected confidence intervals.
- `zero_mixed_bootstrap_fixed.m` — same `B`-replica scheme but with
  fixed active-set counts. Isolates the uncertainty due to severity
  estimation, holding the occurrence structure constant. Useful as a
  diagnostic to see how much of the marginal CI is driven by
  occurrence resampling versus severity resampling.

Reporting:

- `zero_mixed_print_ci_table.m` — pretty-prints a Bonferroni-corrected
  confidence interval table for the calibrated parameters.

## Pipeline

Calibration:

```
RunFINAL_PROJECT.m
  └── zero_mixed_calibration
        ├── zero_mixed_case_data           (active-set metadata)
        └── (per active set with |S| >= 2)
              └── calibration_rho_zero_mixed
```

Simulation:

```
mat_sim.m (Backtest)
  └── zero_mixed_sim
        ├── (sample active set)
        └── copula_sim                     (Gaussian copula draw)
```

Bootstrap:

```
RunFINAL_PROJECT.m
  └── zero_mixed_bootstrap
        ├── (per replica) zero_mixed_sim
        └── zero_mixed_calibration         (refit on the replica)
  └── zero_mixed_print_ci_table
```

## Parameter count

The model has `O(2^d)` parameter blocks. For `d = 3` (this project):
`8` active sets, `7` non-empty marginal blocks, and `4` blocks with
`|S| >= 2` that carry a copula correlation matrix. The full parameter
vector contains roughly 30 numbers, against the 9 of the parsimonious
Comb-Bernoulli specification.

## Conventions

- `X` is `N x d`, non-negative, zeros encode the atom of the marginal
  (not missing values).
- Active sets are encoded as logical row vectors of length `d` or as
  integer index arrays; both representations are used internally,
  `case_data` is the canonical source.
- Output of `zero_mixed_calibration` is a `1 x 2^d` cell of structs.
  The struct format is documented in the function header.
- This folder is fully generic in `d`. The only `d`-specific code in
  the project lives in `Comb_and_Semi/` (Cholesky parametrisation
  hard-coded at 3) and is bypassed here.
