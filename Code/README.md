# Code

MATLAB implementation of Project 6 — copula calibration and VaR
backtest on the Danish multivariate fire-loss dataset.

## Entry point

`RunFINAL_PROJECT.m` runs the full pipeline end-to-end:

1. Loads `danishmulti.csv` and builds `X = [Building Contents Profits]`.
2. Calibrates the three models on the full sample (point estimates).
3. Runs the parametric bootstrap for Comb-Bernoulli and
   Semi-Parametric, and a custom bootstrap for Zero-Mixed.
4. Runs the Static and Rolling-Window VaR backtests.
5. Runs the Christoffersen / Kupiec tests on the exceptions.
6. Re-runs the backtest with Lomax marginals as an extension.

The script adds the model folders to the path at the top
(`addpath('Comb_and_Semi', 'zero_mixed', 'Backtest', 'higher_dim', 'Lomax')`).

## Folder map

| Folder            | Scope                                                                 |
|-------------------|-----------------------------------------------------------------------|
| `Comb_and_Semi/`  | Comb-Bernoulli and Semi-Parametric models (shared copula calibration) |
| `zero_mixed/`     | Zero-Mixed model (per-active-set calibration)                         |
| `Backtest/`       | Out-of-sample VaR backtest, Static and Rolling-Window modes           |
| `Lomax/`          | Extension: Comb-Bernoulli with Lomax marginals                        |
| `higher_dim/`     | Extension: Comb-Bernoulli at arbitrary `d`                            |

Each folder has its own `README.md` with a detailed file-by-file map
and pipeline diagram.

## Data

`danishmulti.csv` — daily insurance claims for three risk classes
(Building, Contents, Profits) over 1980–1990. Loaded by
`readDataset.m` into a MATLAB `timetable` indexed by Date.

## Models in one sentence each

- **Zero-Mixed**: one full calibration per active set, so `2^d`
  calibrations in total. Maximally specification-sensitive; parameter
  count grows exponentially in `d`.
- **Comb-Bernoulli**: shared lognormal marginals across active sets;
  Gaussian copula on the active block. Parameter count `O(d^2)`.
- **Semi-Parametric**: same skeleton as Comb-Bernoulli but with an
  empirical CDF on the positive part instead of a parametric lognormal.
  Better tail fit; noisier 99.5% VaR because of the empirical
  extrapolation at the top.
- **Lomax extension**: Comb-Bernoulli with Lomax marginals
  (power-law tail) instead of lognormal.

## VaR backtest at a glance

Two modes, both implemented by `Backtest/backtest.m`:

- `Fixed`: one calibration on `[window_start, window_end]`, reused for
  every day of the evaluation window.
- `Rolling-window`: at every day, shift both ends of the training
  window by one day and recalibrate.

Christoffersen / Kupiec tests in `Backtest/christoffersen_test.m`
output three p-values per model and per level (POF, Independence,
Conditional Coverage).

## Reproducibility

`rng(762)` is the convention used throughout `RunFINAL_PROJECT.m` for
all the calls that involve Monte Carlo simulation. Always set the seed
immediately before the call if you want to reproduce the exact numbers
in the report.

## Conventions

- `X` is always `N x d`, non-negative, with zeros encoding the atom of
  the marginal (not missing values).
- `U` matrices have the same shape, with entries in `[1 - p_j, 1]` for
  positive observations and exactly `1 - p_j` on the atom.
- `rho` denotes the off-diagonal correlation vector of length
  `d(d-1)/2`; `R` denotes the full `d x d` correlation matrix. Convert
  with `R = squareform(rho) + eye(d)`.
- `exceptions` is a `1 x 3` cell of logical `M x nAlpha` matrices, one
  per model.
- All code is in MATLAB; no external dependencies beyond Statistics
  and Optimization Toolboxes.
