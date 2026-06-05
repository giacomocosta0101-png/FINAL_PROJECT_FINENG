# Backtest

Scope: out-of-sample VaR backtest for the three claim models
(Zero-Mixed, Comb-Bernoulli, Semi-Parametric).

## Files in this folder

- `backtest.m` — main entry point. Splits the dataset into a training
  window and an evaluation window, then computes VaR and exception flags
  in `Fixed` or `Rolling-window` mode.
- `var_calc.m` — slices the timetable on the training window, calls
  `calibr_wrapper`, runs `mat_sim` and returns the `(1 - alpha)`
  quantile of the simulated aggregate losses for each model.
- `calibr_wrapper.m` — fits the three models on the same `X`. Returns a
  `3 x 1` cell of parameter structs. The Semi-Parametric struct also
  carries `X` (the calibration sample), because the simulation needs the
  full empirical CDF.
- `mat_sim.m` — simulates `N` aggregate daily losses for each model.
  Returns a `3 x 1` cell of `N x 1` vectors. Cholesky factors are
  computed once and reused across `N` draws.
- `semi_parametric_losses.m` — Semi-Parametric simulation. Draws
  pseudo-uniforms from the Gaussian copula and inverts them through the
  empirical CDF.
- `data_split.m` — slices a timetable by date range. Used both for the
  training window inside `var_calc` and for the evaluation window inside
  `backtest`.
- `plot_backtest.m` — visual diagnostics. One subplot per model with
  realised losses, VaR levels and exception markers.
- `christoffersen_test.m` — Kupiec POF + Christoffersen Independence +
  Conditional Coverage tests on the exception series.

External dependencies (added to the path by `RunFINAL_PROJECT.m`):

- `../Comb_and_Semi/*` for marginal CDF, calibrate_model and the
  Comb-Bernoulli / Semi-Parametric simulators.
- `../zero_mixed/*` for the Zero-Mixed calibration and simulator.

## Pipeline

```
RunFINAL_PROJECT.m
   └── backtest.m
         ├── data_split.m                 (evaluation window)
         └── var_calc.m
               ├── data_split.m           (training window)
               ├── calibr_wrapper.m       (three model calibrations)
               └── mat_sim.m              (three N-replica simulations)
```

Then:

```
RunFINAL_PROJECT.m
   ├── plot_backtest.m                    (visual diagnostics)
   └── christoffersen_test.m              (formal tests)
```

## Modes

- `Fixed`: calibrate once on `[window_start, window_end]`, reuse the
  same VaR across every day of the evaluation window. Output `VaR` is
  `3 x nAlpha`.
- `Rolling-window`: at every day `i` of the evaluation window, shift
  both ends of the training window by one day and recalibrate. Output
  `VaR` is `3 x nAlpha x M` where `M = height(backtest_window)`.

`plot_backtest` and `christoffersen_test` both branch on `mode` to
handle the two output shapes.

## Output conventions

- `backtest_window`: timetable of the evaluation period (from the day
  after `window_end` to the end of `data`).
- `exceptions`: `1 x 3` cell, one entry per model. Each is `M x nAlpha`
  logical, true where the realised `Total` strictly exceeds the model's
  VaR.
- `VaR`: see the two modes above.

## Things to know about Semi-Parametric

Inversion is done by `interp1` in `semi_parametric_losses.m:34` with
options `'linear', 'extrap'`. There are two consequences worth keeping
in mind:

1. The inverse is piecewise linear between observed positive losses, not
   stepwise. This is intentional: a step inverse would produce a finite
   set of simulated severities, with implausibly thick atoms.
2. With `'extrap'`, simulated uniforms that fall above the maximum
   observed quantile are linearly extrapolated past `max(X(:,j))`. This
   can produce simulated losses larger than anything in the calibration
   window, sometimes by a noticeable margin if the last two observations
   are far apart. This is the dominant source of variance in the
   Semi-Parametric 99.5% VaR.

## Empirical reproducibility notes

- Always set `rng(seed)` immediately before `backtest` if you want
  reproducible comparisons. The two calls in `RunFINAL_PROJECT.m` use
  `rng(762)` for exactly this reason.
- The default `N` in `RunFINAL_PROJECT.m` is `1e6`, which is enough to
  stabilise the 99.5% Semi-Parametric VaR on this dataset. Earlier
  experiments at `N = 1e4` showed Semi-Parametric VaR standard deviation
  roughly 4x that of Zero-Mixed and Comb-Bernoulli at 99.5%: with the
  current `N` the gap narrows substantially, but Semi-Parametric remains
  the noisiest of the three because of the tail extrapolation above.
- If the Semi-Parametric VaR looks odd at a given seed, inspect
  `semi_parametric_losses.m` first (the `interp1` line) and check the
  range of the simulated `U_sim(:,j)` against `max(cdf{j}(X_col))`.

## Christoffersen / Kupiec tests

`christoffersen_test.m` consumes `exceptions` (whatever mode they were
produced under) and prints a table with three p-values per model and
per confidence level:

- `p_POF`: Kupiec unconditional coverage. H0 — observed exception rate
  equals the nominal level.
- `p_ind`: Christoffersen independence. H0 — Markov-1 transition
  probabilities `pi_01 = pi_11 = pi`.
- `p_cc`: joint conditional coverage. H0 — both above hold.

Reject at 5%. The implementation guards `0 * log(0)` via local `xlogy`
and divisions by zero via `safe_div`, so degenerate exception series
(no exceptions, all exceptions, no transitions out of one state) do not
crash the test.
