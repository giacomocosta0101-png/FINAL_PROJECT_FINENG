# Lomax

Scope: extension of the Comb-Bernoulli backtest with **Lomax**
(Pareto type II) positive marginals instead of lognormal.

Motivation: the standard Comb-Bernoulli underestimates the 99% VaR
on this dataset because the lognormal upper tail is too light. The
Lomax has a power-law upper tail `P(X > x) ~ x^{-alpha}`, which
captures the heaviness of the empirical positive part more faithfully
while keeping the parameter count of the marginal unchanged
(two parameters `(alpha, lambda)` plus the atom `p`).

## Files in this folder

Calibration:

- `marginal_parameter_calibration_lomax.m` — per-marginal Lomax fit on
  the positive part. The shape parameter `alpha` is profiled out
  analytically (`alpha = n / sum(log(1 + x/lambda))`) and `lambda` is
  estimated by 1-D MLE via `fminsearch` on the concentrated
  log-likelihood. Returns `(alpha, lambda, p)` per risk class.
- `marginal_cdf_lomax.m` — function handle for the composite CDF
  `F(x) = (1 - p)*1{x >= 0} + p*(1 - (1 + x/lambda)^(-alpha))*1{x > 0}`,
  broadcasted over columns.

Simulation:

- `comb_bern_lomax_sim.m` — Algorithm 1 of Baviera-Manzoni 2026
  adapted to the Lomax marginal: draw correlated Gaussians, map to
  uniforms, set zero where the uniform falls below `1 - p_j`,
  otherwise invert the Lomax in closed form
  `X = lambda * ((1 - u)^(-1/alpha) - 1)`.

Backtest:

- `backtest_lomax.m` — out-of-sample VaR backtest, mirroring the
  structure of `Backtest/backtest.m` but with a single model
  (Comb-Bernoulli + Lomax) instead of three.
- `var_calc_lomax.m` — calibration + simulation step inside the
  backtest loop. Functional analogue of `Backtest/var_calc.m`.

Diagnostics:

- `plot_backtest_lomax.m` — visual diagnostics, same layout as
  `Backtest/plot_backtest.m` but tailored to a single model.
- `christoffersen_test_lomax.m` — Kupiec POF + Christoffersen
  Independence + Conditional Coverage tests on the Lomax exception
  series. Same statistics and conventions as
  `Backtest/christoffersen_test.m`.

## Pipeline

```
RunFINAL_PROJECT.m (Extra block)
  └── backtest_lomax
        ├── data_split                                 (training window)
        ├── var_calc_lomax
        │     ├── marginal_parameter_calibration_lomax
        │     ├── marginal_cdf_lomax                   (build U)
        │     ├── calibrate_model                      (rho from Comb_and_Semi)
        │     └── comb_bern_lomax_sim                  (N daily losses)
        └── (exception flags vs. realised Total)
```

## Calibration mechanics

For each column `j`:

1. Take the positive subset `xp = X(X(:,j) > 0, j)`.
2. Concentrate the Lomax log-likelihood by substituting
   `alpha = n / sum(log(1 + x/lambda))` into `log L(alpha, lambda)`.
3. Minimise the 1-D concentrated negative log-likelihood with
   `fminsearch` starting at `lambda_0 = mean(xp)`.
4. Recover `alpha` in closed form from the MLE relation above.

This is the standard "profile likelihood" trick for Lomax. It is more
robust than a joint 2-D MLE because the joint surface is ill-conditioned
when the tail is close to non-heavy (`alpha >> 2`).

## Things to know

- The CDF is strictly in `(0, 1)` for any finite `x` and any
  `p in (0, 1)`. No risk of `norminv(1)` blowing up downstream.
- The inversion in `comb_bern_lomax_sim` is in closed form, so the
  simulation is as fast as the lognormal variant. No interpolation, no
  extrapolation issues at the upper end.
- `alpha` controls the tail index: `alpha = 2` is the boundary where
  variance becomes infinite. If `alpha` drops too close to `2`, the
  simulated 99.5% VaR becomes very noisy. The calibration routine does
  not clamp `alpha`; if you observe very large bootstrap CIs, inspect
  the calibrated `alpha` first.
- The copula step is unchanged: `calibrate_model` from
  `../Comb_and_Semi/` is used as-is, because it depends only on
  `(U, p)` and not on the marginal family.

## Comparison with the lognormal Comb-Bernoulli

The two pipelines are interchangeable up to the marginal block:

| Step                     | Lognormal                              | Lomax                                  |
|--------------------------|----------------------------------------|----------------------------------------|
| Marginal parameters      | `(p, mu, sigma)`                       | `(p, alpha, lambda)`                   |
| Marginal calibration     | closed-form (`mean(log)`, `std(log)`)  | profile MLE (`fminsearch` on lambda)   |
| Marginal CDF             | `marginal_cdf`                         | `marginal_cdf_lomax`                   |
| Copula calibration       | `calibrate_model`                      | `calibrate_model` (same)               |
| Simulation               | `comb_bern_sim`                        | `comb_bern_lomax_sim`                  |
| Tail behaviour           | sub-exponential (light)                | power-law (heavy)                      |

If you want a direct head-to-head, run both backtests on the same
seed (`rng(762)` is the convention in `RunFINAL_PROJECT.m`) and
compare the Christoffersen output at the 99% level.
