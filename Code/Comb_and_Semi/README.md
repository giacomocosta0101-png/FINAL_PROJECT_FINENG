# Comb-Bernoulli and Semi-Parametric

Scope: shared infrastructure for the two parsimonious copula models.

Both models share the same simulation skeleton (Bernoulli activation +
Gaussian copula on the active block) and the same calibration of the
copula correlation; the only difference is the choice of marginal CDF.

## Files in this folder

Calibration of marginals and the copula:

- `marginal_parameter_calibration.m` — closed-form MLE of
  `(p_j, mu_j, sigma_j)` for each risk class. `p_j` is the empirical
  fraction of positive observations; `(mu_j, sigma_j)` are the
  lognormal MLE on the positive part.
- `marginal_cdf.m` — function handle for the composite CDF
  `F(x) = (1 - p)*1{x >= 0} + p*Phi((log x - mu)/sigma)*1{x > 0}`.
  Returns a single broadcasted handle; can be applied to an `N x d`
  matrix in one shot.
- `cumulative_cdf_semi_parametric_vec.m` — function handle for the
  Semi-Parametric composite CDF: same atom at zero, empirical CDF on
  the positive part. Uses denominator `n + 1` so the empirical CDF
  never reaches 1 on observed data.
- `calibrate_model.m` — IFM second step. Takes pseudo-observations `U`
  and atom probabilities `p`, maps `Y = Phi^{-1}(U)`, maximises the
  Gaussian copula log-likelihood via `fminunc` using the spherical
  Cholesky parametrisation (`d = 3` hard-coded; for general `d` use
  `../higher_dim/calibrate_model_generalized.m`). Returns the
  off-diagonal correlation vector `rho` and a re-estimated atom
  probability `p_hat`.

  **Note**: the routine is model-agnostic in its inputs (`U` and `p`),
  with no dependence on marginal parameters. It is therefore reused
  without modification for both Comb-Bernoulli and Semi-Parametric.

- `log_likelihood_trivariate_copula.m` — Gaussian copula log-likelihood
  for the censored trivariate problem, with `d = 3`. The likelihood is
  assembled by looping over unique active-set patterns rather than over
  rows, via `[C, ia, ic] = unique(state_matrix, 'rows')`. The Gaussian
  copula density on the active block is computed explicitly (no
  `mvnpdf`, no `det`, no `inv`); for `|T| = 2` the bivariate CDF is
  computed via Drezner-Wesolowsky with cached Gauss-Legendre nodes.

Simulation:

- `comb_bern_sim.m` — Algorithm 1 of Baviera-Manzoni 2026. Draws
  `N x d` Gaussian via Cholesky, maps to uniforms with `normcdf`,
  inverts the lognormal marginal where active and sets zero where the
  uniform falls below `1 - p_j`. Takes the Cholesky factor `L`
  precomputed by the caller.
- `semi_parametric_sim.m` — same activation logic, but the output is
  the simulated uniforms `U_sim` (not yet inverted). Inversion to
  losses is left to the caller (typically `Backtest/semi_parametric_losses.m`).

Bootstrap (point-calibration uncertainty):

- `bootstrap.m` — `B` parametric replicas of the full calibration
  pipeline. For each replica: simulate a new dataset, refit marginals
  and copula, store the estimates. At the end, builds componentwise
  bootstrap confidence intervals with Bonferroni correction at level
  `alpha/m`.
- `plot_bootstrap_rho.m` — visual diagnostics for the bootstrap
  distribution of `rho` (one histogram per pair).

## Calibration pipeline (Comb-Bernoulli)

```
RunFINAL_PROJECT.m
  ├── marginal_parameter_calibration   -> (p, mu, sigma)
  ├── marginal_cdf(mu, sigma, p)        -> handle F
  ├── U_CB = F(X)                       -> pseudo-observations
  └── calibrate_model(U_CB, p)
        └── log_likelihood_trivariate_copula  -> rho_CB
```

## Calibration pipeline (Semi-Parametric)

```
RunFINAL_PROJECT.m
  ├── marginal_parameter_calibration   -> (p, mu, sigma)   [only p used here]
  ├── cumulative_cdf_semi_parametric_vec(p, X) -> cell of handles
  ├── U_SP = cellfun(handle, X columns)  -> pseudo-observations
  └── calibrate_model(U_SP, p)
        └── log_likelihood_trivariate_copula  -> rho_SP
```

## Bootstrap pipeline

```
RunFINAL_PROJECT.m
  └── bootstrap(rho_hat, p, mu, sigma, model, N, B, alpha)
        ├── (per replica) comb_bern_sim or semi_parametric_sim
        ├── marginal_parameter_calibration (refit)
        ├── marginal_cdf or cumulative_cdf_semi_parametric_vec (refit)
        └── calibrate_model (refit copula)
```

## Conventions

- All functions assume non-negative observations. Zeros are the atom of
  the Bernoulli component, not missing values.
- `U` matrices have one row per observation, one column per risk class.
- `rho` is an off-diagonal vector of length `d(d-1)/2` ordered as
  `[rho_12, rho_13, rho_23]` for `d = 3`. To recover the full
  correlation matrix use `R = squareform(rho) + eye(d)`.
- All routines except `calibrate_model` and `log_likelihood_trivariate_copula`
  are written generically in `d`. The two exceptions hard-code the
  `d = 3` spherical Cholesky parametrisation; for `d > 3` use the
  generalised counterparts in `../higher_dim/`.
