# Higher-Dimensional Extension

Scope: dimension-agnostic version of the Comb-Bernoulli calibrator,
used to validate the model on synthetic `d > 3` datasets.

Motivation: the production code in `../Comb_and_Semi/` hard-codes the
spherical Cholesky parametrisation for `d = 3`. This folder
re-implements the same calibration pipeline using the unconstrained
Cholesky-based parametrisation of Huang-Ye-Wang (2025), which is well
defined for any `d`.

## Files in this folder

- `Run_Extra.m` — driver script. Loads the dataset, generates a
  synthetic fourth column to obtain `d = 4`, then runs the Zero-Mixed
  calibration + bootstrap and the generalised Comb-Bernoulli
  calibration as a stress test.
- `calibrate_model_generalized.m` — generalised IFM second step.
  Identical in spirit to `../Comb_and_Semi/calibrate_model.m` but with
  no hard-coded `d`. Returns the full `d x d` correlation matrix `R`
  (not the off-diagonal vector).
- `corr_from_cholesky_param.m` — implementation of the Huang-Ye-Wang
  parametrisation. Takes `h in R^{d(d-1)/2}` (unconstrained), maps it
  through `theta = pi / (1 + exp(-h))` and a recursive cosine/sine
  product (paper eq. 5) to a valid correlation matrix `R = L' * L`.
- `log_likelihood_trivariate_copula_generalized.m` — Gaussian copula
  log-likelihood for arbitrary `d`. Same patterns as the trivariate
  version (loop over unique active sets, explicit log density on the
  active block via Cholesky, conditional CDF on the inactive block) but
  parameterised through `corr_from_cholesky_param` rather than the
  hard-coded `d = 3` Cholesky factor.

## Parametrisation

For each pair `(i, j)` with `i < j` there is one free real parameter
`h_{i,j}`. The mapping to the correlation matrix is:

```
h_{i,j}  ->  theta_{i,j} = pi / (1 + exp(-h_{i,j}))  in (0, pi)
         ->  C_{i,j} = cos(theta_{i,j}),  S_{i,j} = sin(theta_{i,j})
         ->  L (upper triangular, paper eq. 5)
         ->  R = L' * L
```

The `h_{i,j}` are *not* correlations. They are unconstrained
optimisation variables; the actual correlations only come out of
`R = L' * L` after optimisation.

The parameter vector is ordered column-by-column:

```
d = 3  ->  h = [h_{1,2}; h_{1,3}; h_{2,3}]
d = 4  ->  h = [h_{1,2}; h_{1,3}; h_{2,3}; h_{1,4}; h_{2,4}; h_{3,4}]
```

## Pipeline

```
Run_Extra.m
  ├── (generate the fourth column of X)
  ├── zero_mixed_calibration                  (sanity check at d = 4)
  ├── zero_mixed_bootstrap
  ├── marginal_parameter_calibration          (lognormal MLE per column)
  ├── marginal_cdf -> U_CB                    (pseudo-observations)
  └── calibrate_model_generalized(U_CB, p)
        ├── log_likelihood_trivariate_copula_generalized
        │     └── corr_from_cholesky_param    (h -> R)
        └── fminunc                            (unconstrained over h)
```

## Initialisation and convergence

Two practical points learned the hard way:

1. **Do not start at `h = 0`**. At `h = 0` all angles equal `pi/2`, all
   cosines vanish, `L = I` and `R = I`. For `d >= 4` this point is a
   saddle: the gradient of the log-likelihood with respect to several
   `h_{i,j}` is exactly zero (compound products of `cos = 0`), and
   `fminunc` either does not move or moves only some entries. The
   resulting `R` will have exact zeros in some off-diagonals.
2. **Do not start at any constant `h * ones(...)`**. By symmetry the
   gradient at any all-equal starting point has small or zero
   components along several directions. The optimiser can declare
   first-order optimality at the initial point and return immediately.

The recommended initialisation is a small random perturbation, e.g.
`h0 = 0.5 * randn(1, d*(d-1)/2)`, combined with a multi-start (5 to 20
restarts) keeping the run with the lowest negative log-likelihood
among those that satisfy `exitflag > 0` and `output.iterations > 5`.

A defensive check after calibration:

```matlab
fprintf('min eig(R) = %.3e\n', min(eig(R)));
fprintf('|R - R''|   = %.2e\n', norm(R - R.', 'fro'));
```

## Relationship to the `d = 3` production code

`calibrate_model_generalized` is functionally a superset of
`../Comb_and_Semi/calibrate_model.m`. For `d = 3` it should produce a
correlation matrix essentially identical to the spherical-Cholesky
output of the production code (up to optimiser tolerances). It is kept
in a separate folder because:

1. The trivariate spherical parametrisation in production is faster and
   numerically slightly different.
2. The generalised version requires multi-start and additional guards
   to be reliable in `d > 3`; folding all of this into the production
   path would have hidden the simpler trivariate logic.

For the report, the recommended convention is: use
`../Comb_and_Semi/calibrate_model.m` for the main `d = 3` results, use
`calibrate_model_generalized` only inside `Run_Extra.m` to demonstrate
that the pipeline extends to higher dimensions.
