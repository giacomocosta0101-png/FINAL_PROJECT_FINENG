# VaR Context

Scope: this note is only for the VaR path.

Files to keep in mind:

- `var_calc.m`
- `calibr_wrapper.m`
- `mat_sim.m`
- `semi_parametric_losses.m`
- `../Comb_and_Semi/*`
- `../zero_mixed/*`

Pipeline:

- `var_calc` slices the timetable and builds `X = [Building Contents Profits]`.
- `calibr_wrapper` fits 3 models on the same `X`.
- `mat_sim` simulates `N` total losses for each model.
- VaR is `quantile(sim_losses{i}, 1 - alpha)`.

Model summary:

- Zero-mixed: draw the active set first, then simulate positive severities inside that set.
- Comb-Bernoulli: Bernoulli zero mass + lognormal positive severity + Gaussian copula.
- Semi-parametric: same zero mass logic, but positive severities are recovered from an empirical CDF through `semi_parametric_losses`.

Important detail for model 3:

- the inverse is not stepwise empirical
- it uses `interp1(..., 'linear', 'extrap')`
- so the top tail depends a lot on the gap between the largest observed losses

Why the semi-parametric VaR moves a lot with the seed:

- with a fixed calibration window, the seed changes only the Monte Carlo draw
- at `alpha = 0.005` and `N = 10000`, the VaR depends on about the worst 50 simulations
- the positive tail is sparse, especially `Profits`
- near the top observed losses, tiny changes in simulated `U` become large changes in simulated loss

Result:

- model 3 is much noisier than models 1 and 2 at the 99.5% level
- in the `1980-01-01` to `1983-12-31` window, across 40 seeds and `N = 10000`, the 99.5% VaR std was about:
- zero-mixed: `1.74`
- comb-bernoulli: `1.43`
- semi-parametric: `6.89`

Practical rules:

- set `rng(seed)` before `var_calc` if you want reproducible comparisons
- `N = 10000` is too small for a stable 99.5% semi-parametric VaR here
- if model 3 looks strange, inspect `semi_parametric_losses.m` first
