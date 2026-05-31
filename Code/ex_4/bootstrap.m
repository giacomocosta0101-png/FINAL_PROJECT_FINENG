function [rho_CI, p_CI, rho_hat, pi_hat]= bootstrap(rho,p,mu,sigma,model,N,B,alpha)
% BOOTSTRAP  Parametric bootstrap for the Comb-Bernoulli copula parameters.
%
%   1. simulate B replicas of length N via Algorithm 1,
%      using the point estimates rho, p, mu, sigma;
%   2. for each replica, re-estimate the marginal parameters (in the
%      Comb-Bernoulli case) and the copula correlation via IFM;
%   3. construct percentile confidence intervals from the empirical
%      quantiles of the bootstrap distribution, with Bonferroni
%      correction.
%
% Marginal-model branches:
%   "Comb-Bernoulli"  : full bootstrap, marginal parameters are re-fitted
%                       on every replica (lognormal MLE), so the CI
%                       captures both marginal and copula uncertainty.
%   "Semi-parametric" : the simulator returns Gaussian-copula uniforms
%                       directly (Algorithm 1, Step 1).
%
% INPUT
%   rho   : (1 x m) point estimate of the off-diagonal correlations,
%           m = d(d-1)/2; ordering as squareform()
%   p     : (1 x d) jump probabilities  p_i = P(X_i > 0)
%   mu    : (1 x d) lognormal location parameters (ignored if SP)
%   sigma : (1 x d) lognormal scale parameters   (ignored if SP)
%   model : "Comb-Bernoulli" or "Semi-parametric"
%   N     : replica length (number of observations per bootstrap sample)
%   B     : number of bootstrap replicas
%   alpha : nominal significance level for the CIs (e.g. 0.05 -> 95% CI)
%
% OUTPUT
%   rho_CI  : (m x 2) Bonferroni-corrected percentile CIs for each rho_ij
%   p_CI    : (m x 2) Bonferroni-corrected percentile CIs for each p_i
%   rho_hat : (B x m) full bootstrap distribution of the rho estimates
%   pi_hat  : (B x m) full bootstrap distribution of the p estimates
%
% NOTE
%   Set the random seed via rng() *before* calling this function; the
%   bootstrap loop relies on randn inside comb_bern_sim / semi_parametric_sim.

if strcmp(model,'Comb-Bernoulli')
    flag = 1;
    generate_replica = @(corr)comb_bern_sim(corr, mu, sigma, p, N);
elseif strcmp(model,'Semi-parametric')
    generate_replica = @(corr)semi_parametric_sim(corr,p,N);
    flag = 2;
end

d = numel(p);
m = d*(d-1)/2;           % number of correlations
rho_hat = zeros(B, m);
pi_hat = zeros(B,m);
p_calibr = p;

rho_CI = zeros(m,2);
p_CI = zeros(m,2);

alpha_correct = alpha / m;

t0 = tic;

parfor i = 1:B

    replica = generate_replica(rho);

    if flag == 1
        % Re-estimation of the parameters from the simulated replica:
        [p_hat,mu_hat,sigma_hat]= marginal_parameter_calibration(replica);
        % Update the cdf with the calibrated parameters:
        cdf = marginal_cdf(mu_hat,sigma_hat,p_hat);
        U = cdf(replica);
        p_calibr = p_hat;
    else
        % As explained in Algorithm 1, we consider the set of Uniform as
        % the 'simulated' ecdf of the replica:

        U = replica;
    end

    [rho_hat(i,:),pi_hat(i,:)] = calibrate_model(U,p_calibr);

    if mod(i, 50) == 0 || i == B   % Prints every 50 iterations
        elapsed = toc(t0);
        eta = elapsed / i * (B - i);
        fprintf('Iter %4d/%d | trascorso %6.1fs | ETA %6.1fs\n', ...
            i, B, elapsed, eta);
    end

end

for i = 1:3
    rho_CI(i,:) = quantile(rho_hat(:,i), [alpha_correct/2, 1-alpha_correct/2]);
    p_CI(i,:) = quantile(pi_hat(:,i), [alpha_correct/2, 1-alpha_correct/2]);
end

end




