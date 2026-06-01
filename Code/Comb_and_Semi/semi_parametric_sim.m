function sim = semi_parametric_sim(L,p, N)
%
% Simulates N observations from a Gaussian copula with
% Bernoulli-lognormal mixed marginals following Algorithm 1 of
% Baviera-Manzoni 2026.
%
% Convention:
%   p_j = P(X_j > 0)
%   P(X_j = 0) = 1 - p_j
%   X_j | X_j > 0 ~ Lognormal(mu_j, sigma_j^2)
%
% INPUT:
%   R      : d x d Gaussian copula correlation matrix
%   mu     : mean vector
%   sigma  : std dev vector
%   p      : positive-claim probability vector
%   N      : number of simulations
%
% OUTPUT:
%   sim    : N x d matrix of simulated observations

%% Input checks

% Force row vectors
p     = p(:)';

% Dimension
d = numel(p);

%% Core

Z_std= randn(N,d);

Z_corr = L*Z_std';
U = (normcdf(Z_corr))';

sim = (1-p).*ones(N,d);

for j = 1:d
    active = U(:,j) > 1 - p(j);

    sim(active,j) = U(active,j);
end

end