function sim = comb_bern_sim(R, mu, sigma, p, N)
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

if ~ismatrix(R) || size(R,1) ~= size(R,2)
    error('R must be a square matrix.');
end

% Force row vectors
mu    = mu(:)';
sigma = sigma(:)';
p     = p(:)';

% Dimension
d = size(R,1);

%% Core

Z = mvnrnd(zeros(1,d), R, N);
U = normcdf(Z);

sim = zeros(N,d);

for j = 1:d
    active = U(:,j) > 1 - p(j);

    U_pos = (U(active,j) - (1 - p(j))) ./ p(j);
    sim(active,j) = exp(mu(j) + sigma(j) .* norminv(U_pos));
end

end