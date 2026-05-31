function sim = copula_sim(R, mu, sigma, N)
%
% Simulates observations from a Gaussian copula with lognormal marginals.
%
% INPUT:
%   R      : d x d Gaussian copula correlation matrix
%   mu     : 1 x d vector of lognormal location parameters
%   sigma  : 1 x d vector of lognormal scale parameters
%   N      : number of simulations
%
% OUTPUT:
%   sim    : N x d matrix of simulated observations

arguments
    R (:,:) double {mustBeNonempty, mustBeReal, mustBeFinite}
    mu double {mustBeNonempty, mustBeReal, mustBeFinite, mustBeVector}
    sigma double {mustBeNonempty, mustBeReal, mustBeFinite, mustBeVector, mustBePositive}
    N (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
end

%% Initial setup
mu = mu(:)';
sigma = sigma(:)';

d = numel(mu);


%% Core simulation
Z = randn(N, d) * chol(R);

sim = exp(mu + sigma .* Z);

end
