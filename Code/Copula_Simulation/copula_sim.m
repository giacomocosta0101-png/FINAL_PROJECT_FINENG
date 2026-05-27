function sim = copula_sim(R, mu, sigma, N)
%
% Simulates N observations from a Gaussian copula with lognormal marginals.
%
% INPUT:
%   R      : d x d Gaussian copula correlation matrix
%            or scalar rho in the bivariate case
%   mu     : mean vector of lognormal marginals
%   sigma  : std dev vector of lognormal marginals
%   N      : number of simulations
%
% OUTPUT:
%   sim    : N x d matrix of simulated observations

%% Input handling and Checks


mu = mu(:)';
sigma = sigma(:)';

d = numel(mu);

if numel(sigma) ~= d
    error('mu and sigma must have the same length.');
end


%I had to add this since bc of tecnical mvnrnd reasons
% Basically if you pass a scalar correlation to mvnrnd it read it as a
% single sigma for 2 independet marginals

% If R is scalar, interpret it as bivariate correlation rho
if isscalar(R)

    rho = R;
    R = [1 rho; rho 1];

else

    if size(R,1) ~= size(R,2)
        error('You passed a non sqaure R matrix ');
    end

    if size(R,1) ~= d
        error('Dimension mismatch between R and mu/sigma.');
    end

end

%% Core simulation


Z = randn(N, d) * chol(R);

sim = exp(mu + sigma .* Z);

end