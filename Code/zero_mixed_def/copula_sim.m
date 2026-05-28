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

    if size(R,1) ~= size(R,2)
        error('You passed a non sqaure R matrix ');
    end

    if size(R,1) ~= d
        error('Dimension mismatch between R and mu/sigma.');
    end



%% Core simulation


Z = randn(N, d) * chol(R);

sim = exp(mu + sigma .* Z);

end
