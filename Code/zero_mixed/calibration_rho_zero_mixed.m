function R = calibration_rho_zero_mixed(X, mu, sigma)
%
% Calibrates the Gaussian copula correlation parameter for lognormal
% positive marginals in the zero-mixed model.
%
% Inputs:
%   X      : N x d matrix of strictly positive observations
%   mu     : 1 x d vector of fitted lognormal location parameters
%   sigma  : 1 x d vector of fitted lognormal scale parameters
%
% Output:
%   R      : if d = 2, 2x2 correlation matrix
%            if d = 3, 3 x 3 correlation matrix
%
% Note: if N<=d return the identiy matrix
%
% The estimator is:
%   Y_ij = (log(X_ij) - mu_j) / sigma_j
%   R_hat = (Y' * Y) / N

%% Input check

[N, d] = size(X);

if any(X <= 0, 'all')
    error('X must contain only positive observations.');
end

mu = mu(:)';
sigma = sigma(:)';

if numel(mu) ~= d || numel(sigma) ~= d
    error('mu and sigma must have one entry per column of X.');
end

if any(sigma <= 0)
    error('All sigma values must be strictly positive.');
end


%% Compute the estimator 

if N > d

    Y = (log(X) - mu) ./ sigma;
    R = (Y'*Y)/N;

else

    R = eye(d);


end
