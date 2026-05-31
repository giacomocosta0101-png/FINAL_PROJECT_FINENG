function R = calibration_rho_zero_mixed(X, mu, sigma)
%
% Calibrates the Gaussian copula correlation matrix for the positive
% marginals of a zero-mixed model.
%
% INPUT:
%   X      : N x d matrix of strictly positive observations
%   mu     : 1 x d vector of fitted lognormal location parameters
%   sigma  : 1 x d vector of fitted lognormal scale parameters
%
% OUTPUT:
%   R      : d x d Gaussian copula correlation matrix
%
% Note:
%   If N <= d, the function returns the identity matrix.
%
% The estimator is:
%   Y_ij = (log(X_ij) - mu_j) / sigma_j
%   R_hat = (Y' * Y) / N

%% Input check

[N, d] = size(X);
if d <2
    error('X must have at least 2 columns');
end

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
