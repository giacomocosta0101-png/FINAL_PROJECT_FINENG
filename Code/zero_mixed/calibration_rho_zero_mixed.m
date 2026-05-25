function rho = calibration_rho_zero_mixed(X, mu, sigma)
%
% Calibrates the Gaussian copula correlation parameter for lognormal
% positive marginals in the zero-mixed model.
%
% X must already contain only the observations belonging to the active set.
%
% Inputs:
%   X      : N x d matrix of positive observations
%   mu     : 1 x d vector of fitted lognormal location parameters
%   sigma  : 1 x d vector of fitted lognormal scale parameters
%
% Output:
%   rho    : if d = 2, scalar rho_12
%            if d = 3, 3 x 3 correlation matrix
%
% The estimator is:
%   Y_ij = (log(X_ij) - mu_j) / sigma_j
%   R_hat = (Y' * Y) / N

%% Input check

[N, d] = size(X);

if d ~= 2 && d ~= 3
    error('X must have either 2 or 3 columns.');
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

%% Trasformation


Y = (log(X) - mu) ./ sigma;

R = (Y'*Y)/N;


%% Output depending on dimension

switch d
    case 2
        rho = R(1,2);

    case 3
        rho = R;
end

end
