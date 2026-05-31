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

arguments
    X (:,:) double {mustBeNonempty, mustBeReal, mustBeFinite, mustBePositive}
    mu double {mustBeNonempty, mustBeReal, mustBeFinite, mustBeVector}
    sigma double {mustBeNonempty, mustBeReal, mustBeFinite, mustBeVector, mustBePositive}
end

%% Initial setup
[N, d] = size(X);

mu = mu(:)';
sigma = sigma(:)';


%% Compute the estimator 

if N > d

    Y = (log(X) - mu) ./ sigma;
    R = (Y'*Y)/N;

else

    R = eye(d);


end
