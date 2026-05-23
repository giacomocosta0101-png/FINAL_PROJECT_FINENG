function rho = calibration_rho_zero_mixed(X)
%
% Calibrates Gaussian copula rho for lognormal positive marginals.
%
% X must already be partitioned outside the function:
%   d = 2 -> returns scalar rho_12
%   d = 3 -> returns 3x3 correlation matrix
%
% This function is precisely correct and follows what Massaria described

%% Input check

[~, d] = size(X);

if d ~= 2 && d ~= 3
    error('X must have either 2 or 3 columns.');
end

if any(X <= 0, 'all')
    error('X must contain only positive observations.');
end

%% Core

Y = log(X);

mu = mean(Y, 1);
sigma = sqrt(mean((Y - mu).^2, 1));   % MLE estimator

Z = (Y - mu) ./ sigma;

R = corr(Z, 'Type', 'Pearson');

%Output depending on d

switch d
    case 2
        rho = R(1,2);     % scalar

    case 3
        rho = R;          % full 3x3 correlation matrix
end

end