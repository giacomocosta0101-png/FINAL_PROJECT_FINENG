function [p,mu,sigma] = marginal_parameter_calibration(X)
% MARGINAL_PARAMETER_CALIBRATION  Calibrate lognormal + prob mass in zero 
% distribution parameters.
%
%   INPUT:
%       X: full data 
%
%   OUTPUT: struct() containing:
%
%   The estimators are the closed-form univariate MLEs:
%       p_i     = fraction of positive observations,
%       mu_i    = mean(log(X_i(X_i > 0))),
%       sigma_i = sqrt(mean((log(X_i(X_i > 0)) - mu_i).^2)).


%% Basic input check

if isempty(X)
    error('Input data X cannot be empty.');
end

%% Core

p = mean(X > 0, 1);

Xpos = X;
Xpos(X == 0) = NaN;
logX = log(Xpos);

mu    = mean(logX, 1, 'omitnan');
sigma = std(logX, 1, 1, 'omitnan');

end
