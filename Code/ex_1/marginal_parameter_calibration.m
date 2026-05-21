function [p, mu, sigma, var_names] = marginal_parameter_calibration(X, var_names)
% Calibrate lognormal + prob mass in zero distribution parameters.
%
%   INPUT:
%       X: full data (timetable)
%       var_names: column we want to calibrate parameters of (array of
%       string)
%
%   The estimators are the closed-form univariate MLEs:
%       p_i     = fraction of positive observations,
%       mu_i    = mean(log(X_i(X_i > 0))),
%       sigma_i = sqrt(mean((log(X_i(X_i > 0)) - mu_i).^2)).


%% Basic input check

if nargin < 2
    var_names = ["Building", "Contents", "Profits"];
end

if isempty(X)
    error('Input data X cannot be empty.');
end

if ~istimetable(X)
    error('Input X must be a timetable.');
end

X = timetable2table(X);


%% Core
X_num = X{:, cellstr(var_names)};

d = size(X_num, 2);
p = mean(X_num > 0, 1);
mu = zeros(1, d);
sigma = zeros(1, d);

for i = 1:d

    positive_obs = X_num(X_num(:, i) > 0, i);
    log_positive_obs = log(positive_obs);

    mu(i) = mean(log_positive_obs);
    sigma(i) = sqrt(mean((log_positive_obs - mu(i)).^2));
end


end
