function marginal_params = marginal_parameter_calibration(X)
% Calibrate lognormal + prob mass in zero distribution parameters.
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

d = size(X, 2);
p = mean(X > 0, 1);
mu = zeros(1, d);
sigma = zeros(1, d);

for i = 1:d

    positive_obs = X(X(:, i) > 0, i);
    log_positive_obs = log(positive_obs);

    mu(i) = mean(log_positive_obs);
    sigma(i) = sqrt(mean((log_positive_obs - mu(i)).^2));
end

marginal_params = struct();
marginal_params.p = p;
marginal_params.mu = mu;
marginal_params.sigma=sigma;



end
