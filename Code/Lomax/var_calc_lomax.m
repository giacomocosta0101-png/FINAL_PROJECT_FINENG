function VaR = var_calc_lomax(data,alpha,window_start,...
    window_end,N)
% VAR_CALC_LOMAX  Calculate Value at Risk (VaR) for the Comb-Bernoulli copula
% model with Lomax marginals.
%
% This function isolates a specific time window from the dataset, calibrates
% a single Comb-Bernoulli copula model with Lomax marginals using that data,
% simulates total losses, and computes the empirical VaR for the specified
% significance levels.
%
% INPUT
%   data         : (timetable) full dataset containing 'Building', 'Contents', 
%                  and 'Profits' variables
%   alpha        : (scalar or vector) significance level(s) for the VaR 
%                  (e.g., 0.05 for a 95% confidence level)
%   window_start : (datetime or string) start of the calibration window
%   window_end   : (datetime or string) end of the calibration window
%   N            : (scalar) number of Monte Carlo simulations to run
%
% OUTPUT
%   VaR          : (1 x m) row vector of VaR estimates, one per alpha level

arguments
    data timetable
    alpha (1,:) double
    window_start datetime
    window_end datetime
    N (1,1) double
end

%% Calibration
data_new = data_split(data, window_start, ...
    window_end);
if isempty(data_new)
    error('var_calc:emptyWindow', ...
        'The calibration window contains no observations.');
end

building = data_new.Building(:);
contents = data_new.Contents(:);
profits = data_new.Profits(:);

X = [building contents profits];


[alpha_lomax, lambda, p]= marginal_parameter_calibration_lomax(X);
cdf_comb_bernoulli_lomax= marginal_cdf_lomax(p, alpha_lomax,lambda);

U_CB = cdf_comb_bernoulli_lomax(X);
[rho_CB,~] = calibrate_model(U_CB,p);

%% Simulation
% define the cholesky factorization at the begninning to avoid doing it each
% iteration

R = squareform(rho_CB)+eye(3);
L = chol(R,'lower');

X_comb_ber = comb_bern_lomax_sim(L, alpha_lomax, lambda, p, N);

sim_losses = sum(X_comb_ber,2);

VaR = zeros(1,2);

for i = 1:length(alpha)
    VaR(i) = quantile(sim_losses,1-alpha(i));
end


end
