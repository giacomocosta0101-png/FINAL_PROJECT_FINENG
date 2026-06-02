function VaR = var_calc(data,alpha,window_start,...
    window_end,N)
% VAR_CALC  Calculate Value at Risk (VaR) across different copula models.
%
% This function isolates a specific time window from the dataset, calibrates 
% three copula models (Zero-Mixed, Comb-Bernoulli, Semi-Parametric) using 
% that data, simulates total losses, and computes the empirical VaR for 
% the specified significance levels.
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
%   VaR          : (3 x m) matrix of VaR estimates, where rows correspond 
%                  to the 3 models and columns correspond to the m alpha levels

%% Calibration
data_new = data_split(data, window_start, ...
    window_end);
building = data_new.Building(:);
contents = data_new.Contents(:);
profits = data_new.Profits(:);

X_new = [building contents profits];

% Pass the calibration window to a wrapper:

calibrated_parameters = calibr_wrapper(X_new);

%% Simulation

% mat_sim takes the 3x1 cell of calibrated-parameter structs and returns a
% 3x1 cell (one entry per model). Each entry is an N x 1 vector of simulated
% TOTAL losses: mat_sim already sums the 3 risk components (Building,
% Contents, Profits) internally.
sim_losses = mat_sim(calibrated_parameters, N);

VaR = zeros(length(sim_losses),length(alpha));

for i = 1:size(VaR,1)
    VaR(i,:) = quantile(sim_losses{i},1-alpha);
end


end
