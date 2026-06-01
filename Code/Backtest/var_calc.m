function VaR = var_calc(data,alpha,window_start,...
    window_end,N)

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
