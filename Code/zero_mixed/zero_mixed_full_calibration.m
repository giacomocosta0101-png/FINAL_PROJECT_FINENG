function calibrated_params = zero_mixed_full_calibration(X)
%
% ZERO_MIXED_FULL_CALIBRATION
%
% Calibrates the zero-mixed model parameters needed for bootstrap CIs.
%
% INPUT:
%   X : N x 3 matrix
%
% OUTPUT:
%   calibrated_params : struct with fields
%
%   .prob      : 1x8 active-set probabilities
%   .rho_cell  : 1x8 cell array
%                rho_cell{5} = rho_12 for Building-Contents
%                rho_cell{6} = rho_13 for Building-Profits
%                rho_cell{7} = rho_23 for Contents-Profits
%                rho_cell{8} = 3x3 trivariate correlation matrix
%
%   .rho_vec   : 1x6 vector
%                [rho12_2, rho13_2, rho23_2, rho12_3, rho13_3, rho23_3]
%

%% Basic setup

N = size(X, 1);

B = X(:,1);
C = X(:,2);
P = X(:,3);

active_sets = {
    [], ...
    1, ...
    2, ...
    3, ...
    [1 2], ...
    [1 3], ...
    [2 3], ...
    [1 2 3]
    };

%% Hard-coded masks

masks = cell(1, 8);

masks{1} = (B == 0) & (C == 0) & (P == 0);

masks{2} = (B > 0) & (C == 0) & (P == 0);
masks{3} = (B == 0) & (C > 0) & (P == 0);
masks{4} = (B == 0) & (C == 0) & (P > 0);

masks{5} = (B > 0) & (C > 0) & (P == 0);
masks{6} = (B > 0) & (C == 0) & (P > 0);
masks{7} = (B == 0) & (C > 0) & (P > 0);

masks{8} = (B > 0) & (C > 0) & (P > 0);

%% Initialize output

calibrated_params = struct();

calibrated_params.prob = NaN(1, 8);
calibrated_params.rho_cell = cell(1, 8);

% [rho12_biv, rho13_biv, rho23_biv, rho12_tri, rho13_tri, rho23_tri]
calibrated_params.rho_vec = NaN(1, 6);

%% Calibration loop

for k = 1:8

    idx_active = active_sets{k};
    mask = masks{k};

    X_case = X(mask, :);
    X_active = X_case(:, idx_active);

    n = sum(mask);

    calibrated_params.prob(k) = n / N;
    

    if numel(idx_active) >= 2

        if n >= 2
        mu    = mean(log(X_active), 1);
        sigma = sqrt(mean((log(X_active) - mu).^2, 1));   % MLE estimator
            rho = calibration_rho_zero_mixed(X_active, mu, sigma);
            %calibration_rho_zero_mixed returns a scalar in the case of
            % d=2
            calibrated_params.rho_cell{k} = rho;

            % Store bivariate correlations
            if k == 5
                calibrated_params.rho_vec(1) = rho;
            elseif k == 6
                calibrated_params.rho_vec(2) = rho;
            elseif k == 7
                calibrated_params.rho_vec(3) = rho;
            elseif k == 8
                R = rho;
                calibrated_params.rho_vec(4) = R(1,2);
                calibrated_params.rho_vec(5) = R(1,3);
                calibrated_params.rho_vec(6) = R(2,3);
            end

        else

            calibrated_params.rho_cell{k} = NaN;

        end

    end

end

end

