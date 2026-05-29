function zero_mixed = zero_mixed_first_calibration(data)
%
% This function has 2 main objective:
%   - Unpack the original dataset into the different subcases
%   - Calibrate the parameters of interest
%
% This function works only with the "Danish fire insurance" dataset turned into a timetable
% since columns names have been hardcoded and the number of cols is fixed to 3.
%
%
% INPUT:
%   dataset [timetable]
%
% OUTPUT:
%   zero_mixed  : 1x8 cell array of struct
%
%          1 -> no jump:      [0 0 0]
%          2 -> Building:     [1 0 0]
%          3 -> Contents:     [0 1 0]
%          4 -> Profits:      [0 0 1]
%          5 -> Building-Contents: [1 1 0]
%          6 -> Building-Profits:  [1 0 1]
%          7 -> Contents-Profits:  [0 1 1]
%          8 -> all three:    [1 1 1]
%
% Each struct contains:
%   .active_idx      active columns
%   .active_names    active variable names
%   .X               full 3-column data for that case
%   .X_active        only active positive columns
%   .n               number of observations in the case
%   .prob            empirical probability of the case
%   .mu              lognormal mu estimates for active marginals
%   .sigma           lognormal sigma estimates for active marginals
%   .rho             calibrated Gaussian copula rho, if dimension >= 2

%% Data and Input check

%here remeber to check input is a timetable, contains the "names" col and
%all entry are positive

X = [data.Building, data.Contents, data.Profits];

N = size(X, 1);

B = X(:,1);
C = X(:,2);
P = X(:,3);

%% Hardcoding active sets

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

% Hard-coded masks

masks = cell(1, 8);

% no jump
masks{1} = (B == 0) & (C == 0) & (P == 0);

% one jump
masks{2} = (B > 0) & (C == 0) & (P == 0);
masks{3} = (B == 0) & (C > 0) & (P == 0);
masks{4} = (B == 0) & (C == 0) & (P > 0);

% two jumps
masks{5} = (B > 0) & (C > 0) & (P == 0);
masks{6} = (B > 0) & (C == 0) & (P > 0);
masks{7} = (B == 0) & (C > 0) & (P > 0);

% three jumps
masks{8} = (B > 0) & (C > 0) & (P > 0);


%% Calibration loop

zero_mixed = cell(1, 8);

for k = 1:8

    idx_active = active_sets{k};
    mask = masks{k};

    X_case = X(mask, :);
    X_active = X_case(:, idx_active);

    out = struct();

    out.active_idx   = idx_active;

    out.X = X_active;

    out.n    = sum(mask);
    out.prob = out.n / N;

    out.mu    = [];
    out.sigma = [];
    out.rho   = [];


    %Calibrate mu, sigma for a lognormal (only in there is something to
    %compute)
    if ~isempty(idx_active) && out.n > 0

        logX = log(X_active);

        out.mu    = mean(logX, 1);
        out.sigma = sqrt(mean((logX - out.mu).^2, 1));   % MLE estimator

    end



    % Calibrate the copula corr matrix

    if numel(idx_active) >= 2

        if out.n >= 2
            out.rho = calibration_rho_zero_mixed(X_active, out.mu, out.sigma);
        else
            out.rho = eye(numel(idx_active));
        end

    end

    zero_mixed{k} = out;

end

end
