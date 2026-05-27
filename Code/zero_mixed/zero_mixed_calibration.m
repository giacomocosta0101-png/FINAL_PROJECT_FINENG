function output = zero_mixed_calibration(data, mode)
%
% This function works only with the "Danish fire insurance" dataset turned into a timetable
% since columns names have been hardcoded and the number of cols is fixed to 3.
%
% INPUT:
%   dataset [timetable]
%
% OUTPUT mode "full":
%   zero_mixed  : 1x8 cell array of struct
%
%   .active_idx      active columns
%   .active_names    active variable names
%   .X               full 3-column data for that case
%   .X_active        only active positive columns
%   .n               number of observations in the case
%   .prob            empirical probability of the case
%   .mu              lognormal mu estimates for active marginals
%   .sigma           lognormal sigma estimates for active marginals
%   .rho             calibrated Gaussian copula rho, if dimension >= 2
%
% OUTPUT mode "rho"
%   zero_mixed : struct with fields
%
%   .prob      : 1x8 active-set probabilities
%   .rho_cell  : 1x8 cell array
%                rho_cell{5} = 2x2 rho matrix for Building-Contents
%                rho_cell{6} = 2x2 rho matrix for Building-Profits
%                rho_cell{7} = 2x2 rho matrix for Contents-Profits
%                rho_cell{8} = 3x3 trivariate correlation matrix
%
%   .rho_vec   : 1x6 vector
%                [rho12_2, rho13_2, rho23_2, rho12_3, rho13_3, rho23_3]

%% Data and Input check

mode = lower(string(mode));
if ~isscalar(mode) || ~ismember(mode, ["full", "rho"])
    error('mode must be either "full" or "rho".');
end

if istimetable(data) || istable(data)
    X = [data.Building, data.Contents, data.Profits];
else
    X = data;
end

if size(X, 2) ~= 3
    error('Input data must have exactly 3 columns.');
end

names = {'Building', 'Contents', 'Profits'};

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

zero_mixed = cell(1,8);
prob = NaN(1, 8);
rho_vec = NaN(1, 6);

for k = 1:8
    %retrive right data for specific case
    idx_active = active_sets{k};
    mask = masks{k};

    %get the right matrix for case k
    X_case = X(mask, :);
    X_active = X_case(:, idx_active);
    n = sum(mask);
    prob(k) = n / N;

    %compute first 2 moments of marginals
    %note that if X_active is empty, returns NaN
    mu = [];
    sigma = [];
    if ~isempty(idx_active) && n > 0
        logX = log(X_active);
        mu = mean(logX, 1);
        sigma = sqrt(mean((logX - mu).^2, 1));   % MLE estimator
    end

    %compute the gaussian copula correlation matrix
    if numel(idx_active) >= 2 && n > numel(idx_active) && all(isfinite(sigma)) && all(sigma > 0)
        R = calibration_rho_zero_mixed(X_active, mu, sigma);
    elseif numel(idx_active) >= 2 && mode == "full"
        % The simulator still needs a valid correlation matrix even when
        % the bootstrap sample has too few observations to estimate one.
        R = eye(numel(idx_active));
    elseif numel(idx_active) >= 2
        R = NaN(numel(idx_active));
    else
        R = NaN;
    end

    if mode == "full"

        out = struct();

        out.active_idx   = idx_active;
        out.active_names = names(idx_active);
        out.X = X_case;
        out.X_active = X_active;
        out.n    = n;
        out.prob = prob(k);
        out.mu = mu ;
        out.sigma = sigma;
        out.rho = R;
        out.R = R;

        zero_mixed{k} = out;
    end

    if mode == "rho"
        zero_mixed{k} = R;

        if k == 5 && isequal(size(R), [2 2])
            rho_vec(1) = R(1,2);
        elseif k == 6 && isequal(size(R), [2 2])
            rho_vec(2) = R(1,2);
        elseif k == 7 && isequal(size(R), [2 2])
            rho_vec(3) = R(1,2);
        elseif k == 8 && isequal(size(R), [3 3])
            rho_vec(4) = R(1,2);
            rho_vec(5) = R(1,3);
            rho_vec(6) = R(2,3);
        end
    end
end

if mode == "full"
    output = zero_mixed;
elseif mode == "rho"
    output.rho_cell = zero_mixed;
    output.prob = prob;
    output.rho_vec = rho_vec;
end

end
