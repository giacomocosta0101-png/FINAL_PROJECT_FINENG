function zero_mixed = zero_mixed_calibration(X)
%
% Calibrates the zero-mixed model on a generic N x d dataset.
%
% INPUT:
%   X          : N x d matrix with non-negative observations
%
% OUTPUT:
%   zero_mixed : 1 x 2^d cell array of structs, one for each active set
%
%   Each struct contains:
%   .n         : number of observations in the active set
%   .idx_active: indices of the active coordinates
%   .prob      : empirical probability of the active set
%   .mu        : lognormal location estimates for the active marginals
%   .sigma     : lognormal scale estimates for the active marginals
%   .R         : Gaussian copula correlation matrix for the active set,
%                when at least two coordinates are active

if ~isnumeric(X) || ~ismatrix(X)
    error('X must be a numeric matrix.');
end

if any(X < 0, 'all')
    error('X must contain only non-negative observations.');
end

%%Input check
[N,d] = size(X);

%% Get active sets and mask

active_sets = get_active_sets(d);
cases = numel(active_sets);

masks = cell(1, cases);

for s = 1:cases

    active = active_sets{s};
    inactive = setdiff(1:d, active);

    mask_active = all(X(:, active) > 0, 2);
    mask_inactive = all(X(:, inactive) == 0, 2);

    masks{s} = mask_active & mask_inactive;

end


%% Calibration loop

zero_mixed = cell(1,numel(active_sets));

for k = 1:cases

    idx_active = active_sets{k};
    mask = masks{k};

    X_case = X(mask, :);
    X_active = X_case(:, idx_active);

    out = struct();

    out.n    = sum(mask);
    out.idx_active = idx_active;
    out.prob = out.n / N;
    
    out.mu    = [];
    out.sigma = [];
    out.R     = [];

    if ~isempty(idx_active) && out.n > 0

        logX = log(X_active);

        out.mu    = mean(logX, 1);
        out.sigma = sqrt(mean((logX - out.mu).^2, 1));   % MLE estimator

    end


    % Calibrate the copula corr matrix
    if numel(idx_active) >= 2
        if out.n > 0 && all(out.sigma > 0)
            out.R = calibration_rho_zero_mixed(X_active, out.mu, out.sigma);
        else
            out.R = eye(numel(idx_active));
        end

    end
    zero_mixed{k} = out;

end



end
