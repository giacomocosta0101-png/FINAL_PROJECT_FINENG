function case_data = zero_mixed_case_data(X)
%
% Precomputes active-set metadata for zero-mixed calibration.
%
% INPUT:
%   X         : N x d matrix with non-negative observations
%
% OUTPUT:
%   case_data : struct containing shared calibration metadata
%
%   Fields:
%   .N          : sample size
%   .d          : data dimension
%   .active_sets: 1 x 2^d cell array of active-coordinate indices
%   .cases      : number of active sets
%   .entries    : 1 x 2^d struct array, one element per active set
%
%   Each entry contains:
%   .idx_active : indices of the active coordinates
%   .rows       : row indices belonging to the active set
%   .n          : number of observations in the active set
%   .prob       : empirical probability of the active set


arguments
    X (:,:) double {mustBeNonempty, mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(X, 0)}
end

[N, d] = size(X);

active_sets = get_active_sets(d);
cases = numel(active_sets);

weights = reshape(2 .^ (0:d-1), [], 1);
signature = (X > 0) * weights;

entries = repmat(struct('idx_active', [], 'rows', [], 'n', 0, 'prob', 0), 1, cases);

for k = 1:cases
    idx_active = active_sets{k};
    code = sum(2 .^ (idx_active - 1));
    rows = find(signature == code);
    n_k = numel(rows);

    entries(k).idx_active = idx_active;
    entries(k).rows = rows;
    entries(k).n = n_k;
    entries(k).prob = n_k / N;
end

case_data = struct();
case_data.N = N;
case_data.d = d;
case_data.active_sets = active_sets;
case_data.cases = cases;
case_data.entries = entries;

end
