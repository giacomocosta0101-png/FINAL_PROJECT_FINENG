function [rho_vec, prob_vec, rho_labels, prob_labels] = zero_mixed_unpacking(zero_mixed)
%
% Unpacks the output of zero_mixed_calibration into flat vectors of
% probabilities and correlations.
%
% INPUT:
%   zero_mixed : 1 x cases cell array returned by zero_mixed_calibration
%
% OUTPUT:
%   rho_vec    : row vector containing the off-diagonal entries of each
%                correlation matrix, ordered by active set and then by
%                pair index within the active set
%   prob_vec   : 1 x cases vector of active-set probabilities
%   rho_labels : 1 x numel(rho_vec) cell array of labels matching rho_vec
%   prob_labels: 1 x cases cell array of labels matching prob_vec
%
% Notes:
%   - Cases with fewer than two active components do not contribute to
%     rho_vec

arguments
    zero_mixed (1,:) cell {mustBeNonempty}
end

cases = numel(zero_mixed);
d = numel(zero_mixed{end}.idx_active);

rho_vec = [];
prob_vec = zeros(1, cases);
rho_labels = {};
prob_labels = cell(1, cases);

for i = 1:cases
    idx_active = zero_mixed{i}.idx_active;

    pattern = repmat('0', 1, d);
    pattern(idx_active) = '1';

    prob_vec(i) = zero_mixed{i}.prob;
    prob_labels{i} = sprintf('p_{%s}', pattern);

    if numel(idx_active) < 2 || ~isfield(zero_mixed{i}, 'R') || isempty(zero_mixed{i}.R)
        continue
    end

    local_pairs = nchoosek(1:numel(idx_active), 2);

    for j = 1:size(local_pairs, 1)
        i1 = local_pairs(j, 1);
        i2 = local_pairs(j, 2);

        rho_vec(end + 1) = zero_mixed{i}.R(i1, i2); %#ok<AGROW>
        rho_labels{end + 1} = sprintf('rho_{%d%d}^{%s}', ... %#ok<AGROW>
            idx_active(i1), idx_active(i2), pattern);
    end
end

end
