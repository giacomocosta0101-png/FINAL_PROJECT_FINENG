function zero_mixed_print_ci_table(ci, alpha)
%
% Prints bootstrap confidence intervals in a compact paper-style format.
%
% INPUT:
%   ci    : one struct, or a cell array of structs returned by
%           zero_mixed_bootstrap
%   alpha : optional significance level, used only if ci is a cell array
%
% OUTPUT:
%   This function does not return values. It prints the table to the
%   command window.

if ~iscell(ci) && ~isstruct(ci)
    error('ci must be a struct or a cell array of structs.');
end

if nargin >= 2 && (~isscalar(alpha) || alpha <= 0 || alpha >= 1)
    error('alpha must be a scalar strictly between 0 and 1.');
end

%% Select the right struct

if iscell(ci)
    if nargin < 2
        ci = ci{1};
    else
        idx = [];
        for i = 1:numel(ci)
            if abs(ci{i}.alpha - alpha) < 1e-12
                idx = i;
                break
            end
        end

        if isempty(idx)
            error('Requested alpha is not available.');
        end

        ci = ci{idx};
    end
end

%% Print

fprintf('\n');
fprintf('Zero-mixed bootstrap confidence intervals\n');
fprintf('alpha = %.4f\n', ci.alpha);

fprintf('\nProbabilities\n');
fprintf('%-18s %12s %12s %12s\n', 'Parameter', 'Estimate', 'Lower', 'Upper');
fprintf('%s\n', repmat('-', 1, 58));
for i = 1:numel(ci.prob_labels)
    fprintf('%-18s %12.3f %12.3f %12.3f\n', ...
        ci.prob_labels{i}, ci.prob_center(i), ci.prob_CI(1,i), ci.prob_CI(2,i));
end

fprintf('\nCorrelations\n');
fprintf('%-18s %12s %12s %12s\n', 'Parameter', 'Estimate', 'Lower', 'Upper');
fprintf('%s\n', repmat('-', 1, 58));
for i = 1:numel(ci.rho_labels)
    fprintf('%-18s %12.3f %12.3f %12.3f\n', ...
        ci.rho_labels{i}, ci.rho_center(i), ci.rho_CI(1,i), ci.rho_CI(2,i));
end

end
