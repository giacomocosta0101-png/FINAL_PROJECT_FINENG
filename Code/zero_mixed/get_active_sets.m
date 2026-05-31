function active_sets = get_active_sets(d)
%
% Builds the list of active sets for a d-dimensional zero-mixed model.
%
% INPUT:
%   d           : model dimension
%
% OUTPUT:
%   active_sets : 1 x 2^d cell array containing all active sets, ordered
%                 by cardinality from the empty set to the full set

arguments
    d (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
end

active_sets = {[]};

for k = 1:d
    C = nchoosek(1:d, k);
    active_sets = [active_sets; num2cell(C, 2)];
end

active_sets = active_sets.';

end
