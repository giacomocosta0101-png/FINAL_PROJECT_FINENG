function sim = zero_mixed_sim_fixed(zero_mixed, N, B)
%
% Simulates replicas from the calibrated zero-mixed model while keeping
% the number of observations in each active set fixed across replicas.
%
% INPUT:
%   zero_mixed : 1 x cases cell array returned by zero_mixed_calibration
%   N          : number of observations per replica
%   B          : number of replicas
%
% OUTPUT:
%   sim        : N x d x B array
%

arguments
    zero_mixed (1,:) cell {mustBeNonempty}
    N (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
    B (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
end

%% Initial setup

K = numel(zero_mixed);
%trick to get the number of cols in the general case
d = size(zero_mixed{K}.R , 2);

active_sets = get_active_sets(d);

prob = zeros(1, K);
for k = 1:K
    prob(k) = zero_mixed{k}.prob;
end

prob = prob / sum(prob);
counts = fixed_case_counts(prob, N);

sim = zeros(N, d, B);

%% Simulation loop

for b = 1:B

    X_b = zeros(N, d);
    i = counts(1);

    for k = 2:K

        n_k = counts(k);
        if n_k == 0
            continue
        end

        rows = (i + 1):(i + n_k);
        i = i + n_k;

        idx_active = active_sets{k};
        s = numel(idx_active);

        mu = zero_mixed{k}.mu(:)';
        sigma = zero_mixed{k}.sigma(:)';

        if s == 1

            X_k = exp(mu + sigma .* randn(n_k, 1));

        else

            R = zero_mixed{k}.R;
            X_k = copula_sim(R, mu, sigma, n_k);

        end

        X_b(rows, idx_active) = X_k;
    end

    sim(:,:,b) = X_b(randperm(N), :);
end

end


function counts = fixed_case_counts(prob, N)
%
% Converts expected active-set frequencies into integer counts that sum
% exactly to N.
%
% INPUT:
%   prob   : vector of active-set probabilities
%   N      : sample size
%
% OUTPUT:
%   counts : integer vector with sum(counts) = N

arguments
    prob (1,:) double {mustBeNonempty, mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(prob, 0)}
    N (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
end

expected = N * prob;
counts = floor(expected);

missing = N - sum(counts);
frac = expected - counts;

if missing > 0
    [~, order] = sort(frac, 'descend');
    counts(order(1:missing)) = counts(order(1:missing)) + 1;
end

end
