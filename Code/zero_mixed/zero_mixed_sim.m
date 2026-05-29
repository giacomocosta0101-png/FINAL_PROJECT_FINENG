function sim = zero_mixed_sim(zero_mixed, B, N)
%
% Simulates B replicas of an N x 3 dataset from the calibrated
% zero-mixed model.
%
% INPUT:
%   zero_mixed : 1x8 cell array returned by zero_mixed_first_calibration
%   B          : number of bootstrap replicas
%   N          : number of observations per replica
%
% OUTPUT:
%   sim        : N x 3 x B array
%                sim(:,:,b) is the b-th simulated dataset
%
% Active-set convention:
%   1 -> no jump:          []
%   2 -> Building:         [1]
%   3 -> Contents:         [2]
%   4 -> Profits:          [3]
%   5 -> Building-Contents:[1 2]
%   6 -> Building-Profits: [1 3]
%   7 -> Contents-Profits: [2 3]
%   8 -> all three:        [1 2 3]

%% Extract active-set probabilities


K = 8;
d = 3;

prob = zeros(1, K);
for k = 1:K
    prob(k) = zero_mixed{k}.prob;
end

cumprob = cumsum(prob);

if (cumprob(end) > (1+1e-4)) || (cumprob(end) < (1-1e-4))
    error("Something off with cases probability")
end


sim = zeros(N, d, B);

%% Bootstrap simulation loop

for b = 1:B

    X_b = zeros(N, d);

    % Step 1: get active set 

    U = rand(N, 1);
    case_id = zeros(N, 1);

    for k = 1:K
        case_id(U <= cumprob(k) & case_id == 0) = k;
    end

    % Step 2: simulate basing on current state
    for k = 2:K   % k = 1 is no-jump, already zero

        rows = find(case_id == k);
        n_k = numel(rows);

        if n_k == 0
            continue
        end

        idx_active = zero_mixed{k}.active_idx;
        s = numel(idx_active);

        mu = zero_mixed{k}.mu(:)';
        sigma = zero_mixed{k}.sigma(:)';

        if s == 1

            % Univariate lognormal
            X_k = exp(mu + sigma .* randn(n_k, 1));

        else

            % Gaussian copula with lognormal marginals
            R = zero_mixed{k}.rho;
            X_k = copula_sim(R, mu, sigma, n_k);

        end

        X_b(rows, idx_active) = X_k;

    end

    sim(:,:,b) = X_b;

end

end


