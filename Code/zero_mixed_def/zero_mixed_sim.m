function sim = zero_mixed_sim(zero_mixed, N, B)
%
% Simulates replicas from the calibrated zero-mixed model.
%
% INPUT:
%   zero_mixed : 1 x cases cell array returned by zero_mixed_calibration
%   N          : number of observations per replica
%   B          : number of replicas
%
% OUTPUT:
%   sim        : N x d x B array
%                sim(:,:,b) is the b-th simulated dataset
%

if ~iscell(zero_mixed) || isempty(zero_mixed)
    error('zero_mixed must be a non-empty cell array.');
end

if ~isscalar(N) || N <= 0 || N ~= floor(N) || ...
        ~isscalar(B) || B <= 0 || B ~= floor(B)
    error('N and B must be positive integers.');
end

%%Input check and some initial work

rng(762)

K = numel(zero_mixed);
%trick to get the number of cols in the general case
d = size(zero_mixed{K}.R , 2);

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

    % Get active set 
    U = rand(N, 1);
    case_id = zeros(N, 1);

    for k = 1:K
        case_id(U <= cumprob(k) & case_id == 0) = k;
    end

    % Simulate basing on current state
    for k = 2:K   % K = 1 is no-jump, already zero

        rows = find(case_id == k);
        n_k = numel(rows);

        if n_k == 0
            continue
        end

        idx_active = zero_mixed{k}.idx_active;
        s = numel(idx_active);

        mu = zero_mixed{k}.mu(:)';
        sigma = zero_mixed{k}.sigma(:)';

        if s == 1

            % Univariate lognormal
            X_k = exp(mu + sigma .* randn(n_k, 1));

        else

            % Gaussian copula with lognormal marginals
            R = zero_mixed{k}.R;
            X_k = copula_sim(R, mu, sigma, n_k);

        end

        X_b(rows, idx_active) = X_k;

    end

    sim(:,:,b) = X_b;

end

end
