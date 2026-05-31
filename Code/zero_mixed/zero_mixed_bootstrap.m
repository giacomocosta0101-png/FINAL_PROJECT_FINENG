function ci = zero_mixed_bootstrap(zero_mixed, alpha, N, B)
%
% Simulates and recalibrates one bootstrap replica at a time, then computes
% Bonferroni-adjusted confidence intervals for probabilities and rho
% parameters.
%
% INPUT:
%   zero_mixed : output of zero_mixed_calibration
%   alpha      : scalar or vector of significance levels
%   N          : number of observations per bootstrap replica
%   B          : number of bootstrap replicas
%
% OUTPUT:
%   ci         : 1 x numel(alpha) cell array of structs
%
%   Each struct contains:
%   .alpha        : significance level
%   .rho_center   : row vector of original rho estimates
%   .prob_center  : row vector of original probability estimates
%   .rho_CI       : 2 x n_rho matrix of confidence intervals
%   .prob_CI      : 2 x n_prob matrix of confidence intervals
%   .rho_labels   : labels matching rho_center and rho_CI
%   .prob_labels  : labels matching prob_center and prob_CI

arguments
    zero_mixed (1,:) cell {mustBeNonempty}
    alpha double {mustBeNonempty, mustBeReal, mustBeFinite, mustBeVector, ...
                  mustBeGreaterThan(alpha, 0), mustBeLessThan(alpha, 1)}
    N (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
    B (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
end

%% Initial setup
alpha = alpha(:).';

%% Extract center and bootstrap matrices

[rho_center, prob_center, rho_labels, prob_labels] = zero_mixed_unpacking(zero_mixed);

n_rho = numel(rho_center);
n_prob = numel(prob_center);

rho_mat = NaN(B, n_rho);
prob_mat = NaN(B, n_prob);

parfor b = 1:B
    X_b = zero_mixed_sim(zero_mixed, N, 1);
    zero_mixed_b = zero_mixed_calibration(X_b);
    [rho_mat(b,:), prob_mat(b,:)] = zero_mixed_unpacking(zero_mixed_b);
end

%% Output construction nightmare

m = n_rho + n_prob; %Bonferroni correction
ci = cell(1, numel(alpha));

for a = 1:numel(alpha)
    alpha_correct = alpha(a) / m;

    rho_CI = NaN(2, n_rho);
    prob_CI = NaN(2, n_prob);

    for j = 1:n_rho
        x = rho_mat(:,j);
        x = x(isfinite(x));
        if isempty(x)
            continue
        end
        rho_CI(:,j) = quantile(x, [alpha_correct/2, 1-alpha_correct/2]);
    end

    for j = 1:n_prob
        x = prob_mat(:,j);
        x = x(isfinite(x));
        if isempty(x)
            continue
        end
        prob_CI(:,j) = quantile(x, [alpha_correct/2, 1-alpha_correct/2]);
    end

    out = struct();
    out.alpha = alpha(a);
    out.rho_center = rho_center;
    out.prob_center = prob_center;
    out.rho_CI = rho_CI;
    out.prob_CI = prob_CI;
    out.rho_labels = rho_labels;
    out.prob_labels = prob_labels;

    ci{a} = out;
end

end
