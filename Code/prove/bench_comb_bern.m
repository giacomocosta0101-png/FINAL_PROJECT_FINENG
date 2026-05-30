function bench_comb_bern()
% BENCH_COMB_BERN  Replicate Table 9 of Baviera et al. 2026 on your code.
%
% Measures the average simulation time of the Comb-Bernoulli generator for
% B replicas of length N, across dimensions d in {2, 3, 20, 100}. The
% reported time is the average over K independent experiments with random
% parameters drawn at each run, mirroring Table 9.
%
% Notes:
%   * The Gaussian copula is used (your comb_bern_sim is Gaussian). The paper
%     uses Student-t for table 9; the scaling story is identical, only the
%     constant differs.
%   * Your comb_bern_sim.m is hardcoded for d=3, so for d>3 the benchmark
%     calls a generalised local twin (simulate_CB_general at the bottom),
%     functionally identical to your simulator but accepting any d.
%   * Set K = 100 for paper-grade averages; K = 10 is enough to see the
%     scaling and runs in a few minutes on a laptop.

    %% --- USER CONFIG ---------------------------------------------------------
    d_list = [2, 3, 20, 100];
    B      = 1e3;       % replicas per experiment
    N      = 1e3;       % length of each replica (paper does not state; Danish set has 4018)
    K      = 10;        % # of experiments to average; bump to 100 for paper-grade
    seed   = 42;
    addpath('../ex_4');  % adapt if your comb_bern_sim lives elsewhere
    % -------------------------------------------------------------------------

    rng(seed);

    %% 1) Sanity check at d=3: your comb_bern_sim VS the generalised twin
    fprintf('=== Sanity check at d = 3 ===\n');
    d = 3;
    R     = gallery('randcorr', d);
    mu    = randn(1, d);
    sigma = 0.5 + rand(1, d);
    p     = 0.1 + 0.5*rand(1, d);
    rho_offdiag = [R(1,2), R(1,3), R(2,3)];

    rng(seed); X_user = comb_bern_sim(rho_offdiag, mu, sigma, p, N);
    rng(seed); X_gen  = simulate_CB_general(R, mu, sigma, p, N);

    delta = max(abs(X_user(:) - X_gen(:)));
    fprintf('  max |X_user - X_gen| = %.3e  -->  %s\n\n', ...
            delta, ternary(delta < 1e-12, 'OK (identical)', 'MISMATCH'));

    %% 2) Scaling benchmark (Table 9)
    fprintf('=== Scaling benchmark: B = %d replicas, N = %d, K = %d experiments ===\n', B, N, K);
    fprintf('| %5s | %14s | %14s | %12s |\n', 'd', 'mean time (s)', 'std time (s)', 'time / d');
    fprintf('|-------|----------------|----------------|--------------|\n');

    results = zeros(numel(d_list), 3);   % [d, mean, std]

    for di = 1:numel(d_list)
        d = d_list(di);
        times = zeros(K, 1);

        for k = 1:K
            % Random parameters for this experiment
            R     = gallery('randcorr', d);
            mu    = -1 + 2*rand(1, d);          % LN location ~ U(-1, 1)
            sigma = 0.5 + 1.5*rand(1, d);       % LN scale ~ U(0.5, 2.0)
            p     = 0.05 + 0.5*rand(1, d);      % jump prob ~ U(0.05, 0.55)

            % Pre-factor R once per experiment (chol is O(d^3) but does not scale with B)
            L = chol(R, 'lower');

            t0 = tic;
            for b = 1:B
                simulate_CB_general_fast(L, mu, sigma, p, N);
            end
            times(k) = toc(t0);
        end

        results(di,:) = [d, mean(times), std(times)];
        fprintf('| %5d | %14.2f | %14.3f | %12.3f |\n', ...
                d, mean(times), std(times), mean(times)/d);
    end

    fprintf('\n');
    fprintf('Linear scaling check: time(d=100)/time(d=20) = %.2f  (expected ~5 for linear)\n', ...
            results(end,2) / results(end-1,2));

    %% 3) Plot
    figure('Color','w');
    errorbar(results(:,1), results(:,2), results(:,3), 'o-', 'LineWidth', 1.6, ...
             'MarkerFaceColor', 'auto');
    set(gca, 'XScale','log', 'YScale','log');
    grid on;
    xlabel('Dimension d');
    ylabel('Mean time (s) over K experiments');
    title(sprintf('Comb-Bernoulli simulation time (B=%d replicas, N=%d)', B, N));

    % Reference line: perfect linear scaling anchored at d=2
    d_ref = d_list;
    t_ref = results(1,2) * d_ref / d_ref(1);
    hold on; plot(d_ref, t_ref, 'k--', 'DisplayName','linear O(d)');
    legend({'Measured', 'O(d) reference'}, 'Location','northwest');
end


% ===== Local helpers ======================================================
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function sim = simulate_CB_general(R, mu, sigma, p, N)
% Generalised d-dimensional Comb-Bernoulli (Algorithm 1, Sec. 2.3).
% Functionally identical to comb_bern_sim.m but accepts the full R matrix
% and any dimension d. Used for the d=3 sanity check.
    d  = size(R, 1);
    L  = chol(R, 'lower');
    Z  = (L * randn(d, N))';                  % N x d correlated normals
    U  = normcdf(Z);
    sim = zeros(N, d);
    for j = 1:d
        active = U(:, j) > (1 - p(j));
        Upos   = (U(active, j) - (1 - p(j))) / p(j);
        sim(active, j) = exp(mu(j) + sigma(j) * norminv(Upos));
    end
end

function sim = simulate_CB_general_fast(L, mu, sigma, p, N)
% Same as simulate_CB_general but takes the pre-computed Cholesky L,
% avoiding O(d^3) factorisation B times per experiment. This is what we use
% in the scaling loop.
    d  = numel(p);
    Z  = (L * randn(d, N))';
    U  = normcdf(Z);
    sim = zeros(N, d);
    for j = 1:d
        active = U(:, j) > (1 - p(j));
        Upos   = (U(active, j) - (1 - p(j))) / p(j);
        sim(active, j) = exp(mu(j) + sigma(j) * norminv(Upos));
    end
end