% Project 6: Copula calibration

filename = "danishmulti.csv";
addpath('Comb_and_Semi', 'zero_mixed', 'Backtest', 'higher_dim', 'Lomax');

data = readDataset(filename);

building = data.Building(:);
contents = data.Contents(:);
profits = data.Profits(:);

X = [building contents profits];

%% Comparison of new algorithm with d = 3

[p,mu,sigma]= marginal_parameter_calibration(X);

fprintf("\nMarginal cdf parameters:\n");
for i = 1:length(mu)
    fprintf("\n mu_%d    = %.2f", i, mu(i));
    fprintf("\n sigma_%d = %.2f\n", i, sigma(i));
end

cdf_comb_bernoulli = marginal_cdf(mu,sigma,p);
U_CB = cdf_comb_bernoulli(X);
R_CB = calibrate_model_generalized(U_CB,p);
toc

fprintf("\n Comb-Bernoulli Correlation matrix:\n");
disp(R_CB);

%% Generate a new column

%For example a "solar panel"
%If there is a fire then the claim on solar panel is more likely and bigger
% If there is no fire, there could have been an hailstorm
rng(762)
N = size(X, 1);

col_pick   = randi(3, N, 1);
linear_idx = sub2ind([N, 3], (1:N).', col_pick);
X4_base    = X(linear_idx);

% perturbazione lognormale solo dove X4_base > 0
sig_noise  = 0.30;                                  % più grande = correlazione più bassa
noise      = exp( sig_noise * randn(N, 1) );
pos        = X4_base > 0;
X4         = zeros(N, 1);
X4(pos)    = X4_base(pos) .* noise(pos);

X(:,4) = X4;

%% General parameters

alpha = 0.05;
B = 1e4;
N = size(X,1);

%% zero mixed calibration

zero_mixed = zero_mixed_calibration(X);

%% Zero-mixed bootstrap

fprintf("\nZero-mixed bootstrap with  %.0f replicas; alpha = %.3f\n", B, alpha);
fprintf("\nZero-mixed bootstrap\n");
rng(762);
ci_zero_mixed = zero_mixed_bootstrap(zero_mixed, alpha, N, B);

zero_mixed_print_ci_table(ci_zero_mixed)


%% Marginals parameters calibration

[p,mu,sigma]= marginal_parameter_calibration(X);

fprintf("\nMarginal cdf parameters:\n");
for i = 1:length(mu)
    fprintf("\n mu_%d    = %.2f", i, mu(i));
    fprintf("\n sigma_%d = %.2f\n", i, sigma(i));
end

%%
fprintf("\nComb. Bernoulli\n");
cdf_comb_bernoulli = marginal_cdf(mu,sigma,p);
U_CB = cdf_comb_bernoulli(X);
R_CB = calibrate_model_generalized(U_CB,p);
toc

fprintf("\n Correlation matrix:\n");
disp(R_CB);


rho_S    = corr(X, 'type', 'Spearman');
R_target = 2 * sin(pi * rho_S / 6)
gap = norm(R_CB - R_target, 'fro');
fprintf('||R_fit - R_target||_F = %.3f\n', gap);
fprintf('max |R_fit - R_target| per entry = %.3f\n', ...
    max(abs(R_CB(:) - R_target(:))));