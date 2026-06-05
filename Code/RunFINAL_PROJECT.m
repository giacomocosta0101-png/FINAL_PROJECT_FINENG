% Project 6: Copula calibration

filename = "danishmulti.csv";
addpath('Comb_and_Semi', 'zero_mixed', 'Backtest', 'higher_dim', 'Lomax');

data = readDataset(filename);

building = data.Building(:);
contents = data.Contents(:);
profits = data.Profits(:);

X = [building contents profits];

%% General parameters

alpha = 0.05;
B = 1e3;
N = size(X,1);

%% Marginals parameters calibration

[p,mu,sigma]= marginal_parameter_calibration(X);

fprintf("Marginal cdf parameters:\n");
fprintf("\n mu_Buildings    = %.2f", mu(1));
fprintf("\n sigma_Buildings = %.2f\n", sigma(1));
fprintf("\n mu_Contents     = %.2f", mu(2));
fprintf("\n sigma_Contents  = %.2f\n", sigma(2));
fprintf("\n mu_Profits      = %.2f", mu(3));
fprintf("\n sigma_Profits   = %.2f\n", sigma(3));

fprintf("\nProbabilities:");
fprintf("\n p_Buildings = %.2f", p(1));
fprintf("\n p_Contents  = %.2f", p(2));
fprintf("\n p_Profits   = %.2f\n", p(3));

%% Zero-mixed calibration
fprintf("\nZero-mixed\n");
zero_mixed = zero_mixed_calibration(X);

%% Zero-mixed bootstrap

fprintf("\nZero-mixed bootstrap with  %.0f replicas; alpha = %.3f\n", B, alpha);
rng(762);
ci_zero_mixed = zero_mixed_bootstrap(zero_mixed, alpha, N, B);

zero_mixed_print_ci_table(ci_zero_mixed)

%% Zero-mixed bootstrap fixed active-set counts

fprintf("\nZero-mixed bootstrap with fixed active-set counts\n");
rng(762);
ci_zero_mixed_fixed = zero_mixed_bootstrap_fixed(zero_mixed, alpha, N, B);

zero_mixed_print_ci_table(ci_zero_mixed_fixed)

%% Comb-Bernoulli calibration

tic
fprintf("\nComb-Bernoulli\n");
cdf_comb_bernoulli = marginal_cdf(mu,sigma, p);
U_CB = cdf_comb_bernoulli(X);
[rho_CB, ~, R_CB] = calibrate_model(U_CB, p);
toc

fprintf("\n Correlation matrix:\n");
disp(R_CB);

%% Comb-Bernoulli bootstrap

fprintf(" Bootstrap:\n");
rng(762);
model2 = 'Comb-Bernoulli';
[rho_CI_CB, p_CI_CB, rho_hat_CB, pi_hat_CB] = bootstrap(rho_CB, p, mu, sigma, model2, N, B, alpha);

fprintf(" \n Confidence intervals:\n\n");
fprintf("  Rho_12: [ %.3f , %.3f ]\n", rho_CI_CB(1,1), rho_CI_CB(1,2));
fprintf("  Rho_13: [ %.3f , %.3f ]\n", rho_CI_CB(2,1), rho_CI_CB(2,2));
fprintf("  Rho_23: [ %.3f , %.3f ]\n", rho_CI_CB(3,1), rho_CI_CB(3,2));
fprintf("\n  p1: [ %.3f , %.3f ]\n", p_CI_CB(1,1), p_CI_CB(1,2));
fprintf("  p2: [ %.3f , %.3f ]\n", p_CI_CB(2,1), p_CI_CB(2,2));
fprintf("  p3: [ %.3f , %.3f ]\n", p_CI_CB(3,1), p_CI_CB(3,2));

%% DA CAPIRE

plot_bootstrap_rho(rho_hat_CB, rho_CB, alpha);          % 3 pannelli pairwise
% plot_bootstrap_rho_3d(rho_hat_CB, rho_CB, alpha);     % opzionale

%% Semi-parametric calibration

fprintf("\nSemi-Parametric\n");

cdf_semiparametric = cumulative_cdf_semi_parametric_vec(p, X);
U_SP = zeros(size(X));

for i = 1:size(X,2)
    U_SP(:,i) = cdf_semiparametric{i}(X(:,i));
end

[rho_SP,~] = calibrate_model(U_SP, p);

R_SP = squareform(rho_SP) + eye(length(rho_SP));
fprintf("\n Correlation matrix:\n");
disp(R_SP);

%% Semi.parametric bootstrap

fprintf(" Bootstrap:\n");
rng(762);
model3 = 'Semi-parametric';

[rho_CI_SP, p_CI_SP] = bootstrap(rho_SP, p, mu, sigma, model3, N, 1000, alpha);

fprintf(" \n Confidence intervals:\n\n");
fprintf("  Rho_12: [ %.3f , %.3f ]\n", rho_CI_SP(1,1), rho_CI_SP(1,2));
fprintf("  Rho_13: [ %.3f , %.3f ]\n", rho_CI_SP(2,1), rho_CI_SP(2,2));
fprintf("  Rho_23: [ %.3f , %.3f ]\n", rho_CI_SP(3,1), rho_CI_SP(3,2));
fprintf("\n  p1: [ %.3f , %.3f ]\n", p_CI_SP(1,1), p_CI_SP(1,2));
fprintf("  p2: [ %.3f , %.3f ]\n", p_CI_SP(2,1), p_CI_SP(2,2));
fprintf("  p3: [ %.3f , %.3f ]\n", p_CI_SP(3,1), p_CI_SP(3,2));

%% 4.BACKTEST
fprintf("\nBacktest\n");

%% General parameters

training_window_start_date = datetime(1980,1,1);
training_window_end_date = datetime(1983,12,31);

N = 1e6;
alpha = [0.05 0.01];

%% Benchmark: VaR via Historical Simulation

losses = data.Total;
losses_sorted = sort(losses, "descend");
n = height(data);

hs_95 = losses_sorted(round(0.05*n), :);
hs_99 = losses_sorted(round(0.01*n), :);

fprintf("95th and 99th loss quantiles from historical sim: %.3f; %.3f\n", hs_95, hs_99);

%% a. Static Calibration

fprintf("\nStatic Calibration\n");
mode = 'Fixed';

rng(762);
[backtest_window,exc_static_calibration,VaR_static_calibration] = backtest(data, alpha, ...
    training_window_start_date, training_window_end_date, N, mode);

%% DA TOGLIERE LA SEPARAZIONE
plot_backtest(backtest_window, exc_static_calibration,...
    VaR_static_calibration, mode, ...
              'ModelNames', {'Zero-mixed','CB','Semi-par'});

%% Chrisoffersen test for static calibration
fprintf("\nChrisoffersen test\n");

res_static = christoffersen_test(exc_static_calibration);

%% b. Rolling-window calibration
fprintf("\nRolling-window calibration\n");

rng(762);
mode = 'Rolling-window';

[~,exc_rolling_window,VaR_rolling_window] = backtest(data, alpha, ...
    training_window_start_date, training_window_end_date, N, mode);

%% DA TOGLIERE LA SEPARAZIONE
plot_backtest(backtest_window, exc_rolling_window,...
    VaR_rolling_window, mode, ...
              'ModelNames', {'Zero-mixed','CB','Semi-par'});

%% Chrisoffersen test for rolling-window calibration
fprintf("\nChrisoffersen test\n");

res_rolling = christoffersen_test(exc_rolling_window);


%% EXTRA: Comb-Bernoulli with Lomax marginals

%data = readDataset(filename);

training_window_start_date = datetime(1980,1,1);
training_window_end_date = datetime(1983,12,31);

N = 1e6;
alpha = [0.05 0.01];

%% a. Static Calibration
rng(732)
mode = 'Fixed';
[backtest_window,exc_lomax_static,VaR_lomax_static] = backtest_lomax(data,alpha,...
    training_window_start_date,training_window_end_date,N,mode);
%%
plot_backtest_lomax(backtest_window, exc_lomax_static,...
    VaR_lomax_static, mode);

%% b. Rolling-window calibration
rng(732)
mode = 'Rolling-window';
[~,exc_lomax_rolling,VaR_lomax_rolling] = backtest_lomax(data,alpha,...
    training_window_start_date,training_window_end_date,N,mode);
%%
plot_backtest_lomax(backtest_window, exc_lomax_rolling,...
    VaR_lomax_rolling, mode);
