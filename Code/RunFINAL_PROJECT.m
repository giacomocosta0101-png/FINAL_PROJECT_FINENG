% Project 6: Copula calibration

filename = "danishmulti.csv";
addpath('utilities','Comb_and_Semi','zero_mixed','Backtest','Extra');

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

fprintf("\nMarginal cdf parameters:\n");
for i = 1:length(mu)
    fprintf("\n mu_%d    = %.2f", i, mu(i));
    fprintf("\n sigma_%d = %.2f\n", i, sigma(i));
end
fprintf("\nProbabilities:");
fprintf("\n p1 = %.2f", p(1));
fprintf("\n p2 = %.2f", p(2));
fprintf("\n p3 = %.2f\n", p(3));

%% Zero-mixed calibration

zero_mixed = zero_mixed_calibration(X);

%% Zero-mixed bootstrap

fprintf("\nZero-mixed bootstrap with  %.0f replicas; alpha = %.3f\n", B, alpha);
fprintf("\nZero-mixed bootstrap\n");
rng(762);
ci_zero_mixed = zero_mixed_bootstrap(zero_mixed, alpha, N, B);

zero_mixed_print_ci_table(ci_zero_mixed)

%% Zero-mixed bootstrap fixed active-set counts

fprintf("\nZero-mixed bootstrap with fixed active-set counts\n");
rng(762);
ci_zero_mixed_fixed = zero_mixed_bootstrap_fixed(zero_mixed, alpha, N, B);


zero_mixed_print_ci_table(ci_zero_mixed_fixed)

%% Comb. Bernoulli
tic
fprintf("\nComb. Bernoulli\n");
cdf_comb_bernoulli = marginal_cdf(mu,sigma,p);
U_CB = cdf_comb_bernoulli(X);
[rho_CB, ~, R_CB] = calibrate_model(U_CB,p);
toc

fprintf("\n Correlation matrix:\n");
disp(R_CB);

%%
fprintf(" Bootstrap:\n");
rng(762);
model2 = 'Comb-Bernoulli';
[rho_CI_CB, p_CI_CB, rho_hat_CB, pi_hat_CB] = bootstrap(rho_CB,p,mu,sigma,model2,N,100,alpha);

fprintf(" \n Confidence intervals:\n\n");
fprintf("  Rho_12: [ %.3f , %.3f ]\n", rho_CI_CB(1,1), rho_CI_CB(1,2));
fprintf("  Rho_13: [ %.3f , %.3f ]\n", rho_CI_CB(2,1), rho_CI_CB(2,2));
fprintf("  Rho_23: [ %.3f , %.3f ]\n", rho_CI_CB(3,1), rho_CI_CB(3,2));
fprintf("\n  p1: [ %.3f , %.3f ]\n", p_CI_CB(1,1), p_CI_CB(1,2));
fprintf("  p2: [ %.3f , %.3f ]\n", p_CI_CB(2,1), p_CI_CB(2,2));
fprintf("  p3: [ %.3f , %.3f ]\n", p_CI_CB(3,1), p_CI_CB(3,2));
%%

plot_bootstrap_rho(rho_hat_CB, rho_CB, alpha);          % 3 pannelli pairwise
% plot_bootstrap_rho_3d(rho_hat_CB, rho_CB, alpha);     % opzionale
%% Semi-parametric 
fprintf("\nSemi-Parametric\n");

cdf_semiparametric = cumulative_cdf_semi_parametric_vec(p,X);
U_SP = zeros(size(X));

for i = 1:size(X,2)
    U_SP(:,i) = cdf_semiparametric{i}(X(:,i));
end

[rho_SP,~] = calibrate_model(U_SP,p);
%%

R_SP = squareform(rho_SP) + eye(length(rho_SP));
fprintf("\n Correlation matrix:\n");
disp(R_SP);
%%
fprintf(" Bootstrap:\n");
rng(762);
model3 = 'Semi-parametric';

[rho_CI_SP, p_CI_SP]= bootstrap(rho_SP,p,mu,sigma,model3,N,1000,alpha);

fprintf(" \n Confidence intervals:\n\n");
fprintf("  Rho_12: [ %.3f , %.3f ]\n", rho_CI_SP(1,1), rho_CI_SP(1,2));
fprintf("  Rho_13: [ %.3f , %.3f ]\n", rho_CI_SP(2,1), rho_CI_SP(2,2));
fprintf("  Rho_23: [ %.3f , %.3f ]\n", rho_CI_SP(3,1), rho_CI_SP(3,2));
fprintf("\n  p1: [ %.3f , %.3f ]\n", p_CI_SP(1,1), p_CI_SP(1,2));
fprintf("  p2: [ %.3f , %.3f ]\n", p_CI_SP(2,1), p_CI_SP(2,2));
fprintf("  p3: [ %.3f , %.3f ]\n", p_CI_SP(3,1), p_CI_SP(3,2));

%% 4.BACKTEST
%% General parameters

training_window_start_date= datetime("01/01/1980");
training_window_end_date = datetime("31/12/1983");

N = size(X,1);
alpha = [0.05 0.01];

%% a. Static calibration

mode = 'Fixed';
[backtest_window,exc_static_calibration,VaR_static_calibration] = backtest(data,alpha,...
    training_window_start_date,training_window_end_date,N,mode);

plot_backtest(backtest_window, exc_static_calibration,...
    VaR_static_calibration, 'Static calibration', ...
              'ModelNames', {'Zero_mixed','CB','Semi_par'});
%% b. Rolling-window calibration

mode = 'Rolling-window';
[~,exc_rolling_window,VaR_rolling_window] = backtest(data,alpha,...
    training_window_start_date,training_window_end_date,N,mode);

plot_backtest(backtest_window, exc_rolling_window,...
    VaR_rolling_window, 'Rolling-window', ...
              'ModelNames', {'Zero_mixed','CB','Semi_par'});
%% Chrisoffersen test:

res_static = christoffersen_test(exc_static_calibration);
res_rolling = christoffersen_test(exc_rolling_window);


%% Extra 1

% Spherical Parameterization
R_new = calibrate_model_generalized(U_CB,p);
disp(R_new)

% Try with a bigger dataset
dim = 5;
U_big = rand(50,dim); 
p_big = 0.3 + 0.3 * rand(1, dim); % between 0.3 and 0.6
for i=1:dim
    U_big(U_big(:,i) < p_big(i),i) = p_big(i);
end

R_big = calibrate_model_generalized(U_big,p_big);
disp(R_big)

%% Extra 2
% Extreme value theory

cdf_semiparam_2 = cumulative_cdf_semi_parametric_pareto(p,X);
U_SP_2 = zeros(size(X));

for i = 1:size(X,2)
    U_SP_2(:,i) = cdf_semiparam_2{i}(X(:,i));
end

[rho_SP_2,~] = calibrate_model(U_SP,p);

rho_SP_2
rho_SP