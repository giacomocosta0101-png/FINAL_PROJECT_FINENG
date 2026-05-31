% Project 6: Copula calibration

filename = "danishmulti.csv";
addpath('utilities','ex_1','ex_2','ex_4');
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

%% Zero-mixed 

tic
zero_mixed = zero_mixed_calibration(X);
toc

%%

sims = zero_mixed_sim(zero_mixed, N, 10000);
ci = zero_mixed_bootstrap(sims, zero_mixed, alpha);

zero_mixed_print_ci_table(ci)

%% Comb. Bernoulli
tic
fprintf("\nComb. Bernoulli\n");
cdf_comb_bernoulli = marginal_cdf(mu,sigma,p);
U_CB = cdf_comb_bernoulli(X);
[rho_CB, ~, R_CB] = calibrate_model(U_CB,p);
toc

%%
tic
bench_comb_bern()
toc
%%

fprintf("\n Correlation matrix:\n");
disp(R_CB);

%%
fprintf(" Bootstrap:\n");
rng default;
model2 = 'Comb-Bernoulli';
[rho_CI_CB, p_CI_CB, rho_hat_CB, pi_hat_CB] = bootstrap(rho_CB,p,mu,sigma,model2,N,B,alpha);

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
%% Semi-parametric - VEDI SOTTO VERSIONE IN UNA SOLA FUNZIONE
fprintf("\nSemi-Parametric\n");

cdf_semiparametric_1 = cumulative_cdf_semi_parametric(p(1),X(:,1));
U_1 = cdf_semiparametric_1(X(:,1));

cdf_semiparametric_2 = cumulative_cdf_semi_parametric(p(2),X(:,2));
U_2 = cdf_semiparametric_2(X(:,2));

cdf_semiparametric_3 = cumulative_cdf_semi_parametric(p(3),X(:,3));
U_3 = cdf_semiparametric_3(X(:,3));

U_SP = [U_1 U_2 U_3];
[rho_SP,~] = calibrate_model(U_SP,p);

%%
scatter(sort(X(:,1)),sort(U_1))
hold on
scatter(sort(X(:,1)),sort(U_CB(:,1)))
%%

sorted = sort(X(:,1));
S_emp = sort(U_1);
S_LN = sort(U_CB(:,1));
semilogx(sorted, S_emp, 'o', 'MarkerSize', 4, 'DisplayName','Empirica'); 
hold on; grid on;
semilogx(sorted, S_LN,  '-', 'LineWidth', 1.6, 'DisplayName','Lognormale');
%%

R_SP = squareform(rho_SP) + eye(length(rho_SP));
fprintf("\n Correlation matrix:\n");
disp(R_SP);
%%
fprintf(" Bootstrap:\n");
rng default;
model3 = 'Semi-parametric';
[rho_CI_SP, p_CI_SP]= bootstrap(rho_SP,p,mu,sigma,model3,N,B,alpha);

fprintf(" \n Confidence intervals:\n\n");
fprintf("  Rho_12: [ %.3f , %.3f ]\n", rho_CI_SP(1,1), rho_CI_SP(1,2));
fprintf("  Rho_13: [ %.3f , %.3f ]\n", rho_CI_SP(2,1), rho_CI_SP(2,2));
fprintf("  Rho_23: [ %.3f , %.3f ]\n", rho_CI_SP(3,1), rho_CI_SP(3,2));
fprintf("\n  p1: [ %.3f , %.3f ]\n", p_CI_SP(1,1), p_CI_SP(1,2));
fprintf("  p2: [ %.3f , %.3f ]\n", p_CI_SP(2,1), p_CI_SP(2,2));
fprintf("  p3: [ %.3f , %.3f ]\n", p_CI_SP(3,1), p_CI_SP(3,2));


%% SEMI-PARAMETRIC IN UNA SOLA FUNZIONE

cdf_semiparametric = cumulative_cdf_semi_parametric_vec(p,X);

U_SP = zeros(size(X));
for i = 1:size(X,2)
    U_SP(:,i) = cdf_semiparametric{i}(X(:,i));
end

[rho_SP2,~] = calibrate_model(U_SP,p);

%% Backtest

%calibrating function has as input data (full, timetable)
% and the calibrating period
start_date = datetime("01/01/1980");
end_date = datetime("31/12/1983");
N = 10000;
alpha = [0.05 0.005];
mode = 'Rolling-window';
data_new = data_split(data,start_date,datetime("31/12/1990"));

tic
[backtest_window,exceptions,VaR] = backtest(data_new,alpha,start_date,end_date,N,mode);
toc
%%
bw = backtest_window;
exc = exceptions;
%%
% default: 5 eccezioni 99% più gravi annotate
plot_backtest(bw, exc, VaR, 'Rolling-window', ...
              'ModelNames', {'Zero_mixed','CB','Semi_par'});
%%

% se vuoi più dettaglio sulle eccezioni
plot_backtest(bw, exc, VaR, 'Rolling-window', 'TopK', 10);

%%

% zero annotazioni: solo marker colorati
plot_backtest(bw, exc, VaR, 'Rolling-window', 'TopK', 0);

%%

% Kupiec POF test (unconditional coverage)
N  = size(backtest_window,1);
x  = sum(exceptions{3}(:,2));        % eccezioni 99%
p  = 0.01;
pi_hat = x/N;
LR_POF = -2*log( ((1-p)^(N-x) * p^x) / ((1-pi_hat)^(N-x) * pi_hat^x) );
pval_POF = 1 - chi2cdf(LR_POF, 1);

% Christoffersen independence test
% (devi contare le transizioni 00, 01, 10, 11 nella serie di eccezioni)
%%
res = christoffersen_test(exceptions);
