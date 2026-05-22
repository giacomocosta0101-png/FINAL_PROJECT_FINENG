filename = "danishmulti.csv";
addpath('utilities','ex_1');
data = readDataset(filename);

building = data.Building(:);
Contents = data.Contents(:);
profits = data.Profits(:);

X = [building Contents profits];

[p, mu, sigma] = marginal_parameter_calibration(data);

cdf = marginal_cdf(p,mu,sigma);
U = cdf(X)
Z = norminv(U);
pdf = marginal_pdf(mu,sigma);
psi = pdf(X);

mask = X>0;

[C,ia,ic] = unique(mask,'rows');

rows = (2 == ic);
X_zeros = X(rows,:);


ss = mask(ia(2),:);

%%

num_comb = size(C,1);

%nll= @(theta)-log_likelihood_CB(mask,theta,p,Z,psi,ia,ic,num_comb);
%%
nll= @(theta)-log_likelihood_semiparametric(mask,theta,Z,ia,ic,num_comb);
%%
opts = optimoptions('fminunc', ...
    'Algorithm',                'quasi-newton', ...   % BFGS, no gradient richiesto
    'Display',                  'iter-detailed', ...
    'OptimalityTolerance',      1e-7, ...
    'StepTolerance',            1e-7, ...
    'FunctionTolerance',        1e-7, ...
    'MaxIterations',            2000, ...
    'MaxFunctionEvaluations',   2e5, ...
    'FiniteDifferenceType',     'central', ...        % +preciso di 'forward'
    'FiniteDifferenceStepSize', 1e-5, ...
    'UseParallel',              false);

tic

theta_opt = fminunc(nll, [pi/2 pi/2 pi/2], opts);

elapsed_time = toc;

fprintf('Tempo totale ottimizzazione: %.2f secondi\n', elapsed_time);

%%
theta = [pi/2; pi/2; pi/2];

profile clear
profile on

log_likelihood_semiparametric(mask, theta, Z, ia, ic, num_comb);

profile off
profile viewer
%%
Z_semi = norminv(U);
nll= @(theta)-log_likelihood_semiparametric(mask,theta,Z_semi,ia,ic,num_comb);
semi_parametric_theta_opt = fminunc(nll, [pi/2 pi/2 pi/2], opts);

semi_rho(1) = cos(semi_parametric_theta_opt(1));
semi_rho(2) = cos(semi_parametric_theta_opt(2));
semi_rho(3) = cos(semi_parametric_theta_opt(1))*cos(semi_parametric_theta_opt(2))+sin(semi_parametric_theta_opt(1))*cos(semi_parametric_theta_opt(3))*sin(semi_parametric_theta_opt(2));
%%
semi_R = squareform(semi_rho)+eye(3);
%%

rho(1) = cos(theta_opt(1));
rho(2) = cos(theta_opt(2));
rho(3) = cos(theta_opt(1))*cos(theta_opt(2))+sin(theta_opt(1))*cos(theta_opt(3))*sin(theta_opt(2));

R = squareform(rho)+eye(3);
R
R([1 0 1],[0 1 0])
theta_opt
%%
cdf_semiparametric_1 = cumulative_cdf_semi_parametric(p(1),X(:,1));
U_1 = cdf_semiparametric_1(X(:,1));
cdf_semiparametric_1(0)

%%
cdf_semiparametric_2 = cumulative_cdf_semi_parametric(p(2),X(:,2));
U_2 = cdf_semiparametric_2(X(:,2));
cdf_semiparametric_3 = cumulative_cdf_semi_parametric(p(3),X(:,3));
U_3 = cdf_semiparametric_3(X(:,3));
U = [U_1 U_2 U_3];

Z = norminv(U);

z1 = norminv(1-p(1));
z2 = norminv(1-p(2));

g = @(theta) g_rho(theta,p([1 2]),U,z1,z2,Z);

rho_hat = fzero(g,pi/3);
%%
sin(rho_hat)