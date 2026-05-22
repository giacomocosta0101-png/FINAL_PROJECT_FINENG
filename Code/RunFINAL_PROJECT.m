% Project 6: Copula calibration

filename = "danishmulti.csv";
addpath('utilities','ex_1');
data = readDataset(filename);

%%just a check

buildings_jump = data.Building>0 & data.Contents == 0 & data.Profits == 0;
contents_jump = data.Building == 0 & data.Contents > 0 & data.Profits == 0;
contents_buildings_jump = data.Building > 0 & data.Contents >0 & data.Profits == 0;
p_i = mean(contents_buildings_jump>0);
p_b = mean(buildings_jump);
p_c = mean(contents_jump);

%% Marginals parameters calibration

[p, mu, sigma] = marginal_parameter_calibration(data);

%%
rho = [1 2 3 4 5 6]';
x = [1 2 3 4];
R = zeros(length(x));
mask = triu(true(length(x)),1)

%%
X = [0 0 0 1;
    4 0 4 0;
    7 4 0 45];

mu = [3 4 5 9];

sigma = [1 3 5 0.3];
p = [0.1 0.3 0.6 0.9];
rho = [0.3 0.2 0.1 0.5 0.4 0.35];
x = X(1,:);

log_likelihood = loglikelihood_Cbernoulli(X,mu,sigma,p,rho);
%FUNGE

%%
%%
X = [0 0 1;
    4 0 4;
    7 4 0];

mu = [3 4 5];

sigma = [1 3 5];
p = [0.1 0.3 0.6];
rho = [0.3 0.2 -0.1];
x = X(1,:);

log_likelihood = loglikelihood_Cbernoulli(X,mu,sigma,p,rho);
%FUNGE
%%

[p, mu, sigma, var_names] = marginal_parameter_calibration(data);

A = [];
b = [];
Aeq = [];
beq = []; 

lb = [-1,-1,-1];
ub = [1,1,1];  

%%
X = [0 0 1;
    4 0 4;
    7 4 0];

%%

building = data.Building(:);
Contents = data.Contents(:);
profits = data.Profits(:);
X = [building Contents profits];

log_likelihood = @(R_chol) -loglikelihood_Cbernoulli(X,mu,sigma,p,R_chol);
%%
n = 3;

% (1) Triangolarità: forziamo a zero gli elementi sopra la diagonale
% R_chol(:) è column-major: gli indici (1,2),(1,3),(2,3) corrispondono a 4,7,8
Aeq = zeros(3, n^2);
Aeq(1,4) = 1;   % R_chol(1,2) = 0
Aeq(2,7) = 1;   % R_chol(1,3) = 0
Aeq(3,8) = 1;   % R_chol(2,3) = 0
beq = zeros(3,1);

% (2) Bound: diagonale > 0, off-diagonal in [-1,1]
lb = -ones(n);    ub = ones(n);
for i = 1:n, lb(i,i) = 1e-6; end

% (3) Norma unitaria delle righe -> ceq
nonlcon = @(R) deal([], sum(R.^2, 2) - 1);

% Punto iniziale ammissibile: l'identità soddisfa tutto
R_chol = fmincon(log_likelihood, eye(n), [],[], Aeq, beq, lb, ub, nonlcon, opts);
%%


opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'EnableFeasibilityMode', false, ...
    'SubproblemAlgorithm','cg', ...
    'Display','iter');

R_chol = fmincon(log_likelihood,eye(3),A,b,Aeq,beq,lb,ub,[],opts); 

%%
filename = "danishmulti.csv";
addpath('utilities','ex_1');
data = readDataset(filename);

building = data.Building(:);
Contents = data.Contents(:);
profits = data.Profits(:);

X = [building Contents profits];

j = 1;
for i = 1:size(X,1)
    if X(i,1)>0 && X(i,3)>0 && X(i,2)== 0
        X_new(j,:) = X(i,[1 3]);
        j = j+1;
    end
end

%%
cdf = marginal_cdf([1 1],mu([1 3]),sigma([1 3]));
Z = cdf(X_new);
Z
%%
rng default  % For reproducibility
Rho= copulafit('Gaussian',Z)

rho_1_3 = Rho(1,2);

sigma = Rho;
L = chol(Rho,'lower');
Rho_calibrated = zeros(1e3,1);
for i = 1:1e3

    Z = randn(5,2);
    Y = L*Z';

    Rho_hat = copulafit("Gaussian",Y);
    Rho_calibrated = 
%%
rho = corr(X_new, X_new, 'Type', 'Spearman')

R = 2*sin((pi/6).*rho)


%%
filename = "danishmulti.csv";
addpath('utilities','ex_1');
data = readDataset(filename);

building = data.Building(:);
contents = data.Contents(:);
profits  = data.Profits(:);

X = [building contents profits];

% Active set {1,3}: Building > 0, Profits > 0, Contents = 0
idx = X(:,1) > 0 & X(:,3) > 0 & X(:,2) == 0;

X_13 = X(idx, [1 3]);

rho_s = corr(X_13, 'Type', 'Spearman');

R = 2*sin((pi/6).*rho_s);

rho_13_zero_mixed = R(1,2)

%%
cdf = marginal_cdf(p,mu,sigma);

z = norminv(cdf(X))




