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

log_likelihood = loglikelihood_Cbernoulli(X,mu,sigma,p);
%%
lb = [-0.99,-0.99,-0.99];
ub = [0.99,0.99,0.99];
A = [];
b = [];
Aeq = [];
beq = [];

function [c,ceq] = semi_positive_def(rho)
c = -eig(squareform(rho)+eye(3));
ceq = [];
end

nonlincon = @semi_positive_def;


rho = fmincon(log_likelihood,[0.5 0.5 0.5],A,b,Aeq,beq,lb,ub,nonlincon); 