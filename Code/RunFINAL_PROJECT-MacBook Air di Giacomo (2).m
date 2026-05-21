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

%%
rho = [1 2 3 4 5 6]';
x = [1 2 3 4];
R = zeros(length(x));
mask = triu(true(length(x)),1)

%%
X = [1 0 3;
    4 5 4;
    7 4 0];
mu = [3 4 5];
sigma = [1 3 5];
p = [0.1 0.3 0.5];
rho = [0.3 0.2 -0.1];
x = X(1,:);

log_likelihood = loglikelihood_Cbernoulli(X,mu,sigma,p,rho);
%FUNGE