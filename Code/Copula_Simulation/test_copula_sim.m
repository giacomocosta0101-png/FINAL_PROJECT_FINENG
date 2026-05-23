
clear all; clc;

%% Test 1
% If sigma = 0, then X_i = exp(mu_i) exactly, regardless of R.

R = [1.0 0.5;
     0.5 1.0];

mu = [0.1, -0.2];
sigma = [0, 0];
N = 1000;

sim = copula_sim(R, mu, sigma, N);

expected = exp(mu);

err = max(abs(sim - expected), [], 'all');

fprintf('Test 1 error: %.3e\n', err);

assert(err < 1e-12, 'Test 1 failed: deterministic lognormal result is wrong.');

disp('Test 1 passed: sigma = 0 gives exact exp(mu).');


%%

R = [1,0.3,0.2; 0.3, 1, 0.4; 0.2, 0.4, 1]; 

mu = [0.1, -0.2, 0.5];
sigma = [1, 0.9, 0.3];
p=[0.5, 0.5, 0.5];
N = 1e8;

sim = comb_bern_sim(R,mu,sigma,p,N);


mu = mean(sim)
sigma = std(sim)

