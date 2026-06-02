function obj = mat_sim(calibrated_parameters, N)
% MAT_SIM  Simulate total losses across three models.
%
% This function takes the calibrated parameters for the Zero-Mixed, 
% Comb-Bernoulli, and Semi-Parametric models, runs N  simulations for each,
% and aggregates the 3-dimensional risk components into a single total loss
% vector per model.
%
% INPUT
%   calibrated_parameters : (3 x 1) cell array of parameter structs generated
%                           by the calibr_wrapper function.
%   N                     : (scalar) number of Monte Carlo simulations to run.
%
% OUTPUT
%   obj                   : (3 x 1) cell array. Each entry contains an 
%                           (N x 1) vector of simulated total aggregate losses.


obj = cell(3, 1); % Initialize the output cell array

%% Zero-Mixed
X_zero_mixed = zero_mixed_sim(calibrated_parameters{1}, N,1);

%% Comb-Bernoulli
mu = calibrated_parameters{2}.mu;
sigma =calibrated_parameters{2}.sigma;
p = calibrated_parameters{2}.p;
rho = calibrated_parameters{2}.rho;

% define the cholesky factorization at the begninning to avoid doing it each
% iteration
R = squareform(rho)+eye(3);
L = chol(R,'lower');

X_comb_ber = comb_bern_sim(L, mu, sigma, p, N);

%% Semi-Parametric
p = calibrated_parameters{3}.p;
rho = calibrated_parameters{3}.rho;
X = calibrated_parameters{3}.X;

% define the cholesky factorization at the begninning to avoid doing it each
% iteration
R = squareform(rho)+eye(3);
L = chol(R,'lower');

X_semi_parametric = semi_parametric_losses(L,p, N,X);


obj{1} = sum(X_zero_mixed,2);
obj{2} = sum(X_comb_ber,2);
obj{3} = sum(X_semi_parametric,2);

end