function calibrated_parameters = calibr_wrapper(X)
% CALIBR_WRAPPER  Main wrapper to calibrate three different copula models 
% for zero-inflated data.
%
% This function runs a full calibration for three distinct 
% approaches: 1) Zero-Mixed, 2) Comb-Bernoulli, and 3) Semi-Parametric. 
% It evaluates the marginals and extracts the copula correlation matrices.
%
% INPUT
%   X : (N x d) matrix of observations (e.g., jump sizes where 0 = no jump)
%
% OUTPUT
%   calibrated_parameters : (3 x 1) cell array of structs containing:
%       {1} zero_mixed : parameters from the zero-mixed calibration
%       {2} comb_ber   : struct with p, mu, sigma, and rho (Comb-Bernoulli)
%       {3} semi_par   : struct with p, mu, sigma, rho, and X (Semi-Parametric)

calibrated_parameters = cell(3,1);


%% Zero-Mixed


zero_mixed = zero_mixed_calibration(X);

calibrated_parameters{1} = zero_mixed; 

%% Comb-Bernoulli
comb_ber = struct();
[p, mu, sigma] = marginal_parameter_calibration(X);
cdf_comb_bernoulli = marginal_cdf(mu,sigma,p);

U_CB = cdf_comb_bernoulli(X);
[rho_CB,~] = calibrate_model(U_CB,p);

comb_ber.p = p;
comb_ber.mu= mu;
comb_ber.sigma = sigma;
comb_ber.rho = rho_CB;


calibrated_parameters{2} = comb_ber; 


%% Semi-Parametric

semi_par = struct();

cdf_semiparametric = cumulative_cdf_semi_parametric_vec(p,X);

U_SP = zeros(size(X));
for i = 1:size(X,2)
    U_SP(:,i) = cdf_semiparametric{i}(X(:,i));
end

[rho_SP2,~] = calibrate_model(U_SP,p);

semi_par.p = p;
semi_par.mu= mu;
semi_par.sigma = sigma;
semi_par.rho = rho_SP2;
semi_par.X = X;


calibrated_parameters{3} = semi_par; 


end
