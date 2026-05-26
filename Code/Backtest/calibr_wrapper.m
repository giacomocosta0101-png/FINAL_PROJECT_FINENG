function calibrated_parameters = calibr_wrapper(data)
%ritorna cell di struct

calibrated_parameters = cell(3,1);


%% zero mixed
zero_mixed = zero_mixed_first_calibration(data);

calibrated_parameters{1} = zero_mixed; 

%% comb ber
comb_ber = struct();
[p, mu, sigma] = marginal_parameter_calibration(data);
cdf_comb_bernoulli = marginal_cdf(p,mu,sigma);

X(:,1)= data.Building;
X(:,2) = data.Contents;
X(:,3) = data.Profits;

U_CB = cdf_comb_bernoulli(X);
[rho_CB,~] = calibrate_model(U_CB,p);

comb_ber.p = p;
comb_ber.mu= mu;
comb_ber.sigma = sigma;
comb_ber.rho = rho_CB;


calibrated_parameters{2} = comb_ber; 


%% semi parametric

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
semi_par.rho = rho_CB;
semi_par.X = X;


calibrated_parameters{3} = semi_par; 


end