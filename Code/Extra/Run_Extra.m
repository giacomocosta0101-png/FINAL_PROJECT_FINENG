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

%%
fprintf("\nComb. Bernoulli\n");
cdf_comb_bernoulli = marginal_cdf(mu,sigma,p);
U_CB = cdf_comb_bernoulli(X);
[rho_CB, ~, R_CB] = calibrate_model(U_CB,p);
toc

fprintf("\n Correlation matrix:\n");
disp(R_CB);


%% Spherical Parameterization
R_new = calibrate_model_generalized(U_CB,p);
disp(R_new)

% Try with a bigger dataset
dim = 4;
U_big = rand(100,dim); 
p_big = 0.3 + 0.3 * rand(1, dim); % between 0.3 and 0.6
for i=1:dim
    U_big(U_big(:,i) < p_big(i),i) = p_big(i);
end

R_big = calibrate_model_generalized(U_big,p_big);
disp(R_big)