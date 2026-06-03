% Project 6: Copula calibration

extra_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(extra_dir);
filename = fullfile(code_dir, "danishmulti.csv");
addpath(code_dir, ...
    fullfile(code_dir, 'Comb_and_Semi'), ...
    fullfile(code_dir, 'zero_mixed'), ...
    fullfile(code_dir, 'Backtest'), ...
    fullfile(code_dir, 'Extra'));

data = readDataset(filename);

building = data.Building(:);
contents = data.Contents(:);
profits = data.Profits(:);

X = [building contents profits];

%% Generate a new column

%For example a "solar panel"
%If there is a fire then the claim on solar panel is more likely and bigger
% If there is no fire, there could have been an hailstorm


for i=1:size(X,1)
    U = rand(1);
    if X(i,1)==0
        if U < 0.8
            X(i,4) = 0;
        else
            noise = exp(randn(1)*0.2 + 0.1);
            X(i,4) = noise;
        end
    else
        if U < 0.3
            X(i,4) = 0;
        else
            noise = exp(randn(1)*0.1 + 0.7);
            X(i,4) = noise;
        end
    end

end

%% General parameters

alpha = 0.05;
B = 1e4;
N = size(X,1);

%% zero mixed calibration

zero_mixed = zero_mixed_calibration(X);

%% Zero-mixed bootstrap

fprintf("\nZero-mixed bootstrap with  %.0f replicas; alpha = %.3f\n", B, alpha);
fprintf("\nZero-mixed bootstrap\n");
rng(762);
ci_zero_mixed = zero_mixed_bootstrap(zero_mixed, alpha, N, B);

zero_mixed_print_ci_table(ci_zero_mixed)


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
R_CB = calibrate_model_generalized(U_CB,p);
toc

fprintf("\n Correlation matrix:\n");
disp(R_CB);
