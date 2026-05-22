%%
filename = "danishmulti.csv";
addpath('utilities','ex_1');
data = readDataset(filename);

building = data.Building(:);
Contents = data.Contents(:);
profits = data.Profits(:);

[~, mu, sigma, ~] = marginal_parameter_calibration(data);

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
U_original = cdf(X_new);
%%
Rho = corr(X_new,X_new, 'Type', 'Pearson');


sigma_1_3 = [sigma(1) sigma(3)];
mu_1_3 = [mu(1) mu(3)];

rho_1_3 = Rho(1,2);

L = chol(Rho,'lower');
Rho_calibrated = zeros(1e3,1);
for i = 1:1e3
    

    Z = randn(5,2);
    Y = (L*Z')';

    U = normcdf(Y);

    for j = 1:2
        X_it(:,j) = exp(sigma_1_3(j).*norminv(U(:,j))+mu_1_3(j));
    end

    X_it = normcdf(X_it);


    Rho_calibr= corr(U,U, 'Type', 'Pearson');
    Rho_calibrated (i) = Rho_calibr(1,2);


end

quantile_1 = quantile(Rho_calibrated,0.025)
quantile_2 = quantile(Rho_calibrated,1-0.025)
%%

filename = "danishmulti.csv";
addpath('utilities','ex_1');

data = readDataset(filename);

building = data.Building(:);
Contents = data.Contents(:);
profits  = data.Profits(:);

[p, mu, sigma, var_names] = ...
    marginal_parameter_calibration(data);

X = [building Contents profits];

idx = X(:,1)>0 & X(:,3)>0 & X(:,2)==0;

X_new = X(idx,[1 3]);

cdf_fun = marginal_cdf([1 1], mu([1 3]), sigma([1 3]));

U_original = cdf_fun(X_new);

Rho = corr( ...
        X_new(:,1), ...
        X_new(:,2), ...
        'Type','Pearson');

sigma_1_3 = sigma([1 3]);
mu_1_3    = mu([1 3]);

L = chol(Rho,'lower');

Nmc = 1e4;
n   = size(X_new,1);

Rho_calibrated = zeros(Nmc,1);

for i = 1:Nmc

    Z = randn(n,2);

    Y = (L*Z')';

    U = normcdf(Y);

    X_it = zeros(n,2);

    for j = 1:2
        X_it(:,j) = exp( ...
            sigma_1_3(j)*norminv(U(:,j)) ...
            + mu_1_3(j));
    end

    Rho_calibrated(i) = corr( ...
        X_it(:,1), ...
        X_it(:,2), ...
        'Type','Pearson');

end

%% Intervallo empirico 95 %
quantile_1 = quantile(Rho_calibrated,0.025)
quantile_2 = quantile(Rho_calibrated,1-0.035)