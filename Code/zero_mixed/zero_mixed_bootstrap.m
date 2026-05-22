
clc; clear all;

addpath("ex_1")

filename = "danishmulti.csv";
addpath('utilities','ex_1');
data = readDataset(filename);





%% Fix consider just the trivariate case

X(:,1) = data.Building;
X(:,2) = data.Contents;
X(:,3) = data.Profits;

j=1;
for i=1:height(data)
    if (X(i,1)>0) && (X(i,2)>0) && (X(i,3)>0)
        X_new(j,:) = X(i,:);
        j=j+1;
    end
end



rho = calibration_rho_zero_mixed(X_new);

rho(1,2)
%result is fine for any subcase! this rho calculations is right
%the only question is if to recopute the mean inside or not


%%

% now let's slowly try to bootstrap

%we need to generate B matrix of N rows and 3 cols;
% each matrix is the simulation of a gaussian copula with the estimated
% parameters
%
% Then, for each realization, we evaluate the rho and take some  

B = 1000;
N = size(X_new,1);

%set initial parameters:
mu = mean(log(X_new), 1);
sigma = std(log(X_new), 1, 1);   % MLE denominator N
for i=1:B
    
    X_sim= copula_sim(rho, mu, sigma, N);
    rho_sim = calibration_rho_zero_mixed(X_sim);
    
    rho_12(i) = rho_sim(1,2);
    rho_23(i) = rho_sim(2,3);
    rho_13(i) = rho_sim(1,3);


end

quantile(rho_12,[0.025 0.975])

