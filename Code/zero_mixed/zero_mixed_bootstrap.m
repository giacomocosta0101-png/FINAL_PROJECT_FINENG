
clc; clear all;

addpath("ex_1")
addpath("zero_mixed")
addpath("Copula_Simulation")

filename = "danishmulti.csv";
addpath('utilities','ex_1');
data = readDataset(filename);

%% let's have fun
 % this is what perfection looks like

zero_mixed = zero_mixed_first_calibration(data);

%now the calibration is perfect, I need to simulate the whole dataset


%simulate the whole dataset

B=1e3;
N=height(data);
zero_sim = zero_mixed_sim(zero_mixed, B, N);


%Ok now we can calibrate again all parameters


for b = 1:B

    X_b = zero_sim(:,:,b);

    params_b = zero_mixed_full_calibration(X_b);

    boot_prob(b,:) = params_b.prob;
    boot_rho(b,:)  = params_b.rho_vec;

end


%% Printing

alpha = 0.05;

% Centers from original calibration

prob_center = zeros(1, 8);

for k = 1:8
    prob_center(k) = zero_mixed{k}.prob;
end

% rho order:
% [rho12_2, rho13_2, rho23_2, rho12_3, rho13_3, rho23_3]

rho_center = NaN(1, 6);

rho_center(1) = zero_mixed{5}.rho;   % Building-Contents bivariate
rho_center(2) = zero_mixed{6}.rho;   % Building-Profits bivariate
rho_center(3) = zero_mixed{7}.rho;   % Contents-Profits bivariate

R3 = zero_mixed{8}.rho;
rho_center(4) = R3(1,2);   % Building-Contents trivariate
rho_center(5) = R3(1,3);   % Building-Profits trivariate
rho_center(6) = R3(2,3);   % Contents-Profits trivariate

% Confidence intervals from bootstrap

prob_CI = NaN(2, 8);
rho_CI  = NaN(2, 6);

for j = 1:8
    x = boot_prob(:,j);
    x = x(~isnan(x));

    prob_CI(:,j) = quantile(x, [alpha/14, 1-alpha/14]);
end

for j = 1:6
    x = boot_rho(:,j);
    x = x(~isnan(x));

    rho_CI(:,j) = quantile(x, [alpha/14, 1-alpha/14]);
end


prob_center   % 1x8 point estimates
prob_CI       % 2x8 confidence intervals

rho_center    % 1x6 point estimates
rho_CI        % 2x6 confidence intervals