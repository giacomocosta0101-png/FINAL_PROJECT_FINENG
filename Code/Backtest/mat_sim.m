function obj = mat_sim(calibrated_parameters, N) 

% prende N e cell 3x1 con tutti i parametri utili alla simulazione
%output -> cell con 3 matrici Nx3

obj = cell(3, 1); % Initialize the output cell array


X_zero_mixed = zero_mixed_sim(calibrated_parameters{1}, N,1);

mu = calibrated_parameters{2}.mu;
sigma =calibrated_parameters{2}.sigma;
p = calibrated_parameters{2}.p;
rho = calibrated_parameters{2}.rho;
X_comb_ber = comb_bern_sim(rho, mu, sigma, p, N);



p = calibrated_parameters{3}.p;
rho = calibrated_parameters{3}.rho;
X = calibrated_parameters{3}.X;
X_semi_parametric = semi_parametric_losses(rho,p, N,X);


obj{1} = sum(X_zero_mixed,2);
obj{2} = sum(X_comb_ber,2);
obj{3} = sum(X_semi_parametric,2);

end