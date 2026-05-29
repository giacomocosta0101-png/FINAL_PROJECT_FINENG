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

q_u  = 0.99;
pts  = cell(1, size(X,2));
U_SP = zeros(size(X));

for i = 1:size(X,2)
    X_pos  = X(X(:,i) > 0, i);
    pts{i} = paretotails(X_pos, 0, q_u, 'ecdf');

    % U totale con massa atomica in 0
    Ui = zeros(size(X(:,i)));
    pos = X(:,i) > 0;
    Ui(~pos) = 1 - p(i);
    Ui( pos) = (1 - p(i)) + p(i) .* cdf(pts{i}, X(pos,i));
    U_SP(:,i) = Ui;
end

[rho_SP2,~] = calibrate_model(U_SP,p);

semi_par.p   = p;
semi_par.rho = rho_SP2;
semi_par.pts = pts;       % salvi gli oggetti paretotails
calibrated_parameters{3} = semi_par;


end