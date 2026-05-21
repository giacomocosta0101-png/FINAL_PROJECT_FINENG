function rho_hat = numerical_mle(loglikelihood, rho)



%define object to minimize and initial guess
x0 = (0.5).*ones(length(rho));
lb = zeros(length(rho));
ub = ones(length(rho));


%optimize 
opt_params = fmincon(-obj, x0, [], [], [], [], lb, ub, @nonlcon);