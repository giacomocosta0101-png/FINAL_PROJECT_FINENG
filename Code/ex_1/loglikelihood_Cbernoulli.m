function log_likelihood = loglikelihood_Cbernoulli(X,mu,sigma,p,R_chol)

R = R_chol*R_chol';

rho = [R(1,2) R(1,3) R(2,3)];

log_likelihood= sum(arrayfun(@(i) likelihood_comb_bernoulli(rho, X(i,:)',...
    p, mu, sigma), 1:size(X,1)));


end
