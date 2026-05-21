function log_likelihood = loglikelihood_Cbernoulli(X,mu,sigma,p)

log_likelihood= @(rho)(sum(arrayfun(@(i) likelihood_comb_bernoulli(rho, X(i,:)',...
    p, mu, sigma), 1:size(X,1))));


end
