function log_likelihood = loglikelihood_Cbernoulli(X,mu,sigma,p,rho)

values = arrayfun(@(i) likelihood_comb_bernoulli(rho, X(i,:)', p, mu, sigma), 1:size(X,1));
log_likelihood= sum(values);
 

end
