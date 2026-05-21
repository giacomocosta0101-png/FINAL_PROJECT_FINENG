function log_likelihood = likelihood_comb_bernoulli(rho,x,p,mu,sigma)

x   = x(:);
p   = p(:);
mu  = mu(:);
sigma = sigma(:);

ss = find(x>0);
tt = find(x==0);

cdf_ss = marginal_cdf(p(ss),mu(ss),sigma(ss));
cdf_tt = marginal_cdf(p(tt),mu(tt),sigma(tt));

zs = norminv(cdf_ss(x(ss)));
zt = norminv(cdf_tt(x(tt)));

pdf     = marginal_pdf(mu(ss), sigma(ss));
pdf_val = pdf(x(ss));
pdf_val = pdf_val(:);

log_likelihood = log(copula_derivative(x, zs, zt, rho));
p = p(:);

if ~isempty(ss)
    log_likelihood = log_likelihood + sum( log( p(ss) .* pdf_val ) );
end

end
