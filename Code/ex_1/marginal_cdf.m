function cdf = marginal_cdf(mu,sigma,p)

cdf = @(x) (1-p).*(x>=0) + p .* normcdf((log(max(x,eps))-mu)./sigma) .* (x>0);

end
