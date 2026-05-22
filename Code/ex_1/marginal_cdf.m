function cdf = marginal_cdf(p,mu,sigma)

cdf = @(x) p.*normcdf((log(x)-mu)./sigma).*(x>0)+(1-p).*(x>=0);

end
