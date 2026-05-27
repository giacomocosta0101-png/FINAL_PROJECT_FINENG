function cdf = marginal_cdf(marginal_params)

mu = marginal_params.mu;
sigma = marginal_params.sigma;
p = marginal_params.p;

cdf = @(x) p.*normcdf((log(x)-mu)./sigma).*(x>0)+(1-p).*(x>=0);

end
