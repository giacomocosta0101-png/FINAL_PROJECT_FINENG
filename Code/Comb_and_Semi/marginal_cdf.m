function cdf = marginal_cdf(mu,sigma,p)
% MARGINAL_CDF  Construct the marginal CDF function handle for a 
% zero-inflated lognormal distribution.
%
% Models a mixed distribution with a point mass at zero (probability 1-p) 
% and a lognormal distribution for strictly positive values (probability p).
%
% INPUT
%   mu    : mean vector of the underlying normal distribution
%   sigma : standard deviation vector of the underlying normal distribution
%   p     : probability vector of positive observations P(X > 0)
%
% OUTPUT
%   cdf   : (1 x d) cell array of function handles. Each handle @(x) evaluates
%           the mixed CDF at given points x

cdf = @(x) (1-p).*(x>=0) + p .* normcdf((log(max(x,eps))-mu)./sigma) .* (x>0);

end
