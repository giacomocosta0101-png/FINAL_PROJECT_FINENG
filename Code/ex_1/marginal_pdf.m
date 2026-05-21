function pdf = marginal_pdf(mu,sigma)

pdf = @(x) normpdf((log(x)-mu)./sigma)./(x.*sigma);

end
