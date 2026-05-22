function derivative = copula_derivative(x,zs,zt,rho)
% copula_derivative is the implementation the formula presented
% in Massaria&Baviera 2026 to compute coupla density appearning in the
% log-likelihood function for rho. In particular this function works for a
% generic dimention of the x vector.
%
%   INPUTS:
%       x:      vector of inputs, x_i, raw of dataset (lenght x > 1)
%       rho:    generic vector of correlations; since this function is
%       mainly used for rho optimization rho is usally parametric
%   
%   Example: x=[1,0,3,19] and rho = [a,b,c,d,e,f]
%
%   OUTPUT:
%       R -> complete correlation matrix
%       Rss, Rtt, Rst, Rts -> as defined in the paper
%
%


%Force row vector
zs = zs(:).';   
zt = zt(:).';

[~,Rss, Rtt, Rst, Rts] = corr_matrix(x,rho);

ss = find(x>0);
tt = find(x==0);

%We check is SPD:
    [~, flag1] = chol(Rss);
    [~, flag2] = chol(Rtt);
    [~, flag3] = chol(Rtt - Rts*(Rss\Rst));

    if flag1 ~= 0 || flag2 ~= 0 || flag3 ~= 0
        derivative = realmin;
        return
    end

if isempty(ss)
    derivative = mvncdf(zt, zeros(size(zt)), Rtt);
elseif isempty(tt)
    derivative = mvnpdf(zs, zeros(size(zs)), Rss) / mvnpdf(zs);
else

    derivative = ( mvnpdf(zs, zeros(size(zs)), Rss) / mvnpdf(zs) ) * ...
                 mvncdf(zt - (Rts * (Rss \ zs.')).', zeros(size(zt)), Rtt - Rts*(Rss\Rst));
end


%%
%Qui bisogna mettere cor_matrix che è solo una helper a questa
