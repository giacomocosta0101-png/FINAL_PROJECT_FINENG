function [R, Rss, Rtt, Rst, Rts] = corr_matrix(x,rho)
% corr_matrix build the 4 matrix needed to implement the formula presented
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


%% Implement a basic input check
if isempty(x) || isempty(rho)
    error('Input vectors x and rho must not be empty.');
end

if length(x)<= 1
    error('Input vector x must have more than one element.');
end

if not(length(x) == 2 && length(rho) == 1) ||not(length(x) == 3 && length(rho) == 3) || not(length(x) == 4 && length(rho) == 6) || not(length(x) == 5 && length(rho) == 10)
    error('rho and x have incompatible dimentions');
end


%% Core

ss = find(x>0);
tt = find(x==0);

R = squareform(rho)+diag(ones(length(x),1));

Rss = R(ss,ss);
Rtt = R(tt,tt);
Rst = R(ss,tt);
Rts = R(tt,ss);


end


