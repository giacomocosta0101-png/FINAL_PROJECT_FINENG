function [R, Rss, Rtt, Rst, Rts] = corr_matrix(x,rho)
% corr_matrix build the 4 matrix needed to implement the formula presented
% in Massaria&Baviera 2026 to compute coupla density appearning in the
% log-likelihood function for rho. In particular this function works for a
% generic dimention of the x vector.
%
%   INPUTS:
%       x:      vector of inputs, x_i, raw of dataset
%       rho:    generic vector of correlations
%       since this function is mainly used for rho optimization, 
%   
%   Example: x=[1,0,3,19] and rho = [a,b,c,d,e,f]
%
%

ss = find(x>0);
tt = find(x==0);

R = squareform(rho)+diag(ones(length(x),1));

Rss = R(ss,ss);
Rtt = R(tt,tt);
Rst = R(ss,tt);
Rts = R(tt,ss);


end


