function [R, Rss, Rtt, Rst, Rts] = corr_matrix(x, rho)
%   corr_matrix build the Gaussian copula correlation matrix R and partitions
%   it according to the active set S(x) = {i : x_i > 0} and its complement T,
%   consistently with the implementation in Baviera-Manzoni 2026
%
%   INPUTS
%       x   Observed claim vector. Entries must be non-negative
%       rho Vector of off-diagonal correlations in MATLAB squareform order:
%           [rho12 rho13 ... rho1d rho23 ... rho(d-1)d]
%
%   OUTPUTS
%       R    Full d-by-d correlation matrix
%       Rss, Rtt, Rst, Rts As defined in the paper


%% Inputs check
%force x col and rho raw
x = x(:);
rho = rho(:).';

d = numel(x);
expected_num_rho = d * (d - 1) / 2; 

if isempty(x) || isempty(rho)
    error('Inputs x and rho must not be empty');
end

if d <= 1
    error('Input x must contain at least two components');
end

if any(x < 0)
    error('Input x must be non-negative');
end

if numel(rho) ~= expected_num_rho
    error('rho and x dimentions are incompatible');
end

if any(abs(rho) > 1)
    error('Each entry of rho must lie in the interval [-1, 1]');
end


%% Core function

%define active set S(x) and T(x) its complementary
ss = find(x > 0);
tt = find(x == 0);

R = squareform(rho);
R(1:d+1:end) = 1;

Rss = R(ss, ss);
Rtt = R(tt, tt);
Rst = R(ss, tt);
Rts = R(tt, ss);

end
