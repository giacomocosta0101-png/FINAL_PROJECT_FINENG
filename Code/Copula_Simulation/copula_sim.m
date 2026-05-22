function sim = copula_sim(R, mu, sigma, N)
%
% Simulates N observations from a Gaussian copula with lognormal marginals.
%
% INPUT:
%   R      : d x d Gaussian copula correlation matrix
%   mu     : mean vector
%   sigma  : std dev vector
%   N      : number of simulations
%
% OUTPUT:
%   sim    : N x d matrix of simulated observations

%% Input checks

if ~ismatrix(R) || size(R,1) ~= size(R,2)
    error('R must be a square matrix.');
end

% Force row vectors
mu = mu(:)';
sigma = sigma(:)';

% Dimension
d = size(R, 1);



%% Core

% Simulate Gaussian copula driver
Z = mvnrnd(zeros(1,d), R, N);

% Transform to lognormal marginals
sim = exp(mu + sigma .* Z);

end



%Unoptimezed old version:

%Z = mvnrnd(zeros(1,d), R, N);
%U = normcdf(Z);
%
%sim = zeros(N,l);
%
%for i=1:l
%    sim(:,i) = exp(norminv(U(:,i)).*sigma(i) + mu(i));
%end
%