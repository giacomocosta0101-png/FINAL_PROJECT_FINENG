function X_new = semi_parametric_losses(L, p, N, X)
% SEMI_PARAMETRIC_LOSSES  Simulate losses using a semi-parametric copula model.
%
% This function generates simulated uniform margins with a specific 
% dependence structure, then maps them back to the loss domain by 
% inverting the semi-parametric empirical CDF via linear interpolation.
%
% INPUT
%   L : (matrix) dependence structure parameter (e.g., Cholesky factor of copula)
%   p : (1 x d) probabilities of positive observation P(X > 0)
%   N : (scalar) number of Monte Carlo simulations to generate
%   X : (M x d) original calibration data matrix
%
% OUTPUT
%   X_new : (N x d) matrix of simulated losses

    U_sim = semi_parametric_sim(L, p, N);
    X_new = zeros(N,size(X,2));
    
    cdf = cumulative_cdf_semi_parametric_vec(p, X);
    
    for i = 1:size(X,2)
        
        X_col = unique(X(:,i));
        U = cdf{i}(X_col);
        X_new(:,i) = interp1(U, X_col, U_sim(:,i), 'linear','extrap');
    end
   
end