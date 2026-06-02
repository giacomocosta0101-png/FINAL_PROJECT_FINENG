function [R] = calibrate_model_generalized(U,p)
% CALIBRATE_MODEL  Calibrate Gaussian copula correlation under the 
% Comb-Bernoulli model via IFM.
%
% INPUT
%   U : (N x d) marginal CDF values F_i(x_i; p_i) for each observation
%   p : (1 x d) jump probabilities p_i = P(X_i > 0)
%
% OUTPUT
%   rho   : (1 x d) off-diagonal correlations
%   p_hat : (d x 1) empirical jump frequencies (sanity check vs. p)
%   R     : (d x d) full Gaussian copula correlation matrix

    % Convert U to a 'state' matrix:
    state_matrix  = U > (1 - p);         % N x d
    
    % We store 3 different variables:
    % C : combination matrix. It stores in one matrix the possible states of
    %     our dataset, without repetitions.
    % ia: vector of indexes such that C(i,:) = state_matrix(ia(i),:)
    % ic: vector, which associate each row of the state matrix to the state 
    %     it belongs

    [C,ia,ic] = unique(state_matrix,'rows');
    num_comb = size(C,1);
    
    Z = norminv(U);
    
    % As far as now, the spherical parametrization works only in 3-d, but
    % it can be extended to d dimensions:

    nll = @(h) -log_likelihood_trivariate_copula_generalized( ...
                  state_matrix, h, Z, ia, ic, num_comb);
    options = optimoptions('fminunc', ...
        'Display',       'off', ...
        'Algorithm',     'quasi-newton', ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance',       1e-10, ...
        'MaxIterations',       500);

    d = size(U,2);
    
    h_opt = fminunc(nll, zeros(1, d*(d-1)/2), options);

    R = corr_from_cholesky_param(h_opt, d);
end
