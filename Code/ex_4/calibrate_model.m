function [rho,p_hat,R] = calibrate_model(U,p)
% Calibrate Gaussian copula correlation under the Comb-Bernoulli model
% via IFM (eq. 9, Sec. 2.4.1 of Baviera et al. 2026). Hardcoded for d=3.
%
% INPUT
%   U : (N x 3) marginal CDF values F_i(x_i; p_i) for each observation
%   p : (1 x 3) jump probabilities p_i = P(X_i > 0)
%
% OUTPUT
%   rho   : (1 x 3) off-diagonal correlations [rho_12, rho_13, rho_23]
%   p_hat : (3 x 1) empirical jump frequencies (sanity check vs. p)
%   R     : (3 x 3) full Gaussian copula correlation matrix

arguments
    U (:,3) double {mustBeReal, mustBeFinite, ...
                    mustBeGreaterThanOrEqual(U, 0), ...
                    mustBeLessThanOrEqual(U, 1)}
    p (1,3) double {mustBeReal, mustBeFinite, ...
                    mustBeGreaterThan(p, 0), ...
                    mustBeLessThanOrEqual(p, 1)}
end

    % Convert U to a 'state' matrix:
    state_matrix  = U > (1 - p);         % N x d
    p_hat = mean(state_matrix, 1).';     % d x 1
    
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

    nll= @(theta)-log_likelihood_trivariate_copula(state_matrix,theta,...
                                                        Z,ia,ic,num_comb);
    options = optimoptions('fminunc', ...
        'Display',       'off', ...
        'Algorithm',     'quasi-newton', ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance',       1e-10, ...
        'MaxIterations',       500);
    
    [theta_opt, ~, exitflag, output] = fminunc(nll, [pi/2 pi/2 pi/2], options);

    if exitflag <= 0
        warning('calibrate_model:NoConvergence', ...
            'fminunc did not converge (exitflag = %d, msg: %s). ', ...
            'Returning best iterate; results may be unreliable.', ...
            exitflag, output.message);
    end

    rho = zeros(size(U,2),1); 
    rho(1) = cos(theta_opt(1));
    rho(2) = cos(theta_opt(2));
    rho(3) = cos(theta_opt(1))*cos(theta_opt(2))+sin(theta_opt(1))*cos(theta_opt(3))*sin(theta_opt(2));
    
    R = squareform(rho) + eye(size(U,2));
end
