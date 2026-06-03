function log_likelihood = log_likelihood_trivariate_copula_generalized(state_matrix,theta,...
    Z,ia,ic,num_comb)
% LOG_LIKELIHOOD_TRIVARIATE_COPULA  Evaluates the Log-likelihood (copula 
% part only) for the trivariate Gaussian copula.
% The marginal contribution is omitted because it is constant in theta
% in the IFM step.
%
% INPUT
%   state_matrix : (N x 3) binary, state_matrix(t,i) = 1 iff x_i(t) > 0
%   theta        : (1 x 3) Cholesky spherical angles
%   Z            : (N x 3) Gaussian quantiles Phi^{-1}(F_i(x_i; p_i))
%   ia, ic       : useful vectors to 'orientate' in the state matrix
%   num_comb     : number of unique active sets observed
%
% OUTPUT
%   log_likelihood : scalar

d = size(Z,2);
[R, ~] = corr_from_cholesky_param(theta, d); 

log_likelihood = 0;

for i = 1:num_comb
    rows = (i == ic);
    s    = state_matrix(ia(i),:);
    ss   = find(s > 0);
    tt   = find(s == 0);

    Rss = R(ss, ss);
    Rtt = R(tt, tt);
    Rst = R(ss, tt);
    Rts = R(tt, ss);

    zs  = Z(rows, ss);
    zt  = Z(rows, tt);

    if isempty(ss)
        % No active component: only the Gaussian CDF on the inactive block.
        % All rows share the same z_T (= norminv(1-p)), evaluate once.
        log_dens_i = log( mvncdf(zt(1,:), zeros(1, numel(tt)), Rtt) );
        log_dens   = log_dens_i * size(zt,1);

    else
        % explicit log Gaussian copula on the active block
        L_R     = chol(Rss, 'lower');
        log_det = 2*sum(log(diag(L_R)));            % log |Rss|
        y       = L_R \ zs.';                        % k x n_rows
        quad_R  = sum(y.^2, 1).';                    % zs * inv(Rss) * zs'
        quad_I  = sum(zs.^2, 2);                     % zs * zs'

        log_copula = -0.5*log_det - 0.5*quad_R + 0.5*quad_I;

        if isempty(tt)
            log_dens = log_copula;

        elseif isscalar(tt)
            sigma2_c = Rtt - Rts*(Rss\Rst);
            mu_c     = (Rts * (Rss \ zs.')).';
            log_cdf  = log( normcdf(zt - mu_c, 0, sqrt(sigma2_c)) );
            log_dens = log_copula + log_cdf;

        else
            Sigma_c  = Rtt - Rts*(Rss\Rst);
            mu_c     = (Rts * (Rss \ zs.')).';
            log_cdf  = log( mvncdf(zt - mu_c, zeros(1, numel(tt)), Sigma_c) );
            log_dens = log_copula + log_cdf;
        end
    end

    log_likelihood = log_likelihood + sum(log_dens);
end

end
