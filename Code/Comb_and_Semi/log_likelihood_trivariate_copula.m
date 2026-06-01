function log_likelihood = log_likelihood_trivariate_copula(state_matrix,theta,...
    Z,ia,ic,num_comb)
% Log-likelihood (copula part only) for the trivariate Gaussian copula.
% The marginal contribution is omitted because it is constant in theta
% in the IFM step.
%
% INPUT
%   state_matrix : (N x 3) binary, state_matrix(t,i) = 1 iff x_i(t) > 0
%   theta        : (1 x 3) Cholesky spherical angles
%   Z            : (N x 3) Gaussian quantiles Phi^{-1}(F_i(x_i; p_i))
%   ia, ic       : useful vectors to 'orientate' in the state matrix
%   num_comb     : number of unique active sets observed
% OUTPUT
%   log_likelihood : scalar

L = [1 0 0;
    cos(theta(1)) sin(theta(1)) 0;
    cos(theta(2)) cos(theta(3))*sin(theta(2)) sin(theta(3))*sin(theta(2))];

R = L*L';

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

    zs      = Z(rows, ss);
    zt      = Z(rows, tt);

    if isempty(ss)

        % NB: in this group all rows share the same z_T = norminv(1-p),
        % so we evaluate mvncdf once and scale the log by the group size.

        log_dens_i = log( mvncdf(zt(1,:), zeros(1, numel(tt)), Rtt) );
        log_dens = log_dens_i*size(zt,1);

    elseif isempty(tt)

        log_dens= log( mvnpdf(zs, zeros(1, numel(ss)), Rss) ) ...
            - log( mvnpdf(zs) );
    elseif isscalar(tt)

        log_copula = log( mvnpdf(zs, zeros(1, numel(ss)), Rss) ) ...
            - log( mvnpdf(zs) );

        % In case of a 'double' jump, the number of non-active positions is
        % one, hence we use 'normcdf' instead of mvncdf, which is less
        % optimized:

        log_cdf    = log( normcdf(zt - (Rts * (Rss \ zs.')).',...
            0, sqrt(Rtt - Rts * (Rss \ Rst))) );

        log_dens   = log_copula + log_cdf;


    else

        log_copula = log( mvnpdf(zs, zeros(1, numel(ss)), Rss) ) ...
            - log( mvnpdf(zs) );


        log_cdf    = log( mvncdf(zt - (Rts * (Rss \ zs.')).',...
            zeros(1, numel(tt)), Rtt - Rts * (Rss \ Rst)) );

        log_dens   = log_copula + log_cdf;
    end

    log_likelihood = log_likelihood + sum(log_dens);
end

end
