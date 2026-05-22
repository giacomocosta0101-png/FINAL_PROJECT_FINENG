function log_likelihood = log_likelihood_2(X,mask,rho,p,Z,psi,ia,ic,num_comb)

R = squareform(rho)+eye(size(X,2));

[~, flag] = chol(R);

if flag ~= 0 
        log_likelihood = -1e10;
        return       
end

log_likelihood = 0;

for i = 1:num_comb
    rows = (i == ic);
    s    = mask(ia(i),:);
    ss   = find(s > 0);
    tt   = find(s == 0);

    Rss = R(ss, ss);
    Rtt = R(tt, tt);
    Rst = R(ss, tt);
    Rts = R(tt, ss);

    zs      = Z(rows, ss);
    zt      = Z(rows, tt);
    psi_ss  = psi(rows, ss);
    p_ss    = reshape(p(ss), 1, []); 

    if isempty(ss)
       
        log_dens = log( mvncdf(zt, zeros(1, numel(tt)), Rtt) );

    elseif isempty(tt)
       
        L = chol(Rss, 'lower');
        y = L \ zs.';
        log_copula = -0.5*sum(y.^2,1).' - sum(log(diag(L))) +0.5*sum(zs.^2,2);
                          
        log_marg   = sum( log(p_ss .* psi_ss), 2 );
        log_dens   = log_copula + log_marg;

    else
        
        L = chol(Rss, 'lower');
        y = L \ zs.';
        log_copula = -0.5*sum(y.^2,1).' - sum(log(diag(L))) +0.5*sum(zs.^2,2);
                     

        log_cdf    = log( mvncdf(zt - (Rts * (Rss \ zs.')).',...
            zeros(1, numel(tt)), Rtt - Rts * (Rss \ Rst)) );
        log_marg   = sum( log(p_ss .* psi_ss), 2 );

        log_dens   = log_copula + log_cdf + log_marg;
    end

    log_likelihood = log_likelihood + sum(log_dens);
end

end
