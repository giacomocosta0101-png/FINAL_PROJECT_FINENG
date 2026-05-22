function log_likelihood = log_likelihood_semiparametric(mask,theta,...
    Z,ia,ic,num_comb)

L = [1 0 0;
    cos(theta(1)) sin(theta(1)) 0;
    cos(theta(2)) cos(theta(3))*sin(theta(2)) sin(theta(3))*sin(theta(2))];

R = L*L';

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

    if isempty(ss)
       
        log_dens_i = log( mvncdf(zt(1,:), zeros(1, numel(tt)), Rtt) );
        log_dens = log_dens_i*size(zt,1);

    elseif isempty(tt)
       
        log_dens= log( mvnpdf(zs, zeros(1, numel(ss)), Rss) ) ...
                   - log( mvnpdf(zs) );
    elseif isscalar(tt)

        log_copula = log( mvnpdf(zs, zeros(1, numel(ss)), Rss) ) ...
                   - log( mvnpdf(zs) );
                     

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
