function X_new = semi_parametric_losses(rho, p, N, X)

    U_sim = semi_parametric_sim(rho, p, N);
    X_new = zeros(N,size(X,2));
    
    for i = 1:size(X,2)
        cdf = cumulative_cdf_semi_parametric(p(i), X(:,i));
        X_col = unique(X(:,i));
        U = cdf(X_col);
        X_new(:,i) = interp1(U, X_col, U_sim(:,i), 'linear','extrap');
    end
   
end