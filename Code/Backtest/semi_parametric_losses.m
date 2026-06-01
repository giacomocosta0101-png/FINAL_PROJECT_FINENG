function X_new = semi_parametric_losses(L, p, N, X)

    U_sim = semi_parametric_sim(L, p, N);
    X_new = zeros(N,size(X,2));
    
    cdf = cumulative_cdf_semi_parametric_vec(p, X);
    
    for i = 1:size(X,2)
        
        X_col = unique(X(:,i));
        U = cdf{i}(X_col);
        X_new(:,i) = interp1(U, X_col, U_sim(:,i), 'linear','extrap');
    end
   
end