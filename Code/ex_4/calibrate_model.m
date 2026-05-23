function [rho,p_hat] = calibrate_model(U,p)

    mask = zeros(size(U));
    p_hat = zeros(3,1);

    for i = 1:size(U,2)
        idx_active = find(U(:,i)>1-p(i));
        mask(idx_active,i) = 1;
        p_hat(i) = nnz(idx_active)/size(U,1);
    end
   
    
    [C,ia,ic] = unique(mask,'rows');
    num_comb = size(C,1);
    
    Z = norminv(U);
    
    nll= @(theta)-log_likelihood_semiparametric(mask,theta,Z,ia,ic,num_comb);
    
    options = optimoptions('fminunc', 'Display', 'off');
    theta_opt = fminunc(nll, [pi/2 pi/2 pi/2],options);
    
    rho(1) = cos(theta_opt(1));
    rho(2) = cos(theta_opt(2));
    rho(3) = cos(theta_opt(1))*cos(theta_opt(2))+sin(theta_opt(1))*cos(theta_opt(3))*sin(theta_opt(2));

end
