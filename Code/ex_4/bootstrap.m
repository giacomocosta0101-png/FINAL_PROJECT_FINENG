function [rho_CI, p_CI]= bootstrap(rho,p,mu,sigma,model,N,B,alpha)

    if strcmp(model,'Comb-Bernoulli')
        flag = 1;
        generate_replica = @(corr)comb_bern_sim(corr, mu, sigma, p, N);
        cdf = marginal_cdf(p,mu,sigma);
    elseif strcmp(model,'Semi-parametric')
        generate_replica = @(corr)semi_parametric_sim(corr,p,N);
        flag = 2;
    end
    
    rho_hat = zeros(B,3);
    pi_hat = zeros(B,3);

    
    rho_CI = zeros(3,2);
    p_CI = zeros(3,2);
    
    t0 = tic;
    for i = 1:B
        
        replica = generate_replica(rho);
        
        if flag == 1
            U = cdf(replica);
        else
            U = replica;
        end
        
        [rho_hat(i,:),pi_hat(i,:)] = calibrate_model(U,p);

        if mod(i, 50) == 0 || i == B   % stampa ogni 50 iterazioni
        elapsed = toc(t0);
        eta = elapsed / i * (B - i);
        fprintf('Iter %4d/%d | trascorso %6.1fs | ETA %6.1fs\n', ...
                i, B, elapsed, eta);
        end

    end
    
    %Dobbiamo fare correzione di bonferroni e capire:
    alpha_correct = alpha/3;
    for i = 1:3
        rho_CI(i,:) = quantile(rho_hat(:,i), [alpha_correct/2, 1-alpha_correct/2]);
        p_CI(i,:) = quantile(pi_hat(:,i), [alpha_correct/2, 1-alpha_correct/2]);
    end

end




