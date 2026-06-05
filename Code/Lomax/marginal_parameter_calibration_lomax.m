function [alpha, lambda, p] = marginal_parameter_calibration_lomax(X)
    d = size(X, 2);
    p     = mean(X > 0, 1);
    alpha = zeros(1, d);
    lambda = zeros(1, d);
    
    for i = 1:d
    xp = X(X(:,i) > 0, i); % track the strictly positive observations
    n  = numel(xp);

    % alpha is profiled outside -> the only parameter to which maximize the
    % likelihood is lambda:

    negll = @(lam) profiled_negll(lam, xp, n);
    lam_hat  = fminsearch(negll, mean(xp));

    alpha(i)  = n / sum(log(1 + xp/lam_hat));
    lambda(i) = lam_hat;
    end
end

function f = profiled_negll(lam, xp, n)
    if lam <= 0, f = Inf; return; end          % lambda > 0
    S = sum(log(1 + xp/lam));
    % loglikelihood concentrated (alpha = n/S from (dl(alpha,lam)/dalpha) = 0)
    % Then alpha into the likelihood:
    %   l(lam) = n*log(n/S) - n*log(lam) - n - S
    
    f = -( n*log(n/S) - n*log(lam) - n - S );
end