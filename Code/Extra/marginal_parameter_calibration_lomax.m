function [alpha, lambda, p] = marginal_parameter_calibration_lomax(X)
    d = size(X, 2);
    p     = mean(X > 0, 1);
    alpha = zeros(1, d);
    lambda = zeros(1, d);

    for i = 1:d
        xp = X(X(:,i) > 0, i);
        m  = mean(xp);
        s2 = var(xp, 1);          % MLE-style denominator n
        CV2 = s2 / m^2;

        if CV2 <= 1
            warning('Column %d: CV^2 = %.3f <= 1, Lomax MoM invalid.', i, CV2);
            alpha(i)  = Inf;       % degenera all'esponenziale
            lambda(i) = m;
        else
            alpha(i)  = 2*CV2 / (CV2 - 1);
            lambda(i) = m * (alpha(i) - 1);
        end
    end

    marginal_params.p      = p;
    marginal_params.alpha  = alpha;
    marginal_params.lambda = lambda;
end