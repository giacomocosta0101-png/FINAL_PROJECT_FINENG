function cdf = cumulative_cdf_semi_parametric_pareto(p, X)

    % Positive observations = jump sizes
    X_sorted = sort(X);
    n = sum(X_sorted>0);
    
    cdf = cell(1, size(X,2));

    for i = 1:size(X,2)
        if n(i) == 0
            cdf{i} = @(x) (1-p(i)) .* (x >= 0);
        else
            Xvec = X(:,i);
            Xpos = Xvec(Xvec>0);
            pv = paretotails(Xpos,0,0.99,'ecdf');

            cdf{i} = @(x) (1-p(i)) .* (x >= 0) + p(i) .* pv.cdf(x) .* (x > 0);
        end
    end
end