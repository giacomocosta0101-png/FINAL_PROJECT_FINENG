function cdf = marginal_cdf_lomax(p, alpha, lambda)
    eps_u = 1e-10;
    cdf = @(X) clamp_unit( eval_lomax_cdf(X, p, alpha, lambda), eps_u );
end

function U = eval_lomax_cdf(X, p, alpha, lambda)
    U = zeros(size(X));
    for j = 1:size(X, 2)
        atom = (X(:, j) == 0);
        pos  = (X(:, j) > 0);
        U(atom, j) = 1 - p(j);
        U(pos , j) = (1 - p(j)) + p(j) .* ...
                     (1 - (1 + X(pos, j)./lambda(j)).^(-alpha(j)));
    end
end

function U = clamp_unit(U, e)
    U = min(max(U, e), 1 - e);
end