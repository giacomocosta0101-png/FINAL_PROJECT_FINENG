function cdf = cumulative_cdf_semi_parametric(p, X)

    % Positive observations = jump sizes
    X_jump = sort(X(X > 0));
    n = numel(X_jump);

    if n == 0
        cdf = @(x) (1-p) .* (x >= 0);
        return
    end

    % Unique jump sizes and empirical CDF evaluated at them
    [x_unique, ia] = unique(X_jump, 'last');

    % ia gives the position of the last occurrence of each unique value
    % Divide by n+1 so the empirical CDF never reaches 1
    F_unique = ia ./ (n + 1);

    cdf = @(x) semi_parametric_eval(x, p, x_unique, F_unique);

end


function u = semi_parametric_eval(x, p, x_unique, F_unique)

    u = zeros(size(x));

    % Mass at zero: P(X = 0) = 1-p
    u(x == 0) = 1 - p;

    positive = x > 0;

    if any(positive)

        F_jump = interp1( ...
            x_unique, ...
            F_unique, ...
            x(positive), ...
            'previous', ...
            'extrap' ...
        );

        % For positive values below the first observed jump
        F_jump(x(positive) < x_unique(1)) = 0;

        % Numerical safety
        F_jump = max(F_jump, 0);
        F_jump = min(F_jump, max(F_unique));

        u(positive) = (1 - p) + p .* F_jump;

    end

end