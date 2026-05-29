function X_new = semi_parametric_losses(rho, p, N, pts)
    U_sim = semi_parametric_sim(rho, p, N);
    d     = numel(pts);
    X_new = zeros(N, d);

    for i = 1:d
        u    = U_sim(:,i);
        mass = u <= (1 - p(i));
        v    = (u - (1 - p(i))) ./ p(i);
        v    = min(max(v, 1e-12), 1 - 1e-12);
        X_new(~mass, i) = icdf(pts{i}, v(~mass));
        % gli `mass` restano 0 (init)
    end
end