function plot_bootstrap_rho(rho_hat, rho_point, alpha)
% PLOT_BOOTSTRAP_RHO  Pairwise scatter of the bootstrap rho distribution
% with the no-bootstrap point estimate and the (1-alpha) confidence ellipse.
%
% INPUT
%   rho_hat   : (B x 3) bootstrap rho estimates, columns = [rho_12, rho_13, rho_23]
%   rho_point : (1 x 3) original point estimate (no bootstrap)
%   alpha     : nominal significance level for the ellipse (e.g. 0.05)
%
% The ellipse is the Mahalanobis 95% region under a bivariate Gaussian fit
% to the bootstrap cloud:
%       { x : (x-mubar)' * Sboot^{-1} * (x-mubar) <= chi2inv(1-alpha, 2) }
% This is a JOINT confidence region for the pair; it is not the same thing
% as the product of the two Bonferroni-corrected marginal CIs returned by
% bootstrap.m (those are conservative rectangles, the ellipse is the
% natural elliptical region).

    if nargin < 3, alpha = 0.05; end

    pairs  = [1 2; 1 3; 2 3];
    labels = {'\rho_{12}','\rho_{13}','\rho_{23}'};
    chi2c  = chi2inv(1 - alpha, 2);

    figure('Color','w','Position',[100 100 1300 420]);
    for k = 1:3
        i = pairs(k,1);  j = pairs(k,2);
        ax = subplot(1,3,k); hold(ax,'on');

        % Bootstrap cloud
        scatter(rho_hat(:,i), rho_hat(:,j), 10, [0.45 0.65 0.95], 'filled', ...
                'MarkerFaceAlpha', 0.30, 'DisplayName','Bootstrap');

        % Mahalanobis ellipse from bootstrap covariance
        mu_b = mean(rho_hat(:, [i j]), 1);
        S_b  = cov( rho_hat(:, [i j]) );
        [ex, ey] = ellipse_xy(mu_b, S_b, chi2c, 200);
        plot(ex, ey, 'r-', 'LineWidth', 1.8, ...
             'DisplayName', sprintf('%.0f%% Mahalanobis ellipse', 100*(1-alpha)));

        % Point estimate (no-bootstrap)
        plot(rho_point(i), rho_point(j), 'p', 'MarkerSize', 16, ...
             'MarkerFaceColor', [1 0.85 0], 'MarkerEdgeColor', 'k', ...
             'LineWidth', 1.2, 'DisplayName', 'Point estimate');

        % Bootstrap mean
        plot(mu_b(1), mu_b(2), '+', 'MarkerSize', 12, 'LineWidth', 2, ...
             'Color', [0.6 0 0], 'DisplayName', 'Bootstrap mean');

        xlabel(labels{i}); ylabel(labels{j});
        title(sprintf('%s vs %s', labels{i}, labels{j}));
        grid on; axis equal;
        if k == 3
            legend('Location','bestoutside');
        end
    end
    sgtitle('Bootstrap distribution of Gaussian-copula correlations');
end

% --- helper ----------------------------------------------------------------
function [x, y] = ellipse_xy(mu, S, chi2c, npts)
% Points on the ellipse (x-mu)' S^{-1} (x-mu) = chi2c
    t = linspace(0, 2*pi, npts);
    L = chol((S+S.')/2, 'lower');     % simmetrizza per sicurezza numerica
    XY = sqrt(chi2c) * L * [cos(t); sin(t)] + mu(:);
    x = XY(1,:);  y = XY(2,:);
end