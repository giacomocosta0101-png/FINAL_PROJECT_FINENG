function plot_backtest(backtest_window, exceptions, VaR, mode, varargin)
% PLOT_BACKTEST  Backtest VaR plot per ogni modello.
%
%   plot_backtest(bw, exc, VaR, mode)
%   plot_backtest(..., 'ModelNames',    {'M1','M2','M3'})
%   plot_backtest(..., 'TopK',          5)      % quante eccezioni annotare
%   plot_backtest(..., 'AnnotateLevel', 99)     % 99 | 95 | 'both'

    p = inputParser;
    p.addParameter('ModelNames',   {'Modello 1','Modello 2','Modello 3'});
    p.addParameter('TopK',         5);
    p.addParameter('AnnotateLevel', 99);
    p.parse(varargin{:});
    model_names = p.Results.ModelNames;
    topK        = p.Results.TopK;
    annLevel    = p.Results.AnnotateLevel;

    dates  = backtest_window.Date;
    losses = backtest_window.Total;
    N      = numel(dates);

    col95 = [0.10 0.55 0.20];
    col99 = [0.80 0.15 0.15];

    for m = 1:3
        figure('Name', model_names{m}, 'Color','w', ...
               'Position',[100 100 1150 520]);
        ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on');

        if strcmpi(mode,'Fixed')
            VaR95 = VaR(m,1) * ones(N,1);
            VaR99 = VaR(m,2) * ones(N,1);
        else
            VaR95 = squeeze(VaR(m,1,:));
            VaR99 = squeeze(VaR(m,2,:));
        end

        plot(dates, VaR95, '-', 'Color', col95, 'LineWidth', 1.6, ...
             'DisplayName','VaR 95%');
        plot(dates, VaR99, '-', 'Color', col99, 'LineWidth', 1.6, ...
             'DisplayName','VaR 99%');

        scatter(dates, losses, 22, [0.45 0.45 0.45], 'filled', ...
                'MarkerFaceAlpha', 0.45, 'DisplayName','Loss');

        exc    = exceptions{m};
        only95 = exc(:,1) & ~exc(:,2);
        idx95  = find(only95);
        idx99  = find(exc(:,2));

        if ~isempty(idx95)
            scatter(dates(idx95), losses(idx95), 50, col95, 'filled', ...
                    'MarkerEdgeColor','k','LineWidth',0.5, ...
                    'DisplayName', sprintf('Eccezione 95%% (%d)', numel(idx95)));
        end
        if ~isempty(idx99)
            scatter(dates(idx99), losses(idx99), 70, col99, 'filled', ...
                    'MarkerEdgeColor','k','LineWidth',0.8, ...
                    'DisplayName', sprintf('Eccezione 99%% (%d)', numel(idx99)));
        end

        title(sprintf('Backtest VaR – %s (%s)', model_names{m}, mode));
        xlabel('Data'); ylabel('Loss / VaR');
        legend('Location','best');

        % --- annoto solo le TopK eccezioni più severe ---
        switch annLevel
            case 99,   cand = idx99;
            case 95,   cand = idx95;
            otherwise, cand = unique([idx95; idx99]);
        end

        if ~isempty(cand) && topK > 0
            [~, ord] = sort(losses(cand), 'descend');
            cand = cand(ord(1:min(topK, numel(cand))));

            ybase = min(ylim(ax));
            for k = 1:numel(cand)
                i = cand(k);
                c = col95;
                if exc(i,2), c = col99; end

                % linea verticale tratteggiata dal marker all'asse x
                plot([dates(i) dates(i)], [ybase losses(i)], '--', ...
                     'Color', c, 'LineWidth', 0.8, 'HandleVisibility','off');

                % label vicino al marker, ruotata 45°
                text(dates(i), losses(i), [' ' datestr(dates(i),'dd/mm/yy')], ...
                     'Rotation', 45, ...
                     'HorizontalAlignment','left', ...
                     'VerticalAlignment','bottom', ...
                     'FontSize', 9, 'Color', c, 'FontWeight','bold');
            end
        end

        hold(ax,'off');
    end
end