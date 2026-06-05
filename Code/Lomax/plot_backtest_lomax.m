function plot_backtest_lomax(backtest_window, exceptions, VaR, mode, varargin)
% PLOT_BACKTEST_LOMAX  Visualize VaR backtesting results for a SINGLE model.
%
% Versione "single-model" di plot_backtest: var_calc_lomax e backtest_lomax
% ora restituiscono il VaR di un solo modello (Comb-Bernoulli Lomax), quindi
% gli output NON sono più impacchettati in una cell 1x3. Questa funzione
% disegna una sola figura con le due soglie di VaR (95% e 99%) e le relative
% eccezioni.
%
% INPUT
%   backtest_window : (timetable) dati di backtest con 'Date' e 'Total'
%   exceptions      : (M x 2) matrice logica/numerica; colonna 1 = eccezioni
%                     vs VaR 95%, colonna 2 = eccezioni vs VaR 99%.
%                     (in 'Rolling-window' accetta anche una 1 x 2 x M).
%   VaR             : stime di VaR per UN modello.
%                       'Fixed'          -> 1 x 2   (le due soglie);
%                       'Rolling-window' -> 1 x 2 x M (3a dim = giorno).
%   mode            : (string) 'Fixed' oppure 'Rolling-window'
%
% OPTIONAL PARAMETERS (Name-Value pairs)
%   'ModelName'     : (string) nome del modello (default: 'Lomax Comb-Bernoulli')
%   'TopK'          : (scalar) numero di eccezioni peggiori da annotare (default: 5)
%   'AnnotateLevel' : (scalar/string) 99, 95, o 'both' (default: 99)
%
% ESEMPI
%   plot_backtest_lomax(bw, exc, VaR, 'Fixed')
%   plot_backtest_lomax(..., 'ModelName', 'Lomax')
%   plot_backtest_lomax(..., 'TopK', 5)
%   plot_backtest_lomax(..., 'AnnotateLevel', 99)   % 99 | 95 | 'both'

arguments
    backtest_window timetable
    exceptions {mustBeNonempty}
    VaR double {mustBeNonempty, mustBeReal, mustBeFinite}
    mode {mustBeTextScalar}
end

arguments (Repeating)
    varargin
end

%% Input checks

mode = string(mode);
if ~isscalar(mode) || ~any(mode == ["Fixed","Rolling-window"])
    error('plot_backtest_lomax:invalidMode', ...
        'mode must be ''Fixed'' or ''Rolling-window''.');
end

N = height(backtest_window);

% Normalizza le eccezioni a una matrice M x 2 (gestisce anche la 1 x 2 x M
% che esce dalla modalità rolling).
if ndims(exceptions) == 3 && isequal(size(exceptions), [1 2 N])
    exceptions = squeeze(exceptions).';   % -> N x 2
end
if ~(isnumeric(exceptions) || islogical(exceptions)) || ~isequal(size(exceptions), [N 2])
    error('plot_backtest_lomax:invalidExceptionShape', ...
        'exceptions must be an %d x 2 logical/numeric matrix.', N);
end
exceptions = logical(exceptions);

% Controllo di forma sul VaR a seconda della modalità.
if mode == "Fixed"
    if ~isequal(size(VaR), [1 2])
        error('plot_backtest_lomax:invalidFixedVaRShape', ...
            'In Fixed mode, VaR must have size 1 x 2.');
    end
else
    if ndims(VaR) ~= 3 || ~isequal(size(VaR), [1 2 N])
        error('plot_backtest_lomax:invalidRollingVaRShape', ...
            'In Rolling-window mode, VaR must have size 1 x 2 x %d.', N);
    end
end

%% Opzioni
p = inputParser;
p.addParameter('ModelName', 'Lomax Comb-Bernoulli');
p.addParameter('TopK', 5);
p.addParameter('AnnotateLevel', 99);
p.parse(varargin{:});
model_name = p.Results.ModelName;
topK       = p.Results.TopK;
annLevel   = p.Results.AnnotateLevel;

dates  = backtest_window.Date;
losses = backtest_window.Total;
N      = numel(dates);

col95 = [0.10 0.55 0.20];
col99 = [0.80 0.15 0.15];

%% Figura (una sola, un solo modello)
figure('Name', model_name, 'Color','w', 'Position',[100 100 1150 520]);
ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on');

if mode == "Fixed"
    VaR95 = VaR(1) * ones(N,1);
    VaR99 = VaR(2) * ones(N,1);
else
    VaR95 = squeeze(VaR(1,1,:));
    VaR99 = squeeze(VaR(1,2,:));
end

plot(dates, VaR95, '-', 'Color', col95, 'LineWidth', 1.6, ...
     'DisplayName','VaR 95%');
plot(dates, VaR99, '-', 'Color', col99, 'LineWidth', 1.6, ...
     'DisplayName','VaR 99%');

scatter(dates, losses, 22, [0.45 0.45 0.45], 'filled', ...
        'MarkerFaceAlpha', 0.45, 'DisplayName','Loss');

only95 = exceptions(:,1) & ~exceptions(:,2);
idx95  = find(only95);
idx99  = find(exceptions(:,2));

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

title(sprintf('Backtest VaR – %s (%s)', model_name, mode));
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
        if exceptions(i,2), c = col99; end

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
