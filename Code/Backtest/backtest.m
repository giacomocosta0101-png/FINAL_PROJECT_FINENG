function [backtest_window,exceptions,VaR] = backtest(data,alpha,...
    window_start,window_end,N,mode)
% BACKTEST  VaR backtest over an out-of-sample window (fixed or rolling).
%
%   Calibrates the claim models on a training window, simulates N loss
%   scenarios, and compares the resulting VaR against the realized total
%   losses on each day of the evaluation window (backtest_window).
%
%   INPUT:
%     data         : timetable indexed by Date, with columns Building,
%                    Contents, Profits and Total (per-day losses).
%     alpha        : 1xnAlpha vector of significance levels; the VaR is the
%                    (1-alpha) quantile of the simulated losses (e.g.
%                    [0.05 0.005] -> VaR 95% and 99.5%).
%     window_start : datetime, first day of the training window
%                    (the initial one, in rolling mode).
%     window_end   : datetime, last day of the training window
%                    (the initial one, in rolling mode).
%     N            : number of Monte Carlo scenarios drawn per calibration.
%     mode         : 'Fixed'          -> calibrate once, reuse the same VaR
%                                        over the whole window;
%                    'Rolling-window' -> re-calibrate each day on a window of
%                                        constant length rolled forward by 1d.
%
%   OUTPUT:
%     backtest_window : timetable of the evaluation period (from the day after
%                       window_end to the end of data).
%     exceptions      : 1x3 cell, one entry per model; each is an MxnAlpha
%                       logical, true where the realized Total exceeds VaR.
%     VaR             : Value-at-Risk estimates.
%                       'Fixed'          -> 3 x nAlpha (one row per model);
%                       'Rolling-window' -> 3 x nAlpha x M (3rd dim = day).
%
%   NOTE:
%     - The leading dimension 3 is the number of models (Zero-mixed,
%       Comb-Bernoulli, Semi-parametric); nAlpha = numel(alpha);
%       M = height(backtest_window).
%     - VaR changes shape between modes (2-D in 'Fixed', 3-D in
%       'Rolling-window'): any consumer must therefore branch on mode
%       (see plot_backtest).

arguments
    data timetable
    alpha (1,:) double {mustBeNonempty, mustBeReal, mustBeFinite}
    window_start {mustBeDatetimeTextScalar}
    window_end {mustBeDatetimeTextScalar}
    N (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
    mode {mustBeTextScalar}
end

%% Input checks

if any(alpha <= 0 | alpha >= 1)
    error('backtest:invalidAlpha', ...
        'alpha must contain values strictly between 0 and 1.');
end

window_start = localDatetimeScalar(window_start, 'window_start');
window_end = localDatetimeScalar(window_end, 'window_end');
mode = string(mode);

if window_start > window_end
    error('backtest:invalidWindow', ...
        'window_start must be earlier than or equal to window_end.');
end

if ~isscalar(mode) || ~any(mode == ["Fixed","Rolling-window"])
    error('backtest:unknownMode', ...
        'mode must be ''Fixed'' or ''Rolling-window''.');
end

if mode == "Rolling-window" && numel(alpha) ~= 2
    error('backtest:rollingAlphaCount', ...
        'Rolling-window mode currently expects exactly two alpha levels.');
end

if window_end >= data.Date(end)
    error('backtest:emptyBacktestWindow', ...
        'window_end must be strictly earlier than the last date in data.');
end

% Consider the backtest window (evaluation) from the next day after 
% the last day of the training window:

backtest_window = data_split(data,window_end+caldays(1));

    if mode == "Fixed"
        % if the training window is 'static', then it is enough to compute
        % one only replica of the dataset and use it to backtest throughout
        % the whole backtest window
        VaR = var_calc(data,alpha,window_start,...
            window_end,N);

        exceptions = cell(1,3);
    for i= 1:3
        exceptions{i} = (backtest_window.Total > VaR(i,:));
    end

    elseif mode == "Rolling-window"
        
        % First/last day of the INITIAL training window, taken from the 
        % inputs. i0 and j are the index in the original dataset associated 
        % to the training start and end date given in input:

        i0 = find(data.Date==window_start);
        j = find(data.Date==window_end);
        VaR = zeros(3,2,size(backtest_window,1));
        exceptions = cell(1,3);
    
    % Set the time to keep track of the backtest
    t0 = tic;

    for i = 1:size(backtest_window,1)
        % Roll the window by one day: both ends shift by +1, so the number 
        % of observations in the window stays constant.
        window_start = data.Date(i0 + i - 1);
        window_end = data.Date(j);
        VaR(:,:,i) = var_calc(data,alpha,window_start,...
            window_end,N);

        j = j+1;

        for k= 1:3
            exceptions{k}(i,:) = (backtest_window.Total(i) > VaR(k,:,i));
        end

        if mod(i, 50) == 0 || i == size(backtest_window,1)   % Prints every 50 iterations
            elapsed = toc(t0);
            eta = elapsed / i * (size(backtest_window,1) - i);
            fprintf('Day %4d/%d | elapsed %6.1fs | ETA %6.1fs\n', ...
                i, size(backtest_window,1), elapsed, eta);
        end
    end

    else
        
        % Guard against an unrecognised mode:
        error('backtest:unknownMode', ...
            "Unknown mode '%s'. Use 'Fixed' or 'Rolling-window'.", mode);
    end                  
end

function x = localDatetimeScalar(x, name)
% LOCALDATETIMESCALAR  Convert a date-like input to a scalar datetime.

    if isdatetime(x)
        if ~isscalar(x)
            error('backtest:dateNotScalar', ...
                '%s must be a scalar datetime.', name);
        end
        return
    end

    if ischar(x) || (isstring(x) && isscalar(x))
        x = localParseDateText(x, name);
        return
    end

    error('backtest:invalidDateType', ...
        '%s must be a scalar datetime or a scalar text value.', name);
end

function x = localParseDateText(x, name)
% LOCALPARSEDATETEXT  Parse supported date text without locale guessing.

    x = string(x);
    formats = ["yyyy-MM-dd", "dd/MM/yyyy"];

    for i = 1:numel(formats)
        try
            candidate = datetime(x, 'InputFormat', formats(i));
        catch
            candidate = NaT;
        end

        if ~isnat(candidate)
            x = candidate;
            return
        end
    end

    error('backtest:invalidDateFormat', ...
        '%s must use yyyy-MM-dd or dd/MM/yyyy.', name);
end

function mustBeDatetimeTextScalar(x)
% MUSTBEDATETIMETEXTSCALAR  Validate a scalar datetime or scalar text input.

    is_text_scalar = ischar(x) || (isstring(x) && isscalar(x));
    is_datetime_scalar = isdatetime(x) && isscalar(x);

    if ~(is_text_scalar || is_datetime_scalar)
        error('backtest:invalidDateInput', ...
            'Date inputs must be scalar datetimes or scalar text values.');
    end
end
