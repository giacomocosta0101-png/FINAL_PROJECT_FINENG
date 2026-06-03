function VaR = var_calc(data,alpha,window_start,...
    window_end,N)
% VAR_CALC  Calculate Value at Risk (VaR) across different copula models.
%
% This function isolates a specific time window from the dataset, calibrates 
% three copula models (Zero-Mixed, Comb-Bernoulli, Semi-Parametric) using 
% that data, simulates total losses, and computes the empirical VaR for 
% the specified significance levels.
%
% INPUT
%   data         : (timetable) full dataset containing 'Building', 'Contents', 
%                  and 'Profits' variables
%   alpha        : (scalar or vector) significance level(s) for the VaR 
%                  (e.g., 0.05 for a 95% confidence level)
%   window_start : (datetime or string) start of the calibration window
%   window_end   : (datetime or string) end of the calibration window
%   N            : (scalar) number of Monte Carlo simulations to run
%
% OUTPUT
%   VaR          : (3 x m) matrix of VaR estimates, where rows correspond 
%                  to the 3 models and columns correspond to the m alpha levels

arguments
    data timetable
    alpha (1,:) double {mustBeNonempty, mustBeReal, mustBeFinite}
    window_start {mustBeDatetimeTextScalar}
    window_end {mustBeDatetimeTextScalar}
    N (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
end

%% Input checks

if any(alpha <= 0 | alpha >= 1)
    error('var_calc:invalidAlpha', ...
        'alpha must contain values strictly between 0 and 1.');
end

window_start = localDatetimeScalar(window_start, 'window_start');
window_end = localDatetimeScalar(window_end, 'window_end');

if window_start > window_end
    error('var_calc:invalidWindow', ...
        'window_start must be earlier than or equal to window_end.');
end

%% Calibration
data_new = data_split(data, window_start, ...
    window_end);
if isempty(data_new)
    error('var_calc:emptyWindow', ...
        'The calibration window contains no observations.');
end
building = data_new.Building(:);
contents = data_new.Contents(:);
profits = data_new.Profits(:);

X_new = [building contents profits];

% Pass the calibration window to a wrapper:

calibrated_parameters = calibr_wrapper(X_new);

%% Simulation

% mat_sim takes the 3x1 cell of calibrated-parameter structs and returns a
% 3x1 cell (one entry per model). Each entry is an N x 1 vector of simulated
% TOTAL losses: mat_sim already sums the 3 risk components (Building,
% Contents, Profits) internally.
sim_losses = mat_sim(calibrated_parameters, N);

VaR = zeros(length(sim_losses),length(alpha));

for i = 1:size(VaR,1)
    VaR(i,:) = quantile(sim_losses{i},1-alpha);
end


end

function x = localDatetimeScalar(x, name)
% LOCALDATETIMESCALAR  Convert a date-like input to a scalar datetime.

    if isdatetime(x)
        if ~isscalar(x)
            error('var_calc:dateNotScalar', ...
                '%s must be a scalar datetime.', name);
        end
        return
    end

    if ischar(x) || (isstring(x) && isscalar(x))
        x = localParseDateText(x, name);
        return
    end

    error('var_calc:invalidDateType', ...
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

    error('var_calc:invalidDateFormat', ...
        '%s must use yyyy-MM-dd or dd/MM/yyyy.', name);
end

function mustBeDatetimeTextScalar(x)
% MUSTBEDATETIMETEXTSCALAR  Validate a scalar datetime or scalar text input.

    is_text_scalar = ischar(x) || (isstring(x) && isscalar(x));
    is_datetime_scalar = isdatetime(x) && isscalar(x);

    if ~(is_text_scalar || is_datetime_scalar)
        error('var_calc:invalidDateInput', ...
            'Date inputs must be scalar datetimes or scalar text values.');
    end
end
