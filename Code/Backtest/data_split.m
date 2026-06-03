function data_new = data_split(data,start, last)
% DATA_SPLIT  Extract a time-based subset of a timetable.
%
%   INPUT
%     data  : input timetable
%     start : scalar datetime or scalar text convertible to datetime
%     last  : scalar datetime or scalar text convertible to datetime
%             (optional, default = last date in data)
%
%   OUTPUT
%     data_new : timetable subset over the inclusive [start, last] range

arguments
    data timetable
    start {mustBeDatetimeTextScalar}
    last {mustBeDatetimeTextScalar} = []
end

%% Input checks

if isempty(last)
    last = data.Date(end);
end

start = localDatetimeScalar(start, 'start');
last = localDatetimeScalar(last, 'last');

if start > last
    error('data_split:invalidRange', ...
        'start must be earlier than or equal to last.');
end

data_new = data(data.Date >= start & data.Date <= last, :);

end

function x = localDatetimeScalar(x, name)
% LOCALDATETIMESCALAR  Convert a date-like input to a scalar datetime.

    if isdatetime(x)
        if ~isscalar(x)
            error('data_split:dateNotScalar', ...
                '%s must be a scalar datetime.', name);
        end
        return
    end

    if ischar(x) || (isstring(x) && isscalar(x))
        x = localParseDateText(x, name);
        return
    end

    error('data_split:invalidDateType', ...
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

    error('data_split:invalidDateFormat', ...
        '%s must use yyyy-MM-dd or dd/MM/yyyy.', name);
end

function mustBeDatetimeTextScalar(x)
% MUSTBEDATETIMETEXTSCALAR  Validate a scalar datetime or scalar text input.

    if isempty(x)
        return
    end

    is_text_scalar = ischar(x) || (isstring(x) && isscalar(x));
    is_datetime_scalar = isdatetime(x) && isscalar(x);

    if ~(is_text_scalar || is_datetime_scalar)
        error('data_split:invalidDateInput', ...
            'Date inputs must be scalar datetimes or scalar text values.');
    end
end
