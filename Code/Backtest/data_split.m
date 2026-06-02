function data_new = data_split(data,start, last)
% DATA_SPLIT  Extract a time-based subset of a timetable.
%
% INPUT
%   data  : (timetable) the input data containing time-stamped rows
%   start : (datetime or string) the start date of the desired subset
%   last  : (datetime or string) the end date of the desired subset (optional)
%
% OUTPUT
%   data_new : (timetable) the filtered subset of data

if nargin < 3
    last = data.Date(end);
end

%add a check last>start

last = datenum(last);
last = last+1;
last = datetime(last, "ConvertFrom", "datenum");

data_new = data(timerange(start,last), :);

end