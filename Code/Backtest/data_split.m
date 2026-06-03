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
    start datetime
    last datetime = data.Date(end)
end

%add a check last>start

last = datenum(last);
last = last+1;
last = datetime(last, "ConvertFrom", "datenum");

data_new = data(timerange(start,last), :);

end
