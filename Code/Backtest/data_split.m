function data_new = data_split(data,start, last)

if nargin < 3
    last = data.Date(end);
end

%add a check last>start

last = datenum(last);
last = last+1;
last = datetime(last, "ConvertFrom", "datenum");

data_new = data(timerange(start,last), :);

end