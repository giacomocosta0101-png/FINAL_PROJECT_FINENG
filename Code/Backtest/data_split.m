function data = data_split(data,start, last)
%datetime

%add a check last>start

last = datenum(last);
last = last+1;
last = datetime(last, "ConvertFrom", "datenum");

data = data(timerange(start,last), :);

end