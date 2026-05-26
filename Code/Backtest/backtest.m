function [backtest_window,exceptions,VaR] = backtest(data,alpha,start_date,end_date,N,mode)


backtest_window = data_split(data,end_date+caldays(1));

if strcmp(mode,'Fixed')
    VaR = var_calc(data,alpha,start_date,end_date,N);

    exceptions = cell(1,3);
    for i= 1:3
        exceptions{i} = (backtest_window.Total > VaR(i,:));
    end

elseif strcmp(mode,'Rolling-window')
    j = find(data.Date==end_date);
    VaR = zeros(3,2,size(backtest_window,1));
    exceptions = cell(1,3);

    for i = 1:size(backtest_window,1)
        start_date = data.Date(i);
        end_date = data.Date(j);
        VaR(:,:,i) = var_calc(data,alpha,start_date,end_date,N);

        j = j+1;

        for k= 1:3
            exceptions{k}(i,:) = (backtest_window.Total(i) > VaR(k,:,i));
        end

    end
end
