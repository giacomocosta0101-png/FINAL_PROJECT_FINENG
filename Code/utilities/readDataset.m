function data = readDataset(filename)

    raw_data = readtable(filename,'ReadRowNames',false);
    

    day_zero = datetime('1980-01-01');     
    C = {day_zero,0,0,0,0};
    day_zero = cell2table(C,'VariableNames',{'Date','Building','Contents','Profits','Total'});

    raw_data = [day_zero; raw_data];
    
    %Group together duplicates:
    data_adj = groupsummary(raw_data,"Date","sum");
    data_adj = removevars(data_adj,{'GroupCount'});
 
    data_adj = table2timetable(data_adj);

    data_filled= retime(data_adj,"daily");

    missingIdx = ~ismember(data_filled.Date,data_adj.Date);

    data_filled(missingIdx,:) = {0};

    data = data_filled;

    data.Properties.VariableNames = ["Building","Contents","Profits","Total"];

end




