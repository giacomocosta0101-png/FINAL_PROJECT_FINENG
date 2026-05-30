function VaR = var_calc(data,alpha,start_date,end_date,N)

%% calibrazione
data_new = data_split(data, start_date, end_date);
building = data_new.Building(:);
contents = data_new.Contents(:);
profits = data_new.Profits(:);

X_new = [building contents profits];
calibrated_parameters = calibr_wrapper(X_new);

% -> mette in output cell di struct con i parametri modello per modello


%% sim
%simulazione prende in input la cell di struct di parametri calibrati

% -> mette in output una cell 3x1, in ciascuna cell una matrice Nx3


sim_losses = mat_sim(calibrated_parameters, N);

VaR = zeros(length(sim_losses),length(alpha));

for i = 1:size(VaR,1)
    VaR(i,:) = quantile(sim_losses{i},1-alpha);
end


end
