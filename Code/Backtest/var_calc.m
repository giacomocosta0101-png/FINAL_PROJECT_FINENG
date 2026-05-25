%backtest prende in input la timetable
%mette in output sei valori (i 2 VAR ai 2 livelli per i 3 modelli)

% caricare i dati



% calibrazione

%calibrating function has as input data (full, timetable)
% and the calibrating period

% -> mette in output cell di struct con i parametri modello per modello



%simulazione prende in input la cell di struct di parametri calibrati

% -> mette in output una cell 3x1, in ciascuna cell una matrice Nx3

N=1000;
sim_losses = mat_sim(calibrated_parameters, N);


% Calcoliamo il var in chiaro per i 3 modelli
% in chiaro -> senza funzione separata






%quantili



%
