%backtest prende in input la timetable
%mette in output sei valori (i 2 VAR ai 2 livelli per i 3 modelli)

%% caricare i dati

clc; clear all;

addpath("ex_1")
addpath("zero_mixed")
addpath("Copula_Simulation")
addpath("Backtest/")
addpath("ex_4")

filename = "danishmulti.csv";
addpath('utilities','ex_1');
data = readDataset(filename);


%% calibrazione

%calibrating function has as input data (full, timetable)
% and the calibrating period
start = datetime("01/01/1980");
last = datetime("31/12/1983");

data_new = data_split(data, start, last);
calibrated_parameters = calibr_wrapper(data_new);

% -> mette in output cell di struct con i parametri modello per modello


%% sim
%simulazione prende in input la cell di struct di parametri calibrati

% -> mette in output una cell 3x1, in ciascuna cell una matrice Nx3

N=1000;
sim_losses = mat_sim(calibrated_parameters, N);


% Calcoliamo il var in chiaro per i 3 modelli
% in chiaro -> senza funzione separata






%quantili



%
