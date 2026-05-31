
%per prima cosa runna il primo blocco di run per avere X
%poi spostati nella cartella zero_mixed per runnare

%% se vuoi verificare che effettivamente sia generico


for i=1:size(X,1)
    U = rand(1);
    if U < 0.7
        X(i,4) = 0;
    else 
        noise = exp(randn(1)*0.5 + 1);
        X(i,4) = X(i,3) + noise;
    end
end

%% test run
zero_mixed = zero_mixed_calibration(X);

N = height(data);
B = 10000;
alpha = 0.05;

ci = zero_mixed_bootstrap(zero_mixed, alpha, N, B);

zero_mixed_print_ci_table(ci)


%% test run fixed

zero_mixed = zero_mixed_calibration(X);

N = height(data);
B = 10000;
alpha = 0.05;

ci = zero_mixed_bootstrap_fixed(zero_mixed, alpha, N, B);

zero_mixed_print_ci_table(ci)
