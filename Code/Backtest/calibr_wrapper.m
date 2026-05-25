function calibrated_parameters = calibr_wrapper(data)
%ritorna cell di struct

calibrated_parameters = cell(3,1);

zero_mixed = struct();
zero_mixed = zero_mixed_first_calibration(data);

comb_ber = struct();
co




calibrated_parameters{1} = zero_mixed; 



end