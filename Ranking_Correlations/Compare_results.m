% Sharing score difference
load('results_inter.mat')
results_i = results;
load('results_plates.mat')
results_p = results;
results=results_i-results_p;