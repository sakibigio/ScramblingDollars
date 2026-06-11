%% Batch driver: set up workspace then run estimate_params.
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1;          % Cobb-Douglas
pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat');
load('data/calibration.mat');
run('functions/params.m');
estimate_params
save('_estimation_result.mat','lambda_opt','eta_opt','iota_ss_opt','fval','exitflag');
