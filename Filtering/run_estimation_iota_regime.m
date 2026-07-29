%% Batch driver: set up workspace then run estimate_params_iota_regime.
% Estimates (lambda, iota_pre, iota_post, eta) with a GFC break in iota at
% t_break = 94 (Oct-2008).
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1;          % Cobb-Douglas
pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat');
load('data/calibration.mat');
run('functions/params.m');
estimate_params_iota_regime
save('_estimation_result_iota_regime.mat', ...
    'lambda_opt', 'eta_opt', 'iota_pre_opt', 'iota_post_opt', 't_break', ...
    'fval', 'exitflag');
