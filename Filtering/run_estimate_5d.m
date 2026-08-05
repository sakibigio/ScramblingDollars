%% Driver: set up workspace and run estimate_params_4d.m
cd(fileparts(mfilename('fullpath')));   % this file's own folder
addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

matching_type = 1;          % Cobb-Douglas (matches paper default)
pi_us_ss = 1;
pi_eu_ss = 1;

load('data/LFX_data.mat');
load('data/calibration.mat');
run('functions/params.m');

% ---- Baseline spec (2026-07-31): penalty-free estimation on LCR-excess mu ----
% mu enters the inversion NET of the LCR requirement, matching main_filter.m
% (mu_minus_lcr = 1). US only; there is no EU requirement series.
mu_us = mu_minus_lcr_level(mu_us, LCR_us);
% No identification penalty: the corr/trend/kpss modes in estimate_params_4d.m
% are retained for reference but the baseline sets their weight to zero.
w_corr_override = 0;
% (To reproduce the pre-2026-07-31 spec, comment the two lines above.)

% Baseline correction is handled inside estimate_params_4d.m via
% use_mu_baseline_estim (defaults to 0 -- no correction, matching main_filter.m).

% Initial guesses for sigma (passed to filter)
sigma_us = 0.20;
sigma_eu = 0.15;

estimate_params_4d;
