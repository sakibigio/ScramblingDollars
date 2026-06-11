%% Driver: set up workspace and run estimate_params_4d.m
cd('/Users/sakibigio/Dropbox/Scrabmling for Dollars/Code 2025 ScramblingDollars/Filtering');
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

% Baseline correction is handled inside estimate_params_4d.m via
% use_mu_baseline_estim (defaults to 0 -- no correction, matching main_filter.m).

% Initial guesses for sigma (passed to filter)
sigma_us = 0.20;
sigma_eu = 0.15;

estimate_params_4d;
