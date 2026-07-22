%% main_robustness.m
%  Driver for the referee-3 robustness exercise (risk aversion).
%
%  Pipeline
%    main_filter.m           (baseline run; produces RW_shock.csv and
%                             data/MS_sigma_us_prob.csv)
%    main_robustness.m       (this file: grids iota and ploss,
%                             re-estimates lambda, runs Julia per pair)
%    plotting/plot_robustness.m   (called internally for both grids)
%
%  Flags below choose which sub-exercises to run.
%
% (c) Saki Bigio

clear; close all;

%% Paths & options
addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

% --- Flags ---
do_estimate_iota   = 0;
do_estimate_ploss  = 0;
do_plot_iota       = 1;
do_plot_ploss      = 1;
run_julia          = 0;    % set 0 to skip per-pair Markov estimation

% --- Setup workspace (matches estimate_params.m prelude) ---
matching_type = 1;          % 1 = Cobb-Douglas, 0 = Leontief
pi_eu_ss = 1;
pi_us_ss = 1;

load('data/LFX_data.mat');
load('data/calibration.mat');
run('functions/params.m');

% Baseline values (frozen before any robustness loop overwrites them)
lambda_baseline = lambda_us;
eta_baseline    = eta;
iota_baseline   = iota_ss;
ploss_baseline  = ploss_us;

% Scale used throughout the filter
abs_scale = 12e4;

% Date variables
T = length(mu_us);
datesperiod = 1:T;
dates = datenum(2001, 1:T, 1);

% Baseline correction (matches main_filter.m)
use_mu_baseline = 0.1;       % keep parity with main_filter.m default
if use_mu_baseline == 1
    [mubar_us_t, mubar_eu_t] = liquidity_baseline(dates);
else
    mubar_us_t = zeros(T,1);
    mubar_eu_t = zeros(T,1);
end

% varrho default (params.m sets to 0; restate for clarity)
if ~exist('varrho', 'var'), varrho = 0; end

% --- Locate Julia (mirrors main_filter.m logic) ---
julia_path = '';
if run_julia == 1
    [status, jp] = system('which julia');
    if status == 0
        julia_path = strtrim(jp);
    else
        candidates = {
            '/usr/local/bin/julia'
            '/opt/homebrew/bin/julia'
            '/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia'
            '/Users/sakibigio/.juliaup/bin/julia'
        };
        for kk = 1:numel(candidates)
            if exist(candidates{kk}, 'file')
                julia_path = candidates{kk}; break
            end
        end
    end
    if isempty(julia_path)
        warning('main_robustness: Julia not found; Markov re-estimation will be skipped.');
    else
        fprintf('main_robustness: Julia at %s\n', julia_path);
    end
end

fprintf('\nBaseline params: lambda=%.4f  eta=%.4f  iota_ss=%.4f  ploss=%.4f\n', ...
    lambda_baseline, eta_baseline, iota_baseline, ploss_baseline);

%% Run sub-exercises
if do_estimate_iota
    run('estimate_robustness_iota.m');
end

if do_estimate_ploss
    run('estimate_robustness_ploss.m');
end

if do_plot_iota
    which_grid = 'iota';
    run('plotting/plot_robustness.m');
end

if do_plot_ploss
    which_grid = 'ploss';
    run('plotting/plot_robustness.m');
end

fprintf('\nmain_robustness complete.\n');
