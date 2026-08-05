%% Liquidity Exchange Rate Model
%
% Version 1: Los Angeles, July, 2019
% Version 9: Los Angeles, February, 2025
% Version 10: Refactored into subscripts, February 2026
%
% Subscripts (in functions/):
%   solve_steady_state.m  — Calibration, SS solution, tests, report
%   setup_markov.m        — Regime parameters, Tauchen grids, transition matrix
%   solve_global.m        — Exogenous paths, fsolve global equilibrium, save
%   simulate_model.m      — MC simulation, regressions, ergodic distribution
%   compute_moments.m     — Theoretical + simulation moments, LaTeX tables
% Plotting (in plotting/):
%   plot_paper_figures.m  — All paper figures + phase diagrams
% Archive:
%   archive/main_LFX_commented.m — Inactive code removed during cleanup

% TO DO LIST:
% [1] Solves for Steady State of the model..
% [2] check interior solutions

%% Preserve matching_type (+ optional test-driver flags) through clear
if ~exist('matching_type', 'var')
    matching_type = 1;  % 0 = Leontief, 1 = Cobb-Douglas (default)
end
% Optional flags a test driver may set (see run_lfx_config.m). Defaults
% reproduce the original behavior exactly.
if ~exist('lfx_printit','var'),    lfx_printit    = 1;  end
if ~exist('lfx_foldername','var'), lfx_foldername = ''; end  % '' = Overleaf default
if ~exist('lfx_tag','var'),        lfx_tag        = ''; end
if ~exist('lfx_target_ebp','var'), lfx_target_ebp = 200; end % SS bond-premium target (bps)
save('temp_matching_type.mat', 'matching_type', 'lfx_printit', ...
     'lfx_foldername', 'lfx_tag', 'lfx_target_ebp');

clear; close all;

%% Load matching_type
load('temp_matching_type.mat');
delete('temp_matching_type.mat');

if matching_type == 0
    mt_suffix = '_l';
else
    mt_suffix = '_cd';
end

addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

foldername = overleaf_dir();
% Test-driver overrides: redirect ALL file output away from Overleaf.
% (compute_moments writes .tex tables to foldername UNGUARDED by printit,
%  so redirecting foldername is the only safe way to run tests.)
if ~isempty(lfx_foldername)
    foldername = lfx_foldername;
    if ~endsWith(foldername, filesep), foldername = [foldername filesep]; end
end
if ~exist(foldername, 'dir'); mkdir(foldername); end

printit = lfx_printit;
plotit= 0;

tic;

% Plot Preferences
T = 5;
LFX_plotprefs;
freq = 12;
rate_scale=1e4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'Z_im_us.mat';
load 'Z_im_eu.mat';

%load 'im_parameter.mat';
%load 'inf_mean.mat';
% LFX_empty_mats;

% Parameters
% pi_eu_ss = 1;
pi_eu_ss = 1;
pi_us_ss = 1;
load dynare_calibration_param.mat;
params;

% Override sigma initial guesses for calibration solver
% (feqm_calibrate solves for sigma — these are just starting points)
sigma_eu = 0.15;
sigma_us = 0.2;

%% Name Scenario
nameplot='calibration_dynare_sigma_us';
xperiment='$\epsilon^{\lambda^{*}}$';

%% Load Eq Equations
LFX_nt_0e_eqs_2;

%% [1] Steady State: solve, test, report
solve_steady_state;

%% [2] Markov Regime Setup
setup_markov;

%% [3] Global Solution
solve_global;

%% [4] Simulation
simulate_model;

%% [5] Paper Figures
plot_paper_figures;

%% [6] Moments and Tables
compute_moments;
