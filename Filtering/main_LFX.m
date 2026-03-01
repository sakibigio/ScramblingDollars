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

%% Preserve matching_type through clear
if ~exist('matching_type', 'var')
    matching_type = 0;  % 0 = Leontief, 1 = Cobb-Douglas (default)
end
save('temp_matching_type.mat', 'matching_type');

clear; close all;

%% Load matching_type
load('temp_matching_type.mat');
delete('temp_matching_type.mat');

addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

if strcmp(getenv('HOME'),'/Users/sakiclaudia')
    foldername='/Users/sakiclaudia/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/quantfigs/';
else
    foldername='/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/quantfigs/';
end

printit=1;
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
sigma_us = 0.20;

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
