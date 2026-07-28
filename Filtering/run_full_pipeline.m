%% run_full_pipeline.m — one-shot driver for the full paper pipeline
%
% Runs, in order:
%   [1] load_data              xlsx -> data/LFX_data.mat  (ALWAYS first: main_filter
%                              reads the .mat, NOT the workbook — skipping this
%                              step silently reuses stale data)
%   [2] main_filter            TED-inversion filter -> RW_shock.csv, then
%                              markov_estimation.jl (regime MLE) and the
%                              regime/fit figures (F_*.pdf)
%   [3] main_LFX               steady state, Markov setup, global solution,
%                              simulation, paper figures (nlmod_*.pdf) + moment tables
%   [4] compute_data_moments   data-side conditional-moment tables
%                              (Data_*.tex — NOT run by any other stage)
%
% Re-estimation of (lambda, iota, eta) — estimate_params_4d via run_estimate_5d —
% is NOT part of this driver; run it only when the data basis changes, then
% re-run this pipeline (it reads _calibration_override.mat).
%
% Output destination: with printit = 1 the figure/table scripts write STRAIGHT
% INTO THE OVERLEAF quantfigs/ FOLDER (path hardcoded by username in
% plotting/plot_regimes.m, plot_baseline.m, and main_LFX.m). Set printit = 0
% below for a dry run that touches nothing outside this folder.

clear; close all;
set(groot, 'DefaultFigureVisible', 'off');   % headless-friendly; comment out to watch

printit_all = 1;    % 1 = export figures/tables (Overleaf!), 0 = dry run
save('temp_pipeline_flags.mat', 'printit_all');  % survives the clears between stages

%% [1] Data
load_data;

%% [2] Filter + Markov + fit figures
clear; close all;
load('temp_pipeline_flags.mat', 'printit_all');
matching_type    = 1;            % 1 = Cobb-Douglas (paper), 0 = Leontief
printit          = printit_all;
run_julia        = 1;
do_regimes       = 1;
do_plot_baseline = 1;
main_filter

%% [3] Model stage
clear; close all;
load('temp_pipeline_flags.mat', 'printit_all');
matching_type = 1;
lfx_printit   = printit_all;
main_LFX

%% [4] Data-side moment tables
clear; close all;
compute_data_moments

if exist('temp_pipeline_flags.mat', 'file'), delete('temp_pipeline_flags.mat'); end
fprintf('\nFULL PIPELINE DONE\n');
