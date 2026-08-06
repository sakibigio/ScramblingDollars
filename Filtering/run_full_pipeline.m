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

% Every stage below writes relative paths (temp_*.mat, RW_shock.csv) and calls
% addpath('data'), so the whole chain assumes pwd == this folder. Anchor it to
% the script's own location instead of the caller's cwd -- otherwise running
% with a stale or read-only current folder dies at the first save with a bare
% "No such file or directory" that says nothing about the real cause.
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = which('run_full_pipeline');   % pasted into the console
end
if ~isempty(this_file)
    cd(fileparts(this_file));
    % Also put this folder on the MATLAB path, not just in pwd. Stages call
    % run('plotting/plot_baseline.m'), and run() cds into plotting/ for the
    % duration -- at which point helpers living here (overleaf_dir, acf_local)
    % stop resolving, because only pwd and the path are searched. addpath
    % survives both the run() shifts and the clear calls between stages.
    addpath(fileparts(this_file));
end
clear this_file

% Prove the current folder is actually writable before the pipeline commits to
% it, and say exactly what is wrong if it is not.
[fid, msg] = fopen('.pipeline_write_test', 'w');
if fid < 0
    error(['run_full_pipeline: cannot write to the current folder.\n' ...
           '  pwd  : %s\n' ...
           '  error: %s\n' ...
           'Fix: cd to the Filtering folder of your clone and re-run. If the ' ...
           'folder was moved, deleted or re-cloned while MATLAB was open, run ' ...
           '"cd(''~'')" then "rehash toolboxcache" first.'], pwd, msg);
end
fclose(fid);
delete('.pipeline_write_test');
fprintf('run_full_pipeline: running in %s\n', pwd);

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

%% [3b] In-sample CIP row (Table 4's in-sample line). Audit High-4 fix
% (2026-08-06): this row's only producer was previously a manual out-of-band
% call, so the documented driver never regenerated Mod_CIP_Moments_insample.tex.
clear; close all;
load('temp_pipeline_flags.mat', 'printit_all');
if printit_all == 1
    make_insample_cip_row
end

%% [4] Data-side moment tables
clear; close all;
compute_data_moments

if exist('temp_pipeline_flags.mat', 'file'), delete('temp_pipeline_flags.mat'); end
fprintf('\nFULL PIPELINE DONE\n');
