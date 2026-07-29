%% run_main_filter_iota_tv.m
% Run the FULL filtering pipeline (main_filter.m + plot_baseline.m) using the
% TIME-VARYING iota calibration:   iota_t = (DW - Tbill)_t + premium.
%
% Calibration comes from Filtering_iota_tv/_calibration_iota_tv.mat
% (lambda, eta, premium) and the monthly corridor in data/iota_corridor_monthly.mat.
%
% SAFE BY DESIGN:
%   printit    = 0  -> does NOT write into the Overleaf quantfigs/ folder,
%                      so the paper's existing figures are untouched.
%   run_julia  = 0  -> skips the Julia Markov step.
%   do_regimes = 0  -> skips regime plots (they consume the Julia output).
% All figures are instead captured to output/paper_iota_tv/.

clear; close all;

%% --- Pull the time-varying calibration ---
L = load('../Filtering_iota_tv/_calibration_iota_tv.mat');
fprintf('TV calibration: lambda=%.4f  eta=%.4f  premium=%.5f (%.1f bps)\n', ...
    L.lambda_opt, L.eta_opt, L.premium_opt, L.premium_opt*1e4);

% Install it as the params.m override (backing up whatever is there now).
lambda_us_override_val = L.lambda_opt;
lambda_eu_override_val = L.lambda_opt;
eta_override_val       = L.eta_opt;
iota_ss_override_val   = L.premium_opt;   % level placeholder; tv path set in main_filter
if exist('_calibration_override.mat','file')
    copyfile('_calibration_override.mat','_calibration_override_PREV.mat');
end
save('_calibration_override.mat', 'lambda_us_override_val', ...
     'lambda_eu_override_val', 'eta_override_val', 'iota_ss_override_val');

%% --- Flags consumed by main_filter (survive its internal `clear`) ---
matching_type = 1;          % Cobb-Douglas
use_iota_tv   = 1;          % <-- time-varying iota ON
use_new_liq   = 1;          % <-- must match the estimation (log LiqRatio H_18)
iota_premium  = L.premium_opt;
printit       = 0;          % do NOT overwrite the paper's Overleaf figures
run_julia     = 1;          % run markov_estimation.jl on the NEW sigma
do_regimes    = 1;          % produce the paper figures (plot_regimes.m)

% Back up the existing Markov output before Julia regenerates it.
for fbk = {'data/MS_sigma_us_prob.csv','data/MS_sigma_us_params.csv'}
    if exist(fbk{1},'file') && ~exist([fbk{1} '.PREV'],'file')
        copyfile(fbk{1}, [fbk{1} '.PREV']);
    end
end

%% --- Run the pipeline ---
main_filter

%% --- Capture every figure locally (workspace was cleared; use literals) ---
outdir = fullfile('output','paper_iota_tv');
if ~exist(outdir,'dir'), mkdir(outdir); end
figs = findobj('Type','figure');
[~, ord] = sort([figs.Number]);
figs = figs(ord);
for k = 1:numel(figs)
    f  = figs(k);
    nm = get(f,'Name');
    if isempty(nm), nm = sprintf('fig%d', f.Number); end
    nm = regexprep(strtrim(nm), '[^\w]+', '_');
    fn = fullfile(outdir, sprintf('%02d_%s.png', f.Number, nm));
    try
        exportgraphics(f, fn, 'Resolution', 150);
    catch ME
        fprintf('  [skip fig %d: %s]\n', f.Number, ME.message);
    end
end
fprintf('\n[figures] captured %d figures to %s\n', numel(figs), outdir);

%% --- Restore the previous calibration override ---
if exist('_calibration_override_PREV.mat','file')
    copyfile('_calibration_override_PREV.mat','_calibration_override.mat');
    delete('_calibration_override_PREV.mat');
    fprintf('[restore] previous _calibration_override.mat restored\n');
end
