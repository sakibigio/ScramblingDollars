%% run_main_filter_cfgB.m
% Config B end-to-end: time-varying iota with the ORIGINAL (LCR-adjusted) mu.
% Runs main_filter -> markov_estimation.jl -> plot_regimes so we get B's own
% Markov regime parameters and regime-conditional moments.
% printit = 0, so the Overleaf figures are untouched.

clear; close all;

L = load('../Filtering_iota_tv/_calibration_tv_origliq.mat');
fprintf('cfg B calibration: lambda=%.4f eta=%.4f premium=%.5f (%.0f bps)\n', ...
    L.lambda_opt, L.eta_opt, L.premium_opt, L.premium_opt*1e4);

lambda_us_override_val = L.lambda_opt;
lambda_eu_override_val = L.lambda_opt;
eta_override_val       = L.eta_opt;
iota_ss_override_val   = L.premium_opt;
if exist('_calibration_override.mat','file') && ~exist('_calibration_override_BBK.mat','file')
    copyfile('_calibration_override.mat','_calibration_override_BBK.mat');
end
save('_calibration_override.mat','lambda_us_override_val','lambda_eu_override_val', ...
     'eta_override_val','iota_ss_override_val');

matching_type = 1;
use_iota_tv   = 1;
use_new_liq   = 0;      % <-- ORIGINAL mu (main_filter's LCR-adjusted series)
iota_premium  = L.premium_opt;
printit       = 0;
run_julia     = 1;      % re-estimate the Markov process on B's sigma
do_regimes    = 1;
cfg_tag       = 'B';

main_filter

% capture figures
outdir = fullfile('output','paper_cfgB');
if ~exist(outdir,'dir'), mkdir(outdir); end
figs = findobj('Type','figure'); [~,o]=sort([figs.Number]); figs=figs(o);
for k=1:numel(figs)
    f=figs(k); nm=get(f,'Name'); if isempty(nm), nm=sprintf('fig%d',f.Number); end
    nm=regexprep(strtrim(nm),'[^\w]+','_');
    try, exportgraphics(f,fullfile(outdir,sprintf('%02d_%s.png',f.Number,nm)),'Resolution',150); catch, end
end
fprintf('[figures] %d -> %s\n', numel(figs), outdir);

if exist('_calibration_override_BBK.mat','file')
    copyfile('_calibration_override_BBK.mat','_calibration_override.mat');
    delete('_calibration_override_BBK.mat');
end
