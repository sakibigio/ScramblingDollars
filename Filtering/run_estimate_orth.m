%% Driver: parameter estimation with a selectable identification penalty.
%
% Runs estimate_params_4d.m with orth_mode taken from the ORTH_MODE
% environment variable (default 'trend'). Modes:
%   'corr'  : corr(log sigma_us, mu_us)^2   (original orthogonality condition)
%   'trend' : corr(log sigma_us, t)^2       (no secular trend in filtered stress)
%   'kpss'  : KPSS-type partial-sum statistic (level stationarity)
%
% Non-'corr' modes snapshot to _calibration_override_<mode>.mat, so the
% committed baseline _calibration_override.mat is untouched.
%
% Usage (headless):
%   ORTH_MODE=trend /Applications/MATLAB_R2025b.app/bin/matlab -batch run_estimate_orth

em = getenv('ORTH_MODE');
if isempty(em), em = 'trend'; end
orth_mode = em;
% Optional: ETA_FIX pins eta (2-D search); ORTH_W overrides the penalty weight.
ef = getenv('ETA_FIX');
if isempty(ef), eta_fix = []; else, eta_fix = str2double(ef); end
ow = getenv('ORTH_W');
if isempty(ow), w_corr_override = []; else, w_corr_override = str2double(ow); end
cw = getenv('COND_W');
if isempty(cw), w_cond = 0; else, w_cond = str2double(cw); end
fprintf('run_estimate_orth: orth_mode = %s, eta_fix = %s, w_override = %s\n', ...
    orth_mode, mat2str(eta_fix), mat2str(w_corr_override));

addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1;          % Cobb-Douglas
pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat');
load('data/calibration.mat');
run('functions/params.m');

% MU_LCR=1: estimate on the SAME mu the live pipeline filter uses --
% excess liquidity above the LCR requirement (main_filter.m: mu_minus_lcr=1).
% US only; there is no EU LCR series (matches main_filter).
ml = getenv('MU_LCR');
if ~isempty(ml) && str2double(ml) == 1
    mu_us = mu_minus_lcr_level(mu_us, LCR_us);
    override_suffix = '_lcr';
    fprintf('MU_LCR ON: mu_us = log(exp(mu) - LCR/100), exp(mu_us) range=[%.3f, %.3f]\n', ...
        min(exp(mu_us)), max(exp(mu_us)));
else
    override_suffix = '';
end

estimate_params_4d

if isempty(eta_fix)
    res_file = sprintf('_estimation_result_%s%s.mat', orth_mode, override_suffix);
else
    res_file = sprintf('_estimation_result_eta%02.0f%s.mat', 100 * eta_fix, override_suffix);
end
save(res_file, 'lambda_opt', 'eta_opt', 'iota_ss_opt', 'fval', 'exitflag', 'orth_mode');
