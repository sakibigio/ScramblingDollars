%% Driver: time-varying penalty-spread estimation (audit/Refine item #1 redo).
%   iota_t = (DW - 1M Tbill corridor)_t + kappa, applied to BOTH currencies.
% Same criterion, weights, and data basis as the committed baseline
% (run_estimate_orth with ORTH_MODE=corr, MU_LCR=1: symmetric LCR netting).
% Estimates x = (lambda, kappa, eta) via estimate_params_tvk (generated from
% estimate_params_4d by scratchpad/make_tvk.py).
%
% Outputs: _calibration_override_tvk.mat and _estimation_result_tvk.mat.
% NEVER clobbers the committed baseline _calibration_override.mat.
%
% Headless: /Applications/MATLAB_R2025b.app/bin/matlab -batch run_estimate_tvk

orth_mode = 'corr';
eta_fix   = [];
w_cond    = 0;
% Optional: ORTH_W overrides the orthogonality-penalty weight (default: file's
% w_corr). Result files get a _w<val> suffix so runs never clobber each other.
ow = getenv('ORTH_W');
if ~isempty(ow)
    w_corr_override = str2double(ow);
    tvk_wtag = sprintf('_w%g', w_corr_override);
else
    w_corr_override = [];
    tvk_wtag = '';
end

addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1;          % Cobb-Douglas
pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat');
load('data/calibration.mat');
run('functions/params.m');

% Symmetric LCR netting: SAME mu the live pipeline filter uses
% (main_filter.m: mu_minus_lcr = 1, netted from BOTH currencies).
mu_us = mu_minus_lcr_level(mu_us, LCR_us);
mu_eu = mu_minus_lcr_level(mu_eu, LCR_us);
override_suffix = ['_tvk' tvk_wtag];
fprintf('TVK driver: symmetric LCR netting ON; exp(mu_us) range=[%.3f, %.3f]\n', ...
    min(exp(mu_us)), max(exp(mu_us)));

% DW - 1M T-bill corridor (annual decimal), padded to the model grid:
% hold-first for the lead, hold-last for the tail, linear fill inside.
S = load('data/iota_corridor_tbill.mat');
iota_corr_dec = S.iota_sprd_dec(:);
cov = logical(S.cover_mask(:));
fc = find(cov, 1); lc = find(cov, 1, 'last');
iota_corr_dec(1:fc-1)   = iota_corr_dec(fc);
iota_corr_dec(lc+1:end) = iota_corr_dec(lc);
bad = ~isfinite(iota_corr_dec);
if any(bad)
    gi = find(~bad);
    iota_corr_dec(bad) = interp1(gi, iota_corr_dec(gi), find(bad), 'linear', 'extrap');
end
if length(iota_corr_dec) ~= length(mu_us)
    error('run_estimate_tvk: corridor length %d ~= sample length %d', ...
        length(iota_corr_dec), length(mu_us));
end

estimate_params_tvk

save(['_estimation_result_tvk' tvk_wtag '.mat'], 'lambda_opt', 'eta_opt', 'kappa_opt', ...
     'iota_ss_opt', 'fval', 'exitflag', 'orth_mode');
fprintf('\nTVK done: lambda=%.4f  eta=%.4f  kappa=%.4f  (mean iota_t %.4f)\n', ...
    lambda_opt, eta_opt, kappa_opt, iota_ss_opt);
