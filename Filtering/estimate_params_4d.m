%% Estimate (lambda_us, lambda_eu, iota_ss, eta, varrho) by mixing
%   price-variable time-series targets (squared distances) with two scalar
%   volume targets (average DW, average FF).
%
%   Price targets        (squared distance over valid sample, normalized by 1/var(data)):
%     P1: CIP deviation         (model CIP_t            vs cip)            [bps]
%     P2: US Bond Premium       (BP_us_t demeaned       vs Rb_Rm demeaned) [bps]
%     P3: EU Bond Premium       (BP_eu_t demeaned       vs Rb_Rm_eu demeaned) [bps]
%
%   Volume targets       (squared log-mean distance):
%     V1: mean DW (US)          ( log mean DW_us_t      vs log mean DW_n )
%     V2: mean FF (US)          ( log mean FF_us_t      vs log mean FF_n )
%
%   Objective = w_cip   * mean( r_cip^2  ) / var(cip_d * abs_scale)
%             + w_bpus  * mean( r_bpus^2 ) / var(bpus_dm_d * abs_scale)
%             + w_bpeu  * mean( r_bpeu^2 ) / var(bpeu_dm_d * abs_scale)
%             + w_meanDW * ( log mean DW_m - log mean DW_d )^2
%             + w_meanFF * ( log mean FF_m - log mean FF_d )^2
%
%   Each price term is unit-free after the 1/var normalization, and the
%   volume terms are already small log-deviations, so the default
%   user weights (all = 1) yield a roughly balanced objective.
%
%   Filter pipeline mirrors main_filter.m (TED-inversion in both US and EU,
%   Cobb-Douglas backstop with find_sigma_min, baseline-corrected
%   mu = exp(mu) - mubar). lambda_us and lambda_eu are estimated separately.
%
%   Run after main_filter.m setup, i.e.:
%       load('data/LFX_data.mat'); load('data/calibration.mat');
%       pi_us_ss = 1; pi_eu_ss = 1;
%       matching_type = 1;
%       run('functions/params.m');

fprintf('\n====== Parameter Estimation: (lambda_us, lambda_eu, iota_ss, eta, varrho) ======\n');
fprintf('  Price targets:   CIP, US BP, EU BP  (time series, squared distance)\n');
fprintf('  Volume targets:  mean DW, mean FF   (scalar log-mean moments)\n');

%% ====================== USER SETTINGS ======================
% Target weights (each multiplies a unit-free term; set 0 to drop a target).
w_cip      = 0.0;
w_bpus     = 1.0;
w_bpeu     = 1.0;
w_meanDW   = 0.0;   % DW FLOW level (mismatched: model flow vs data stock)
w_meanFF   = 1.0;
w_ratioDF  = 0.0;
w_meanDWS  = 1.0;   % DW STOCK level (model stock vs data stock -- like-with-like)

% Orthogonality penalty:  corr(log sigma_us_t, mu_us_t)^2.
% At the η=0.5 baseline, |corr|≈0.90 in levels (very strong geometric coupling).
% w_corr=10 produced corr=0.88 at optimum (almost no decoupling).
% w_corr=100 forces the optimizer to accept worse moment fit if necessary.
% A corr=0.5 fit then contributes 25 to the objective -- larger than any
% reasonable moment-fit cost.
w_corr     = 100.0;

% Identification-penalty mode (all unit-free, so w_corr is comparable across modes):
%   'corr'  : corr(log sigma_us, mu_us)^2      -- original orthogonality condition
%   'trend' : corr(log sigma_us, t)^2          -- R^2 of a linear time trend in
%             log sigma; imposes that the filtered stress process has no secular
%             trend (stationarity of the structural shock) without referencing mu.
%   'kpss'  : KPSS-type partial-sum statistic  -- sum(S_t^2)/(T^2*var), S_t the
%             partial sums of demeaned log sigma; penalizes ANY low-frequency
%             drift (level nonstationarity), not just a linear trend.
% A driver script may predefine orth_mode before running this file.
if ~exist('orth_mode', 'var'), orth_mode = 'corr'; end
% A driver may also override the penalty weight (e.g. 0 to drop the penalty).
if exist('w_corr_override', 'var') && ~isempty(w_corr_override)
    w_corr = w_corr_override;
end
fprintf('Identification penalty mode: %s (w=%.1f)\n', orth_mode, w_corr);

% Fixed-eta mode: a driver may predefine eta_fix to pin eta at a constant
% (e.g. 0.5 = symmetric Nash bargaining) and estimate only (lambda, iota).
if ~exist('eta_fix', 'var'), eta_fix = []; end
if ~isempty(eta_fix)
    fprintf('Fixed-eta mode: eta = %.4f (2-D search over lambda, iota)\n', eta_fix);
end

% Conditional-amplification target (default OFF). When w_cond > 0, adds
%   w_cond * ((mean(CIP_m|scr) - mean(CIP_m|nor)) - cond_diff_data)^2 / cond_diff_data^2
% where months are classified scrambling/normal by P(scr) > 0.5 from the
% committed baseline regime probabilities (data/MS_sigma_us_prob_cbase.csv,
% rows start 2001-02 = series index 2). This feeds the paper's headline
% UNTARGETED moment into the estimation transparently -- it identifies the
% amplification role of lambda that the unconditional moments leave flat.
if ~exist('w_cond', 'var'), w_cond = 0; end

% Stock construction:  DWS_t = DW_t + (1 - delta_DWS_estim) * DWS_{t-1}
delta_DWS_estim = 0.2;   % monthly depreciation rate (must match main_filter.m)

% Parameter bounds (closed intervals; penalty pushes search back inside)
b_lambda    = [0.5, 12.0];   % common lambda for US and EU (paper simplification)
b_iota      = [0.005, 0.12];
b_eta       = [0.05, 0.95];  % bargaining-power bounds (avoid degenerate extremes)

% Asymptote-feasibility constraint: at sigma -> infty the model TED
% approaches (1-eta)*iota_loc, so we require
%   (1-eta)*iota_ss/pi_ss * 1e4  >  ted_asymp_margin * max(TED_bps)
ted_asymp_margin = 1.05;

% Initial guess (3-D: lambda, iota_ss, eta).  lambda_us = lambda_eu (paper
% simplification); varrho=0.  Start eta at the previous fixed value (0.4)
% so the search is anchored near the recently-estimated point.
lambda_init = 0.5 * (lambda_us + lambda_eu);
eta_init    = 0.4;
x0 = [ ...
    min(max(lambda_init, b_lambda(1) + 1e-3), b_lambda(2) - 1e-3), ...
    min(max(iota_ss,     b_iota(1)   + 1e-4), b_iota(2)   - 1e-4), ...
    min(max(eta_init,    b_eta(1)    + 1e-3), b_eta(2)    - 1e-3)];
if ~isempty(eta_fix)
    x0 = x0(1:2);   % 2-D search; eta pinned at eta_fix throughout
end
% ===========================================================

%% Baseline correction (recompute if not already set up)
T_obs = length(mu_us);
% Match main_filter.m: use_mu_baseline = 0 means mubar = 0 (no correction).
% Set to 1 to subtract the annual regulatory baseline.
use_mu_baseline_estim = 0;
if use_mu_baseline_estim == 1
    dates_obs = datenum(2001, 1:T_obs, 1);
    [mubar_us_t, mubar_eu_t] = liquidity_baseline(dates_obs);
    fprintf('Baseline correction ON: mubar_us=[%.3f,%.3f], mubar_eu=[%.3f,%.3f]\n', ...
        min(mubar_us_t), max(mubar_us_t), min(mubar_eu_t), max(mubar_eu_t));
else
    mubar_us_t = zeros(T_obs, 1);
    mubar_eu_t = zeros(T_obs, 1);
    fprintf('Baseline correction OFF: mubar = 0 (matches main_filter.m default)\n');
end

%% Constants & exogenous time series
abs_sc    = 12e4;            % freq * 1e4 (bps scaling)
Rm_us_vec = exp(im_us);
Rm_eu_vec = exp(im_eu);

%% Valid-sample indices (mirror plot_regimes truncation rules)
NaN_cip   = (cip == 0)      | isnan(cip);
NaN_bpus  = (Rb_Rm == 0)    | isnan(Rb_Rm);
NaN_bpeu  = (Rb_Rm_eu == 0) | isnan(Rb_Rm_eu);
dw_valid  = ~isnan(DW_n) & DW_n > 0;
ff_valid  = ~isnan(FF_n) & FF_n > 0;

idx_cip  = find(~NaN_cip);
idx_bpus = find(~NaN_bpus);
idx_bpeu = find(~NaN_bpeu);
idx_dw   = find(dw_valid);
idx_ff   = find(ff_valid);

fprintf('Valid periods: CIP=%d, BPus=%d, BPeu=%d, DW=%d, FF=%d\n', ...
    numel(idx_cip), numel(idx_bpus), numel(idx_bpeu), numel(idx_dw), numel(idx_ff));

%% Pre-compute data targets and normalizations
cip_d           = cip(idx_cip);
bpus_d_demeaned = Rb_Rm(idx_bpus)    - mean(Rb_Rm(idx_bpus));
bpeu_d_demeaned = Rb_Rm_eu(idx_bpeu) - mean(Rb_Rm_eu(idx_bpeu));
logmean_DW_d    = log(mean(DW_n(idx_dw)));
logmean_FF_d    = log(mean(FF_n(idx_ff)));

% Data log-ratio DW/FF: use overlapping sample for both
idx_both     = intersect(idx_dw, idx_ff);
logratio_DF_d = log(mean(DW_n(idx_both))) - log(mean(FF_n(idx_both)));

v_cip  = max(var(cip_d           * abs_sc), 1e-12);
v_bpus = max(var(bpus_d_demeaned * abs_sc), 1e-12);
v_bpeu = max(var(bpeu_d_demeaned * abs_sc), 1e-12);

fprintf('Data std: CIP=%.2f bps, BPus(dm)=%.2f bps, BPeu(dm)=%.2f bps\n', ...
    sqrt(v_cip), sqrt(v_bpus), sqrt(v_bpeu));
fprintf('Data:     log mean DW = %.4f, log mean FF = %.4f, log(DW/FF) = %.4f\n', ...
    logmean_DW_d, logmean_FF_d, logratio_DF_d);

%% Conditional-target setup (only when w_cond > 0)
idx_scr = []; idx_nor = []; cond_diff_data_bps = NaN;
if w_cond > 0
    prob_tbl = readtable('data/MS_sigma_us_prob_cbase.csv');
    % prob rows start 2001-02 -> series index offset +1
    scr_mask_series = false(T_obs, 1);
    scr_mask_series(1 + find(prob_tbl.prob_scr > 0.5)) = true;
    idx_scr = intersect(idx_cip, find(scr_mask_series));
    idx_nor = intersect(idx_cip, find(~scr_mask_series));
    cond_diff_data_bps = (mean(cip(idx_scr)) - mean(cip(idx_nor))) * abs_sc;
    fprintf('Conditional CIP target ON (w_cond=%.1f): %d scr / %d nor months, data diff = %.1f bps\n', ...
        w_cond, numel(idx_scr), numel(idx_nor), cond_diff_data_bps);
end

%% TED-asymptote feasibility (varrho = 0 case)
%   Model TED at sigma -> infty  =  (1-eta) * iota_ss / pi_ss     [per-period frac]
%                                =  (1-eta) * iota_ss * 1e4 / pi_ss   [bps, since abs_sc=12e4]
%   Require this to exceed margin * max(TED_data_bps) on both regions.
max_ted_us_bps = max(TED_s_us_t) * abs_sc;
max_ted_eu_bps = max(TED_s_eu_t) * abs_sc;
max_ted_bps    = max(max_ted_us_bps, max_ted_eu_bps);
fprintf('Max TED data: US=%.2f bps, EU=%.2f bps  (binding = %.2f bps)\n', ...
    max_ted_us_bps, max_ted_eu_bps, max_ted_bps);
fprintf('Asymptote constraint:  (1-eta)*iota_ss*1e4 > %.4f * %.2f = %.2f bps  (varrho=0)\n', ...
    ted_asymp_margin, max_ted_bps, ted_asymp_margin * max_ted_bps);

% Quick check that the initial guess satisfies the constraint.
% Now eta is free too -- bump iota if needed to satisfy at the starting eta.
if isempty(eta_fix), eta0 = x0(3); else, eta0 = eta_fix; end
asy0 = (1 - eta0) * x0(2) * 1e4 / pi_us_ss;
if asy0 < ted_asymp_margin * max_ted_bps
    iota_min_feasible = ted_asymp_margin * max_ted_bps / ((1 - eta0) * 1e4 / pi_us_ss);
    warning(['Initial guess violates asymptote constraint: ' ...
        '(1-eta)*iota*1e4=%.2f bps < %.2f bps. Bumping x0(iota) to %.4f.'], ...
        asy0, ted_asymp_margin * max_ted_bps, iota_min_feasible * 1.02);
    x0(2) = min(b_iota(2) - 1e-4, iota_min_feasible * 1.02);
    fprintf('Adjusted x0(iota) -> %.4f\n', x0(2));
end

%% Anonymous objective (3D: lambda, iota, eta)
obj_fun = @(x) estimate_obj_2d(x, ...
    mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
    TED_s_us_t, TED_s_eu_t, Rm_us_vec, Rm_eu_vec, ...
    Rb_Rm, Rb_Rm_eu, cip, ...
    pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us, sigma_eu, ...
    idx_cip, idx_bpus, idx_bpeu, idx_dw, idx_ff, idx_both, ...
    v_cip, v_bpus, v_bpeu, logmean_DW_d, logmean_FF_d, logratio_DF_d, ...
    w_cip, w_bpus, w_bpeu, w_meanDW, w_meanFF, w_ratioDF, w_meanDWS, w_corr, ...
    delta_DWS_estim, abs_sc, ...
    b_lambda, b_iota, b_eta, ...
    max_ted_bps, ted_asymp_margin, orth_mode, eta_fix, ...
    w_cond, idx_scr, idx_nor, cond_diff_data_bps);

%% Run fminsearch
fprintf('Starting fminsearch from x0:\n');
if isempty(eta_fix)
    fprintf('  lambda=%.3f, iota=%.4f, eta=%.3f  (3-D search; lambda_us = lambda_eu)\n', ...
        x0(1), x0(2), x0(3));
else
    fprintf('  lambda=%.3f, iota=%.4f  (2-D search; eta fixed at %.3f)\n', ...
        x0(1), x0(2), eta_fix);
end
fprintf('Weights: w_cip=%.2f w_bpus=%.2f w_bpeu=%.2f w_meanDW=%.2f w_meanFF=%.2f w_ratioDF=%.2f w_meanDWS=%.2f w_corr=%.2f (delta=%.2f)\n', ...
    w_cip, w_bpus, w_bpeu, w_meanDW, w_meanFF, w_ratioDF, w_meanDWS, w_corr, delta_DWS_estim);

options = optimset('Display', 'iter', 'MaxFunEvals', 12000, 'MaxIter', 6000, ...
                   'TolX', 1e-7, 'TolFun', 1e-10);

%% Multi-start: screen several initial points, then polish from the best.
% Each starting point must satisfy the asymptote constraint
%   (1-eta) * iota * 1e4 > 1.05 * max_TED_bps
% so we choose (iota, eta) pairs whose asymptote sits ~10-20% above max TED.
x0_list = {
    x0(:)',                       % the script's default
    [0.70, 0.060, 0.35],          % low lambda, low eta
    [0.70, 0.090, 0.55],          % low lambda, high eta
    [1.50, 0.070, 0.30],          % mid lambda, low eta
    [1.50, 0.060, 0.45],          % mid lambda, mid eta (≈ Nash)
    [1.50, 0.090, 0.60],          % mid lambda, high eta
    [3.00, 0.060, 0.40],          % high lambda
    [5.00, 0.070, 0.50],          % very high lambda
};
if ~isempty(eta_fix)
    % 2-D mode: drop the eta component of every starting point.
    x0_list = cellfun(@(v) v(1:2), x0_list, 'UniformOutput', false);
end

options_screen = optimset('Display', 'off', 'MaxFunEvals', 1500, 'MaxIter', 400, ...
                          'TolX', 1e-5, 'TolFun', 1e-6);

best_fval = Inf;
best_x    = x0_list{1};
fprintf('\n====== Multi-start screening (%d initial points) ======\n', numel(x0_list));
fprintf('  #   x0 [lambda, iota, eta]        ->  fval        exit  \n');
for ii = 1:numel(x0_list)
    x0_try = x0_list{ii};
    try
        [x_try, fval_try, ef_try] = fminsearch(obj_fun, x0_try, options_screen);
        fprintf('  %d   [%s]  ->  %.4f      %d\n', ...
            ii, sprintf('%.4f ', x0_try), fval_try, ef_try);
        if fval_try < best_fval
            best_fval = fval_try;
            best_x    = x_try;
        end
    catch ME
        fprintf('  %d   [%s]  FAILED: %s\n', ...
            ii, sprintf('%.4f ', x0_try), ME.message);
    end
end
fprintf('Best screening fval = %.4f at x = [%s]\n', ...
    best_fval, sprintf('%.5f ', best_x));
fprintf('Polishing from best with full-precision options...\n');
[x_opt, fval, exitflag] = fminsearch(obj_fun, best_x, options);

%% Results
lambda_opt    = x_opt(1);
lambda_us_opt = lambda_opt;     % common lambda (paper simplification)
lambda_eu_opt = lambda_opt;
iota_ss_opt   = x_opt(2);
if isempty(eta_fix), eta_opt = x_opt(3); else, eta_opt = eta_fix; end
varrho_opt    = 0;
asy_opt_bps   = (1 - eta_opt) * iota_ss_opt * 1e4 / pi_us_ss;

fprintf('\n====== Estimation Results ======\n');
fprintf('  lambda    = %.6f   (lambda_us = lambda_eu)\n', lambda_opt);
fprintf('  iota_ss   = %.6f\n', iota_ss_opt);
if isempty(eta_fix)
    fprintf('  eta       = %.6f   (estimated; init at %.3f)\n', eta_opt, eta_init);
else
    fprintf('  eta       = %.6f   (FIXED)\n', eta_opt);
end
fprintf('  varrho    = %.6f   (removed; hardcoded to 0)\n', varrho_opt);
fprintf('  Asymptote (1-eta)*iota*1e4 = %.2f bps  (vs max TED %.2f bps)\n', ...
    asy_opt_bps, max_ted_bps);
fprintf('  Objective = %.10e\n', fval);
fprintf('  Exit flag = %d\n', exitflag);

%% Implied fit at the optimum (request diagnostic outputs too)
[sigma_us_opt, sigma_eu_opt, BP_us_opt, BP_eu_opt, CIP_opt, DW_us_opt, FF_us_opt, ...
    flag_us_t, flag_eu_t, sig_min_us_t, sig_min_eu_t, ted_err_us_t, ted_err_eu_t] = ...
    run_filter_5d(lambda_us_opt, lambda_eu_opt, iota_ss_opt, eta_opt, varrho_opt, ...
        mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
        TED_s_us_t, TED_s_eu_t, ...
        pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us, sigma_eu, ...
        Rm_us_vec, Rm_eu_vec);

% Price RMSEs in bps
fit_cip   = sqrt(mean(((CIP_opt(idx_cip)   - cip(idx_cip)) * abs_sc).^2));
bpus_m_dm = BP_us_opt(idx_bpus) - mean(BP_us_opt(idx_bpus));
bpeu_m_dm = BP_eu_opt(idx_bpeu) - mean(BP_eu_opt(idx_bpeu));
fit_bpus  = sqrt(mean(((bpus_m_dm - bpus_d_demeaned) * abs_sc).^2));
fit_bpeu  = sqrt(mean(((bpeu_m_dm - bpeu_d_demeaned) * abs_sc).^2));

% Volume moment fits (log scale)
logmean_DW_m = log(mean(DW_us_opt(idx_dw)));
logmean_FF_m = log(mean(FF_us_opt(idx_ff)));

if w_cond > 0
    diff_model_opt = (mean(CIP_opt(idx_scr)) - mean(CIP_opt(idx_nor))) * abs_sc;
    fprintf('\nConditional CIP diff (scr - nor): model = %.1f bps, data = %.1f bps\n', ...
        diff_model_opt, cond_diff_data_bps);
end

fprintf('\n--- Price fit (RMSE, model vs data) ---\n');
fprintf('  CIP        : %8.2f bps    (data std = %8.2f bps)\n', fit_cip,  sqrt(v_cip));
fprintf('  BP_us (dm) : %8.2f bps    (data std = %8.2f bps)\n', fit_bpus, sqrt(v_bpus));
fprintf('  BP_eu (dm) : %8.2f bps    (data std = %8.2f bps)\n', fit_bpeu, sqrt(v_bpeu));

% Build DW STOCK at the optimum for reporting
Tf_opt = length(DW_us_opt);
DWS_us_opt = zeros(Tf_opt, 1);
for tt = 1:Tf_opt
    if tt == 1, prev = 0; else, prev = DWS_us_opt(tt-1); end
    DWS_us_opt(tt) = DW_us_opt(tt) + (1 - delta_DWS_estim) * prev;
end
logmean_DWS_m = log(mean(DWS_us_opt(idx_dw)));
logmean_DWS_d = logmean_DW_d;   % data DW_n is already a stock

logratio_DF_m = log(mean(DW_us_opt(idx_both))) - log(mean(FF_us_opt(idx_both)));
fprintf('\n--- Volume moments (logs) ---\n');
fprintf('                 model       data       diff\n');
fprintf('  mean(DW)   :  %8.4f   %8.4f   %+8.4f    (flow; mismatched units)\n', ...
    logmean_DW_m,  logmean_DW_d,  logmean_DW_m - logmean_DW_d);
fprintf('  mean(DWS)  :  %8.4f   %8.4f   %+8.4f    (stock; delta=%.2f)\n', ...
    logmean_DWS_m, logmean_DWS_d, logmean_DWS_m - logmean_DWS_d, delta_DWS_estim);
fprintf('  mean(FF)   :  %8.4f   %8.4f   %+8.4f\n', logmean_FF_m,  logmean_FF_d,  logmean_FF_m - logmean_FF_d);
fprintf('  log(DW/FF) :  %8.4f   %8.4f   %+8.4f\n', logratio_DF_m, logratio_DF_d, logratio_DF_m - logratio_DF_d);
fprintf('\n--- Levels ---\n');
fprintf('  mean(DW) :  %8.4f   %8.4f    (model flow)\n', mean(DW_us_opt(idx_dw)),  mean(DW_n(idx_dw)));
fprintf('  mean(DWS):  %8.4f   %8.4f    (model stock)\n', mean(DWS_us_opt(idx_dw)), mean(DW_n(idx_dw)));
fprintf('  mean(FF) :  %8.4f   %8.4f\n', mean(FF_us_opt(idx_ff)), mean(FF_n(idx_ff)));
fprintf('================================================\n');

%% ----- Sigma-solver diagnostics at the optimum -----
fprintf('\n--- Solver diagnostics at the optimum ---\n');

T_f = length(sigma_us_opt);

% NaN / failure counts
n_nan_us  = sum(~isfinite(sigma_us_opt));
n_nan_eu  = sum(~isfinite(sigma_eu_opt));
n_fail_us = sum(flag_us_t == -1);
n_fail_eu = sum(flag_eu_t == -1);

% Solver-flag breakdown
n_fsolve_us  = sum(flag_us_t == 1);
n_corner_us  = sum(flag_us_t == 2);    % CD-only: within 1e-3 of sig_min
n_fminbnd_us = sum(flag_us_t == 10);
n_fsolve_eu  = sum(flag_eu_t == 1);
n_corner_eu  = sum(flag_eu_t == 2);
n_fminbnd_eu = sum(flag_eu_t == 10);

fprintf('Total periods:  %d\n', T_f);
fprintf('US solver:   fsolve=%d (%.1f%%)  fminbnd=%d (%.1f%%)  near-corner=%d  fail/NaN=%d\n', ...
    n_fsolve_us, n_fsolve_us/T_f*100, n_fminbnd_us, n_fminbnd_us/T_f*100, n_corner_us, n_fail_us + n_nan_us);
fprintf('EU solver:   fsolve=%d (%.1f%%)  fminbnd=%d (%.1f%%)  near-corner=%d  fail/NaN=%d\n', ...
    n_fsolve_eu, n_fsolve_eu/T_f*100, n_fminbnd_eu, n_fminbnd_eu/T_f*100, n_corner_eu, n_fail_eu + n_nan_eu);

% TED inversion residuals (bps)
resid_us = abs(ted_err_us_t(isfinite(ted_err_us_t)));
resid_eu = abs(ted_err_eu_t(isfinite(ted_err_eu_t)));
fprintf('TED resid US (bps): max=%.3f  mean=%.3f  median=%.3f  pct<0.1=%.1f%%\n', ...
    max(resid_us), mean(resid_us), median(resid_us), sum(resid_us<0.1)/numel(resid_us)*100);
fprintf('TED resid EU (bps): max=%.3f  mean=%.3f  median=%.3f  pct<0.1=%.1f%%\n', ...
    max(resid_eu), mean(resid_eu), median(resid_eu), sum(resid_eu<0.1)/numel(resid_eu)*100);

% Distance to sig_min (CD only) — flag suspiciously close periods
if matching_type == 1
    gap_us = sigma_us_opt - sig_min_us_t;
    gap_eu = sigma_eu_opt - sig_min_eu_t;
    n_tight_us = sum(gap_us < 1e-3);
    n_tight_eu = sum(gap_eu < 1e-3);
    fprintf('sigma_us - sig_min: min=%.5f  median=%.4f   periods <1e-3: %d\n', ...
        min(gap_us), median(gap_us), n_tight_us);
    fprintf('sigma_eu - sig_min: min=%.5f  median=%.4f   periods <1e-3: %d\n', ...
        min(gap_eu), median(gap_eu), n_tight_eu);
end

fprintf('sigma_us: range [%.4f, %.4f]  mean=%.4f  std=%.4f\n', ...
    min(sigma_us_opt), max(sigma_us_opt), mean(sigma_us_opt), std(sigma_us_opt));
fprintf('sigma_eu: range [%.4f, %.4f]  mean=%.4f  std=%.4f\n', ...
    min(sigma_eu_opt), max(sigma_eu_opt), mean(sigma_eu_opt), std(sigma_eu_opt));

% Identification-penalty diagnostics at the optimum (all three measures,
% regardless of which one was penalized, for cross-mode comparison)
valid_corr = isfinite(sigma_us_opt) & sigma_us_opt > 0 & isfinite(mu_us);
if sum(valid_corr) > 10
    lsig_opt = log(sigma_us_opt(valid_corr));
    lsig_opt = lsig_opt(:);
    corr_lvl = corr(lsig_opt, mu_us(valid_corr));
    dlsig    = diff(lsig_opt);
    dmu      = diff(mu_us(valid_corr));
    corr_dff = corr(dlsig, dmu(:));
    tvec_opt  = (1:numel(lsig_opt))';
    corr_trnd = corr(lsig_opt, tvec_opt);
    kpss_opt  = kpss_stat_lrv(lsig_opt);
    fprintf('Orthogonality:  corr(log σ, μ) levels = %+.4f   diffs = %+.4f\n', ...
        corr_lvl, corr_dff);
    fprintf('Trend:          corr(log σ, t)        = %+.4f   (R^2 = %.4f)\n', ...
        corr_trnd, corr_trnd^2);
    fprintf('Stationarity:   KPSS-type statistic   = %.4f\n', kpss_opt);
    switch orth_mode
        case 'corr',  pen_opt = corr_lvl^2;
        case 'trend', pen_opt = corr_trnd^2;
        case 'kpss',  pen_opt = kpss_opt;
    end
    fprintf('                penalty (%s) contribution to objective = %.3f (w=%.1f)\n', ...
        orth_mode, w_corr * pen_opt, w_corr);
end

% m_eff diagnostics (negative-reserve regimes)
m_eff_us = exp(mu_us) - mubar_us_t - varrho_opt;
m_eff_eu = exp(mu_eu) - mubar_eu_t - varrho_opt;
fprintf('m_eff_us: min=%.4f  mean=%.4f  max=%.4f  (negative in %d/%d periods)\n', ...
    min(m_eff_us), mean(m_eff_us), max(m_eff_us), sum(m_eff_us<0), T_f);
fprintf('m_eff_eu: min=%.4f  mean=%.4f  max=%.4f  (negative in %d/%d periods)\n', ...
    min(m_eff_eu), mean(m_eff_eu), max(m_eff_eu), sum(m_eff_eu<0), T_f);

% List problem periods if any
problem_us = find(flag_us_t ~= 1);
problem_eu = find(flag_eu_t ~= 1);
if ~isempty(problem_us) && numel(problem_us) <= 20
    fprintf('US non-fsolve periods (t, flag, sigma, sig_min, ted_err_bps):\n');
    for k = 1:numel(problem_us)
        tt = problem_us(k);
        fprintf('  t=%3d  flag=%2d  sigma=%.4f  sig_min=%.4f  ted_err=%.3f\n', ...
            tt, flag_us_t(tt), sigma_us_opt(tt), sig_min_us_t(tt), ted_err_us_t(tt));
    end
elseif ~isempty(problem_us)
    fprintf('US non-fsolve periods: %d total (suppressed; > 20)\n', numel(problem_us));
end
if ~isempty(problem_eu) && numel(problem_eu) <= 20
    fprintf('EU non-fsolve periods (t, flag, sigma, sig_min, ted_err_bps):\n');
    for k = 1:numel(problem_eu)
        tt = problem_eu(k);
        fprintf('  t=%3d  flag=%2d  sigma=%.4f  sig_min=%.4f  ted_err=%.3f\n', ...
            tt, flag_eu_t(tt), sigma_eu_opt(tt), sig_min_eu_t(tt), ted_err_eu_t(tt));
    end
elseif ~isempty(problem_eu)
    fprintf('EU non-fsolve periods: %d total (suppressed; > 20)\n', numel(problem_eu));
end

fprintf('================================================\n');

%% Push estimated values back into the workspace for downstream scripts
lambda_us = lambda_us_opt;
lambda_eu = lambda_eu_opt;
eta       = eta_opt;
iota_ss   = iota_ss_opt;
varrho    = varrho_opt;
iota_us   = iota_ss / freq / pi_us_ss;
iota_eu   = iota_ss / freq / pi_eu_ss;

%% Snapshot the calibration to disk for downstream pipeline use
% File schema (loaded by params.m if present):
%   lambda_us_override_val, lambda_eu_override_val,
%   eta_override_val, iota_ss_override_val
lambda_us_override_val = lambda_us_opt;
lambda_eu_override_val = lambda_eu_opt;
eta_override_val       = eta_opt;
iota_ss_override_val   = iota_ss_opt;
% Non-default penalty modes save to a suffixed file so the committed baseline
% calibration (_calibration_override.mat, orth_mode='corr') is never clobbered.
% Fixed-eta runs always get their own suffix.
if ~isempty(eta_fix) && w_cond > 0
    override_file = sprintf('_calibration_override_eta%02.0f_cond%g.mat', 100 * eta_fix, w_cond);
elseif ~isempty(eta_fix)
    override_file = sprintf('_calibration_override_eta%02.0f.mat', 100 * eta_fix);
elseif strcmp(orth_mode, 'corr')
    override_file = '_calibration_override.mat';
else
    override_file = sprintf('_calibration_override_%s.mat', orth_mode);
end
% A driver may set override_suffix (e.g. '_lcr' for the mu-minus-LCR variant)
% so suffixed runs never clobber the unsuffixed files (incl. the committed
% baseline: corr + suffix -> _calibration_override_lcr.mat).
if ~exist('override_suffix', 'var'), override_suffix = ''; end
if ~isempty(override_suffix)
    override_file = strrep(override_file, '.mat', [override_suffix '.mat']);
end
save(override_file, ...
    'lambda_us_override_val', 'lambda_eu_override_val', ...
    'eta_override_val', 'iota_ss_override_val', ...
    'fval', 'fit_cip', 'fit_bpus', 'fit_bpeu', 'orth_mode');
fprintf('\n[snapshot] Saved calibration to %s\n', override_file);

%% ===== Local functions =====

function obj = estimate_obj_2d(x, ...
        mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
        TED_s_us_t, TED_s_eu_t, Rm_us_vec, Rm_eu_vec, ...
        Rb_Rm, Rb_Rm_eu, cip, ...
        pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us_init, sigma_eu_init, ...
        idx_cip, idx_bpus, idx_bpeu, idx_dw, idx_ff, idx_both, ...
        v_cip, v_bpus, v_bpeu, logmean_DW_d, logmean_FF_d, logratio_DF_d, ...
        w_cip, w_bpus, w_bpeu, w_meanDW, w_meanFF, w_ratioDF, w_meanDWS, w_corr, ...
        delta_DWS_estim, abs_sc, ...
        b_lambda, b_iota, b_eta, ...
        max_ted_bps, ted_asymp_margin, orth_mode, eta_fix, ...
        w_cond, idx_scr, idx_nor, cond_diff_data_bps)

    lam    = x(1);          % common lambda (paper simplification)
    lam_us = lam;
    lam_eu = lam;
    iot    = x(2);
    if isempty(eta_fix)
        et = x(3);          % eta estimated (3-D mode)
    else
        et = eta_fix;       % eta pinned (2-D mode)
    end
    vrh    = 0;             % varrho removed; hardcoded to 0

    % Box-bound penalty (continuous, scaled by distance outside the box)
    pen = 0;
    pen = pen + bound_pen(lam, b_lambda);
    pen = pen + bound_pen(iot, b_iota);
    pen = pen + bound_pen(et,  b_eta);

    % Asymptote-feasibility penalty: (1-eta)*iota_ss*1e4/pi_ss must exceed
    % the margin-scaled max TED in bps.  Holds independent of varrho.
    asy_bps = (1 - et) * iot * 1e4 / pi_us_ss;
    target_bps = ted_asymp_margin * max_ted_bps;
    if asy_bps < target_bps
        pen = pen + 1e6 * (1 + (target_bps - asy_bps) + (target_bps - asy_bps)^2);
    end

    if pen > 0
        obj = pen;
        return;
    end

    try
        [sigma_us_t, ~, BP_us_t, BP_eu_t, CIP_t, DW_us_t, FF_us_t] = run_filter_5d( ...
            lam_us, lam_eu, iot, et, vrh, ...
            mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
            TED_s_us_t, TED_s_eu_t, ...
            pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us_init, sigma_eu_init, ...
            Rm_us_vec, Rm_eu_vec);
    catch
        obj = 1e8;
        return;
    end

    % Reject if ANY period (not just indexed) is bad.  The filter is supposed
    % to invert TED at every t; if it can't for some, the parameter set is
    % infeasible and must not be rewarded by partial-sample averaging.
    if any(~isreal(BP_us_t))  || any(~isreal(BP_eu_t))  || ...
       any(~isreal(CIP_t))    || any(~isreal(DW_us_t))  || any(~isreal(FF_us_t)) || ...
       any(~isfinite(BP_us_t))|| any(~isfinite(BP_eu_t))|| ...
       any(~isfinite(CIP_t))  || any(~isfinite(DW_us_t))|| any(~isfinite(FF_us_t)) || ...
       any(DW_us_t <= 0)      || any(FF_us_t <= 0)
        obj = 1e8;
        return;
    end

    % --- Price residuals (bps) ---
    r_cip = (CIP_t(idx_cip) - cip(idx_cip)) * abs_sc;

    bpus_m_dm = BP_us_t(idx_bpus) - mean(BP_us_t(idx_bpus));
    bpus_d_dm = Rb_Rm(idx_bpus)   - mean(Rb_Rm(idx_bpus));
    r_bpus = (bpus_m_dm - bpus_d_dm) * abs_sc;

    bpeu_m_dm = BP_eu_t(idx_bpeu)  - mean(BP_eu_t(idx_bpeu));
    bpeu_d_dm = Rb_Rm_eu(idx_bpeu) - mean(Rb_Rm_eu(idx_bpeu));
    r_bpeu = (bpeu_m_dm - bpeu_d_dm) * abs_sc;

    % --- Volume log-mean moments ---
    m_DW = log(mean(DW_us_t(idx_dw))) - logmean_DW_d;
    m_FF = log(mean(FF_us_t(idx_ff))) - logmean_FF_d;
    % Log-ratio DW/FF: unit-free, immune to deposit-base normalization gap.
    m_ratioDF = (log(mean(DW_us_t(idx_both))) - log(mean(FF_us_t(idx_both)))) ...
                - logratio_DF_d;

    % DW STOCK moment: build DWS_t = DW_t + (1 - delta) * DWS_{t-1} and
    % compare its log mean to the data DW_n (which is already a stock).
    Tf = length(DW_us_t);
    DWS_us_t = zeros(Tf, 1);
    for tt = 1:Tf
        if tt == 1, prev = 0; else, prev = DWS_us_t(tt-1); end
        DWS_us_t(tt) = DW_us_t(tt) + (1 - delta_DWS_estim) * prev;
    end
    if any(DWS_us_t(idx_dw) <= 0)
        obj = 1e8;
        return;
    end
    m_DWS = log(mean(DWS_us_t(idx_dw))) - logmean_DW_d;   % data DW_n is the stock

    % --- Identification penalty (mode-dependent, unit-free) ---
    %   'corr'  : corr(log sigma, mu_us)^2  -- original orthogonality condition
    %   'trend' : corr(log sigma, t)^2      -- no secular trend in filtered stress
    %   'kpss'  : partial-sum stationarity statistic on demeaned log sigma
    valid_sm = isfinite(sigma_us_t) & sigma_us_t > 0 & isfinite(mu_us);
    if sum(valid_sm) > 10
        lsig = log(sigma_us_t(valid_sm));
        lsig = lsig(:);
        switch orth_mode
            case 'corr'
                mvec = mu_us(valid_sm);
                c = corr(lsig, mvec(:));
                if ~isfinite(c), c = 1; end   % penalize hard if undefined
                corr_pen = c^2;
            case 'trend'
                tvec = (1:numel(lsig))';
                c = corr(lsig, tvec);
                if ~isfinite(c), c = 1; end
                corr_pen = c^2;
            case 'kpss'
                % KPSS level-stationarity statistic with Bartlett long-run
                % variance (bandwidth l12). The LRV normalization is essential:
                % sigma is persistent by construction (stress regimes), and a
                % plain-variance denominator would punish persistence itself
                % rather than nonstationary drift.
                corr_pen = kpss_stat_lrv(lsig);
                if ~isfinite(corr_pen), corr_pen = 1; end
            otherwise
                error('Unknown orth_mode: %s', orth_mode);
        end
    else
        corr_pen = 1;   % too few valid periods -> hard penalty
    end

    % --- Conditional-amplification target (0 when w_cond = 0) ---
    cond_term = 0;
    if w_cond > 0
        diff_model_bps = (mean(CIP_t(idx_scr)) - mean(CIP_t(idx_nor))) * abs_sc;
        cond_term = ((diff_model_bps - cond_diff_data_bps) / cond_diff_data_bps)^2;
    end

    obj = w_cip    * mean(r_cip.^2)  / v_cip  ...
        + w_bpus   * mean(r_bpus.^2) / v_bpus ...
        + w_bpeu   * mean(r_bpeu.^2) / v_bpeu ...
        + w_meanDW  * m_DW^2 ...
        + w_meanFF  * m_FF^2 ...
        + w_ratioDF * m_ratioDF^2 ...
        + w_meanDWS * m_DWS^2 ...
        + w_corr    * corr_pen ...
        + w_cond    * cond_term;
end


function stat = kpss_stat_lrv(y)
    % KPSS level-stationarity statistic: sum(S_t^2) / (T^2 * s2_LRV), where
    % S_t are partial sums of the demeaned series and s2_LRV is the Bartlett
    % (Newey-West) long-run variance with the standard l12 bandwidth,
    % l = floor(12*(T/100)^(1/4)).  Unit-free; 5%% critical value = 0.463.
    y  = y(:);
    e  = y - mean(y);
    T  = numel(e);
    S  = cumsum(e);
    l  = floor(12 * (T / 100)^(1/4));
    s2 = (e' * e) / T;
    for j = 1:l
        w  = 1 - j / (l + 1);
        s2 = s2 + 2 * w * (e(1:T-j)' * e(1+j:T)) / T;
    end
    if s2 <= 0 || ~isfinite(s2)
        stat = 1;
    else
        stat = sum(S.^2) / (T^2 * s2);
    end
end


function p = bound_pen(v, b)
    % Quadratic penalty for v outside [b(1), b(2)].
    p = 0;
    if v < b(1)
        p = 1e6 * (1 + (b(1) - v) + (b(1) - v)^2);
    elseif v > b(2)
        p = 1e6 * (1 + (v - b(2)) + (v - b(2))^2);
    end
end


function [sigma_us_t, sigma_eu_t, BP_us_t, BP_eu_t, CIP_t, DW_us_t, FF_us_t, ...
          flag_us_t, flag_eu_t, sig_min_us_t, sig_min_eu_t, ted_err_us_t, ted_err_eu_t] = ...
    run_filter_5d(lam_us, lam_eu, iot, et, vrh, ...
        mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
        TED_s_us_t, TED_s_eu_t, ...
        pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us_init, sigma_eu_init, ...
        Rm_us_vec, Rm_eu_vec)
    % Run TED-inversion filter with separate lambda_us, lambda_eu, then
    % build the derived series: BP, DW, FF (US only), CIP.
    %
    % Optional diagnostic outputs (8-13):
    %   flag_us_t, flag_eu_t   -- per-period solver flag (see solve_sigma_*)
    %   sig_min_us_t, sig_min_eu_t -- per-period sig_min (CD only; NaN for Leontief)
    %   ted_err_us_t, ted_err_eu_t -- per-period TED inversion residual in bps

    iota_us_loc = iot / freq / pi_us_ss;
    iota_eu_loc = iot / freq / pi_eu_ss;
    abs_sc = 12e4;
    T_f = length(mu_us);

    sigma_us_t   = zeros(T_f, 1);
    sigma_eu_t   = zeros(T_f, 1);
    BP_us_t      = zeros(T_f, 1);
    BP_eu_t      = zeros(T_f, 1);
    DW_us_t      = zeros(T_f, 1);
    FF_us_t      = zeros(T_f, 1);
    flag_us_t    = zeros(T_f, 1);
    flag_eu_t    = zeros(T_f, 1);
    sig_min_us_t = NaN(T_f, 1);
    sig_min_eu_t = NaN(T_f, 1);
    ted_err_us_t = NaN(T_f, 1);
    ted_err_eu_t = NaN(T_f, 1);

    if matching_type == 1
        theta_plus_us = ((exp(lam_us) - 1) / (exp(lam_us) + 1))^2;
        theta_plus_eu = ((exp(lam_eu) - 1) / (exp(lam_eu) + 1))^2;
    end

    sig_us_guess = sigma_us_init;
    sig_eu_guess = sigma_eu_init;
    fsolve_opt = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12);

    for tt = 1:T_f
        mu_us_yt = exp(mu_us(tt)) - mubar_us_t(tt);
        mu_eu_yt = exp(mu_eu(tt)) - mubar_eu_t(tt);
        TED_us_tgt = TED_s_us_t(tt) * abs_sc;
        TED_eu_tgt = TED_s_eu_t(tt) * abs_sc;

        % --- US sigma ---
        if matching_type == 1
            sig_min_us = find_sigma_min(mu_us_yt, ploss_us, theta_plus_us, vrh);
            sig_min_us_t(tt) = sig_min_us;
            [sig_us, flag_us, ted_err_us] = solve_sigma_cd(mu_us_yt, ploss_us, iota_us_loc, lam_us, et, vrh, ...
                                    sig_us_guess, sig_min_us, TED_us_tgt, abs_sc, fsolve_opt);
        else
            [sig_us, flag_us, ted_err_us] = solve_sigma_leon(mu_us_yt, ploss_us, iota_us_loc, lam_us, et, vrh, ...
                                      sig_us_guess, TED_us_tgt, abs_sc, fsolve_opt);
        end
        flag_us_t(tt)    = flag_us;
        ted_err_us_t(tt) = ted_err_us;
        if ~isfinite(sig_us)
            sigma_us_t(:) = NaN; sigma_eu_t(:) = NaN;
            BP_us_t(:) = NaN; BP_eu_t(:) = NaN;
            DW_us_t(:) = NaN; FF_us_t(:) = NaN;
            CIP_t = NaN(T_f, 1);
            return;
        end
        sig_us_guess   = sig_us;
        sigma_us_t(tt) = sig_us;

        % --- EU sigma ---
        if matching_type == 1
            sig_min_eu = find_sigma_min(mu_eu_yt, ploss_eu, theta_plus_eu, vrh);
            sig_min_eu_t(tt) = sig_min_eu;
            [sig_eu, flag_eu, ted_err_eu] = solve_sigma_cd(mu_eu_yt, ploss_eu, iota_eu_loc, lam_eu, et, vrh, ...
                                    sig_eu_guess, sig_min_eu, TED_eu_tgt, abs_sc, fsolve_opt);
        else
            [sig_eu, flag_eu, ted_err_eu] = solve_sigma_leon(mu_eu_yt, ploss_eu, iota_eu_loc, lam_eu, et, vrh, ...
                                      sig_eu_guess, TED_eu_tgt, abs_sc, fsolve_opt);
        end
        flag_eu_t(tt)    = flag_eu;
        ted_err_eu_t(tt) = ted_err_eu;
        if ~isfinite(sig_eu)
            sigma_us_t(:) = NaN; sigma_eu_t(:) = NaN;
            BP_us_t(:) = NaN; BP_eu_t(:) = NaN;
            DW_us_t(:) = NaN; FF_us_t(:) = NaN;
            CIP_t = NaN(T_f, 1);
            return;
        end
        sig_eu_guess   = sig_eu;
        sigma_eu_t(tt) = sig_eu;

        % --- Derived series ---
        BP_us_t(tt) = Echi_m(mu_us_yt, ploss_us, sig_us, iota_us_loc, lam_us, et, matching_type, vrh);
        BP_eu_t(tt) = Echi_m(mu_eu_yt, ploss_eu, sig_eu, iota_eu_loc, lam_eu, et, matching_type, vrh);
        [~, ~, ~, ~, DW_us_t(tt), FF_us_t(tt)] = Chi_sys(mu_us_yt, ploss_us, sig_us, iota_us_loc, lam_us, et, matching_type, vrh);

        % Reject this parameter set entirely if any output is imaginary,
        % non-finite, or non-positive (DW/FF must be > 0 for log()).  Without
        % this, the optimizer is rewarded for choosing parameter regions
        % where the filter breaks down on a subset of periods.
        if ~isreal(BP_us_t(tt)) || ~isreal(BP_eu_t(tt)) || ...
           ~isreal(DW_us_t(tt)) || ~isreal(FF_us_t(tt)) || ...
           ~isfinite(BP_us_t(tt)) || ~isfinite(BP_eu_t(tt)) || ...
           ~isfinite(DW_us_t(tt)) || ~isfinite(FF_us_t(tt)) || ...
           DW_us_t(tt) <= 0 || FF_us_t(tt) <= 0
            sigma_us_t(:) = NaN; sigma_eu_t(:) = NaN;
            BP_us_t(:) = NaN; BP_eu_t(:) = NaN;
            DW_us_t(:) = NaN; FF_us_t(:) = NaN;
            CIP_t = NaN(T_f, 1);
            return;
        end
    end

    % --- CIP from main_filter.m definition ---
    Rb_us_t   = Rm_us_vec + BP_us_t;
    Rb_eu_t   = Rm_eu_vec + BP_eu_t;
    riskprm_t = Rb_us_t ./ Rb_eu_t - 1;
    UIP_t     = Rm_eu_vec - Rm_us_vec;
    CIP_t     = UIP_t + Rm_eu_vec .* riskprm_t;
end


function [sig, flag, ted_err] = solve_sigma_cd(mu_yt, ploss, iota_loc, lam, et, vrh, ...
                              sig_guess, sig_min, TED_tgt, abs_sc, fsolve_opt)
    % Robust Cobb-Douglas sigma solver.
    %   flag:   1 = fsolve success
    %           2 = fsolve success but sigma within 1e-3 of sig_min (CD corner)
    %          10 = fminbnd fallback used
    %          -1 = no finite residual found / both solvers failed
    %   ted_err = signed residual in bps (model TED - target TED) at returned sigma

    res_fun = @(s) Chi_p_psi(mu_yt, ploss, s, iota_loc, lam, et, 1, vrh) * abs_sc - TED_tgt;

    sig_try = sig_guess;
    if sig_try < sig_min + 1e-4
        sig_try = sig_min + max(1e-4, 2 * sig_min);
    end

    for k = 1:6
        if isfinite(safe_eval(res_fun, sig_try)), break; end
        sig_try = sig_try * 1.5 + 0.1;
    end
    if ~isfinite(safe_eval(res_fun, sig_try))
        sig = NaN; flag = -1; ted_err = NaN;
        return;
    end

    sig = NaN;
    ef  = -1;
    try
        [sig_out, ~, ef] = fsolve(res_fun, sig_try, fsolve_opt);
        sig = sig_out;
    catch
        ef = -1;
    end

    if isfinite(sig) && ef > 0 && sig >= sig_min
        flag = 1;
    else
        sq_fun = @(s) safe_sq(res_fun, s);
        try
            sig = fminbnd(sq_fun, sig_min + 1e-6, 100);
            flag = 10;
        catch
            sig = NaN; flag = -1; ted_err = NaN;
            return;
        end
    end

    if isfinite(sig) && (sig - sig_min) < 1e-3
        flag = 2;       % near-corner solution
    end

    ted_err = safe_eval(res_fun, sig);
end


function [sig, flag, ted_err] = solve_sigma_leon(mu_yt, ploss, iota_loc, lam, et, vrh, ...
                                sig_guess, TED_tgt, abs_sc, fsolve_opt)
    res_fun = @(s) Chi_p_psi(mu_yt, ploss, s, iota_loc, lam, et, 0, vrh) * abs_sc - TED_tgt;
    sig_try = max(sig_guess, 1e-3);
    for k = 1:6
        if isfinite(safe_eval(res_fun, sig_try)), break; end
        sig_try = sig_try * 1.5 + 0.1;
    end
    if ~isfinite(safe_eval(res_fun, sig_try))
        sig = NaN; flag = -1; ted_err = NaN;
        return;
    end
    sig = NaN;
    ef  = -1;
    try
        [sig_out, ~, ef] = fsolve(res_fun, sig_try, fsolve_opt);
        if ef > 0, sig = sig_out; end
    catch
    end
    if isfinite(sig) && ef > 0
        flag = 1;
    else
        sq_fun = @(s) safe_sq(res_fun, s);
        try
            sig = fminbnd(sq_fun, 1e-4, 100);
            flag = 10;
        catch
            sig = NaN; flag = -1; ted_err = NaN;
            return;
        end
    end
    ted_err = safe_eval(res_fun, sig);
end


function v = safe_eval(f, x)
    try
        v = f(x);
        if ~isfinite(v), v = NaN; end
    catch
        v = NaN;
    end
end


function v = safe_sq(f, x)
    r = safe_eval(f, x);
    if ~isfinite(r)
        v = 1e30;
    else
        v = r^2;
    end
end
