%% Estimate (lambda, iota_pre, iota_post, eta) with a GFC break in iota.
%
%   Same objective and filter pipeline as estimate_params_4d.m, but the
%   intermediation-cost parameter iota is allowed to differ across two
%   regimes split at the Global Financial Crisis:
%
%       Regime 1 (PRE) : months 1 .. t_break        (Jan-2001 .. Oct-2008)
%       Regime 2 (POST): months t_break+1 .. T       (Nov-2008 .. end)
%
%   with t_break = 94  (Oct-2008; sample starts Jan-2001, monthly).
%
%   Search vector x = [lambda, iota_pre, iota_post, eta]   (4-D).
%   lambda_us = lambda_eu and eta are common across regimes (paper
%   simplification); ONLY iota breaks.  varrho = 0.
%
%   The asymptote-feasibility constraint is applied PER REGIME:
%     (1-eta)*iota_r*1e4/pi_ss > ted_asymp_margin * max(TED over regime r).
%   The post-GFC TED spike is much larger, so iota_post generally faces a
%   tighter floor than iota_pre.
%
%   Run after the standard workspace setup, i.e.:
%       load('data/LFX_data.mat'); load('data/calibration.mat');
%       pi_us_ss = 1; pi_eu_ss = 1;
%       matching_type = 1;
%       run('functions/params.m');
%   (see run_estimation_iota_regime.m for the batch driver).

fprintf('\n====== Parameter Estimation: (lambda, iota_pre, iota_post, eta) [GFC break] ======\n');
fprintf('  Price targets:   CIP, US BP, EU BP  (time series, squared distance)\n');
fprintf('  Volume targets:  mean DW, mean FF   (scalar log-mean moments)\n');

%% ====================== USER SETTINGS ======================
% GFC breakpoint (month index; sample starts Jan-2001 monthly).
t_break = 94;               % Oct-2008 belongs to the PRE regime (t <= t_break)

% Target weights (each multiplies a unit-free term; set 0 to drop a target).
w_cip      = 0.0;
w_bpus     = 1.0;
w_bpeu     = 1.0;
w_meanDW   = 0.0;   % DW FLOW level (mismatched: model flow vs data stock)
w_meanFF   = 1.0;
w_ratioDF  = 0.0;
w_meanDWS  = 1.0;   % DW STOCK level (model stock vs data stock -- like-with-like)

% Orthogonality penalty:  corr(log sigma_us_t, mu_us_t)^2.
w_corr     = 100.0;

% Stock construction:  DWS_t = DW_t + (1 - delta_DWS_estim) * DWS_{t-1}
delta_DWS_estim = 0.2;   % monthly depreciation rate (must match main_filter.m)

% Parameter bounds (closed intervals; penalty pushes search back inside)
b_lambda    = [0.5, 12.0];   % common lambda for US and EU (paper simplification)
b_iota      = [0.005, 0.12];  % applied to BOTH iota_pre and iota_post
b_eta       = [0.05, 0.95];  % bargaining-power bounds (avoid degenerate extremes)

% Asymptote-feasibility margin (per regime).
ted_asymp_margin = 1.05;

% Initial guess.  Start both iota regimes at the pooled value and eta at 0.4.
lambda_init = 0.5 * (lambda_us + lambda_eu);
eta_init    = 0.4;
iota_init   = min(max(iota_ss, b_iota(1) + 1e-4), b_iota(2) - 1e-4);
x0 = [ ...
    min(max(lambda_init, b_lambda(1) + 1e-3), b_lambda(2) - 1e-3), ...
    iota_init, ...            % iota_pre
    iota_init, ...            % iota_post
    min(max(eta_init,    b_eta(1)    + 1e-3), b_eta(2)    - 1e-3)];
% ===========================================================

%% Baseline correction (recompute if not already set up)
T_obs = length(mu_us);
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

fprintf('GFC break at t_break = %d  ->  PRE = %s..%s (%d obs),  POST = %s..%s (%d obs)\n', ...
    t_break, datestr(datenum(2001,1,1),'mmm-yyyy'), datestr(datenum(2001,t_break,1),'mmm-yyyy'), ...
    t_break, datestr(datenum(2001,t_break+1,1),'mmm-yyyy'), datestr(datenum(2001,T_obs,1),'mmm-yyyy'), ...
    T_obs - t_break);

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

idx_both     = intersect(idx_dw, idx_ff);
logratio_DF_d = log(mean(DW_n(idx_both))) - log(mean(FF_n(idx_both)));

v_cip  = max(var(cip_d           * abs_sc), 1e-12);
v_bpus = max(var(bpus_d_demeaned * abs_sc), 1e-12);
v_bpeu = max(var(bpeu_d_demeaned * abs_sc), 1e-12);

fprintf('Data std: CIP=%.2f bps, BPus(dm)=%.2f bps, BPeu(dm)=%.2f bps\n', ...
    sqrt(v_cip), sqrt(v_bpus), sqrt(v_bpeu));
fprintf('Data:     log mean DW = %.4f, log mean FF = %.4f, log(DW/FF) = %.4f\n', ...
    logmean_DW_d, logmean_FF_d, logratio_DF_d);

%% TED-asymptote feasibility -- computed SEPARATELY for each regime.
pre_idx  = 1:t_break;
post_idx = (t_break+1):T_obs;

max_ted_pre  = max([TED_s_us_t(pre_idx);  TED_s_eu_t(pre_idx)])  * abs_sc;
max_ted_post = max([TED_s_us_t(post_idx); TED_s_eu_t(post_idx)]) * abs_sc;
fprintf('Max TED data:  PRE=%.2f bps,  POST=%.2f bps\n', max_ted_pre, max_ted_post);
fprintf('Asymptote floors (varrho=0):  iota_pre*(1-eta)*1e4 > %.2f,  iota_post*(1-eta)*1e4 > %.2f\n', ...
    ted_asymp_margin * max_ted_pre, ted_asymp_margin * max_ted_post);

% Bump initial iotas if they violate their regime floor at the starting eta.
for r = 1:2
    if r == 1, mted = max_ted_pre; else, mted = max_ted_post; end
    asy0 = (1 - x0(4)) * x0(1+r) * 1e4 / pi_us_ss;
    if asy0 < ted_asymp_margin * mted
        iota_min_feasible = ted_asymp_margin * mted / ((1 - x0(4)) * 1e4 / pi_us_ss);
        x0(1+r) = min(b_iota(2) - 1e-4, iota_min_feasible * 1.02);
        fprintf('Adjusted x0(iota_%s) -> %.4f (regime floor)\n', ...
            ternary(r==1,'pre','post'), x0(1+r));
    end
end

%% Anonymous objective (4D: lambda, iota_pre, iota_post, eta)
obj_fun = @(x) estimate_obj_regime(x, ...
    mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
    TED_s_us_t, TED_s_eu_t, Rm_us_vec, Rm_eu_vec, ...
    Rb_Rm, Rb_Rm_eu, cip, ...
    pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us, sigma_eu, ...
    idx_cip, idx_bpus, idx_bpeu, idx_dw, idx_ff, idx_both, ...
    v_cip, v_bpus, v_bpeu, logmean_DW_d, logmean_FF_d, logratio_DF_d, ...
    w_cip, w_bpus, w_bpeu, w_meanDW, w_meanFF, w_ratioDF, w_meanDWS, w_corr, ...
    delta_DWS_estim, abs_sc, ...
    b_lambda, b_iota, b_eta, ...
    max_ted_pre, max_ted_post, ted_asymp_margin, t_break);

%% Run fminsearch
fprintf('Starting fminsearch from x0:\n');
fprintf('  lambda=%.3f, iota_pre=%.4f, iota_post=%.4f, eta=%.3f  (4-D; lambda_us=lambda_eu)\n', ...
    x0(1), x0(2), x0(3), x0(4));
fprintf('Weights: w_cip=%.2f w_bpus=%.2f w_bpeu=%.2f w_meanDW=%.2f w_meanFF=%.2f w_ratioDF=%.2f w_meanDWS=%.2f w_corr=%.2f (delta=%.2f)\n', ...
    w_cip, w_bpus, w_bpeu, w_meanDW, w_meanFF, w_ratioDF, w_meanDWS, w_corr, delta_DWS_estim);

options = optimset('Display', 'iter', 'MaxFunEvals', 16000, 'MaxIter', 8000, ...
                   'TolX', 1e-7, 'TolFun', 1e-10);

%% Multi-start: screen several initial points, then polish from the best.
% Each start must satisfy BOTH regime asymptote constraints. Include points
% where iota_post > iota_pre (the natural GFC hypothesis) and a few equal.
x0_list = {
    x0(:)',                              % script default (equal iotas)
    [0.70, 0.030, 0.090, 0.35],          % low lambda, cost rises post-GFC
    [0.70, 0.020, 0.070, 0.55],          % low lambda, high eta
    [1.50, 0.030, 0.080, 0.30],          % mid lambda, low eta
    [1.50, 0.025, 0.060, 0.45],          % mid lambda, mid eta (~Nash)
    [1.50, 0.040, 0.090, 0.60],          % mid lambda, high eta
    [1.50, 0.060, 0.060, 0.45],          % mid lambda, EQUAL iotas (nested test)
    [3.00, 0.030, 0.080, 0.40],          % high lambda
    [5.00, 0.035, 0.090, 0.50],          % very high lambda
};

options_screen = optimset('Display', 'off', 'MaxFunEvals', 2000, 'MaxIter', 600, ...
                          'TolX', 1e-5, 'TolFun', 1e-6);

best_fval = Inf;
best_x    = x0_list{1};
fprintf('\n====== Multi-start screening (%d initial points) ======\n', numel(x0_list));
fprintf('  #   x0 [lambda, iota_pre, iota_post, eta]     ->  fval        exit  \n');
for ii = 1:numel(x0_list)
    x0_try = x0_list{ii};
    try
        [x_try, fval_try, ef_try] = fminsearch(obj_fun, x0_try, options_screen);
        fprintf('  %d   [%.3f, %.4f, %.4f, %.3f]  ->  %.4f      %d\n', ...
            ii, x0_try(1), x0_try(2), x0_try(3), x0_try(4), fval_try, ef_try);
        if fval_try < best_fval
            best_fval = fval_try;
            best_x    = x_try;
        end
    catch ME
        fprintf('  %d   [%.3f, %.4f, %.4f, %.3f]  FAILED: %s\n', ...
            ii, x0_try(1), x0_try(2), x0_try(3), x0_try(4), ME.message);
    end
end
fprintf('Best screening fval = %.4f at x = [%.4f, %.5f, %.5f, %.4f]\n', ...
    best_fval, best_x(1), best_x(2), best_x(3), best_x(4));
fprintf('Polishing from best with full-precision options...\n');
[x_opt, fval, exitflag] = fminsearch(obj_fun, best_x, options);

%% Results
lambda_opt     = x_opt(1);
lambda_us_opt  = lambda_opt;
lambda_eu_opt  = lambda_opt;
iota_pre_opt   = x_opt(2);
iota_post_opt  = x_opt(3);
eta_opt        = x_opt(4);
varrho_opt     = 0;
asy_pre_bps    = (1 - eta_opt) * iota_pre_opt  * 1e4 / pi_us_ss;
asy_post_bps   = (1 - eta_opt) * iota_post_opt * 1e4 / pi_us_ss;

fprintf('\n====== Estimation Results (iota GFC break at t=%d) ======\n', t_break);
fprintf('  lambda      = %.6f   (lambda_us = lambda_eu)\n', lambda_opt);
fprintf('  iota_pre    = %.6f   (Jan-2001 .. Oct-2008)\n', iota_pre_opt);
fprintf('  iota_post   = %.6f   (Nov-2008 .. end)\n', iota_post_opt);
fprintf('  iota ratio  = %.4f   (post / pre)\n', iota_post_opt / iota_pre_opt);
fprintf('  eta         = %.6f   (estimated; init at %.3f; common across regimes)\n', eta_opt, eta_init);
fprintf('  varrho      = %.6f   (removed; hardcoded to 0)\n', varrho_opt);
fprintf('  Asymptote PRE  (1-eta)*iota_pre*1e4  = %.2f bps  (vs max TED_pre  %.2f bps)\n', ...
    asy_pre_bps, max_ted_pre);
fprintf('  Asymptote POST (1-eta)*iota_post*1e4 = %.2f bps  (vs max TED_post %.2f bps)\n', ...
    asy_post_bps, max_ted_post);
fprintf('  Objective   = %.10e\n', fval);
fprintf('  Exit flag   = %d\n', exitflag);

%% Implied fit at the optimum (request diagnostic outputs too)
[sigma_us_opt, sigma_eu_opt, BP_us_opt, BP_eu_opt, CIP_opt, DW_us_opt, FF_us_opt, ...
    flag_us_t, flag_eu_t, sig_min_us_t, sig_min_eu_t, ted_err_us_t, ted_err_eu_t] = ...
    run_filter_regime(lambda_us_opt, lambda_eu_opt, iota_pre_opt, iota_post_opt, eta_opt, varrho_opt, ...
        mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
        TED_s_us_t, TED_s_eu_t, ...
        pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us, sigma_eu, ...
        Rm_us_vec, Rm_eu_vec, t_break);

% Price RMSEs in bps
fit_cip   = sqrt(mean(((CIP_opt(idx_cip)   - cip(idx_cip)) * abs_sc).^2));
bpus_m_dm = BP_us_opt(idx_bpus) - mean(BP_us_opt(idx_bpus));
bpeu_m_dm = BP_eu_opt(idx_bpeu) - mean(BP_eu_opt(idx_bpeu));
fit_bpus  = sqrt(mean(((bpus_m_dm - bpus_d_demeaned) * abs_sc).^2));
fit_bpeu  = sqrt(mean(((bpeu_m_dm - bpeu_d_demeaned) * abs_sc).^2));

logmean_DW_m = log(mean(DW_us_opt(idx_dw)));
logmean_FF_m = log(mean(FF_us_opt(idx_ff)));

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

logratio_DF_m = log(mean(DW_us_opt(idx_both))) - log(mean(FF_us_opt(idx_both)));
fprintf('\n--- Volume moments (logs) ---\n');
fprintf('                 model       data       diff\n');
fprintf('  mean(DWS)  :  %8.4f   %8.4f   %+8.4f    (stock; delta=%.2f)\n', ...
    logmean_DWS_m, logmean_DW_d, logmean_DWS_m - logmean_DW_d, delta_DWS_estim);
fprintf('  mean(FF)   :  %8.4f   %8.4f   %+8.4f\n', logmean_FF_m,  logmean_FF_d,  logmean_FF_m - logmean_FF_d);
fprintf('  log(DW/FF) :  %8.4f   %8.4f   %+8.4f\n', logratio_DF_m, logratio_DF_d, logratio_DF_m - logratio_DF_d);

% Orthogonality diagnostic at the optimum
valid_corr = isfinite(sigma_us_opt) & sigma_us_opt > 0 & isfinite(mu_us);
if sum(valid_corr) > 10
    corr_lvl = corr(log(sigma_us_opt(valid_corr)), mu_us(valid_corr));
    fprintf('\nOrthogonality:  corr(log sigma, mu) levels = %+.4f  (corr^2 obj contribution = %.3f, w_corr=%.1f)\n', ...
        corr_lvl, w_corr * corr_lvl^2, w_corr);
end

% Solver health
n_ok_us = sum(flag_us_t == 1); n_ok_eu = sum(flag_eu_t == 1);
fprintf('Solver OK (fsolve): US=%d/%d, EU=%d/%d\n', n_ok_us, Tf_opt, n_ok_eu, Tf_opt);
fprintf('sigma_us: range [%.4f, %.4f]  mean=%.4f\n', min(sigma_us_opt), max(sigma_us_opt), mean(sigma_us_opt));
fprintf('================================================\n');

%% Push estimated values back into the workspace for downstream scripts
lambda_us = lambda_us_opt;
lambda_eu = lambda_eu_opt;
eta       = eta_opt;
varrho    = varrho_opt;
% Downstream (single-iota) scripts read iota_ss; expose the POST value as the
% default and keep both regime values as separate variables.
iota_ss   = iota_post_opt;
iota_us   = iota_ss / freq / pi_us_ss;
iota_eu   = iota_ss / freq / pi_eu_ss;

%% Snapshot the calibration to disk
lambda_us_override_val = lambda_us_opt;
lambda_eu_override_val = lambda_eu_opt;
eta_override_val       = eta_opt;
iota_ss_override_val   = iota_post_opt;   % post-GFC as the headline single value
save('_calibration_override_iota_regime.mat', ...
    'lambda_us_override_val', 'lambda_eu_override_val', ...
    'eta_override_val', 'iota_ss_override_val', ...
    'iota_pre_opt', 'iota_post_opt', 't_break', ...
    'fval', 'fit_cip', 'fit_bpus', 'fit_bpeu');
fprintf('\n[snapshot] Saved calibration to _calibration_override_iota_regime.mat\n');

%% ===== Local functions =====

function obj = estimate_obj_regime(x, ...
        mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
        TED_s_us_t, TED_s_eu_t, Rm_us_vec, Rm_eu_vec, ...
        Rb_Rm, Rb_Rm_eu, cip, ...
        pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us_init, sigma_eu_init, ...
        idx_cip, idx_bpus, idx_bpeu, idx_dw, idx_ff, idx_both, ...
        v_cip, v_bpus, v_bpeu, logmean_DW_d, logmean_FF_d, logratio_DF_d, ...
        w_cip, w_bpus, w_bpeu, w_meanDW, w_meanFF, w_ratioDF, w_meanDWS, w_corr, ...
        delta_DWS_estim, abs_sc, ...
        b_lambda, b_iota, b_eta, ...
        max_ted_pre, max_ted_post, ted_asymp_margin, t_break)

    lam      = x(1);
    lam_us   = lam;
    lam_eu   = lam;
    iot_pre  = x(2);
    iot_post = x(3);
    et       = x(4);
    vrh      = 0;

    % Box-bound penalty
    pen = 0;
    pen = pen + bound_pen(lam,      b_lambda);
    pen = pen + bound_pen(iot_pre,  b_iota);
    pen = pen + bound_pen(iot_post, b_iota);
    pen = pen + bound_pen(et,       b_eta);

    % Asymptote-feasibility penalty PER REGIME.
    asy_pre  = (1 - et) * iot_pre  * 1e4 / pi_us_ss;
    asy_post = (1 - et) * iot_post * 1e4 / pi_us_ss;
    tgt_pre  = ted_asymp_margin * max_ted_pre;
    tgt_post = ted_asymp_margin * max_ted_post;
    if asy_pre < tgt_pre
        pen = pen + 1e6 * (1 + (tgt_pre - asy_pre) + (tgt_pre - asy_pre)^2);
    end
    if asy_post < tgt_post
        pen = pen + 1e6 * (1 + (tgt_post - asy_post) + (tgt_post - asy_post)^2);
    end

    if pen > 0
        obj = pen;
        return;
    end

    try
        [sigma_us_t, ~, BP_us_t, BP_eu_t, CIP_t, DW_us_t, FF_us_t] = run_filter_regime( ...
            lam_us, lam_eu, iot_pre, iot_post, et, vrh, ...
            mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
            TED_s_us_t, TED_s_eu_t, ...
            pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us_init, sigma_eu_init, ...
            Rm_us_vec, Rm_eu_vec, t_break);
    catch
        obj = 1e8;
        return;
    end

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
    m_ratioDF = (log(mean(DW_us_t(idx_both))) - log(mean(FF_us_t(idx_both)))) ...
                - logratio_DF_d;

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
    m_DWS = log(mean(DWS_us_t(idx_dw))) - logmean_DW_d;

    % --- Orthogonality penalty ---
    valid_sm = isfinite(sigma_us_t) & sigma_us_t > 0 & isfinite(mu_us);
    if sum(valid_sm) > 10
        lsig = log(sigma_us_t(valid_sm));
        mvec = mu_us(valid_sm);
        c = corr(lsig(:), mvec(:));
        if ~isfinite(c), c = 1; end
        corr_pen = c^2;
    else
        corr_pen = 1;
    end

    obj = w_cip    * mean(r_cip.^2)  / v_cip  ...
        + w_bpus   * mean(r_bpus.^2) / v_bpus ...
        + w_bpeu   * mean(r_bpeu.^2) / v_bpeu ...
        + w_meanDW  * m_DW^2 ...
        + w_meanFF  * m_FF^2 ...
        + w_ratioDF * m_ratioDF^2 ...
        + w_meanDWS * m_DWS^2 ...
        + w_corr    * corr_pen;
end


function p = bound_pen(v, b)
    p = 0;
    if v < b(1)
        p = 1e6 * (1 + (b(1) - v) + (b(1) - v)^2);
    elseif v > b(2)
        p = 1e6 * (1 + (v - b(2)) + (v - b(2))^2);
    end
end


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end


function [sigma_us_t, sigma_eu_t, BP_us_t, BP_eu_t, CIP_t, DW_us_t, FF_us_t, ...
          flag_us_t, flag_eu_t, sig_min_us_t, sig_min_eu_t, ted_err_us_t, ted_err_eu_t] = ...
    run_filter_regime(lam_us, lam_eu, iot_pre, iot_post, et, vrh, ...
        mu_us, mu_eu, mubar_us_t, mubar_eu_t, ploss_us, ploss_eu, ...
        TED_s_us_t, TED_s_eu_t, ...
        pi_us_ss, pi_eu_ss, freq, matching_type, sigma_us_init, sigma_eu_init, ...
        Rm_us_vec, Rm_eu_vec, t_break)
    % Same TED-inversion filter as estimate_params_4d.m's run_filter_5d, but
    % iota switches from iot_pre (t <= t_break) to iot_post (t > t_break).

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
        % --- Regime-dependent iota ---
        if tt <= t_break
            iot = iot_pre;
        else
            iot = iot_post;
        end
        iota_us_loc = iot / freq / pi_us_ss;
        iota_eu_loc = iot / freq / pi_eu_ss;

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
        flag = 2;
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
