%% Main Filter Script — joint (sigma, lambda) filter targeting TED + DW
% Code for Scrambling for Dollars
% (c) Saki Bigio
%
% Variant of main_filter.m: on the US side, when DW_n(tt) is available,
% search for (sigma_us, lambda_us) jointly by minimizing a weighted
% squared-residual objective against (TED, DW). For periods where DW_n
% is NaN (pre-2003), lambda_us is held at its calibrated value and only
% sigma_us is solved from TED, exactly as in main_filter.m.
%
% EU side and other currencies are filtered exactly as in main_filter.m.

%% Setup - preserve matching_type through clear
if ~exist('matching_type', 'var')
    matching_type = 1;  % 0 = Leontief, 1 = Cobb-Douglas (default)
end
save('temp_matching_type.mat', 'matching_type');

clear; close all;
load('temp_matching_type.mat');
delete('temp_matching_type.mat');

% This folder itself, so helpers here (overleaf_dir, acf_local) still resolve
% inside run('plotting/...'), which cds into plotting/ for the duration.
addpath(fileparts(mfilename('fullpath')));
addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

% Print/plot options
do_plot_baseline    = 1;
do_plot_diagnostics = 0;
do_counterfactual   = 0;
do_sensitivity      = 0;
plot_baseline_curr  = 0;

run_julia = 0;  % Skip Julia for now while we evaluate the joint fit

%% Load data
load('data/LFX_data.mat');

%% Parameters
load('data/calibration.mat');
abs_scale = 12e4;

curlist = {'au','ca','jp','nz','no','sw','ch','uk'};

share_us = 1; share_eu = 1;
pi_eu_ss = 1; pi_us_ss = 1;

run('functions/params.m');

epsilon_b = -0.5;
zeta_us = -epsilon_b;
zeta_eu = -epsilon_b;
Theta_b = 3;

T = length(mu_us);
datesperiod = 1:T;
dates = datenum(2001, 1:T, 1);

FSize = 20;
formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize, 'Box', 'Off', 'PlotBoxAspectRatio', [1 0.75 1]);
label_x = @(x) xlabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize, 'Interpreter', 'latex');
label_y = @(x) ylabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize, 'interpreter', 'latex');
desiredDecimalPlaces = 1;

% Re-load LFX_data.mat to restore the time-series im_us, im_eu (which
% params.m overwrites with steady-state scalars). This mirrors the
% structure of main_filter.m, which also reloads the data after params.m.
load('data/LFX_data.mat');

%% Identify DW coverage
first_dw = find(~isnan(DW_n) & DW_n > 0, 1, 'first');
if isempty(first_dw)
    error('No DW_n data found.');
end
fprintf('DW_n first available at tt=%d (%s)\n', first_dw, datestr(dates(first_dw),'mmm yyyy'));
fprintf('Periods with DW_n available: %d of %d\n', sum(~isnan(DW_n)), T);

%% Pre-allocate
sigma_us_t  = zeros(T,1);
sigma_eu_t  = zeros(T,1);
lambda_us_t = zeros(T,1);   % time-varying lambda on US side (filtered)
lambda_eu_t = zeros(T,1);   % EU lambda mirrors US (same shock per period)

Echi_m_us_t = zeros(T,1); Echi_d_us_t = zeros(T,1); Chi_p_psi_us_t = zeros(T,1);
Echi_m_eu_t = zeros(T,1); Echi_d_eu_t = zeros(T,1); Chi_p_psi_eu_t = zeros(T,1);

theta_us_t = zeros(T,1); psi_us_t = zeros(T,1); Smin_us_t = zeros(T,1);
DW_us_t    = zeros(T,1); FF_us_t  = zeros(T,1); Q_us_t    = zeros(T,1);
theta_eu_t = zeros(T,1); psi_eu_t = zeros(T,1); Smin_eu_t = zeros(T,1);
DW_eu_t    = zeros(T,1); FF_eu_t  = zeros(T,1); Q_eu_t    = zeros(T,1);

sigma_us_flag = zeros(T,1);
sigma_eu_flag = zeros(T,1);
joint_flag    = zeros(T,1);   % NEW: status of joint (sigma,lambda) solve
joint_resid_TED = zeros(T,1); % NEW: TED residual at solution (bps)
joint_resid_DW  = zeros(T,1); % NEW: DW residual at solution (ratio)

%% Filter loop
if matching_type == 0
    matching_name = 'Leontief';
else
    matching_name = 'Cobb-Douglas';
end
fprintf('Starting joint TED+DW filter with matching_type = %d (%s)\n', matching_type, matching_name);

% Weights for joint residual: scale by sample standard deviation so the
% two residuals are commensurate.
TED_scale = std(TED_s_us_t * abs_scale, 'omitnan');   % bps
DW_scale  = std(DW_n,                    'omitnan');  % ratio
fprintf('Residual scales: TED_scale=%.2f bps, DW_scale=%.4f\n', TED_scale, DW_scale);

% fminsearch options (handles NaN by treating as Inf)
opts_search = optimset('Display','off','TolX',1e-7,'TolFun',1e-9, ...
    'MaxIter',2000,'MaxFunEvals',2000);
opts_fsolve = optimoptions('fsolve','Display','off','TolFun',1e-12, ...
    'MaxFunEval',1e9,'MaxIter',1e6);

sigma_us_guess = sigma_us;
lambda_us_guess = lambda_us;
sigma_eu_guess = sigma_eu;

for tt = 1:T
    mu_us_yt = exp(mu_us(tt));
    mu_eu_yt = exp(mu_eu(tt));

    TED_us_target = TED_s_us_t(tt) * abs_scale;   % bps
    TED_eu_target = TED_s_eu_t(tt) * abs_scale;
    DW_us_target  = DW_n(tt);                     % ratio (may be NaN)

    %=========================================================
    % US side: joint (sigma, lambda) if DW available,
    %           otherwise single-equation sigma from TED.
    %=========================================================
    if isnan(DW_us_target)
        % --- Single-equation TED solve, lambda held at calibrated value ---
        lambda_use = lambda_us;
        [sig_sol, exitflag] = solve_sigma_from_TED(mu_us_yt, ploss_us, ...
            iota_us, lambda_use, eta, matching_type, TED_us_target, ...
            abs_scale, sigma_us_guess, opts_fsolve);
        sigma_us_t(tt)    = sig_sol;
        lambda_us_t(tt)   = lambda_use;
        sigma_us_flag(tt) = exitflag;
        joint_flag(tt)    = 0;
        sigma_us_guess    = sig_sol;
    else
        % --- Two-stage: outer 1D search over lambda; inner exact TED match ---
        %
        % For each candidate lambda, solve sigma(lambda) such that
        % TED_model(sigma, lambda) = TED_data exactly (1D fsolve — robust).
        % Then minimize squared DW gap, with a small smoothness penalty
        % anchoring lambda to the previous-period value. This stabilizes
        % the path in periods where DW is near zero across many lambdas.
        sigma_warm = sigma_us_guess;
        if tt > 1 && lambda_us_t(tt-1) > 0
            lambda_anchor = lambda_us_t(tt-1);
        else
            lambda_anchor = lambda_us;  % calibrated
        end
        tau = 1e-4;   % smoothness weight (tiny — only binds when DW gap is flat)
        outer_obj = @(lam) dw_gap_along_TED_iso(lam, mu_us_yt, ploss_us, ...
            iota_us, eta, matching_type, TED_us_target, DW_us_target, ...
            abs_scale, sigma_warm, opts_fsolve) ...
            + tau * (lam - lambda_anchor)^2;

        % Bracketed 1D search. Widened range to let the filter explore.
        % Infeasible (sigma_min > SIG_MAX) lambdas auto-return 1e8 from
        % the helper, so fminbnd will steer away.
        lam_lo = 0.05; lam_hi = 6.0;
        [lam_sol, fval, exitflag] = fminbnd(outer_obj, lam_lo, lam_hi, ...
            optimset('Display','off','TolX',1e-5,'MaxIter',200));

        % Recompute sigma at the chosen lambda
        [sig_sol, ~] = solve_sigma_from_TED(mu_us_yt, ploss_us, iota_us, ...
            lam_sol, eta, matching_type, TED_us_target, abs_scale, ...
            sigma_warm, opts_fsolve);

        sigma_us_t(tt)    = sig_sol;
        lambda_us_t(tt)   = lam_sol;
        joint_flag(tt)    = exitflag;
        sigma_us_guess    = sig_sol;
        lambda_us_guess   = lam_sol;
    end

    % Record fitted moments and residuals for US
    ted_mod = Chi_p_psi(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, ...
        lambda_us_t(tt), eta, matching_type) * abs_scale;
    [~, theta_us_t(tt), psi_us_t(tt), Smin_us_t(tt), DW_us_t(tt), ...
        FF_us_t(tt), Q_us_t(tt)] = Chi_sys(mu_us_yt, ploss_us, ...
        sigma_us_t(tt), iota_us, lambda_us_t(tt), eta, matching_type);
    joint_resid_TED(tt) = ted_mod - TED_us_target;
    if ~isnan(DW_us_target)
        joint_resid_DW(tt) = DW_us_t(tt) - DW_us_target;
    else
        joint_resid_DW(tt) = NaN;
    end

    Echi_m_us_t(tt)   = Echi_m(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us_t(tt), eta, matching_type);
    Echi_d_us_t(tt)   = Echi_d(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us_t(tt), eta, matching_type);
    Chi_p_psi_us_t(tt)= Chi_p_psi(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us_t(tt), eta, matching_type);

    %=========================================================
    % EU side: single-equation TED solve, using the US-filtered
    % lambda_us_t(tt) as the EU matching efficiency (per user request).
    %=========================================================
    lam_eu_use = lambda_us_t(tt);
    lambda_eu_t(tt) = lam_eu_use;
    [sig_eu_sol, eflag] = solve_sigma_from_TED(mu_eu_yt, ploss_eu, ...
        iota_eu, lam_eu_use, eta, matching_type, TED_eu_target, abs_scale, ...
        sigma_eu_guess, opts_fsolve);
    if isnan(sig_eu_sol) || ~isscalar(sig_eu_sol)
        % Infeasible at lam_us(tt): fall back to previous-period sigma_eu.
        if tt > 1, sig_eu_sol = sigma_eu_t(tt-1); else, sig_eu_sol = sigma_eu_guess; end
        eflag = -50;
    end
    sigma_eu_t(tt) = sig_eu_sol;
    sigma_eu_flag(tt) = eflag;
    sigma_eu_guess = sig_eu_sol;

    Echi_m_eu_t(tt)   = Echi_m(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lam_eu_use, eta, matching_type);
    Echi_d_eu_t(tt)   = Echi_d(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lam_eu_use, eta, matching_type);
    Chi_p_psi_eu_t(tt)= Chi_p_psi(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lam_eu_use, eta, matching_type);
    [~, theta_eu_t(tt), psi_eu_t(tt), Smin_eu_t(tt), DW_eu_t(tt), FF_eu_t(tt), Q_eu_t(tt)] = Chi_sys(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lam_eu_use, eta, matching_type);

    if mod(tt, 50) == 0
        fprintf('  Period %d/%d  sigma_us=%.3f  lambda_us=%.3f  TED_res=%.1f  DW_res=%.4f\n', ...
            tt, T, sigma_us_t(tt), lambda_us_t(tt), joint_resid_TED(tt), joint_resid_DW(tt));
    end
end

%% Construct rates
Rm_us = exp(im_us);
Rm_eu = exp(im_eu);
inv_e_yt = exp(inv_e);
Rb_Rm_yt = Rb_Rm;

%% Derive other variables
BP_us_t = Echi_m_us_t;
TED_us_t = Chi_p_psi_us_t;
Rd_us_t = Rm_us + Echi_m_us_t + Echi_d_us_t;
Rb_us_t = Rm_us + Echi_m_us_t;
BP_eu_t = Echi_m_eu_t;
TED_eu_t = Chi_p_psi_eu_t;
Rd_eu_t = Rm_eu + Echi_m_eu_t + Echi_d_eu_t;
Rb_eu_t = Rm_eu + Echi_m_eu_t;

riskprm_t = Rb_us_t ./ Rb_eu_t - 1;
UIP_t = Rm_eu - Rm_us;
CIP_t = UIP_t + Rm_eu .* riskprm_t;

f_t = (1 + riskprm_t) ./ exp(inv_e);

%% Quantity variables
inv_e_t  = exp(inv_e);
mu_us_t  = exp(mu_us);
mu_eu_t  = exp(mu_eu);
M_us_t   = exp(M_us);
M_eu_t   = exp(M_eu);
e_euus_t = 1 ./ inv_e_t;
M_rat_t  = M_us_t  ./ M_eu_t;
mu_rat_t = mu_us_t ./ mu_eu_t;

nu_t = (M_eu_t ./ M_us_t) .* inv_e_t .* mu_us_t ./ mu_eu_t;
b_t = Theta_b * (Rm_us + Echi_m_us_t).^(1/epsilon_b);
d_us_t = b_t ./ ((1 - mu_us_t) + nu_t .* (1 - mu_eu_t));
d_eu_t = nu_t .* d_us_t;

p_us_t = M_us_t ./ (d_us_t .* mu_us_t);
p_eu_t = M_eu_t ./ (d_eu_t .* mu_eu_t);

N_pre = 72;
MBS_us_ss = mean(d_us_t(1:N_pre));
MBS_eu_ss = mean(d_eu_t(1:N_pre));
Theta_d_us_t = (d_us_t - share_us * (mu_us_t .* d_us_t) + MBS_us_ss) ./ (Rd_us_t).^(1/zeta_us);
Theta_d_eu_t = (d_eu_t - share_eu * (mu_eu_t .* d_eu_t) + MBS_eu_ss) ./ (Rd_eu_t).^(1/zeta_eu);

riskprm_shock_t = riskprm_t(datesperiod) - (exp(im_us(datesperiod)) - exp(im_eu(datesperiod)));

%% Aliases for compatibility with existing plotting scripts
% (plot_baseline / plot_diagnostics expect *_TED_t and *_TED_flag variants)
sigma_us_TED_t    = sigma_us_t;
sigma_eu_TED_t    = sigma_eu_t;
sigma_us_TED_flag = joint_flag;        % flag from joint solve (or 0 = TED-only fallback)
sigma_eu_TED_flag = sigma_eu_flag;

%% Save output
sigma_mat = [sigma_us_t sigma_eu_t lambda_us_t];
RW_shock = array2table(sigma_mat, 'VariableNames', {'sigma_us','sigma_eu','lambda_us'});
writetable(RW_shock, 'RW_shock_DW.csv');
fprintf('Filter complete. Output saved to RW_shock_DW.csv\n');

%% Diagnostics
fprintf('\n=== Joint TED+DW Filter Diagnostics ===\n');
fprintf('matching_type=%d, lambda_us(calibrated)=%.2f, eta=%.2f, ploss_us=%.2f, iota_us=%.6f\n', ...
    matching_type, lambda_us, eta, ploss_us, iota_us);

% Sigma and lambda summary
fprintf('\n--- Estimated paths ---\n');
fprintf('sigma_us: mean=%.4f, std=%.4f, min=%.4f, max=%.4f\n', ...
    mean(sigma_us_t), std(sigma_us_t), min(sigma_us_t), max(sigma_us_t));
mask_joint = ~isnan(DW_n);
fprintf('lambda_us (over DW-available periods): mean=%.3f, std=%.3f, min=%.3f, max=%.3f\n', ...
    mean(lambda_us_t(mask_joint)), std(lambda_us_t(mask_joint)), min(lambda_us_t(mask_joint)), max(lambda_us_t(mask_joint)));

% Residuals
fprintf('\n--- Joint residuals (DW-available periods) ---\n');
resid_TED = joint_resid_TED(mask_joint);
resid_DW  = joint_resid_DW(mask_joint);
fprintf('TED residual (bps): mean=%.2f, median=%.2f, max|.|=%.2f, p95|.|=%.2f\n', ...
    mean(resid_TED), median(resid_TED), max(abs(resid_TED)), prctile(abs(resid_TED),95));
fprintf('DW residual:        mean=%.4f, median=%.4f, max|.|=%.4f, p95|.|=%.4f\n', ...
    mean(resid_DW), median(resid_DW), max(abs(resid_DW)), prctile(abs(resid_DW),95));

% Fit comparison
fprintf('\n--- Model vs Data correlations (full sample) ---\n');
fprintf('corr(TED_model, TED_data) = %.3f\n', corr(TED_us_t, TED_s_us_t));
mask = ~isnan(DW_n);
fprintf('corr(DW_model,  DW_data)  = %.3f   (over %d valid periods)\n', ...
    corr(DW_us_t(mask), DW_n(mask)), sum(mask));

% Volume comparison
fprintf('\n--- Means over DW-available periods ---\n');
fprintf('DW: model=%.4f   data=%.4f\n', mean(DW_us_t(mask)), mean(DW_n(mask)));
fprintf('FF: model=%.4f   data=%.4f\n', mean(FF_us_t(mask)), mean(FF_n(mask),'omitnan'));

%% Plotting
% Time-series plot of the three key series
figure('Name','Joint filter outputs','Position',[80 80 1200 700]);
subplot(3,1,1);
plot(dates, sigma_us_t,'LineWidth',1.5); hold on;
plot(dates, sigma_eu_t,'LineWidth',1.0,'Color',[0.5 0.5 0.5]);
datetick('x','yyyy'); grid on; axis tight;
ylabel('\sigma'); legend('US','EU','Location','best');
title('Estimated \sigma');

subplot(3,1,2);
plot(dates, lambda_us_t,'LineWidth',1.5); hold on;
yline(lambda_us,'r--','calibrated');
datetick('x','yyyy'); grid on; axis tight;
ylabel('\lambda_{us}');
title('Estimated \lambda_{us} (joint TED+DW solve where DW available)');

subplot(3,1,3);
yyaxis left;
plot(dates, TED_us_t*abs_scale,'LineWidth',1.0); hold on;
plot(dates, TED_s_us_t*abs_scale,'k:','LineWidth',1.0);
ylabel('TED (bps)');
yyaxis right;
plot(dates, DW_us_t,'LineWidth',1.0); hold on;
plot(dates, DW_n,'k--','LineWidth',1.0);
ylabel('DW (ratio)');
datetick('x','yyyy'); grid on; axis tight;
legend('TED model','TED data','DW model','DW data','Location','best');
title('Fit: TED and DW (model vs data)');

saveas(gcf, 'joint_filter_DW.png');
fprintf('Saved joint_filter_DW.png\n');

%% Existing plotting scripts (run with time-varying lambda already baked into series)
if do_plot_baseline == 1
    run('plotting/plot_baseline.m');
end
if do_plot_diagnostics == 1
    run('plotting/plot_diagnostics.m');
end
if do_counterfactual == 1
    run('plotting/plot_counterfactual.m');
end
if do_sensitivity == 1
    run('plotting/plot_sensitivity.m');
end

%% Local helpers

% Solve sigma such that model TED matches target, holding lambda fixed.
% Returns NaN sigma + exitflag<0 when the model is infeasible (e.g. lambda
% so high that the matching cliff sigma_min exceeds any reasonable scale).
function [sig, exitflag] = solve_sigma_from_TED(mu_yt, ploss, iota, lam, eta, mt, TED_target, abs_scale, sig_warm, opts_fsolve)
    res = @(sig) Chi_p_psi(mu_yt, ploss, sig, iota, lam, eta, mt) * abs_scale - TED_target;
    SIG_MAX = 15;
    if mt == 1
        % Cobb-Douglas: respect the sigma_min cliff for this lambda
        theta_plus = ((exp(lam) - 1) / (exp(lam) + 1))^2;
        sigma_min = find_sigma_min(mu_yt, ploss, theta_plus);
        if ~isfinite(sigma_min) || sigma_min >= SIG_MAX - 1e-3
            % Infeasible: matching cliff above search range. Return NaN.
            sig = NaN; exitflag = -1; return;
        end
        guess = max(sig_warm, sigma_min + max(1e-4, 2*sigma_min));
        try
            [sig, ~, exitflag] = fsolve(res, guess, opts_fsolve);
        catch
            sig = guess; exitflag = -99;
        end
        if exitflag <= 0 || sig < sigma_min
            try
                sig = fminbnd(@(s) (Chi_p_psi(mu_yt, ploss, s, iota, lam, eta, mt) * abs_scale - TED_target)^2, ...
                    sigma_min + 1e-6, SIG_MAX);
                exitflag = 10;
            catch
                sig = NaN; exitflag = -2;
            end
        end
    else
        % Leontief
        try
            [sig, ~, exitflag] = fsolve(res, max(sig_warm, 0.01), opts_fsolve);
        catch
            sig = max(sig_warm, 0.01); exitflag = -99;
        end
    end
end

% For a candidate lambda, recover sigma(lambda) along the TED iso and
% return the squared DW gap.
function gap2 = dw_gap_along_TED_iso(lam, mu_yt, ploss, iota, eta, mt, TED_target, DW_target, abs_scale, sig_warm, opts_fsolve)
    if lam <= 0
        gap2 = 1e8; return;
    end
    try
        [sig, ef] = solve_sigma_from_TED(mu_yt, ploss, iota, lam, eta, mt, TED_target, abs_scale, sig_warm, opts_fsolve);
        if isnan(sig) || ~isscalar(sig) || (ef <= 0 && ef ~= 10)
            gap2 = 1e8; return;
        end
        [~, ~, ~, ~, DW_v, ~, ~] = Chi_sys(mu_yt, ploss, sig, iota, lam, eta, mt);
        if ~isfinite(DW_v) || ~isscalar(DW_v)
            gap2 = 1e8; return;
        end
        gap2 = (DW_v - DW_target)^2;
        if ~isfinite(gap2) || ~isscalar(gap2)
            gap2 = 1e8;
        end
    catch
        gap2 = 1e8;
    end
end
