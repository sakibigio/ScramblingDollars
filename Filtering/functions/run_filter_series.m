function out = run_filter_series(lambda, eta, iota_ss, ploss, in)
% RUN_FILTER_SERIES  Run the TED-inversion filter for a single
% parameter vector and return the full set of time series consumed
% downstream by plot_regimes.m / plot_robustness.m.
%
% Inputs
% ------
%   lambda, eta, iota_ss, ploss : scalars (iota_ss is annual, decimal)
%   in : struct with fields
%        mu_us, mu_eu, im_us, im_eu, TED_s_us_t, TED_s_eu_t   (column vecs)
%        matching_type, varrho, freq, pi_us_ss, pi_eu_ss      (scalars)
%        use_mu_baseline                                       (0/1)
%        mubar_us_t, mubar_eu_t                                (T x 1, zeros if baseline off)
%        sigma_us_init, sigma_eu_init                          (initial guesses)
%        abs_scale                                             (e.g. 12e4)
%
% Output
% ------
%   out : struct with fields used by plot_regimes
%        sigma_us_t, sigma_eu_t, sigma_us_TED_t, sigma_eu_TED_t
%        BP_us_t, BP_eu_t, TED_us_t, TED_eu_t
%        DW_us_t, FF_us_t, DWS_us_t, DW_eu_t, FF_eu_t
%        theta_us_t, psi_us_t, Smin_us_t, Q_us_t
%        Rm_us, Rm_eu, Rb_us_t, Rb_eu_t, Rd_us_t, Rd_eu_t
%        riskprm_t, UIP_t, CIP_t
%        Echi_m_us_t, Echi_d_us_t, Echi_m_eu_t, Echi_d_eu_t
%        sigma_us_TED_flag, sigma_eu_TED_flag
%        solver_summary (struct: n_fsolve, n_fminbnd, n_fail)
%
% This mirrors the inner loop of main_filter.m but is encapsulated so it
% can be called for arbitrary (lambda, eta, iota_ss, ploss) without
% touching globals.  Both currencies are filtered because downstream
% Markov estimation expects both sigma series.

% --- Unpack ---
mu_us = in.mu_us; mu_eu = in.mu_eu;
im_us = in.im_us; im_eu = in.im_eu;
TED_s_us_t = in.TED_s_us_t; TED_s_eu_t = in.TED_s_eu_t;
matching_type = in.matching_type; varrho = in.varrho; freq = in.freq;
pi_us_ss = in.pi_us_ss; pi_eu_ss = in.pi_eu_ss;
use_mu_baseline = in.use_mu_baseline;
mubar_us_t = in.mubar_us_t; mubar_eu_t = in.mubar_eu_t;
sigma_us_init = in.sigma_us_init; sigma_eu_init = in.sigma_eu_init;
abs_scale = in.abs_scale;

iota_us_loc = iota_ss / freq / pi_us_ss;
iota_eu_loc = iota_ss / freq / pi_eu_ss;

T = length(mu_us);

% --- Pre-allocate ---
sigma_us_t = zeros(T,1); sigma_eu_t = zeros(T,1);
sigma_us_TED_flag = zeros(T,1); sigma_eu_TED_flag = zeros(T,1);
Echi_m_us_t = zeros(T,1); Echi_d_us_t = zeros(T,1); Chi_p_psi_us_t = zeros(T,1);
Echi_m_eu_t = zeros(T,1); Echi_d_eu_t = zeros(T,1); Chi_p_psi_eu_t = zeros(T,1);
theta_us_t = zeros(T,1); psi_us_t = zeros(T,1); Smin_us_t = zeros(T,1);
DW_us_t = zeros(T,1); FF_us_t = zeros(T,1); Q_us_t = zeros(T,1);
theta_eu_t = zeros(T,1); psi_eu_t = zeros(T,1); Smin_eu_t = zeros(T,1);
DW_eu_t = zeros(T,1); FF_eu_t = zeros(T,1); Q_eu_t = zeros(T,1);

% --- Cobb-Douglas Walrasian thresholds ---
if matching_type == 1
    theta_plus_us = ((exp(lambda) - 1) / (exp(lambda) + 1))^2;
    theta_plus_eu = theta_plus_us;
end

sigma_us_TED_guess = sigma_us_init;
sigma_eu_TED_guess = sigma_eu_init;

for tt = 1:T
    if use_mu_baseline == 1
        mu_us_yt = exp(mu_us(tt)) - mubar_us_t(tt);
        mu_eu_yt = exp(mu_eu(tt)) - mubar_eu_t(tt);
    else
        mu_us_yt = exp(mu_us(tt));
        mu_eu_yt = exp(mu_eu(tt));
    end

    TED_us_target = TED_s_us_t(tt) * abs_scale;
    TED_eu_target = TED_s_eu_t(tt) * abs_scale;

    % ---- US sigma ----
    res_us = @(s) Chi_p_psi(mu_us_yt, ploss, s, iota_us_loc, lambda, eta, ...
                            matching_type, varrho) * abs_scale - TED_us_target;
    [sigma_us_t(tt), sigma_us_TED_flag(tt), sigma_us_TED_guess] = ...
        solve_sigma(res_us, mu_us_yt, ploss, iota_us_loc, lambda, eta, ...
                    matching_type, varrho, sigma_us_TED_guess, TED_us_target, abs_scale);

    % ---- EU sigma ----
    res_eu = @(s) Chi_p_psi(mu_eu_yt, ploss, s, iota_eu_loc, lambda, eta, ...
                            matching_type, varrho) * abs_scale - TED_eu_target;
    [sigma_eu_t(tt), sigma_eu_TED_flag(tt), sigma_eu_TED_guess] = ...
        solve_sigma(res_eu, mu_eu_yt, ploss, iota_eu_loc, lambda, eta, ...
                    matching_type, varrho, sigma_eu_TED_guess, TED_eu_target, abs_scale);

    % ---- Auxiliary functions of (mu, sigma) ----
    Echi_m_us_t(tt)    = Echi_m(mu_us_yt, ploss, sigma_us_t(tt), iota_us_loc, lambda, eta, matching_type, varrho);
    Echi_d_us_t(tt)    = Echi_d(mu_us_yt, ploss, sigma_us_t(tt), iota_us_loc, lambda, eta, matching_type, varrho);
    Chi_p_psi_us_t(tt) = Chi_p_psi(mu_us_yt, ploss, sigma_us_t(tt), iota_us_loc, lambda, eta, matching_type, varrho);
    Echi_m_eu_t(tt)    = Echi_m(mu_eu_yt, ploss, sigma_eu_t(tt), iota_eu_loc, lambda, eta, matching_type, varrho);
    Echi_d_eu_t(tt)    = Echi_d(mu_eu_yt, ploss, sigma_eu_t(tt), iota_eu_loc, lambda, eta, matching_type, varrho);
    Chi_p_psi_eu_t(tt) = Chi_p_psi(mu_eu_yt, ploss, sigma_eu_t(tt), iota_eu_loc, lambda, eta, matching_type, varrho);

    [~, theta_us_t(tt), psi_us_t(tt), Smin_us_t(tt), DW_us_t(tt), FF_us_t(tt), Q_us_t(tt)] = ...
        Chi_sys(mu_us_yt, ploss, sigma_us_t(tt), iota_us_loc, lambda, eta, matching_type, varrho);
    [~, theta_eu_t(tt), psi_eu_t(tt), Smin_eu_t(tt), DW_eu_t(tt), FF_eu_t(tt), Q_eu_t(tt)] = ...
        Chi_sys(mu_eu_yt, ploss, sigma_eu_t(tt), iota_eu_loc, lambda, eta, matching_type, varrho);
end

% --- Sanity check (mirror main_filter.m) ---
bad = ~isreal(DW_us_t) | ~isreal(FF_us_t) | ~isfinite(DW_us_t) | ~isfinite(FF_us_t) | ...
      DW_us_t <= 0 | FF_us_t <= 0 | ...
      ~isreal(DW_eu_t) | ~isreal(FF_eu_t) | ~isfinite(DW_eu_t) | ~isfinite(FF_eu_t) | ...
      DW_eu_t <= 0 | FF_eu_t <= 0;
if any(bad)
    warning('run_filter_series: %d periods with invalid DW/FF (lambda=%.4f, iota=%.4f, ploss=%.4f).', ...
        sum(bad), lambda, iota_ss, ploss);
end

% --- Derived rates / spreads (match main_filter.m) ---
Rm_us = exp(im_us);
Rm_eu = exp(im_eu);

BP_us_t  = Echi_m_us_t;
BP_eu_t  = Echi_m_eu_t;
TED_us_t = Chi_p_psi_us_t;
TED_eu_t = Chi_p_psi_eu_t;
Rd_us_t  = Rm_us + Echi_m_us_t + Echi_d_us_t;
Rd_eu_t  = Rm_eu + Echi_m_eu_t + Echi_d_eu_t;
Rb_us_t  = Rm_us + Echi_m_us_t;
Rb_eu_t  = Rm_eu + Echi_m_eu_t;
riskprm_t = Rb_us_t ./ Rb_eu_t - 1;
UIP_t     = Rm_eu - Rm_us;
CIP_t     = UIP_t + Rm_eu .* riskprm_t;

% --- DW loan stock (delta = 0.2 matches main_filter.m) ---
delta_DWS = 0.2;
DWS_us_t = zeros(T,1);
for tt = 1:T
    if tt == 1
        prev = 0;
    else
        prev = DWS_us_t(tt-1);
    end
    DWS_us_t(tt) = DW_us_t(tt) + (1 - delta_DWS) * prev;
end

% --- Pack output ---
out = struct();
out.sigma_us_t     = sigma_us_t;
out.sigma_eu_t     = sigma_eu_t;
out.sigma_us_TED_t = sigma_us_t;   % alias to match plot_regimes conventions
out.sigma_eu_TED_t = sigma_eu_t;
out.BP_us_t        = BP_us_t;
out.BP_eu_t        = BP_eu_t;
out.TED_us_t       = TED_us_t;
out.TED_eu_t       = TED_eu_t;
out.DW_us_t        = DW_us_t;
out.FF_us_t        = FF_us_t;
out.DWS_us_t       = DWS_us_t;
out.DW_eu_t        = DW_eu_t;
out.FF_eu_t        = FF_eu_t;
out.theta_us_t     = theta_us_t;
out.psi_us_t       = psi_us_t;
out.Smin_us_t      = Smin_us_t;
out.Q_us_t         = Q_us_t;
out.Rm_us          = Rm_us;
out.Rm_eu          = Rm_eu;
out.Rb_us_t        = Rb_us_t;
out.Rb_eu_t        = Rb_eu_t;
out.Rd_us_t        = Rd_us_t;
out.Rd_eu_t        = Rd_eu_t;
out.riskprm_t      = riskprm_t;
out.UIP_t          = UIP_t;
out.CIP_t          = CIP_t;
out.Echi_m_us_t    = Echi_m_us_t;
out.Echi_d_us_t    = Echi_d_us_t;
out.Echi_m_eu_t    = Echi_m_eu_t;
out.Echi_d_eu_t    = Echi_d_eu_t;
out.sigma_us_TED_flag = sigma_us_TED_flag;
out.sigma_eu_TED_flag = sigma_eu_TED_flag;

out.solver_summary = struct( ...
    'n_fsolve',  sum(sigma_us_TED_flag > 0 & sigma_us_TED_flag ~= 10), ...
    'n_fminbnd', sum(sigma_us_TED_flag == 10), ...
    'n_fail',    sum(sigma_us_TED_flag <= 0));
end


function [sig_out, ef_out, sig_guess_next] = solve_sigma( ...
        res_fun, mu_yt, ploss, iota_loc, lambda, eta, ...
        matching_type, varrho, sig_guess, TED_target, abs_scale)
% Solve TED residual for sigma, mirroring main_filter.m's per-period logic.
opt = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, ...
                   'MaxFunEval', 1e9, 'MaxIter', 1e6);

if matching_type == 1
    theta_plus = ((exp(lambda) - 1) / (exp(lambda) + 1))^2;
    sig_min = find_sigma_min(mu_yt, ploss, theta_plus, varrho);
    if sig_guess < sig_min + 1e-4
        sig_guess = sig_min + max(0.0001, 2 * sig_min);
    end
    % Walk guess up if residual is NaN there
    for kk = 1:6
        if isfinite(res_fun(sig_guess)), break; end
        sig_guess = sig_guess * 1.5 + 0.1;
    end

    try
        [sig_out, ~, ef, ~] = fsolve(res_fun, sig_guess, opt);
    catch
        sig_out = NaN; ef = -1;
    end

    if ef <= 0 || ~isfinite(sig_out) || sig_out < sig_min
        sq_fun = @(s) (Chi_p_psi(mu_yt, ploss, s, iota_loc, lambda, eta, 1, varrho) ...
                       * abs_scale - TED_target)^2;
        sig_hi = max(sig_min + 1e-3, 200);
        sig_out = fminbnd(sq_fun, sig_min + 1e-6, sig_hi);
        if isempty(sig_out) || ~isfinite(sig_out)
            sig_out = sig_min + 1e-3;
        end
        ef = 10;
    end
else
    [sig_out, ~, ef, ~] = fsolve(res_fun, sig_guess, opt);
end

ef_out = ef;
sig_guess_next = sig_out;
end
