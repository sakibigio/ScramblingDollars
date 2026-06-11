%% Estimate (lambda, eta, iota_ss) by matching 3 aggregate moments:
%    M1: log mean(FF_model)    vs log mean(FF_data)   -- average FF volume
%    M2: log mean(DW_model)    vs log mean(DW_data)   -- average DW volume
%    M3: log max (DW_model)    vs log max (DW_data)   -- DW crisis peak
%
%  Objective = M1^2 + M2^2 + M3^2  (just-identified: 3 moments, 3 params).
%
%  Run after main_filter.m lines 1-76 have set up the workspace, i.e.:
%    load('data/LFX_data.mat');  load('data/calibration.mat');
%    run('functions/params.m');
%
%  Requires in workspace:
%    mu_us, ploss_us, TED_s_us_t, pi_us_ss, freq,
%    matching_type, sigma_us (initial guess), FF_n, DW_n
%  Requires on path: Chi_p_psi, Chi_sys, find_sigma_min, liquidity_baseline

fprintf('\n====== Parameter Estimation: (lambda, eta, iota_ss) ======\n');
fprintf('  Targets: mean(FF), mean(DW), max(DW) -- all in logs.\n');

%% Baseline correction (matches main_filter.m: mu_yt = exp(mu) - mubar)
T_obs = length(mu_us);
dates_obs = datenum(2001, 1:T_obs, 1);
[mubar_us_t, ~] = liquidity_baseline(dates_obs);
fprintf('Baseline correction ON: mubar_us range=[%.3f,%.3f]\n', ...
    min(mubar_us_t), max(mubar_us_t));

%% Moment-specific valid samples
valid_FF = ~isnan(FF_n) & FF_n > 0;
valid_DW = ~isnan(DW_n) & DW_n > 0;
idx_FF   = find(valid_FF);
idx_DW   = find(valid_DW);
fprintf('FF sample: %d periods.\n', numel(idx_FF));
fprintf('DW sample: %d periods.\n', numel(idx_DW));

%% Precompute data targets in log space
logmean_FF_d = log(mean(FF_n(idx_FF)));
logmean_DW_d = log(mean(DW_n(idx_DW)));
logmax_DW_d  = log(max(DW_n(idx_DW)));
fprintf('Data: log mean FF = %.4f, log mean DW = %.4f, log max DW = %.4f\n', ...
    logmean_FF_d, logmean_DW_d, logmax_DW_d);

%% Build anonymous objective
obj_fun = @(x) estimate_obj_fn(x, mu_us, mubar_us_t, ploss_us, TED_s_us_t, ...
    pi_us_ss, freq, matching_type, sigma_us, idx_FF, idx_DW, ...
    logmean_FF_d, logmean_DW_d, logmax_DW_d);

%% Run fminsearch
x0 = [1.0, 0.65, 0.1];  % [lambda, eta, iota_ss]
options = optimset('Display', 'iter', 'MaxFunEvals', 5000, 'MaxIter', 3000, ...
                   'TolX', 1e-7, 'TolFun', 1e-12);

fprintf('Starting fminsearch from x0 = [%.2f, %.2f, %.2f] ...\n', x0);
[x_opt, fval, exitflag] = fminsearch(obj_fun, x0, options);

%% Results
lambda_opt  = x_opt(1);
eta_opt     = x_opt(2);
iota_ss_opt = x_opt(3);

fprintf('\n====== Estimation Results ======\n');
fprintf('  lambda    = %.6f\n', lambda_opt);
fprintf('  eta       = %.6f\n', eta_opt);
fprintf('  iota_ss   = %.6f\n', iota_ss_opt);
fprintf('  Objective = %.10e\n', fval);
fprintf('  Exit flag = %d\n', exitflag);

%% Implied moment fit at optimum
[~, FF_opt, DW_opt] = run_filter_for_params(lambda_opt, eta_opt, iota_ss_opt, ...
    mu_us, mubar_us_t, ploss_us, TED_s_us_t, pi_us_ss, freq, matching_type, sigma_us);

logmean_FF_m = log(mean(FF_opt(idx_FF)));
logmean_DW_m = log(mean(DW_opt(idx_DW)));
logmax_DW_m  = log(max(DW_opt(idx_DW)));

fprintf('\n--- Moment fit ---\n');
fprintf('              model       data       diff (logs)\n');
fprintf('  mean(FF):  %8.4f   %8.4f   %+8.4f\n', logmean_FF_m, logmean_FF_d, logmean_FF_m - logmean_FF_d);
fprintf('  mean(DW):  %8.4f   %8.4f   %+8.4f\n', logmean_DW_m, logmean_DW_d, logmean_DW_m - logmean_DW_d);
fprintf('  max (DW):  %8.4f   %8.4f   %+8.4f\n', logmax_DW_m,  logmax_DW_d,  logmax_DW_m  - logmax_DW_d);
fprintf('\n  Levels:\n');
fprintf('  mean(FF):  %8.4f   %8.4f\n', mean(FF_opt(idx_FF)), mean(FF_n(idx_FF)));
fprintf('  mean(DW):  %8.4f   %8.4f\n', mean(DW_opt(idx_DW)), mean(DW_n(idx_DW)));
fprintf('  max (DW):  %8.4f   %8.4f\n', max(DW_opt(idx_DW)),  max(DW_n(idx_DW)));
fprintf('================================================\n');

%% ===== Local functions =====

function obj = estimate_obj_fn(x, mu_us, mubar_us_t, ploss_us, TED_s_us_t, ...
        pi_us_ss, freq, matching_type, sigma_us_init, idx_FF, idx_DW, ...
        logmean_FF_d, logmean_DW_d, logmax_DW_d)
    % Unpack
    lam  = x(1);
    et   = x(2);
    iot  = x(3);

    % Penalty for out-of-bounds. eta < 0.8 per numerical-stability note.
    pen = 0;
    if lam <= 0,    pen = pen + 1e6 * (1 + abs(lam));     end
    if et  <= 0.1,  pen = pen + 1e6 * (1 + abs(0.1-et));  end
    if et  >= 0.8,  pen = pen + 1e6 * (1 + abs(et-0.8));  end
    if iot <= 0,    pen = pen + 1e6 * (1 + abs(iot));     end
    if iot >= 0.15, pen = pen + 1e6 * (1 + abs(iot-0.15));end
    if pen > 0
        obj = pen;
        return;
    end

    [~, FF_m, DW_m] = run_filter_for_params(lam, et, iot, ...
        mu_us, mubar_us_t, ploss_us, TED_s_us_t, pi_us_ss, freq, matching_type, sigma_us_init);

    % Guard against non-positive model values on the relevant samples
    if any(FF_m(idx_FF) <= 0) || any(DW_m(idx_DW) <= 0)
        obj = 1e8;
        return;
    end

    % Three aggregate moments, log space, squared
    m1 = log(mean(FF_m(idx_FF))) - logmean_FF_d;
    m2 = log(mean(DW_m(idx_DW))) - logmean_DW_d;
    m3 = log(max (DW_m(idx_DW))) - logmax_DW_d;
    obj = m1^2 + m2^2 + m3^2;
end

function [sigma_t, FF_t, DW_t] = run_filter_for_params(lam, et, iot, ...
        mu_us, mubar_us_t, ploss_us, TED_s_us_t, pi_us_ss, freq, matching_type, sigma_us_init)
    % Run the TED-inversion filter for given (lambda, eta, iota_ss).
    % Applies baseline correction: mu_t = exp(mu_us(tt)) - mubar_us_t(tt).
    % varrho defaults to 0 inside the chi functions (not threaded here).

    iota_us_loc = iot / freq / pi_us_ss;
    abs_sc = 12e4;
    T_f = length(mu_us);

    sigma_t = zeros(T_f, 1);
    FF_t    = zeros(T_f, 1);
    DW_t    = zeros(T_f, 1);

    % Cobb-Douglas threshold
    if matching_type == 1
        theta_plus = ((exp(lam) - 1) / (exp(lam) + 1))^2;
    end

    sig_guess = sigma_us_init;

    for tt = 1:T_f
        mu_t = exp(mu_us(tt)) - mubar_us_t(tt);
        TED_tgt = TED_s_us_t(tt) * abs_sc;

        % Invert Chi_p_psi to find sigma
        res_fun = @(sig) Chi_p_psi(mu_t, ploss_us, sig, iota_us_loc, ...
                                    lam, et, matching_type) * abs_sc - TED_tgt;

        if matching_type == 1
            sig_min = find_sigma_min(mu_t, ploss_us, theta_plus);
            if sig_guess < sig_min + 1e-4
                sig_guess = sig_min + max(0.0001, 2 * sig_min);
            end
            [sig_out, ~, ef] = fsolve(res_fun, sig_guess, ...
                optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12));
            if ef <= 0 || sig_out < sig_min
                sq_fun = @(s) (Chi_p_psi(mu_t, ploss_us, s, iota_us_loc, ...
                               lam, et, 1) * abs_sc - TED_tgt)^2;
                sig_out = fminbnd(sq_fun, sig_min + 1e-6, 15);
            end
        else
            [sig_out, ~, ~] = fsolve(res_fun, sig_guess, ...
                optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12));
        end
        sig_guess = sig_out;
        sigma_t(tt) = sig_out;

        [~, ~, ~, ~, DW_t(tt), FF_t(tt)] = ...
            Chi_sys(mu_t, ploss_us, sig_out, iota_us_loc, lam, et, matching_type);
    end
end
