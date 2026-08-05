%% Compare current calibration vs estimated parameters
cd(fileparts(mfilename('fullpath')));   % this file's own folder
addpath('functions'); addpath('functions/chi'); addpath('data');
load('data/LFX_data.mat'); load('data/calibration.mat');
pi_eu_ss = 1; pi_us_ss = 1; matching_type = 1;
run('functions/params.m');

% Valid periods
valid = ~isnan(FF_n) & ~isnan(DW_n) & FF_n > 0 & DW_n > 0;
idx = find(valid);
abs_sc = 12e4;
T_f = length(mu_us);

%% --- Helper: run filter for given params ---
run_filt = @(lam, et, iot) local_filter(lam, et, iot, mu_us, ploss_us, TED_s_us_t, pi_us_ss, freq, matching_type, sigma_us, T_f, abs_sc);

%% Current calibration
[FF_c, DW_c] = run_filt(lambda_us, eta, iota_ss);
obj_c = sum((log(FF_c(idx)) - log(FF_n(idx))).^2) + sum((log(DW_c(idx)) - log(DW_n(idx))).^2);

fprintf('\n=== Current Calibration: lambda=%.2f, eta=%.2f, iota_ss=%.4f ===\n', lambda_us, eta, iota_ss);
fprintf('  Objective (log SSE)  = %.4f\n', obj_c);
fprintf('  Mean FF  (model)     = %.6f    (data) = %.6f\n', mean(FF_c(idx)), mean(FF_n(idx)));
fprintf('  Mean DW  (model)     = %.6f    (data) = %.6f\n', mean(DW_c(idx)), mean(DW_n(idx)));
fprintf('  Mean log(FF) model   = %.4f    data   = %.4f\n', mean(log(FF_c(idx))), mean(log(FF_n(idx))));
fprintf('  Mean log(DW) model   = %.4f    data   = %.4f\n', mean(log(DW_c(idx))), mean(log(DW_n(idx))));
fprintf('  Std  log(FF) model   = %.4f    data   = %.4f\n', std(log(FF_c(idx))), std(log(FF_n(idx))));
fprintf('  Std  log(DW) model   = %.4f    data   = %.4f\n', std(log(DW_c(idx))), std(log(DW_n(idx))));
fprintf('  Corr log(FF) m vs d  = %.4f\n', corr(log(FF_c(idx)), log(FF_n(idx))));
fprintf('  Corr log(DW) m vs d  = %.4f\n', corr(log(DW_c(idx)), log(DW_n(idx))));

%% Estimated
lam_e = 0.917868; eta_e = 0.715252; iot_e = 0.100000;
[FF_e, DW_e] = run_filt(lam_e, eta_e, iot_e);
obj_e = sum((log(FF_e(idx)) - log(FF_n(idx))).^2) + sum((log(DW_e(idx)) - log(DW_n(idx))).^2);

fprintf('\n=== Estimated: lambda=%.4f, eta=%.4f, iota_ss=%.4f ===\n', lam_e, eta_e, iot_e);
fprintf('  Objective (log SSE)  = %.4f\n', obj_e);
fprintf('  Mean FF  (model)     = %.6f    (data) = %.6f\n', mean(FF_e(idx)), mean(FF_n(idx)));
fprintf('  Mean DW  (model)     = %.6f    (data) = %.6f\n', mean(DW_e(idx)), mean(DW_n(idx)));
fprintf('  Mean log(FF) model   = %.4f    data   = %.4f\n', mean(log(FF_e(idx))), mean(log(FF_n(idx))));
fprintf('  Mean log(DW) model   = %.4f    data   = %.4f\n', mean(log(DW_e(idx))), mean(log(DW_n(idx))));
fprintf('  Std  log(FF) model   = %.4f    data   = %.4f\n', std(log(FF_e(idx))), std(log(FF_n(idx))));
fprintf('  Std  log(DW) model   = %.4f    data   = %.4f\n', std(log(DW_e(idx))), std(log(DW_n(idx))));
fprintf('  Corr log(FF) m vs d  = %.4f\n', corr(log(FF_e(idx)), log(FF_n(idx))));
fprintf('  Corr log(DW) m vs d  = %.4f\n', corr(log(DW_e(idx)), log(DW_n(idx))));

fprintf('\n=== Improvement ===\n');
fprintf('  Objective: %.4f -> %.4f  (%.1f%% reduction)\n', obj_c, obj_e, (1 - obj_e/obj_c)*100);

%% Local filter function
function [FF_t, DW_t] = local_filter(lam, et, iot, mu_us, ploss_us, TED_s_us_t, pi_us_ss, freq, matching_type, sigma_us_init, T_f, abs_sc)
    iota_us_loc = iot / freq / pi_us_ss;
    FF_t = zeros(T_f, 1);
    DW_t = zeros(T_f, 1);
    if matching_type == 1
        theta_plus = ((exp(lam) - 1) / (exp(lam) + 1))^2;
    end
    sig_guess = sigma_us_init;
    for tt = 1:T_f
        mu_t = exp(mu_us(tt));
        TED_tgt = TED_s_us_t(tt) * abs_sc;
        res_fun = @(sig) Chi_p_psi(mu_t, ploss_us, sig, iota_us_loc, lam, et, matching_type) * abs_sc - TED_tgt;
        if matching_type == 1
            sig_min = find_sigma_min(mu_t, ploss_us, theta_plus);
            if sig_guess < sig_min + 1e-4
                sig_guess = sig_min + max(0.0001, 2 * sig_min);
            end
            [sig_out, ~, ef] = fsolve(res_fun, sig_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12));
            if ef <= 0 || sig_out < sig_min
                sq_fun = @(s) (Chi_p_psi(mu_t, ploss_us, s, iota_us_loc, lam, et, 1) * abs_sc - TED_tgt)^2;
                sig_out = fminbnd(sq_fun, sig_min + 1e-6, 15);
            end
        else
            [sig_out, ~, ~] = fsolve(res_fun, sig_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12));
        end
        sig_guess = sig_out;
        [~, ~, ~, ~, DW_t(tt), FF_t(tt)] = Chi_sys(mu_t, ploss_us, sig_out, iota_us_loc, lam, et, matching_type);
    end
end
