%% Robustness: optimize (lambda, eta) for fixed iota_ss values
%  Run after workspace is set up (load_data, calibration, params).
%
%  Saves results to data/estimation_robustness.mat

fprintf('\n====== Robustness: (lambda, eta) | fixed iota_ss ======\n');

addpath('functions'); addpath('functions/chi'); addpath('data');
load('data/LFX_data.mat'); load('data/calibration.mat');
pi_eu_ss = 1; pi_us_ss = 1; matching_type = 1;
run('functions/params.m');

% Valid periods
valid = ~isnan(FF_n) & ~isnan(DW_n) & FF_n > 0 & DW_n > 0;
idx = find(valid);
N_obs = length(idx);
abs_sc = 12e4;
T_f = length(mu_us);

fprintf('Using %d valid periods.\n', N_obs);

%% Grid of iota_ss values
iota_grid = 0.03:0.01:0.12;
N_grid = length(iota_grid);

% Pre-allocate result arrays
res_lambda  = zeros(N_grid, 1);
res_eta     = zeros(N_grid, 1);
res_iota    = iota_grid(:);
res_obj     = zeros(N_grid, 1);
res_exitflag = zeros(N_grid, 1);
res_mean_FF_model = zeros(N_grid, 1);
res_mean_DW_model = zeros(N_grid, 1);
res_mean_FF_data  = mean(FF_n(idx));
res_mean_DW_data  = mean(DW_n(idx));
res_mean_logFF_model = zeros(N_grid, 1);
res_mean_logDW_model = zeros(N_grid, 1);
res_mean_logFF_data  = mean(log(FF_n(idx)));
res_mean_logDW_data  = mean(log(DW_n(idx)));
% Sigma moments
res_sigma_mean = zeros(N_grid, 1);
res_sigma_std  = zeros(N_grid, 1);
res_sigma_min  = zeros(N_grid, 1);
res_sigma_max  = zeros(N_grid, 1);
res_sigma_p05  = zeros(N_grid, 1);
res_sigma_p95  = zeros(N_grid, 1);

options = optimset('Display', 'iter', 'MaxFunEvals', 5000, 'MaxIter', 3000, ...
                   'TolX', 1e-6, 'TolFun', 1e-10);

for ii = 1:N_grid
    iot_fixed = iota_grid(ii);
    fprintf('\n--- iota_ss = %.4f (%d/%d) ---\n', iot_fixed, ii, N_grid);

    % Objective: optimize over [lambda, eta] only
    obj_fun = @(x) obj_2d(x, iot_fixed, mu_us, ploss_us, TED_s_us_t, ...
        pi_us_ss, freq, matching_type, sigma_us, FF_n, DW_n, idx, T_f, abs_sc);

    x0 = [1.0, 0.65];  % [lambda, eta]
    [x_opt, fval, ef] = fminsearch(obj_fun, x0, options);

    res_lambda(ii)  = x_opt(1);
    res_eta(ii)     = x_opt(2);
    res_obj(ii)     = fval;
    res_exitflag(ii) = ef;

    % Re-run filter at optimum to get moments
    [sigma_t, FF_t, DW_t] = run_filter_local(x_opt(1), x_opt(2), iot_fixed, ...
        mu_us, ploss_us, TED_s_us_t, pi_us_ss, freq, matching_type, sigma_us, T_f, abs_sc);

    res_mean_FF_model(ii) = mean(FF_t(idx));
    res_mean_DW_model(ii) = mean(DW_t(idx));
    res_mean_logFF_model(ii) = mean(log(FF_t(idx)));
    res_mean_logDW_model(ii) = mean(log(DW_t(idx)));

    % Sigma moments (full sample)
    res_sigma_mean(ii) = mean(sigma_t);
    res_sigma_std(ii)  = std(sigma_t);
    res_sigma_min(ii)  = min(sigma_t);
    res_sigma_max(ii)  = max(sigma_t);
    res_sigma_p05(ii)  = prctile(sigma_t, 5);
    res_sigma_p95(ii)  = prctile(sigma_t, 95);

    fprintf('  lambda=%.4f, eta=%.4f, obj=%.2f, exit=%d\n', x_opt(1), x_opt(2), fval, ef);
    fprintf('  sigma: mean=%.4f, std=%.4f, [%.4f, %.4f]\n', ...
        res_sigma_mean(ii), res_sigma_std(ii), res_sigma_min(ii), res_sigma_max(ii));
end

%% Print summary table
fprintf('\n====== Summary Table ======\n');
fprintf('%8s %8s %8s %10s %4s | %8s %8s %8s %8s | %8s %8s %8s\n', ...
    'iota_ss', 'lambda', 'eta', 'Objective', 'ef', ...
    'logFF_m', 'logFF_d', 'logDW_m', 'logDW_d', ...
    'sig_mean', 'sig_std', 'sig_p95');
fprintf('%s\n', repmat('-', 1, 115));
for ii = 1:N_grid
    fprintf('%8.4f %8.4f %8.4f %10.2f %4d | %8.4f %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f\n', ...
        res_iota(ii), res_lambda(ii), res_eta(ii), res_obj(ii), res_exitflag(ii), ...
        res_mean_logFF_model(ii), res_mean_logFF_data, ...
        res_mean_logDW_model(ii), res_mean_logDW_data, ...
        res_sigma_mean(ii), res_sigma_std(ii), res_sigma_p95(ii));
end
fprintf('================================================\n');

%% Save to .mat
save('data/estimation_robustness.mat', ...
    'iota_grid', 'res_lambda', 'res_eta', 'res_iota', 'res_obj', 'res_exitflag', ...
    'res_mean_FF_model', 'res_mean_DW_model', 'res_mean_FF_data', 'res_mean_DW_data', ...
    'res_mean_logFF_model', 'res_mean_logDW_model', 'res_mean_logFF_data', 'res_mean_logDW_data', ...
    'res_sigma_mean', 'res_sigma_std', 'res_sigma_min', 'res_sigma_max', 'res_sigma_p05', 'res_sigma_p95');
fprintf('Results saved to data/estimation_robustness.mat\n');

%% ===== Local functions =====

function obj = obj_2d(x, iot, mu_us, ploss_us, TED_s_us_t, ...
        pi_us_ss, freq, matching_type, sigma_us_init, FF_n, DW_n, idx, T_f, abs_sc)
    lam = x(1);
    et  = x(2);

    % Bounds: lambda > 0, 0.1 < eta < 0.9
    pen = 0;
    if lam <= 0,   pen = pen + 1e6 * (1 + abs(lam));     end
    if et  <= 0.1, pen = pen + 1e6 * (1 + abs(0.1 - et)); end
    if et  >= 0.9, pen = pen + 1e6 * (1 + abs(et - 0.9)); end
    if pen > 0, obj = pen; return; end

    [~, FF_m, DW_m] = run_filter_local(lam, et, iot, ...
        mu_us, ploss_us, TED_s_us_t, pi_us_ss, freq, matching_type, sigma_us_init, T_f, abs_sc);

    if any(FF_m(idx) <= 0) || any(DW_m(idx) <= 0)
        obj = 1e8;
        return;
    end

    obj = sum((log(FF_m(idx)) - log(FF_n(idx))).^2) + ...
          sum((log(DW_m(idx)) - log(DW_n(idx))).^2);
end

function [sigma_t, FF_t, DW_t] = run_filter_local(lam, et, iot, ...
        mu_us, ploss_us, TED_s_us_t, pi_us_ss, freq, matching_type, sigma_us_init, T_f, abs_sc)
    iota_us_loc = iot / freq / pi_us_ss;
    sigma_t = zeros(T_f, 1);
    FF_t    = zeros(T_f, 1);
    DW_t    = zeros(T_f, 1);
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
        sigma_t(tt) = sig_out;
        [~, ~, ~, ~, DW_t(tt), FF_t(tt)] = Chi_sys(mu_t, ploss_us, sig_out, iota_us_loc, lam, et, matching_type);
    end
end
