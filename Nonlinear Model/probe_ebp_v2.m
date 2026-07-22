%% probe_ebp_v2.m
% Sweep EBP target with the NEW fixed-sigma setup (sigma fixed at Markov
% ergodic mean, Theta_d free) at the current lambda = 1.2778.
% Warm-start continuation to keep on the real branch.

clear; close all;

T = 5;
freq = 12;
rate_scale = 1e4;
LFX_plotprefs;

load 'Z_im_us.mat';
load 'Z_im_eu.mat';

pi_eu_ss = 1;
pi_us_ss = 1;
load dynare_calibration_param.mat;
LFX_params_cd_v3;

Theta_b   = 1;
Theta_d_eu = 1;
Theta_d_us = 1;
epsilon_b = -0.001;
zeta_us   = 1000;
zeta_eu   = 1000;
im_eu_adj = 0.0006;
share     = 1;

LFX_nt_0e_eqs_2;

% Baseline SS (needed for x0)
[mu_us_ss, ~, ~, ~] = mu_us_star_f();
mu_eu_ss = mu_eu_ame(mu_us_ss);
Rd_us_ss = Rd_us_f(mu_us_ss);
Rd_eu_ss = Rd_eu_f(mu_us_ss);
Rb_us_ss = Rb_us_f(mu_us_ss);

load exchange_rate_data.mat ln_eu_us_ss;
load LFX_data2.mat mu_us;
mu_us_data_mean = mean(mu_us);
clear mu_us;

% Markov ergodic mean of sigma_us
mu_r1_erg   = -0.7370; rho_r1_erg = 0.9723; Sig_r1_erg = 0.0445;
mu_r2_erg   = -0.4914; rho_r2_erg = 0.1477; Sig_r2_erg = 0.6918;
P_regime    = [0.9925 0.0075; 0.0979 0.9021];
E_sig_r1    = exp(mu_r1_erg + 0.5 * (Sig_r1_erg / sqrt(1 - rho_r1_erg^2))^2);
E_sig_r2    = exp(mu_r2_erg + 0.5 * (Sig_r2_erg / sqrt(1 - rho_r2_erg^2))^2);
pi_r1       = P_regime(2,1) / (P_regime(2,1) + P_regime(1,2));
pi_r2       = P_regime(1,2) / (P_regime(2,1) + P_regime(1,2));
sigma_us_erg = pi_r1 * E_sig_r1 + pi_r2 * E_sig_r2;
sigma_eu_erg = sigma_us_erg;

fprintf('\nMarkov ergodic mean sigma_us = %.4f\n\n', sigma_us_erg);

% Sweep from 200 (known feasible) downward
sweep = [200, 150, 100, 50, 24];
nS = numel(sweep);
results = struct();
results.target_bps = sweep(:);
results.is_real    = false(nS,1);
results.exitflag   = nan(nS,1);
results.fnorm      = nan(nS,1);
results.Theta_d    = nan(nS,1);
results.mu_us      = nan(nS,1);
results.mu_eu      = nan(nS,1);
results.EBP_bps    = nan(nS,1);
results.LP_bps     = nan(nS,1);
results.log_e      = nan(nS,1);

x_warm = [mu_eu_ss; mu_us_ss; Rd_eu_ss; Rd_us_ss; Rb_us_ss; bard_us; bard_eu/bard_us; Theta_d_us];

opts = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, ...
                    'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e4);

for k = 1:nS
    tgt_bps = sweep(k);
    fprintf('===== Target EBP = %g bps =====\n', tgt_bps);
    [x, fval, exitflag, ~] = fsolve(@(x) ...
        feqm_calibrate_fixed_sigma(x, Echi_d, Echi_m, Rm_eu, Rm_us, ...
            ploss_eu, ploss_us, iota_eu, iota_us, lambda_eu, lambda_us, ...
            Theta_b, epsilon_b, zeta_eu, zeta_us, M_eu, M_us, ...
            sigma_us_erg, sigma_eu_erg, tgt_bps, share), ...
        x_warm, opts);

    is_real = isreal(x);
    if is_real && exitflag > 0
        x_warm = x;   % advance warm start only on clean success
    end

    Theta_d_s = x(8);
    mu_eu_s   = x(1);
    mu_us_s   = x(2);
    Rb_us_s   = x(5);
    d_us_s    = x(6);
    d_eu_s    = Theta_d_s * (x(3))^(1/zeta_eu);
    p_eu_s    = M_eu / (d_eu_s * mu_eu_s);
    inv_e_s   = M_us / (d_us_s * mu_us_s) / p_eu_s;
    e_euus_s  = 1 / inv_e_s;

    Rm_us_eff = Rm_us;
    Rm_eu_eff = Rm_eu - im_eu_adj;
    EBP_bps   = (Rb_us_s - Rm_us_eff) * 1e4 * 12;
    LP_bps    = (Rm_eu_eff - Rm_us_eff) * 1e4 * 12;

    results.is_real(k)  = is_real;
    results.exitflag(k) = exitflag;
    results.fnorm(k)    = norm(fval);
    results.Theta_d(k)  = real(Theta_d_s);
    results.mu_us(k)    = real(mu_us_s);
    results.mu_eu(k)    = real(mu_eu_s);
    results.EBP_bps(k)  = real(EBP_bps);
    results.LP_bps(k)   = real(LP_bps);
    results.log_e(k)    = real(log(e_euus_s));

    fprintf('  is_real=%d, exitflag=%d, |fval|=%.2e\n', is_real, exitflag, norm(fval));
    fprintf('  Theta_d=%.4f, mu_us=%.4f, mu_eu=%.4f\n', real(Theta_d_s), real(mu_us_s), real(mu_eu_s));
    fprintf('  EBP=%.2f bps [target %g]  |  log(e)=%+.4f (data %+.4f)\n\n', ...
        real(EBP_bps), tgt_bps, real(log(e_euus_s)), ln_eu_us_ss);
end

fprintf('\n=========== SS Calibration Sweep (fixed sigma, calibrate Theta_d) ===========\n');
fprintf(' Target | Real? | Exit | |F|       | Theta_d  | mu_us   | mu_eu   | EBP    | LP    | log(e)\n');
fprintf('--------+-------+------+-----------+----------+---------+---------+--------+-------+--------\n');
for k = 1:nS
    fprintf(' %5d  |  %d    |  %2d  | %.3e | %8.4f | %7.4f | %7.4f | %6.2f | %5.2f | %+.4f\n', ...
        sweep(k), results.is_real(k), results.exitflag(k), results.fnorm(k), ...
        results.Theta_d(k), results.mu_us(k), results.mu_eu(k), ...
        results.EBP_bps(k), results.LP_bps(k), results.log_e(k));
end
fprintf('===============================================================================\n');
fprintf('log(e_euus) data reference = %+.4f\n', ln_eu_us_ss);
fprintf('log(mu_us)  data reference = %+.4f\n', mu_us_data_mean);
fprintf('=== probe_ebp_v2 complete ===\n');
