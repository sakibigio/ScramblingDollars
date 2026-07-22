%% test_ss_only.m
% Runs SS calibration only (skips global Markov-switching solve).
% Uses feqm_calibrate_fixed_sigma with BOTH Theta_b and Theta_d calibrated:
%   Theta_b -> beta-fix (Rb-Rm = 200 bps)
%   Theta_d -> log(mu_us) = data mean

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

Theta_b   = 1;         % starting guess only; will be calibrated
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

% ---- Markov ergodic mean of sigma_us ----
mu_r1_erg   = -0.7370; rho_r1_erg = 0.9723; Sig_r1_erg = 0.0445;
mu_r2_erg   = -0.4914; rho_r2_erg = 0.1477; Sig_r2_erg = 0.6918;
P_regime    = [0.9925 0.0075; 0.0979 0.9021];
E_sig_r1    = exp(mu_r1_erg + 0.5 * (Sig_r1_erg / sqrt(1 - rho_r1_erg^2))^2);
E_sig_r2    = exp(mu_r2_erg + 0.5 * (Sig_r2_erg / sqrt(1 - rho_r2_erg^2))^2);
pi_r1       = P_regime(2,1) / (P_regime(2,1) + P_regime(1,2));
pi_r2       = P_regime(1,2) / (P_regime(2,1) + P_regime(1,2));
sigma_us_erg = pi_r1 * E_sig_r1 + pi_r2 * E_sig_r2;
sigma_eu_erg = sigma_us_erg;

fprintf('\n=== Markov ergodic decomposition ===\n');
fprintf('  E[sigma | r1] = %.4f, E[sigma | r2] = %.4f\n', E_sig_r1, E_sig_r2);
fprintf('  pi_r1 = %.4f, pi_r2 = %.4f\n', pi_r1, pi_r2);
fprintf('  Markov ergodic mean sigma_us = %.4f\n\n', sigma_us_erg);

target_ebp_bps    = 200;

% x0: 8 elements — sigma_us and sigma_eu are fixed inputs, not unknowns
x0 = [mu_eu_ss; mu_us_ss; Rd_eu_ss; Rd_us_ss; Rb_us_ss; bard_us; bard_eu/bard_us; Theta_d_us];

fprintf('=== SS calibration (sigma_us=sigma_eu fixed at Markov mean, calibrate Theta_d) ===\n');
[x_ss, fval, exitflag, ~] = fsolve(@(x) ...
    feqm_calibrate_fixed_sigma(x, Echi_d, Echi_m, Rm_eu, Rm_us, ...
        ploss_eu, ploss_us, iota_eu, iota_us, lambda_eu, lambda_us, ...
        Theta_b, epsilon_b, zeta_eu, zeta_us, M_eu, M_us, ...
        sigma_us_erg, sigma_eu_erg, target_ebp_bps, share), ...
    x0, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, ...
                     'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e4));

fprintf('  exitflag = %d, |fval| = %.3e, isreal = %d\n\n', exitflag, norm(fval), isreal(x_ss));

% Unpack (8 elements)
mu_eu_s    = x_ss(1);
mu_us_s    = x_ss(2);
Rd_eu_s    = x_ss(3);
Rd_us_s    = x_ss(4);
Rb_us_s    = x_ss(5);
d_us_s     = x_ss(6);
nu_s       = x_ss(7);
Theta_d_s  = x_ss(8);
sigma_eu_s = sigma_eu_erg;

d_eu_s   = Theta_d_s * (Rd_eu_s)^(1/zeta_eu);
p_eu_s   = M_eu / (d_eu_s * mu_eu_s);
inv_e_s  = M_us / (d_us_s * mu_us_s) / p_eu_s;
e_euus_s = 1 / inv_e_s;

Rm_us_eff = Rm_us;
Rm_eu_eff = Rm_eu - im_eu_adj;

EBP_bps = (Rb_us_s - Rm_us_eff) * 1e4 * 12;
LP_bps  = (Rm_eu_eff - Rm_us_eff) * 1e4 * 12;
CIP_bps = (Rm_eu_eff^12 - Rm_us_eff^12) * 1e4;

fprintf('=== SS Calibration Results ===\n');
fprintf('  sigma_us    = %.4f  (fixed at Markov ergodic mean)\n', sigma_us_erg);
fprintf('  sigma_eu    = %.4f  (calibrated)\n', sigma_eu_s);
fprintf('  Theta_d     = %.4f  (calibrated)\n', Theta_d_s);
fprintf('  mu_us       = %.4f\n', mu_us_s);
fprintf('  mu_eu       = %.4f\n', mu_eu_s);

fprintf('\n=== Moment fit ===\n');
fprintf('  EBP  (target 200 bps)                : %.2f bps  [diff %+.4f]\n', EBP_bps, EBP_bps - target_ebp_bps);
fprintf('  LP   (rate-implied, not targeted)    : %.2f bps\n', LP_bps);
fprintf('  CIP  (rate-implied, not targeted)    : %.2f bps\n', CIP_bps);
fprintf('  log(e_euus) (DIAGNOSTIC): model = %+.4f  |  data = %+.4f  |  diff = %+.4f\n', log(e_euus_s), ln_eu_us_ss, log(e_euus_s) - ln_eu_us_ss);
fprintf('  log(mu_us)  (DIAGNOSTIC): model = %+.4f  |  data = %+.4f  |  diff = %+.4f\n', log(mu_us_s), mu_us_data_mean, log(mu_us_s) - mu_us_data_mean);

fprintf('\n=== test_ss_only complete ===\n');
