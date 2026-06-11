%% probe_ebp.m
% Sweep the EBP calibration target over [200, 150, 100, 50, 25] bps and
% report SS calibration moments. Uses warm-start continuation: solve at
% 200 first, then use that x as the starting point for 150, then 100, ...
% This keeps fsolve on the real-valued branch as long as a real solution
% exists.
%
% Skips the expensive global Markov-switching solve.

clear; close all;

% --- Mimic main_LFX.m setup (lines 1-65) ---
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

% --- Solve baseline SS for starting point ---
[mu_us_ss, ~, ~, ~] = mu_us_star_f();
mu_eu_ss = mu_eu_ame(mu_us_ss);
Rd_us_ss = Rd_us_f(mu_us_ss);
Rd_eu_ss = Rd_eu_f(mu_us_ss);
Rb_us_ss = Rb_us_f(mu_us_ss);

load exchange_rate_data.mat ln_eu_us_ss;
load LFX_data2.mat mu_us;
mu_us_data = mu_us;  % keep separate to avoid name collision with x(2)

% --- Sweep targets, warm-start continuation ---
sweep = [200, 150, 100, 50, 25];  % descend from real-valued baseline
nS = numel(sweep);

results = struct();
results.target_bps = sweep(:);
results.is_real    = false(nS,1);
results.exitflag   = nan(nS,1);
results.fnorm      = nan(nS,1);
results.sigma_us   = nan(nS,1);
results.sigma_eu   = nan(nS,1);
results.Theta_d    = nan(nS,1);
results.mu_us      = nan(nS,1);
results.mu_eu      = nan(nS,1);
results.EBP_bps    = nan(nS,1);
results.LP_bps     = nan(nS,1);
results.CIP_bps    = nan(nS,1);
results.log_e      = nan(nS,1);

x0 = [mu_eu_ss; mu_us_ss; Rd_eu_ss; Rd_us_ss; Rb_us_ss; ...
      bard_us; bard_eu/bard_us; sigma_eu; sigma_us; Theta_d_us];

opts = optimoptions('fsolve', 'Display', 'off', ...
                    'TolFun', 1e-12, ...
                    'MaxFunctionEvaluations', 1e6, ...
                    'MaxIterations', 1e4);

x_warm = x0;  % start at the LFX_params_cd_v3 defaults

for k = 1:nS
    tgt_bps = sweep(k);
    target  = [tgt_bps; ln_eu_us_ss; mean(mu_us_data)];

    fprintf('\n===== Target EBP = %g bps =====\n', tgt_bps);
    fprintf('  Starting fsolve from %s start...\n', ...
            iif(k==1, 'cold (params defaults)', 'warm (prev solution)'));

    [x, fval, exitflag, ~] = fsolve(...
        @(x) feqm_calibrate(x, Echi_d, Echi_m, Rm_eu, Rm_us, ...
                            ploss_eu, ploss_us, iota_eu, iota_us, ...
                            lambda_eu, lambda_us, Theta_b, epsilon_b, ...
                            Theta_d_eu, Theta_d_us, zeta_eu, zeta_us, ...
                            M_eu, M_us, target, share), ...
        x_warm, opts);

    is_real = isreal(x);
    if is_real
        x_warm = x;  % advance warm start
    end

    % Unpack
    mu_eu_sol    = x(1);
    mu_us_sol    = x(2);
    Rd_eu_sol    = x(3);
    Rd_us_sol    = x(4);
    Rb_us_sol    = x(5);
    d_us_sol     = x(6);
    sigma_eu_sol = x(8);
    sigma_us_sol = x(9);
    Theta_d_sol  = x(10);

    % Reconstruct equilibrium prices
    d_eu_sol  = Theta_d_sol * (Rd_eu_sol)^(1/zeta_eu);
    p_eu_sol  = M_eu / (d_eu_sol * mu_eu_sol);
    inv_e_sol = M_us / (d_us_sol * mu_us_sol) / p_eu_sol;
    e_euus    = 1/inv_e_sol;

    Rm_us_eff = Rm_us;
    Rm_eu_eff = Rm_eu - im_eu_adj;

    EBP_bps_mdl = (Rb_us_sol - Rm_us_eff) * 1e4 * 12;
    LP_bps_mdl  = (Rm_eu_eff - Rm_us_eff) * 1e4 * 12;
    CIP_bps_mdl = (Rm_eu_eff^12 - Rm_us_eff^12) * 1e4;

    % Store (take real part for printing; flag if complex)
    results.is_real(k)  = is_real;
    results.exitflag(k) = exitflag;
    results.fnorm(k)    = norm(fval);
    results.sigma_us(k) = real(sigma_us_sol);
    results.sigma_eu(k) = real(sigma_eu_sol);
    results.Theta_d(k)  = real(Theta_d_sol);
    results.mu_us(k)    = real(mu_us_sol);
    results.mu_eu(k)    = real(mu_eu_sol);
    results.EBP_bps(k)  = real(EBP_bps_mdl);
    results.LP_bps(k)   = real(LP_bps_mdl);
    results.CIP_bps(k)  = real(CIP_bps_mdl);
    results.log_e(k)    = real(log(e_euus));

    fprintf('  is_real=%d, exitflag=%d, |fval|=%.2e\n', is_real, exitflag, norm(fval));
    fprintf('  sigma_us=%.4f, sigma_eu=%.4f, Theta_d=%.4f\n', ...
            real(sigma_us_sol), real(sigma_eu_sol), real(Theta_d_sol));
    fprintf('  EBP=%.2f bps, LP=%.2f bps, CIP=%.2f bps, log(e)=%+.4f\n', ...
            real(EBP_bps_mdl), real(LP_bps_mdl), real(CIP_bps_mdl), real(log(e_euus)));
    if ~is_real
        fprintf('  ** WARNING: complex solution -- model has no real-valued SS at this target\n');
        fprintf('  ** max(|imag|): mu_us=%.2e  Rb_us=%.2e  sigma_us=%.2e  EBP=%.2e\n', ...
                abs(imag(mu_us_sol)), abs(imag(Rb_us_sol)), ...
                abs(imag(sigma_us_sol)), abs(imag(EBP_bps_mdl)));
    end
end

%% Print summary table
fprintf('\n\n=========== SS Calibration Sweep: EBP Target ===========\n');
fprintf(' Target |Real?| Exit | |F|       | sigma_us | sigma_eu | Theta_d  | mu_us   | mu_eu   | EBP_mdl| LP    | CIP   | log_e\n');
fprintf('--------+-----+------+-----------+----------+----------+----------+---------+---------+--------+-------+-------+--------\n');
for k = 1:nS
    fprintf(' %5d  |  %d  |  %2d  | %.3e | %8.4f | %8.4f | %8.4f | %7.4f | %7.4f | %6.2f | %5.2f | %5.2f | %+.4f\n', ...
        sweep(k), results.is_real(k), results.exitflag(k), results.fnorm(k), ...
        results.sigma_us(k), results.sigma_eu(k), results.Theta_d(k), ...
        results.mu_us(k), results.mu_eu(k), ...
        results.EBP_bps(k), results.LP_bps(k), results.CIP_bps(k), results.log_e(k));
end
fprintf('=========================================================\n');

save('probe_ebp_results.mat', 'results');
disp('=== probe_ebp complete ===');

% Inline helper since older MATLAB doesn't have a built-in
function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end
