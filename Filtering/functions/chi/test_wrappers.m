%% Test Wrapper Functions
% Tests Echi_m, Echi_d, Chi_p_psi, Chi_sys against expected behavior
clear; clc;

%% Parameters (typical filter values)
mu = 0.15;          % Reserve ratio
ploss = 0.75;       % Probability of loss
sigma = 2.0;        % Volatility
iota = 0.0012;      % 12 bps spread
lambda = 3.5;       % Matching efficiency
eta = 0.5;          % Bargaining power

fprintf('=== Test Parameters ===\n');
fprintf('mu = %.2f, ploss = %.2f, sigma = %.2f\n', mu, ploss, sigma);
fprintf('iota = %.4f, lambda = %.1f, eta = %.1f\n\n', iota, lambda, eta);

%% Test 1: Compare Leontief (matching_type = 0)
fprintf('=== LEONTIEF (matching_type = 0) ===\n\n');

% Call wrapper functions
echi_m_val = Echi_m(mu, ploss, sigma, iota, lambda, eta, 0);
echi_d_val = Echi_d(mu, ploss, sigma, iota, lambda, eta, 0);
[ted_val, chi_p_val, psi_p_val] = Chi_p_psi(mu, ploss, sigma, iota, lambda, eta, 0);
[sys_echi_m, theta_val, psi_val, Smin_val, DW_val, FF_val, Q_val] = Chi_sys(mu, ploss, sigma, iota, lambda, eta, 0);

fprintf('Echi_m (bond premium):    %.6f (%.2f bps annualized)\n', echi_m_val, echi_m_val * 1e4 * 12);
fprintf('Echi_d (deposit premium): %.6f (%.2f bps annualized)\n', echi_d_val, echi_d_val * 1e4 * 12);
fprintf('Chi_p_psi (TED spread):   %.6f (%.2f bps annualized)\n', ted_val, ted_val * 1e4 * 12);
fprintf('\n');
fprintf('Chi_sys outputs:\n');
fprintf('  theta (tightness):      %.4f\n', theta_val);
fprintf('  psi (match prob):       %.4f\n', psi_val);
fprintf('  Smin (deficit):         %.4f\n', Smin_val);
fprintf('  DW (discount window):   %.4f\n', DW_val);
fprintf('  FF (fed funds vol):     %.4f\n', FF_val);
fprintf('\n');

% Consistency check
fprintf('Consistency check: Echi_m from wrapper vs Chi_sys: %s\n', ...
    string(abs(echi_m_val - sys_echi_m) < 1e-12));

%% Test 2: Compare Cobb-Douglas (matching_type = 1)
fprintf('\n=== COBB-DOUGLAS (matching_type = 1) ===\n\n');

echi_m_cd = Echi_m(mu, ploss, sigma, iota, lambda, eta, 1);
echi_d_cd = Echi_d(mu, ploss, sigma, iota, lambda, eta, 1);
[ted_cd, chi_p_cd, psi_p_cd] = Chi_p_psi(mu, ploss, sigma, iota, lambda, eta, 1);
[sys_echi_m_cd, theta_cd, psi_cd, Smin_cd, DW_cd, FF_cd, Q_cd] = Chi_sys(mu, ploss, sigma, iota, lambda, eta, 1);

fprintf('Echi_m (bond premium):    %.6f (%.2f bps annualized)\n', echi_m_cd, echi_m_cd * 1e4 * 12);
fprintf('Echi_d (deposit premium): %.6f (%.2f bps annualized)\n', echi_d_cd, echi_d_cd * 1e4 * 12);
fprintf('Chi_p_psi (TED spread):   %.6f (%.2f bps annualized)\n', ted_cd, ted_cd * 1e4 * 12);
fprintf('\n');
fprintf('Chi_sys outputs:\n');
fprintf('  theta (tightness):      %.4f\n', theta_cd);
fprintf('  psi (match prob):       %.4f\n', psi_cd);
fprintf('  Smin (deficit):         %.4f\n', Smin_cd);
fprintf('  DW (discount window):   %.4f\n', DW_cd);
fprintf('  FF (fed funds vol):     %.4f\n', FF_cd);
fprintf('\n');

%% Test 3: Comparison table
fprintf('\n=== COMPARISON: Leontief vs Cobb-Douglas ===\n\n');
fprintf('Variable              Leontief    Cobb-Douglas    Ratio\n');
fprintf('--------              --------    ------------    -----\n');
fprintf('Echi_m (bps)          %8.2f    %12.2f    %.3f\n', echi_m_val*1e4*12, echi_m_cd*1e4*12, echi_m_cd/echi_m_val);
fprintf('Echi_d (bps)          %8.2f    %12.2f    %.3f\n', echi_d_val*1e4*12, echi_d_cd*1e4*12, echi_d_cd/echi_d_val);
fprintf('TED (bps)             %8.2f    %12.2f    %.3f\n', ted_val*1e4*12, ted_cd*1e4*12, ted_cd/ted_val);
fprintf('Psi (match prob)      %8.4f    %12.4f    %.3f\n', psi_val, psi_cd, psi_cd/psi_val);
fprintf('DW usage              %8.4f    %12.4f    %.3f\n', DW_val, DW_cd, DW_cd/DW_val);
fprintf('FF volume             %8.4f    %12.4f    %.3f\n', FF_val, FF_cd, FF_cd/FF_val);

%% Test 4: Default arguments (should use Leontief, eta=0.5)
fprintf('\n=== TESTING DEFAULTS ===\n');
echi_m_default = Echi_m(mu, ploss, sigma, iota, lambda);
echi_m_explicit = Echi_m(mu, ploss, sigma, iota, lambda, 0.5, 0);
fprintf('Default call:   %.8f\n', echi_m_default);
fprintf('Explicit call:  %.8f\n', echi_m_explicit);
fprintf('Match: %s\n', string(abs(echi_m_default - echi_m_explicit) < 1e-12));

%% Test 5: Vector inputs
fprintf('\n=== TESTING VECTOR INPUTS ===\n');
mu_vec = [0.10, 0.15, 0.20, 0.25];
echi_m_vec = Echi_m(mu_vec, ploss, sigma, iota, lambda, eta, 0);
fprintf('mu values:    ');
fprintf('%.2f  ', mu_vec);
fprintf('\n');
fprintf('Echi_m (bps): ');
fprintf('%.2f  ', echi_m_vec * 1e4 * 12);
fprintf('\n');

fprintf('\nAll tests complete.\n');
