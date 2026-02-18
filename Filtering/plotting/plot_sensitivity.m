%% PLOT_SENSITIVITY.M
% Sensitivity analysis for Chi functions and equilibrium conditions
% Part of: Scrambling for Dollars filter pipeline
%
% Purpose:
%   Generates 3D surface plots showing how:
%   - Echi_m, Echi_d, Chi_p_psi vary with (sigma, mu)
%   - Arbitrage and budget constraint conditions interact
%   - Equilibrium responds to parameter perturbations
%
% Requires:
%   - Run main_filter.m first (populates workspace with parameters and paths)
%   - Chi functions in path: Echi_m, Echi_d, Chi_p_psi, Chi_sys
%
% Outputs:
%   - ~12 figures showing sensitivity surfaces and cross-conditions
%
% Authors: Bianchi, Bigio, Engel
% Last updated: February 2025
%==========================================================================

%% Settings
% close all;  % Uncomment to close existing figures

% Define formatting functions if not in workspace
if ~exist('formataxis', 'var')
    FSize = 20;
    formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', FSize, 'Box', 'Off', 'PlotBoxAspectRatio', [1 0.75 1]);
end
if ~exist('label_x', 'var')
    label_x = @(x) xlabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', 16, 'Interpreter', 'latex');
end
if ~exist('label_y', 'var')
    label_y = @(x) ylabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', 16, 'Interpreter', 'latex');
end

% Scale for BPS
if ~exist('abs_scale', 'var')
    abs_scale = 12e4;
end

%% Grid Parameters
t_shock_i = 93;  % Time period for shock analysis (adjust as needed)
N_sigma = 60;    % Grid points for sigma
N_mu_us = 80;    % Grid points for mu_us
N_mu_eu = 82;    % Grid points for mu_eu

% Create grids
vec_sigma_us = linspace(min(sigma_us_t), max(sigma_us_t), N_sigma);
vec_mu_us = linspace(0.1, 0.9, N_mu_us);
vec_mu_eu = linspace(0.1, 0.9, N_mu_eu);

%% =========================================================================
%  SECTION 1: Chi Function Surfaces (sigma x mu)
%  =========================================================================
fprintf('Computing Chi function surfaces...\n');

% Pre-allocate
Echi_m_mat = NaN(N_sigma, N_mu_us);
Echi_d_mat = NaN(N_sigma, N_mu_us);
Echi_p_psi_mat = NaN(N_sigma, N_mu_us);

% Compute surfaces
for ss = 1:N_sigma
    for mm = 1:N_mu_us
        mu_us_yt = vec_mu_us(mm);
        sigma_us_yt = vec_sigma_us(ss);
        Echi_m_mat(ss, mm) = Echi_m(mu_us_yt, ploss_us, sigma_us_yt, iota_us, lambda_us, eta, matching_type) * abs_scale;
        Echi_d_mat(ss, mm) = Echi_d(mu_us_yt, ploss_us, sigma_us_yt, iota_us, lambda_us, eta, matching_type) * abs_scale;
        Echi_p_psi_mat(ss, mm) = Chi_p_psi(mu_us_yt, ploss_us, sigma_us_yt, iota_us, lambda_us, eta, matching_type) * abs_scale;
    end
end

% Plot: Echi_m surface
figure('Name', 'Echi_m Sensitivity', 'NumberTitle', 'off');
surf(vec_mu_us, vec_sigma_us, Echi_m_mat);
xlabel('$\mu$ (US)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\sigma$ (US)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('BPS', 'FontSize', 12);
title('$E[\chi^m]$ Sensitivity', 'Interpreter', 'latex', 'FontSize', 18);
colorbar;

% Plot: Echi_d surface
figure('Name', 'Echi_d Sensitivity', 'NumberTitle', 'off');
surf(vec_mu_us, vec_sigma_us, Echi_d_mat);
xlabel('$\mu$ (US)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\sigma$ (US)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('BPS', 'FontSize', 12);
title('$E[\chi^d]$ Sensitivity', 'Interpreter', 'latex', 'FontSize', 18);
colorbar;

% Plot: Chi_p_psi surface
figure('Name', 'Chi_p_psi Sensitivity', 'NumberTitle', 'off');
surf(vec_mu_us, vec_sigma_us, Echi_p_psi_mat);
xlabel('$\mu$ (US)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\sigma$ (US)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('BPS', 'FontSize', 12);
title('$\chi^f$ (Fed Funds) Sensitivity', 'Interpreter', 'latex', 'FontSize', 18);
colorbar;

% Plot: Overlay comparison
figure('Name', 'Echi_m vs Chi_p_psi', 'NumberTitle', 'off');
surf(vec_mu_us, vec_sigma_us, Echi_m_mat, 'FaceAlpha', 0.5, 'FaceColor', 'b'); hold on;
surf(vec_mu_us, vec_sigma_us, Echi_p_psi_mat, 'FaceAlpha', 0.5, 'FaceColor', 'r');
xlabel('$\mu$ (US)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\sigma$ (US)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('BPS', 'FontSize', 12);
title('$E[\chi^m]$ vs $\chi^f$ Comparison', 'Interpreter', 'latex', 'FontSize', 18);
legend('$E[\chi^m]$', '$\chi^f$', 'Interpreter', 'latex', 'Location', 'best');

%% =========================================================================
%  SECTION 2: Cross-Condition Analysis at Fixed Time
%  =========================================================================
fprintf('Computing cross-conditions at t=%d...\n', t_shock_i);

% Check required variables exist
if ~exist('Rm_us', 'var') || ~exist('Theta_d_us_t', 'var')
    warning('Required variables not found. Run main_filter.m first.');
    fprintf('Skipping sensitivity analysis.\n');
    return;
end

% Extract shock values at reference time
sigma_us_yt = sigma_us_t(t_shock_i);
sigma_eu_yt = sigma_eu_t(t_shock_i);
riskprm_yt = riskprm_t(t_shock_i);

% Define rate functions
Echi_m_us_f = @(mu) Echi_m(mu, ploss_us, sigma_us_yt, iota_us, lambda_us, eta, matching_type);
Echi_d_us_f = @(mu) Echi_d(mu, ploss_us, sigma_us_yt, iota_us, lambda_us, eta, matching_type);
Echi_m_eu_f = @(mu) Echi_m(mu, ploss_eu, sigma_eu_yt, iota_eu, lambda_eu, eta, matching_type);
Echi_d_eu_f = @(mu) Echi_d(mu, ploss_eu, sigma_eu_yt, iota_eu, lambda_eu, eta, matching_type);

Rb_us_f = @(mu) Rm_us(t_shock_i) + Echi_m_us_f(mu);
Rd_us_f = @(mu) Rm_us(t_shock_i) + Echi_m_us_f(mu) + Echi_d_us_f(mu);
Rb_eu_f = @(mu) Rm_eu(t_shock_i) + Echi_m_eu_f(mu);
Rd_eu_f = @(mu) Rm_eu(t_shock_i) + Echi_m_eu_f(mu) + Echi_d_eu_f(mu);

% Quantity functions
d_us_f = @(mu) (Theta_d_us_t(t_shock_i) * Rd_us_f(mu).^(1/zeta_us) - MBS_us_ss) / (1 - share_us*mu);
d_eu_f = @(mu) (Theta_d_eu_t(t_shock_i) * Rd_eu_f(mu).^(1/zeta_eu) - MBS_eu_ss) / (1 - share_eu*mu);
b_f = @(mu) Theta_b * (Rb_us_f(mu))^(1/epsilon_b);

% Residual functions
nu_f = @(mu_us, mu_eu) d_eu_f(mu_eu) / d_us_f(mu_us);
res_arb = @(mu_us, mu_eu) (1 + riskprm_yt) * Rb_eu_f(mu_eu) - Rb_us_f(mu_us);
res_bc = @(mu_us, mu_eu) b_f(mu_us) - (nu_f(mu_us, mu_eu)*(1-mu_eu) + 1-mu_us) * d_us_f(mu_us);

%% Solve for equilibrium curves
fprintf('Solving for equilibrium curves...\n');

opts = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, ...
                    'MaxFunEvals', 1e9, 'MaxIter', 1e6);

mu_eu_arb_vec = NaN(N_mu_us, 1);
mu_eu_bc_vec = NaN(N_mu_us, 1);

for mm = 1:N_mu_us
    mu_us_yt = vec_mu_us(mm);
    
    % Arbitrage condition
    [mu_eu_out, ~, exitflag, ~] = fsolve(@(mu_eu) res_arb(mu_us_yt, exp(mu_eu)), log(mu_us_yt), opts);
    if exitflag > 0
        mu_eu_arb_vec(mm) = exp(mu_eu_out);
    end
    
    % Budget constraint
    [mu_eu_out, ~, exitflag, ~] = fsolve(@(mu_eu) res_bc(mu_us_yt, exp(mu_eu)), log(mu_us_yt), opts);
    if exitflag > 0
        mu_eu_bc_vec(mm) = exp(mu_eu_out);
    end
end

%% Compute residual surfaces
fprintf('Computing residual surfaces...\n');

arb_mat = NaN(N_mu_eu, N_mu_us);
bc_mat = NaN(N_mu_eu, N_mu_us);
nu_mat = NaN(N_mu_eu, N_mu_us);
b_mat = NaN(N_mu_eu, N_mu_us);
d_us_mat = NaN(N_mu_eu, N_mu_us);
d_eu_mat = NaN(N_mu_eu, N_mu_us);
f_mat = NaN(N_mu_eu, N_mu_us);

for nn = 1:N_mu_eu
    for mm = 1:N_mu_us
        mu_us_yt = vec_mu_us(mm);
        mu_eu_yt = vec_mu_eu(nn);
        arb_mat(nn, mm) = res_arb(mu_us_yt, mu_eu_yt);
        bc_mat(nn, mm) = res_bc(mu_us_yt, mu_eu_yt);
        nu_mat(nn, mm) = nu_f(mu_us_yt, mu_eu_yt);
        b_mat(nn, mm) = b_f(mu_us_yt);
        d_us_mat(nn, mm) = d_us_f(mu_us_yt);
        d_eu_mat(nn, mm) = d_eu_f(mu_eu_yt);
        f_mat(nn, mm) = (nu_mat(nn, mm)*(1-mu_eu_yt) + 1-mu_us_yt) * d_us_mat(nn, mm);
    end
end

%% Plot: Residual surfaces
figure('Name', 'Arbitrage Residual', 'NumberTitle', 'off');
surf(vec_mu_us, vec_mu_eu, arb_mat, 'FaceAlpha', 0.5);
xlabel('$\mu^{us}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu^{eu}$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Excess', 'FontSize', 12);
title('Arbitrage Condition Residual', 'Interpreter', 'latex', 'FontSize', 18);

figure('Name', 'Budget Constraint Residual', 'NumberTitle', 'off');
surf(vec_mu_us, vec_mu_eu, bc_mat, 'FaceAlpha', 0.5);
xlabel('$\mu^{us}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu^{eu}$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Excess', 'FontSize', 12);
title('Budget Constraint Residual', 'Interpreter', 'latex', 'FontSize', 18);

figure('Name', 'Deposits', 'NumberTitle', 'off');
surf(vec_mu_us, vec_mu_eu, d_us_mat, 'FaceAlpha', 0.5, 'FaceColor', 'b'); hold on;
surf(vec_mu_us, vec_mu_eu, d_eu_mat, 'FaceAlpha', 0.5, 'FaceColor', 'r');
xlabel('$\mu^{us}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu^{eu}$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Value', 'FontSize', 12);
title('Deposit Demands', 'Interpreter', 'latex', 'FontSize', 18);
legend('$d^{us}$', '$d^{eu}$', 'Interpreter', 'latex');

figure('Name', 'Funding vs Bonds', 'NumberTitle', 'off');
surf(vec_mu_us, vec_mu_eu, f_mat, 'FaceAlpha', 0.5, 'FaceColor', 'b'); hold on;
surf(vec_mu_us, vec_mu_eu, b_mat, 'FaceAlpha', 0.5, 'FaceColor', 'r');
xlabel('$\mu^{us}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu^{eu}$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Value', 'FontSize', 12);
title('Funding Conditions', 'Interpreter', 'latex', 'FontSize', 18);
legend('Funding', 'Bonds');

% Contour plot
figure('Name', 'Equilibrium Contours', 'NumberTitle', 'off');
contour(vec_mu_us, vec_mu_eu, arb_mat, 10, 'Color', 'k', 'ShowText', true); hold on;
contour(vec_mu_us, vec_mu_eu, bc_mat, 10, 'Color', 'r', 'ShowText', true);
xlabel('$\mu^{us}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu^{eu}$', 'Interpreter', 'latex', 'FontSize', 14);
title('Equilibrium Conditions (Contours)', 'Interpreter', 'latex', 'FontSize', 18);
legend('Arbitrage', 'Budget Constraint', 'Location', 'best');

%% =========================================================================
%  SECTION 3: Shock Comparative Statics
%  =========================================================================
fprintf('Computing shock comparative statics...\n');

% Perturb sigma downward (less volatility)
sigma_us_yt_shock = sigma_us_yt * 0.8;

% Redefine rate functions with shocked sigma
Echi_m_us_f = @(mu) Echi_m(mu, ploss_us, sigma_us_yt_shock, iota_us, lambda_us, eta, matching_type);
Echi_d_us_f = @(mu) Echi_d(mu, ploss_us, sigma_us_yt_shock, iota_us, lambda_us, eta, matching_type);

Rb_us_f = @(mu) Rm_us(t_shock_i) + Echi_m_us_f(mu);
Rd_us_f = @(mu) Rm_us(t_shock_i) + Echi_m_us_f(mu) + Echi_d_us_f(mu);

d_us_f = @(mu) (Theta_d_us_t(t_shock_i) * Rd_us_f(mu).^(1/zeta_us) - MBS_us_ss) / (1 - share_us*mu);
b_f = @(mu) Theta_b * (Rb_us_f(mu))^(1/epsilon_b);

nu_f = @(mu_us, mu_eu) d_eu_f(mu_eu) / d_us_f(mu_us);
res_arb = @(mu_us, mu_eu) (1 + riskprm_t(t_shock_i)) * Rb_eu_f(mu_eu) - Rb_us_f(mu_us);
res_bc = @(mu_us, mu_eu) b_f(mu_us) - (nu_f(mu_us, mu_eu)*(1-mu_eu) + 1-mu_us) * d_us_f(mu_us);

% Solve for new equilibrium curves
mu_eu_arb_vec2 = NaN(N_mu_us, 1);
mu_eu_bc_vec2 = NaN(N_mu_us, 1);

for mm = 1:N_mu_us
    mu_us_yt = vec_mu_us(mm);
    
    [mu_eu_out, ~, exitflag, ~] = fsolve(@(mu_eu) res_arb(mu_us_yt, exp(mu_eu)), log(0.1), opts);
    mu_eu_arb_vec2(mm) = exp(mu_eu_out);
    
    [mu_eu_out, ~, exitflag, ~] = fsolve(@(mu_eu) res_bc(mu_us_yt, exp(mu_eu)), log(0.1), opts);
    mu_eu_bc_vec2(mm) = exp(mu_eu_out);
end

%% Plot: Comparative Statics
figure('Name', 'Cross-Analysis: Arbitrage', 'NumberTitle', 'off');
plot(vec_mu_us, mu_eu_arb_vec, 'k-', 'LineWidth', 2); hold on;
plot(vec_mu_us, mu_eu_arb_vec2, 'r-', 'LineWidth', 2);
grid on; axis tight;
label_x('$\mu^{us}$');
label_y('$\mu^{eu}$');
formataxis(gca);
legend('Baseline', 'Low $\sigma^{us}$', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
title('Arbitrage Condition', 'Interpreter', 'latex', 'FontSize', 18);

figure('Name', 'Cross-Analysis: Budget Constraint', 'NumberTitle', 'off');
plot(vec_mu_us, mu_eu_bc_vec, 'k-.', 'LineWidth', 2); hold on;
plot(vec_mu_us, mu_eu_bc_vec2, 'r-.', 'LineWidth', 2);
grid on; axis tight;
label_x('$\mu^{us}$');
label_y('$\mu^{eu}$');
formataxis(gca);
legend('Baseline', 'Low $\sigma^{us}$', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
title('Budget Constraint', 'Interpreter', 'latex', 'FontSize', 18);

figure('Name', 'Cross-Analysis: Full', 'NumberTitle', 'off');
plot(vec_mu_us, mu_eu_arb_vec, 'k-', 'LineWidth', 2); hold on;
plot(vec_mu_us, mu_eu_arb_vec2, 'r-', 'LineWidth', 2);
plot(vec_mu_us, mu_eu_bc_vec, 'k-.', 'LineWidth', 2);
plot(vec_mu_us, mu_eu_bc_vec2, 'r-.', 'LineWidth', 2);
grid on; axis tight;
label_x('$\mu^{us}$');
label_y('$\mu^{eu}$');
formataxis(gca);
legend('Arb (baseline)', 'Arb (low $\sigma$)', 'BC (baseline)', 'BC (low $\sigma$)', ...
       'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
title('Equilibrium Cross-Conditions', 'Interpreter', 'latex', 'FontSize', 18);

fprintf('Sensitivity analysis complete.\n');
