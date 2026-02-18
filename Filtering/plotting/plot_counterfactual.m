%% PLOT_COUNTERFACTUAL.M
% Counterfactual analysis: What if there were no scrambling episodes?
% Part of: Scrambling for Dollars filter pipeline
%
% Purpose:
%   Simulates the model with sigma held at regime-1 mean (no crisis spikes)
%   and compares actual vs counterfactual paths for key variables.
%
% Requires:
%   - Run main_filter.m first (populates workspace)
%   - Run markov_estimation.jl first (creates MS_sigma_*_counterfactuals.csv)
%   - Files: MS_sigma_us_counterfactuals.csv, MS_sigma_eu_counterfactuals.csv
%
% Outputs:
%   - ~10 figures comparing actual vs counterfactual paths
%   - Saved results in LFX_rwfilter.mat
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

%% Load Data
load('data/LFX_data.mat');

%% Store Baseline Simulation Values
sigma_us_sim = sigma_us_t;
sigma_eu_sim = sigma_eu_t;
Theta_d_us_sim = Theta_d_us_t;
Theta_d_eu_sim = Theta_d_eu_t;
riskprm_sim = riskprm_t;
Rm_us_sim = exp(im_us);
Rm_eu_sim = exp(im_eu);
M_us_sim = M_us_t;
M_eu_sim = M_eu_t;

%% =========================================================================
%  SECTION 1: Load Counterfactual Shocks from Julia
%  =========================================================================
fprintf('Loading counterfactual shocks from Julia output...\n');

% Check if Julia output files exist
us_file = 'data/MS_sigma_us_counterfactuals.csv';
eu_file = 'data/MS_sigma_eu_counterfactuals.csv';

if ~exist(us_file, 'file') || ~exist(eu_file, 'file')
    warning('Counterfactual CSV files not found. Run markov_estimation.jl first.');
    fprintf('Expected files:\n  %s\n  %s\n', us_file, eu_file);
    fprintf('Skipping counterfactual analysis.\n');
    return;
end

% US counterfactual - use readmatrix for better compatibility
try
    sigma_us_r1sim = readmatrix(us_file, 'NumHeaderLines', 1);
catch
    % Fallback for older MATLAB versions
    sigma_us_r1sim = csvread(us_file, 1, 0);
end
sigma_us_r1sim = [sigma_us_t(1); sigma_us_r1sim];  % Prepend initial value

% EU counterfactual
try
    sigma_eu_r1sim = readmatrix(eu_file, 'NumHeaderLines', 1);
catch
    sigma_eu_r1sim = csvread(eu_file, 1, 0);
end
sigma_eu_r1sim = [sigma_eu_t(1); sigma_eu_r1sim];

%% Plot: Counterfactual vs Actual Shocks
figure('Name', 'Sigma Counterfactual (US)', 'NumberTitle', 'off');
plot(dates(datesperiod), sigma_us_t(datesperiod), 'LineWidth', 3); hold on;
plot(dates(datesperiod), sigma_us_r1sim(datesperiod), 'LineWidth', 2);
plot(dates(datesperiod), sigma_us_sim(datesperiod), 'LineWidth', 1, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Actual', 'Regime 1 Only', 'Simulation', 'Location', 'best', 'Box', 'off');
title('$\sigma^{us}$: Actual vs Counterfactual', 'Interpreter', 'latex', 'FontSize', 18);

figure('Name', 'Sigma Counterfactual (EU)', 'NumberTitle', 'off');
plot(dates(datesperiod), sigma_eu_t(datesperiod), 'LineWidth', 3); hold on;
plot(dates(datesperiod), sigma_eu_r1sim(datesperiod), 'LineWidth', 2);
plot(dates(datesperiod), sigma_eu_sim(datesperiod), 'LineWidth', 1, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Actual', 'Regime 1 Only', 'Simulation', 'Location', 'best', 'Box', 'off');
title('$\sigma^{eu}$: Actual vs Counterfactual', 'Interpreter', 'latex', 'FontSize', 18);

%% Set Counterfactual Scenario
% Use mean of regime-1 shock (no crisis spikes)
sigma_us_sim = ones(length(sigma_us_r1sim), 1) * mean(sigma_us_r1sim);

%% =========================================================================
%  SECTION 2: Pre-allocate Counterfactual Simulation Vectors
%  =========================================================================
T = length(endopath);

% Price variables
Echi_m_us_sim = zeros(T, 1);
Echi_d_us_sim = zeros(T, 1);
Chi_p_psi_us_sim = zeros(T, 1);
Echi_m_eu_sim = zeros(T, 1);
Echi_d_eu_sim = zeros(T, 1);
Chi_p_psi_eu_sim = zeros(T, 1);

% Interbank variables
mu_us_sim = zeros(T, 1);
theta_us_sim = zeros(T, 1);
psi_us_sim = zeros(T, 1);
Smin_us_sim = zeros(T, 1);
DW_us_sim = zeros(T, 1);
FF_us_sim = zeros(T, 1);
Q_us_sim = zeros(T, 1);

mu_eu_sim = zeros(T, 1);
theta_eu_sim = zeros(T, 1);
psi_eu_sim = zeros(T, 1);
Smin_eu_sim = zeros(T, 1);
DW_eu_sim = zeros(T, 1);
FF_eu_sim = zeros(T, 1);
Q_eu_sim = zeros(T, 1);

d_us_sim = zeros(T, 1);
d_eu_sim = zeros(T, 1);

% Diagnostics
rest_flag = zeros(T, 1);

%% =========================================================================
%  SECTION 3: Main Counterfactual Simulation Loop
%  =========================================================================
fprintf('Running counterfactual simulation...\n');

for tt = 1:T
    % Update rate functions with counterfactual sigma
    Echi_m_us_f = @(mu) Echi_m(mu, ploss_us, sigma_us_sim(tt), iota_us, lambda_us, eta, matching_type);
    Echi_d_us_f = @(mu) Echi_d(mu, ploss_us, sigma_us_sim(tt), iota_us, lambda_us, eta, matching_type);
    Echi_m_eu_f = @(mu) Echi_m(mu, ploss_eu, sigma_eu_sim(tt), iota_eu, lambda_eu, eta, matching_type);
    Echi_d_eu_f = @(mu) Echi_d(mu, ploss_eu, sigma_eu_sim(tt), iota_eu, lambda_eu, eta, matching_type);

    Rb_us_f = @(mu) Rm_us_sim(tt) + Echi_m_us_f(mu);
    Rd_us_f = @(mu) Rm_us_sim(tt) + Echi_m_us_f(mu) + Echi_d_us_f(mu);
    Rb_eu_f = @(mu) Rm_eu_sim(tt) + Echi_m_eu_f(mu);
    Rd_eu_f = @(mu) Rm_eu_sim(tt) + Echi_m_eu_f(mu) + Echi_d_eu_f(mu);

    % Quantity functions
    d_us_f = @(mu) (Theta_d_us_sim(tt) * Rd_us_f(mu).^(1/zeta_us) - MBS_us_ss) / (1 - share_us*mu);
    d_eu_f = @(mu) (Theta_d_eu_sim(tt) * Rd_eu_f(mu).^(1/zeta_eu) - MBS_eu_ss) / (1 - share_eu*mu);
    b_f = @(mu) Theta_b * (Rb_us_f(mu))^(1/epsilon_b);

    % Residual functions for equilibrium
    nu_f = @(mu_us, mu_eu) d_eu_f(mu_eu) / d_us_f(mu_us);
    res_arb = @(mu_us, mu_eu) (1 + riskprm_sim(tt)) * Rb_eu_f(mu_eu) - Rb_us_f(mu_us);
    res_bc = @(mu_us, mu_eu) b_f(mu_us) - (nu_f(mu_us, mu_eu)*(1-mu_eu) + 1-mu_us) * d_us_f(mu_us);

    % Total residual
    res = @(mu_vec) [res_arb(mu_vec(1), mu_vec(2)); res_bc(mu_vec(1), mu_vec(2))];

    % Solve for optimal mu's
    opts = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, ...
                        'MaxFunEvals', 1e9, 'MaxIter', 1e6);
    [mu_out, ~, exitflag, ~] = fsolve(@(mu_vec) res(exp(mu_vec)), [0.5; 0.5], opts);
    rest_flag(tt) = exitflag;

    % Store results
    mu_us_sim(tt) = exp(mu_out(1));
    mu_eu_sim(tt) = exp(mu_out(2));
    d_us_sim(tt) = d_us_f(mu_us_sim(tt));
    d_eu_sim(tt) = d_eu_f(mu_eu_sim(tt));

    % Chi values at counterfactual equilibrium
    Echi_m_us_sim(tt) = Echi_m_us_f(mu_us_sim(tt));
    Echi_d_us_sim(tt) = Echi_d_us_f(mu_us_sim(tt));
    Chi_p_psi_us_sim(tt) = Chi_p_psi(mu_us_sim(tt), ploss_us, sigma_us_sim(tt), iota_us, lambda_us, eta, matching_type);
    Echi_m_eu_sim(tt) = Echi_m_eu_f(mu_eu_sim(tt));
    Echi_d_eu_sim(tt) = Echi_d_eu_f(mu_eu_sim(tt));
    Chi_p_psi_eu_sim(tt) = Chi_p_psi(mu_eu_sim(tt), ploss_eu, sigma_eu_sim(tt), iota_eu, lambda_eu, eta, matching_type);

    % Interbank variables
    [~, theta_us_sim(tt), psi_us_sim(tt), Smin_us_sim(tt), DW_us_sim(tt), FF_us_sim(tt), Q_us_sim(tt)] = ...
        Chi_sys(mu_us_sim(tt), ploss_us, sigma_us_sim(tt), iota_us, lambda_us, eta, matching_type);
    [~, theta_eu_sim(tt), psi_eu_sim(tt), Smin_eu_sim(tt), DW_eu_sim(tt), FF_eu_sim(tt), Q_eu_sim(tt)] = ...
        Chi_sys(mu_eu_sim(tt), ploss_eu, sigma_eu_sim(tt), iota_eu, lambda_eu, eta, matching_type);
end

fprintf('Counterfactual simulation complete. Solver flags: %d successes\n', sum(rest_flag > 0));

%% =========================================================================
%  SECTION 4: Construct Counterfactual Variables
%  =========================================================================

% Rates
Rb_us_sim = Rm_us_sim + Echi_m_us_sim;
Rd_us_sim = Rm_us_sim + Echi_m_us_sim + Echi_d_us_sim;
Rb_eu_sim = Rm_eu_sim + Echi_m_eu_sim;
Rd_eu_sim = Rm_eu_sim + Echi_m_eu_sim + Echi_d_eu_sim;

% Price system
nu_sim = d_eu_sim ./ d_us_sim;
p_us_sim = M_us_sim ./ (d_us_sim .* mu_us_sim);
p_eu_sim = M_eu_sim ./ (d_eu_sim .* mu_eu_sim);
e_euus_sim = p_us_sim ./ p_eu_sim;

% Ratios for decomposition
M_rat_sim = M_us_sim ./ M_eu_sim;
M_rat_t = M_us_t ./ M_eu_t;
mu_rat_sim = mu_us_sim ./ mu_eu_sim;
mu_rat_t = mu_us_t ./ mu_eu_t;

%% =========================================================================
%  SECTION 5: Counterfactual Plots
%  =========================================================================

% Liquidity ratio (US)
figure('Name', 'Mu Counterfactual (US)', 'NumberTitle', 'off');
plot(dates(datesperiod), ones(1, length(datesperiod)) * mean(mu_us_t), 'k--', 'LineWidth', 1); hold on;
plot(dates(datesperiod), mu_us_t(datesperiod), 'LineWidth', 3);
plot(dates(datesperiod), mu_us_sim(datesperiod), 'r:', 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Mean', 'Actual', 'Counterfactual', 'Location', 'best', 'Box', 'off');
title('$\mu^{us}$ Counterfactual', 'Interpreter', 'latex', 'FontSize', 18);

% Liquidity ratio (EU)
figure('Name', 'Mu Counterfactual (EU)', 'NumberTitle', 'off');
plot(dates(datesperiod), ones(1, length(datesperiod)) * mean(mu_eu_t), 'k--', 'LineWidth', 1); hold on;
plot(dates(datesperiod), mu_eu_t(datesperiod), 'LineWidth', 3);
plot(dates(datesperiod), mu_eu_sim(datesperiod), 'r:', 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Mean', 'Actual', 'Counterfactual', 'Location', 'best', 'Box', 'off');
title('$\mu^{eu}$ Counterfactual', 'Interpreter', 'latex', 'FontSize', 18);

% Bond premium (US)
figure('Name', 'Echi_m Counterfactual (US)', 'NumberTitle', 'off');
plot(dates(datesperiod), ones(1, length(datesperiod)) * mean(Echi_m_us_t) * abs_scale, 'k--', 'LineWidth', 1); hold on;
plot(dates(datesperiod), Echi_m_us_t(datesperiod) * abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), Echi_m_us_sim(datesperiod) * abs_scale, 'r:', 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
label_y('BPS');
formataxis(gca);
legend('Mean', 'Actual', 'Counterfactual', 'Location', 'best', 'Box', 'off');
title('$E[\chi^m]^{us}$ Counterfactual', 'Interpreter', 'latex', 'FontSize', 18);

% Deposit premium (US)
figure('Name', 'Echi_d Counterfactual (US)', 'NumberTitle', 'off');
plot(dates(datesperiod), ones(1, length(datesperiod)) * mean(Echi_d_us_t) * abs_scale, 'k--', 'LineWidth', 1); hold on;
plot(dates(datesperiod), Echi_d_us_t(datesperiod) * abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), Echi_d_us_sim(datesperiod) * abs_scale, 'r:', 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
label_y('BPS');
formataxis(gca);
legend('Mean', 'Actual', 'Counterfactual', 'Location', 'best', 'Box', 'off');
title('$E[\chi^d]^{us}$ Counterfactual', 'Interpreter', 'latex', 'FontSize', 18);

% Deposit premium (EU)
figure('Name', 'Echi_d Counterfactual (EU)', 'NumberTitle', 'off');
plot(dates(datesperiod), ones(1, length(datesperiod)) * mean(Echi_d_eu_t) * abs_scale, 'k--', 'LineWidth', 1); hold on;
plot(dates(datesperiod), Echi_d_eu_t(datesperiod) * abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), Echi_d_eu_sim(datesperiod) * abs_scale, 'r:', 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
label_y('BPS');
formataxis(gca);
legend('Mean', 'Actual', 'Counterfactual', 'Location', 'best', 'Box', 'off');
title('$E[\chi^d]^{eu}$ Counterfactual', 'Interpreter', 'latex', 'FontSize', 18);

% FX counterfactual (log)
figure('Name', 'FX Counterfactual (log)', 'NumberTitle', 'off');
plot(dates(datesperiod), ones(1, length(datesperiod)) * mean(log(e_euus_t)), 'k--', 'LineWidth', 1); hold on;
plot(dates(datesperiod), log(e_euus_t(datesperiod)), 'LineWidth', 3);
plot(dates(datesperiod), log(e_euus_sim(datesperiod)), 'Color', [0.8 0.6 0.6], 'LineWidth', 2, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Mean', 'Actual', 'Counterfactual', 'Location', 'best', 'Box', 'off');
title('Log FX Counterfactual', 'Interpreter', 'latex', 'FontSize', 18);

% FX decomposition
figure('Name', 'FX Decomposition', 'NumberTitle', 'off');
plot(dates(datesperiod), zeros(1, length(datesperiod)), 'k--', 'LineWidth', 1); hold on;
plot(dates(datesperiod), log(e_euus_sim(datesperiod)) - log(e_euus_t(datesperiod)), 'b', 'LineWidth', 2);
plot(dates(datesperiod), log(M_rat_sim(datesperiod)) - log(M_rat_t(datesperiod)), 'Color', [0.5 0.6 0], 'LineWidth', 2, 'LineStyle', '-.');
plot(dates(datesperiod), log(mu_rat_sim(datesperiod)) - log(mu_rat_t(datesperiod)), 'Color', [0.5 0.6 0.7], 'LineWidth', 2, 'LineStyle', '--');
plot(dates(datesperiod), log(nu_sim(datesperiod)) - log(nu_t(datesperiod)), 'm--', 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('', 'FX', 'M diff', '$\mu$ diff', '$\nu$ diff', 'Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
title('FX Decomposition (Counterfactual - Actual)', 'Interpreter', 'latex', 'FontSize', 18);

%% =========================================================================
%  SECTION 6: Save Results
%  =========================================================================
fprintf('Saving counterfactual results...\n');

% Store baseline shocks for comparison
sigma_us_rw = sigma_us_t;
sigma_eu_rw = sigma_eu_t;
Theta_d_eu_rw = Theta_d_eu_t;
Theta_d_us_rw = Theta_d_us_t;
riskprm_rw = riskprm_t;

save('LFX_rwfilter.mat', 'sigma_us_rw', 'sigma_eu_rw', 'Theta_d_eu_rw', 'Theta_d_us_rw', 'riskprm_rw');

fprintf('Counterfactual analysis complete.\n');
