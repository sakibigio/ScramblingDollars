%% Compare Leontief vs Cobb-Douglas Filtering
% Runs the filter with both matching functions and compares results
%
% (c) Saki Bigio

clear; close all;

%% Setup paths
addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

%% Load data and parameters (shared across both runs)
load('data/LFX_data3.mat');
load('data/calibration.mat');

% Scale
abs_scale = 12e4;

% Parameters
share_us = 1; 
share_eu = 1;
pi_eu_ss = 1;
pi_us_ss = 1;

run('functions/params.m');

epsilon_b = -0.5;
zeta_us = -epsilon_b; 
zeta_eu = -epsilon_b;
Theta_b = 3;

% Date variables
datesperiod = 1:234;
dates = datenum(2001, 1:234, 1);

%% Run filter for LEONTIEF (matching_type = 0)
fprintf('\n========================================\n');
fprintf('Running LEONTIEF filter (matching_type = 0)\n');
fprintf('========================================\n');

matching_type = 0;
[sigma_us_leontief, sigma_eu_leontief, results_leontief] = run_filter(matching_type);

%% Run filter for COBB-DOUGLAS (matching_type = 1)
fprintf('\n========================================\n');
fprintf('Running COBB-DOUGLAS filter (matching_type = 1)\n');
fprintf('========================================\n');

matching_type = 1;
[sigma_us_cd, sigma_eu_cd, results_cd] = run_filter(matching_type);

%% Summary Statistics
fprintf('\n========================================\n');
fprintf('COMPARISON: Leontief vs Cobb-Douglas\n');
fprintf('========================================\n\n');

fprintf('                        Leontief    Cobb-Douglas    Ratio\n');
fprintf('                        --------    ------------    -----\n');
fprintf('Mean sigma_us:          %8.4f    %12.4f    %5.1f\n', ...
    mean(sigma_us_leontief), mean(sigma_us_cd), mean(sigma_us_leontief)/mean(sigma_us_cd));
fprintf('Mean sigma_eu:          %8.4f    %12.4f    %5.1f\n', ...
    mean(sigma_eu_leontief), mean(sigma_eu_cd), mean(sigma_eu_leontief)/mean(sigma_eu_cd));
fprintf('Std sigma_us:           %8.4f    %12.4f    %5.1f\n', ...
    std(sigma_us_leontief), std(sigma_us_cd), std(sigma_us_leontief)/std(sigma_us_cd));
fprintf('Std sigma_eu:           %8.4f    %12.4f    %5.1f\n', ...
    std(sigma_eu_leontief), std(sigma_eu_cd), std(sigma_eu_leontief)/std(sigma_eu_cd));

fprintf('\nImplied Laplace std (sqrt(2)*sigma):\n');
fprintf('  Leontief US:     %.3f\n', sqrt(2)*mean(sigma_us_leontief));
fprintf('  Cobb-Douglas US: %.3f\n', sqrt(2)*mean(sigma_us_cd));
fprintf('\nReference values:\n');
fprintf('  D''Erasmo (daily):     0.04\n');
fprintf('  Bigio (quarterly):    0.06\n');

fprintf('\nRisk Premium (bps):\n');
fprintf('  Leontief:      %.2f\n', results_leontief.riskprm_mean * 1e4);
fprintf('  Cobb-Douglas:  %.2f\n', results_cd.riskprm_mean * 1e4);

%% Plot Comparison
FSize = 14;
formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize, 'Box', 'Off');

% Figure 1: Sigma time series comparison
figure('Name', 'Sigma Comparison', 'Position', [100 100 1200 800]);

subplot(2,2,1);
plot(dates(datesperiod), sigma_us_leontief(datesperiod), 'b-', 'LineWidth', 1.5);
datetick('x', 'yyyy', 'keeplimits');
title('US $\sigma$ - Leontief', 'Interpreter', 'latex', 'FontSize', FSize);
ylabel('$\sigma$', 'Interpreter', 'latex');
formataxis(gca); grid on;

subplot(2,2,2);
plot(dates(datesperiod), sigma_us_cd(datesperiod), 'r-', 'LineWidth', 1.5);
datetick('x', 'yyyy', 'keeplimits');
title('US $\sigma$ - Cobb-Douglas', 'Interpreter', 'latex', 'FontSize', FSize);
ylabel('$\sigma$', 'Interpreter', 'latex');
formataxis(gca); grid on;

subplot(2,2,3);
plot(dates(datesperiod), sigma_eu_leontief(datesperiod), 'b-', 'LineWidth', 1.5);
datetick('x', 'yyyy', 'keeplimits');
title('EU $\sigma$ - Leontief', 'Interpreter', 'latex', 'FontSize', FSize);
ylabel('$\sigma$', 'Interpreter', 'latex');
xlabel('Year');
formataxis(gca); grid on;

subplot(2,2,4);
plot(dates(datesperiod), sigma_eu_cd(datesperiod), 'r-', 'LineWidth', 1.5);
datetick('x', 'yyyy', 'keeplimits');
title('EU $\sigma$ - Cobb-Douglas', 'Interpreter', 'latex', 'FontSize', FSize);
ylabel('$\sigma$', 'Interpreter', 'latex');
xlabel('Year');
formataxis(gca); grid on;

sgtitle('Filtered $\sigma$ Comparison: Leontief vs Cobb-Douglas', 'Interpreter', 'latex', 'FontSize', 16);

% Figure 2: Overlay comparison (normalized)
figure('Name', 'Normalized Comparison', 'Position', [100 100 1000 500]);

subplot(1,2,1);
% Normalize to compare dynamics
sigma_us_leon_norm = (sigma_us_leontief - mean(sigma_us_leontief)) / std(sigma_us_leontief);
sigma_us_cd_norm = (sigma_us_cd - mean(sigma_us_cd)) / std(sigma_us_cd);
plot(dates(datesperiod), sigma_us_leon_norm(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_us_cd_norm(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yyyy', 'keeplimits');
title('US $\sigma$ (normalized)', 'Interpreter', 'latex', 'FontSize', FSize);
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
ylabel('Std deviations from mean');
formataxis(gca); grid on;

% Correlation
corr_us = corr(sigma_us_leontief(datesperiod), sigma_us_cd(datesperiod));
text(0.05, 0.95, sprintf('Corr = %.3f', corr_us), 'Units', 'normalized', ...
    'FontSize', 12, 'VerticalAlignment', 'top');

subplot(1,2,2);
sigma_eu_leon_norm = (sigma_eu_leontief - mean(sigma_eu_leontief)) / std(sigma_eu_leontief);
sigma_eu_cd_norm = (sigma_eu_cd - mean(sigma_eu_cd)) / std(sigma_eu_cd);
plot(dates(datesperiod), sigma_eu_leon_norm(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_eu_cd_norm(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yyyy', 'keeplimits');
title('EU $\sigma$ (normalized)', 'Interpreter', 'latex', 'FontSize', FSize);
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
ylabel('Std deviations from mean');
formataxis(gca); grid on;

corr_eu = corr(sigma_eu_leontief(datesperiod), sigma_eu_cd(datesperiod));
text(0.05, 0.95, sprintf('Corr = %.3f', corr_eu), 'Units', 'normalized', ...
    'FontSize', 12, 'VerticalAlignment', 'top');

sgtitle('Normalized $\sigma$: Do dynamics match?', 'Interpreter', 'latex', 'FontSize', 16);

% Figure 3: Scatter plot
figure('Name', 'Scatter Comparison', 'Position', [100 100 800 400]);

subplot(1,2,1);
scatter(sigma_us_leontief(datesperiod), sigma_us_cd(datesperiod), 20, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('$\sigma_{US}$ Leontief', 'Interpreter', 'latex');
ylabel('$\sigma_{US}$ Cobb-Douglas', 'Interpreter', 'latex');
title(sprintf('US: Corr = %.3f', corr_us), 'Interpreter', 'latex');
formataxis(gca); grid on;
% Add reference line
hold on;
xl = xlim; 
p = polyfit(sigma_us_leontief(datesperiod), sigma_us_cd(datesperiod), 1);
plot(xl, polyval(p, xl), 'r-', 'LineWidth', 1);

subplot(1,2,2);
scatter(sigma_eu_leontief(datesperiod), sigma_eu_cd(datesperiod), 20, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('$\sigma_{EU}$ Leontief', 'Interpreter', 'latex');
ylabel('$\sigma_{EU}$ Cobb-Douglas', 'Interpreter', 'latex');
title(sprintf('EU: Corr = %.3f', corr_eu), 'Interpreter', 'latex');
formataxis(gca); grid on;
hold on;
xl = xlim;
p = polyfit(sigma_eu_leontief(datesperiod), sigma_eu_cd(datesperiod), 1);
plot(xl, polyval(p, xl), 'r-', 'LineWidth', 1);

sgtitle('Leontief vs Cobb-Douglas: Level Relationship', 'Interpreter', 'latex', 'FontSize', 16);

%% Save comparison results
save('data/sigma_comparison.mat', ...
    'sigma_us_leontief', 'sigma_eu_leontief', ...
    'sigma_us_cd', 'sigma_eu_cd', ...
    'results_leontief', 'results_cd', ...
    'dates', 'datesperiod');

fprintf('\nResults saved to data/sigma_comparison.mat\n');

%% ========================================================================
%  FILTER FUNCTION (embedded)
%  ========================================================================
function [sigma_us_t, sigma_eu_t, results] = run_filter(matching_type)
    % Load data
    load('data/LFX_data3.mat');
    load('data/calibration.mat');
    
    % Required before params.m
    pi_eu_ss = 1;
    pi_us_ss = 1;
    
    % Parameters
    run('functions/params.m');
    abs_scale = 12e4;
    epsilon_b = -0.5;
    zeta_us = -epsilon_b; 
    zeta_eu = -epsilon_b;
    Theta_b = 3;
    share_us = 1;
    share_eu = 1;
    
    % Matching type name
    if matching_type == 0
        matching_name = 'Leontief';
    else
        matching_name = 'Cobb-Douglas';
    end
    
    % Min test values
    sigma_vec = (0.01:0.01:10);
    min_test_us = Chi_p_psi(min(exp(mu_us)), ploss_us, min(sigma_vec), iota_us, lambda_us, eta, matching_type) * abs_scale + 0.00012 * abs_scale;
    min_test_eu = Chi_p_psi(min(exp(mu_eu)), ploss_eu, min(sigma_vec), iota_eu, lambda_eu, eta, matching_type) * abs_scale + 0.00012 * abs_scale;
    
    % Initialize
    T = length(mu_us);
    sigma_us_t = zeros(T, 1);
    sigma_eu_t = zeros(T, 1);
    sigma_us_TED_t = zeros(T, 1);
    sigma_eu_TED_t = zeros(T, 1);
    Echi_m_us_t = zeros(T, 1);
    Echi_d_us_t = zeros(T, 1);
    Echi_m_eu_t = zeros(T, 1);
    Echi_d_eu_t = zeros(T, 1);
    theta_us_t = zeros(T, 1);
    theta_eu_t = zeros(T, 1);
    DW_us_t = zeros(T, 1);
    FF_us_t = zeros(T, 1);
    
    % Initial guesses
    sigma_us_TED_guess = sigma_us;
    sigma_eu_TED_guess = sigma_eu;
    
    % Main loop
    for tt = 1:T
        % Targets
        TED_us_target = min_test_us + (TED_s_us_t(tt) - min(TED_s_us_t)) * abs_scale; 
        TED_eu_target = min_test_eu + (TED_s_eu_t(tt) - min(TED_s_eu_t)) * abs_scale;
        
        mu_us_yt = exp(mu_us(tt));
        mu_eu_yt = exp(mu_eu(tt));
        
        % US Sigma - TED Target
        % For Cobb-Douglas, use mu-dependent initial guess to avoid cliff
        if matching_type == 1
            sigma_us_TED_guess = max(sigma_us_TED_guess, mu_us_yt + 0.15);
        end
        sigma_us_res = @(sigma) Chi_p_psi(mu_us_yt, ploss_us, sigma, iota_us, lambda_us, eta, matching_type) * 1e4 * 12 - TED_us_target;
        [sigma_out, ~, ~, ~] = fsolve(@(sigma) sigma_us_res(sigma), sigma_us_TED_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12));
        sigma_us_TED_t(tt) = sigma_out;
        sigma_us_TED_guess = sigma_out;
        
        % EU Sigma - TED target
        % For Cobb-Douglas, use mu-dependent initial guess to avoid cliff
        if matching_type == 1
            sigma_eu_TED_guess = max(sigma_eu_TED_guess, mu_eu_yt + 0.15);
        end
        sigma_eu_res = @(sigma) Chi_p_psi(mu_eu_yt, ploss_eu, sigma, iota_eu, lambda_eu, eta, matching_type) * 1e4 * 12 - TED_eu_target;
        [sigma_out, ~, ~, ~] = fsolve(@(sigma) sigma_eu_res(sigma), sigma_eu_TED_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12));
        sigma_eu_TED_t(tt) = sigma_out;
        sigma_eu_TED_guess = sigma_out;
        
        % Store
        sigma_us_t(tt) = sigma_us_TED_t(tt);
        sigma_eu_t(tt) = sigma_eu_TED_t(tt);
        
        % Interbank variables
        Echi_m_us_t(tt) = Echi_m(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type);
        Echi_d_us_t(tt) = Echi_d(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type);
        Echi_m_eu_t(tt) = Echi_m(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type);
        Echi_d_eu_t(tt) = Echi_d(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type);
        
        [~, theta_us_t(tt), ~, ~, DW_us_t(tt), FF_us_t(tt), ~] = Chi_sys(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type);
        [~, theta_eu_t(tt), ~, ~, ~, ~, ~] = Chi_sys(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type);
        
        if mod(tt, 50) == 0
            fprintf('  %s: Period %d/%d\n', matching_name, tt, T);
        end
    end
    
    % Derived variables
    Rm_us = exp(im_us);
    Rm_eu = exp(im_eu);
    Rb_us_t = Rm_us + Echi_m_us_t;
    Rb_eu_t = Rm_eu + Echi_m_eu_t;
    riskprm_t = Rb_us_t ./ Rb_eu_t - 1;
    
    % Pack results
    results.Echi_m_us = Echi_m_us_t;
    results.Echi_m_eu = Echi_m_eu_t;
    results.Echi_d_us = Echi_d_us_t;
    results.Echi_d_eu = Echi_d_eu_t;
    results.theta_us = theta_us_t;
    results.theta_eu = theta_eu_t;
    results.DW_us = DW_us_t;
    results.FF_us = FF_us_t;
    results.riskprm = riskprm_t;
    results.riskprm_mean = mean(riskprm_t);
    results.matching_type = matching_type;
    results.matching_name = matching_name;
end
