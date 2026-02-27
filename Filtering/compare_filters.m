%% Compare Leontief vs Cobb-Douglas Filtering
% Runs main_filter.m with both matching functions and compares results
%
% (c) Saki Bigio
% 
% This script calls main_filter.m directly to ensure all solver improvements
% (e.g., find_sigma_min, fminbnd backup) are used.

clear; close all;

%% ========================================================================
%  RUN LEONTIEF FILTER
%  ========================================================================
fprintf('\n');
fprintf('========================================================\n');
fprintf('  RUNNING LEONTIEF FILTER (matching_type = 0)\n');
fprintf('========================================================\n');

matching_type = 0;
run('main_filter.m');

% Save Leontief results
sigma_us_L = sigma_us_t;
sigma_eu_L = sigma_eu_t;
BP_us_L = BP_us_t;
BP_eu_L = BP_eu_t;
CIP_L = CIP_t;
TED_us_L = TED_us_t;
TED_eu_L = TED_eu_t;
theta_us_L = theta_us_t;
theta_eu_L = theta_eu_t;
psi_us_L = psi_us_t;
psi_eu_L = psi_eu_t;
riskprm_L = riskprm_t;
DW_us_L = DW_us_t;
FF_us_L = FF_us_t;
Echi_m_us_L = Echi_m_us_t;
Echi_m_eu_L = Echi_m_eu_t;
sigma_us_TED_flag_L = sigma_us_TED_flag;
sigma_eu_TED_flag_L = sigma_eu_TED_flag;
dates_L = dates;
datesperiod_L = datesperiod;

% Save Leontief results to disk
save('data/filter_comparison_leontief.mat', ...
    'sigma_us_L', 'sigma_eu_L', 'BP_us_L', 'BP_eu_L', ...
    'CIP_L', 'TED_us_L', 'TED_eu_L', 'theta_us_L', 'theta_eu_L', ...
    'psi_us_L', 'psi_eu_L', 'riskprm_L', 'DW_us_L', 'FF_us_L', ...
    'Echi_m_us_L', 'Echi_m_eu_L', 'sigma_us_TED_flag_L', 'sigma_eu_TED_flag_L', ...
    'dates_L', 'datesperiod_L');
fprintf('Leontief results saved to data/filter_comparison_leontief.mat\n');

close all;  % Close Leontief plots

%% ========================================================================
%  RUN COBB-DOUGLAS FILTER
%  ========================================================================
fprintf('\n');
fprintf('========================================================\n');
fprintf('  RUNNING COBB-DOUGLAS FILTER (matching_type = 1)\n');
fprintf('========================================================\n');

matching_type = 1;
run('main_filter.m');

% Save Cobb-Douglas results
sigma_us_CD = sigma_us_t;
sigma_eu_CD = sigma_eu_t;
BP_us_CD = BP_us_t;
BP_eu_CD = BP_eu_t;
CIP_CD = CIP_t;
TED_us_CD = TED_us_t;
TED_eu_CD = TED_eu_t;
theta_us_CD = theta_us_t;
theta_eu_CD = theta_eu_t;
psi_us_CD = psi_us_t;
psi_eu_CD = psi_eu_t;
riskprm_CD = riskprm_t;
DW_us_CD = DW_us_t;
FF_us_CD = FF_us_t;
Echi_m_us_CD = Echi_m_us_t;
Echi_m_eu_CD = Echi_m_eu_t;
sigma_us_TED_flag_CD = sigma_us_TED_flag;
sigma_eu_TED_flag_CD = sigma_eu_TED_flag;
dates_CD = dates;
datesperiod_CD = datesperiod;

% Save CD results to disk
save('data/filter_comparison_cd.mat', ...
    'sigma_us_CD', 'sigma_eu_CD', 'BP_us_CD', 'BP_eu_CD', ...
    'CIP_CD', 'TED_us_CD', 'TED_eu_CD', 'theta_us_CD', 'theta_eu_CD', ...
    'psi_us_CD', 'psi_eu_CD', 'riskprm_CD', 'DW_us_CD', 'FF_us_CD', ...
    'Echi_m_us_CD', 'Echi_m_eu_CD', 'sigma_us_TED_flag_CD', 'sigma_eu_TED_flag_CD', ...
    'dates_CD', 'datesperiod_CD');
fprintf('Cobb-Douglas results saved to data/filter_comparison_cd.mat\n');

close all;  % Close CD plots

%% ========================================================================
%  LOAD RESULTS FOR COMPARISON
%  ========================================================================
load('data/filter_comparison_leontief.mat');
load('data/filter_comparison_cd.mat');
dates = dates_L;
datesperiod = datesperiod_L;

%% ========================================================================
%  SUMMARY STATISTICS
%  ========================================================================
fprintf('\n');
fprintf('========================================================\n');
fprintf('  COMPARISON: LEONTIEF vs COBB-DOUGLAS\n');
fprintf('========================================================\n\n');

fprintf('--- Solver Convergence ---\n');
fprintf('Leontief US:     %d/%d converged (%.1f%%)\n', sum(sigma_us_TED_flag_L > 0), length(sigma_us_TED_flag_L), 100*mean(sigma_us_TED_flag_L > 0));
fprintf('Leontief EU:     %d/%d converged (%.1f%%)\n', sum(sigma_eu_TED_flag_L > 0), length(sigma_eu_TED_flag_L), 100*mean(sigma_eu_TED_flag_L > 0));
fprintf('Cobb-Douglas US: %d/%d converged (%.1f%%), %d via fminbnd\n', ...
    sum(sigma_us_TED_flag_CD > 0), length(sigma_us_TED_flag_CD), 100*mean(sigma_us_TED_flag_CD > 0), sum(sigma_us_TED_flag_CD == 10));
fprintf('Cobb-Douglas EU: %d/%d converged (%.1f%%), %d via fminbnd\n', ...
    sum(sigma_eu_TED_flag_CD > 0), length(sigma_eu_TED_flag_CD), 100*mean(sigma_eu_TED_flag_CD > 0), sum(sigma_eu_TED_flag_CD == 10));

fprintf('\n--- Sigma Statistics ---\n');
fprintf('                        Leontief    Cobb-Douglas    Ratio\n');
fprintf('                        --------    ------------    -----\n');
fprintf('sigma_us Mean:          %8.4f    %12.4f    %5.2f\n', ...
    mean(sigma_us_L), mean(sigma_us_CD), mean(sigma_us_L)/mean(sigma_us_CD));
fprintf('sigma_us Std:           %8.4f    %12.4f    %5.2f\n', ...
    std(sigma_us_L), std(sigma_us_CD), std(sigma_us_L)/std(sigma_us_CD));
fprintf('sigma_eu Mean:          %8.4f    %12.4f    %5.2f\n', ...
    mean(sigma_eu_L), mean(sigma_eu_CD), mean(sigma_eu_L)/mean(sigma_eu_CD));
fprintf('sigma_eu Std:           %8.4f    %12.4f    %5.2f\n', ...
    std(sigma_eu_L), std(sigma_eu_CD), std(sigma_eu_L)/std(sigma_eu_CD));

fprintf('\n--- Spreads (bps annualized) ---\n');
fprintf('TED_us Mean:            %8.2f    %12.2f\n', mean(TED_us_L)*12*1e4, mean(TED_us_CD)*12*1e4);
fprintf('TED_us Std:             %8.2f    %12.2f\n', std(TED_us_L)*12*1e4, std(TED_us_CD)*12*1e4);
fprintf('BP_us Mean:             %8.2f    %12.2f\n', mean(BP_us_L)*12*1e4, mean(BP_us_CD)*12*1e4);
fprintf('CIP Mean:               %8.2f    %12.2f\n', mean(CIP_L)*12*1e4, mean(CIP_CD)*12*1e4);
fprintf('Risk Prm Mean:          %8.2f    %12.2f\n', mean(riskprm_L)*12*1e4, mean(riskprm_CD)*12*1e4);

fprintf('\n--- Correlations (Leontief vs Cobb-Douglas) ---\n');
fprintf('sigma_us:    %.4f\n', corr(sigma_us_L, sigma_us_CD));
fprintf('sigma_eu:    %.4f\n', corr(sigma_eu_L, sigma_eu_CD));
fprintf('theta_us:    %.4f\n', corr(theta_us_L, theta_us_CD));
fprintf('BP_us:       %.4f\n', corr(BP_us_L, BP_us_CD));
fprintf('CIP:         %.4f\n', corr(CIP_L, CIP_CD));
fprintf('Risk Prm:    %.4f\n', corr(riskprm_L, riskprm_CD));

%% ========================================================================
%  COMPARISON PLOTS
%  ========================================================================
FSize = 14;
formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
    'Fontsize', FSize, 'Box', 'Off');

%% Figure 1: Sigma time series overlay
figure('Name', 'Sigma Comparison', 'Position', [100 100 1200 600]);

subplot(2,2,1);
plot(dates(datesperiod), sigma_us_L(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_us_CD(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title('$\sigma_{US}$', 'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
formataxis(gca);

subplot(2,2,2);
plot(dates(datesperiod), sigma_eu_L(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_eu_CD(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title('$\sigma_{EU}$', 'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
formataxis(gca);

subplot(2,2,3);
plot(dates(datesperiod), theta_us_L(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), theta_us_CD(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title('$\theta_{US}$ (Market Tightness)', 'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
formataxis(gca);

subplot(2,2,4);
plot(dates(datesperiod), BP_us_L(datesperiod)*12*1e4, 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), BP_us_CD(datesperiod)*12*1e4, 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title('Bond Premium US (bps)', 'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
formataxis(gca);

sgtitle('Filter Comparison: Leontief vs Cobb-Douglas', 'Interpreter', 'latex', 'FontSize', 16);

%% Figure 2: Normalized dynamics comparison
figure('Name', 'Normalized Dynamics', 'Position', [100 100 1000 500]);

subplot(1,2,1);
sigma_us_L_norm = (sigma_us_L - mean(sigma_us_L)) / std(sigma_us_L);
sigma_us_CD_norm = (sigma_us_CD - mean(sigma_us_CD)) / std(sigma_us_CD);
plot(dates(datesperiod), sigma_us_L_norm(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_us_CD_norm(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title(sprintf('$\\sigma_{US}$ Normalized (corr = %.3f)', corr(sigma_us_L, sigma_us_CD)), ...
    'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
ylabel('Std. deviations from mean');
formataxis(gca);

subplot(1,2,2);
sigma_eu_L_norm = (sigma_eu_L - mean(sigma_eu_L)) / std(sigma_eu_L);
sigma_eu_CD_norm = (sigma_eu_CD - mean(sigma_eu_CD)) / std(sigma_eu_CD);
plot(dates(datesperiod), sigma_eu_L_norm(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_eu_CD_norm(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title(sprintf('$\\sigma_{EU}$ Normalized (corr = %.3f)', corr(sigma_eu_L, sigma_eu_CD)), ...
    'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
ylabel('Std. deviations from mean');
formataxis(gca);

sgtitle('Do Dynamics Match Across Matching Functions?', 'Interpreter', 'latex', 'FontSize', 16);

%% Figure 3: Scatter plots
figure('Name', 'Scatter Comparison', 'Position', [100 100 1000 400]);

subplot(1,3,1);
scatter(sigma_us_L(datesperiod), sigma_us_CD(datesperiod), 25, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
p = polyfit(sigma_us_L(datesperiod), sigma_us_CD(datesperiod), 1);
xl = xlim; plot(xl, polyval(p, xl), 'r-', 'LineWidth', 1.5);
xlabel('Leontief', 'Interpreter', 'latex');
ylabel('Cobb-Douglas', 'Interpreter', 'latex');
title(sprintf('$\\sigma_{US}$ (corr=%.3f)', corr(sigma_us_L, sigma_us_CD)), 'Interpreter', 'latex');
formataxis(gca); grid on;

subplot(1,3,2);
scatter(sigma_eu_L(datesperiod), sigma_eu_CD(datesperiod), 25, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
p = polyfit(sigma_eu_L(datesperiod), sigma_eu_CD(datesperiod), 1);
xl = xlim; plot(xl, polyval(p, xl), 'r-', 'LineWidth', 1.5);
xlabel('Leontief', 'Interpreter', 'latex');
ylabel('Cobb-Douglas', 'Interpreter', 'latex');
title(sprintf('$\\sigma_{EU}$ (corr=%.3f)', corr(sigma_eu_L, sigma_eu_CD)), 'Interpreter', 'latex');
formataxis(gca); grid on;

subplot(1,3,3);
scatter(BP_us_L(datesperiod)*12*1e4, BP_us_CD(datesperiod)*12*1e4, 25, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
p = polyfit(BP_us_L(datesperiod)*12*1e4, BP_us_CD(datesperiod)*12*1e4, 1);
xl = xlim; plot(xl, polyval(p, xl), 'r-', 'LineWidth', 1.5);
xlabel('Leontief (bps)', 'Interpreter', 'latex');
ylabel('Cobb-Douglas (bps)', 'Interpreter', 'latex');
title(sprintf('Bond Premium (corr=%.3f)', corr(BP_us_L, BP_us_CD)), 'Interpreter', 'latex');
formataxis(gca); grid on;

sgtitle('Level Relationship: Leontief vs Cobb-Douglas', 'Interpreter', 'latex', 'FontSize', 16);

%% Figure 4: CIP and Risk Premium
figure('Name', 'CIP and Risk Premium', 'Position', [100 100 1000 400]);

subplot(1,2,1);
plot(dates(datesperiod), CIP_L(datesperiod)*12*1e4, 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), CIP_CD(datesperiod)*12*1e4, 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title(sprintf('CIP Deviation (corr = %.3f)', corr(CIP_L, CIP_CD)), 'Interpreter', 'latex', 'FontSize', FSize+2);
ylabel('bps (annualized)');
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
formataxis(gca);

subplot(1,2,2);
plot(dates(datesperiod), riskprm_L(datesperiod)*12*1e4, 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), riskprm_CD(datesperiod)*12*1e4, 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title(sprintf('Risk Premium (corr = %.3f)', corr(riskprm_L, riskprm_CD)), 'Interpreter', 'latex', 'FontSize', FSize+2);
ylabel('bps (annualized)');
legend('Leontief', 'Cobb-Douglas', 'Location', 'best');
formataxis(gca);

sgtitle('International Variables: Leontief vs Cobb-Douglas', 'Interpreter', 'latex', 'FontSize', 16);

%% ========================================================================
%  SAVE RESULTS
%  ========================================================================
save('data/filter_comparison.mat', ...
    'sigma_us_L', 'sigma_eu_L', 'sigma_us_CD', 'sigma_eu_CD', ...
    'BP_us_L', 'BP_eu_L', 'BP_us_CD', 'BP_eu_CD', ...
    'CIP_L', 'CIP_CD', 'TED_us_L', 'TED_eu_L', 'TED_us_CD', 'TED_eu_CD', ...
    'theta_us_L', 'theta_eu_L', 'theta_us_CD', 'theta_eu_CD', ...
    'riskprm_L', 'riskprm_CD', ...
    'dates', 'datesperiod');

fprintf('\n========================================================\n');
fprintf('Results saved to data/filter_comparison.mat\n');
fprintf('========================================================\n');
