%% PLOT_DIAGNOSTICS.M
% Diagnostic plots and time series analysis for filtered shocks
% Part of: Scrambling for Dollars filter pipeline
%
% Purpose:
%   Analyzes the statistical properties of filtered shocks:
%   - ARMA estimation for rate differentials and funding shocks
%   - Money supply rule estimation
%   - Risk premium VAR analysis
%   - FX return predictability tests
%
% Requires:
%   - Run main_filter.m first to populate workspace with:
%     sigma_us_t, sigma_eu_t, Theta_d_us_t, Theta_d_eu_t, riskprm_t
%     dates, datesperiod, im_us, im_eu, M_us_t, M_eu_t, inv_e, f_t, etc.
%
% Outputs:
%   - ~15 diagnostic figures
%   - Console output with ARMA coefficients and white noise tests
%
% Authors: Bianchi, Bigio, Engel
% Last updated: February 2025
%==========================================================================

%% Settings
close all;

% Date range for estimation
year_target = (1:225);
estimate_range = year_target;

%% Handle Dynare path collision (Dynare's arima shadows MATLAB's)
% Save current path and remove Dynare directories in one operation
original_path = path;
path_cells = strsplit(path, pathsep);
dynare_idx = contains(path_cells, 'Dynare', 'IgnoreCase', true);
if any(dynare_idx)
    clean_path = strjoin(path_cells(~dynare_idx), pathsep);
    path(clean_path);
    restore_dynare = true;
    fprintf('Temporarily removed Dynare from path (arima conflict)\n');
else
    restore_dynare = false;
end

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

%% Rename Shocks and Compute Logs
sigma_us_rw = sigma_us_t;
sigma_eu_rw = sigma_eu_t;
Theta_d_us_rw = Theta_d_us_t;
Theta_d_eu_rw = Theta_d_eu_t;
riskprm_rw = riskprm_t;

% Log transformations
sigma_us_lt = log(sigma_us_t);
sigma_us_lag_lt = lagmatrix(sigma_us_lt, 1);
sigma_eu_lt = log(sigma_eu_t);
sigma_eu_lag_lt = lagmatrix(sigma_eu_lt, 1);
dsigma_us_lt = diff(sigma_us_lt);
dsigma_us_lag_lt = lagmatrix(dsigma_us_lt, 1);
dsigma_eu_lt = diff(sigma_eu_lt);
dsigma_eu_lag_lt = lagmatrix(dsigma_eu_lt, 1);
Theta_d_us_lt = log(Theta_d_us_t);
Theta_d_eu_lt = log(Theta_d_eu_t);

% Rate differential
Delta_im = exp(im_us(datesperiod)) - exp(im_eu(datesperiod));
Delta_lsigma = sigma_us_lt - sigma_eu_lt;
Delta_mu = mu_us - mu_eu;
Delta_M = M_us - M_eu;

% Money growth
dlogM_us = [NaN; diff(log(M_us_t))];
dlogM_eu = [NaN; diff(log(M_eu_t))];

% FX returns
fx_ret_t = (exp(inv_e(datesperiod(1:end-1))) ./ exp(inv_e(datesperiod(2:end))) - 1) * abs_scale;
fx_ret_lt = lagmatrix(fx_ret_t, 1);

%% =========================================================================
%  SECTION 1: ARMA for Rate Differential (Delta i_m)
%  =========================================================================

% Check for Econometrics Toolbox
has_econometrics = license('test', 'Econometrics_Toolbox');
if ~has_econometrics
    warning('Econometrics Toolbox not available. Skipping ARIMA sections.');
end

if has_econometrics
try
    fprintf('\n=== ARMA(2,1) for Rate Differential ===\n');
    
    Mdl = arima('ARLags', 1:2, 'MALags', 1, 'Constant', NaN);
    EstMdl = estimate(Mdl, Delta_im);

fprintf('AR coefficients: [%.4f, %.4f]\n', EstMdl.AR{1}, EstMdl.AR{2});
fprintf('Constant: %.6f\n', EstMdl.Constant);

residuals = infer(EstMdl, Delta_im);
fittedValues = Delta_im - residuals;

% Plot: Fitted vs Actual
figure('Name', 'Rate Differential ARMA', 'NumberTitle', 'off');
plot(dates(datesperiod), Delta_im, '-o', 'MarkerSize', 3); hold on;
plot(dates(datesperiod), fittedValues, '-r', 'LineWidth', 1.5);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Actual', 'Fitted', 'Location', 'best', 'Box', 'off');
title('ARMA(2,1) - $\Delta i_m$ Process', 'Interpreter', 'latex', 'FontSize', 18);

% Plot: Residuals
figure('Name', 'Rate Differential Residuals', 'NumberTitle', 'off');
plot(dates(datesperiod), residuals, '-o', 'MarkerSize', 3);
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
title('ARMA(2,1) Residuals - $\Delta i_m$', 'Interpreter', 'latex', 'FontSize', 18);
grid on;

% White noise test
[h, pValue] = lbqtest(residuals);
if h == 0
    fprintf('Residuals ARE white noise (p=%.3f)\n', pValue);
else
    fprintf('Residuals are NOT white noise (p=%.3f)\n', pValue);
end

%% =========================================================================
%  SECTION 2: ARMA for Theta_d (US)
%  =========================================================================
fprintf('\n=== ARMA(2,1) for Theta_d (US) ===\n');

Mdl = arima('ARLags', 1:2, 'MALags', 1, 'Constant', NaN);
EstMdl = estimate(Mdl, Theta_d_us_lt);

fprintf('AR coefficients: [%.4f, %.4f]\n', EstMdl.AR{1}, EstMdl.AR{2});
fprintf('Constant: %.6f\n', EstMdl.Constant);

residuals = infer(EstMdl, Theta_d_us_lt);
fittedValues = Theta_d_us_lt - residuals;

% Plot: Fitted vs Actual
figure('Name', 'Theta_d US ARMA', 'NumberTitle', 'off');
plot(dates(datesperiod), Theta_d_us_lt, '-o', 'MarkerSize', 3); hold on;
plot(dates(datesperiod), fittedValues, '-r', 'LineWidth', 1.5);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Actual', 'Fitted', 'Location', 'best', 'Box', 'off');
title('ARMA(2,1) - $\Theta_{us}$ Process', 'Interpreter', 'latex', 'FontSize', 18);

% Plot: Residuals
figure('Name', 'Theta_d US Residuals', 'NumberTitle', 'off');
plot(dates(datesperiod), residuals, '-o', 'MarkerSize', 3);
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
title('$\Theta_{us}$ Residuals', 'Interpreter', 'latex', 'FontSize', 18);
grid on;

[h, pValue] = lbqtest(residuals);
if h == 0
    fprintf('Residuals ARE white noise (p=%.3f)\n', pValue);
else
    fprintf('Residuals are NOT white noise (p=%.3f)\n', pValue);
end

%% =========================================================================
%  SECTION 3: ARMA for Theta_d (EU)
%  =========================================================================
fprintf('\n=== ARMA(2,2) for Theta_d (EU) ===\n');

Mdl = arima('ARLags', 1:2, 'MALags', 1:2, 'Constant', NaN);
EstMdl = estimate(Mdl, Theta_d_eu_lt);

fprintf('AR coefficients: [%.4f, %.4f]\n', EstMdl.AR{1}, EstMdl.AR{2});
fprintf('Constant: %.6f\n', EstMdl.Constant);

residuals = infer(EstMdl, Theta_d_eu_lt);
fittedValues = Theta_d_eu_lt - residuals;

% Plot: Fitted vs Actual
figure('Name', 'Theta_d EU ARMA', 'NumberTitle', 'off');
plot(dates(datesperiod), Theta_d_eu_lt, '-o', 'MarkerSize', 3); hold on;
plot(dates(datesperiod), fittedValues, '-r', 'LineWidth', 1.5);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Actual', 'Fitted', 'Location', 'best', 'Box', 'off');
title('ARMA(2,2) - $\Theta_{eu}$ Process', 'Interpreter', 'latex', 'FontSize', 18);

% Plot: Residuals
figure('Name', 'Theta_d EU Residuals', 'NumberTitle', 'off');
plot(dates(datesperiod), residuals, '-o', 'MarkerSize', 3);
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
title('$\Theta_{eu}$ Residuals', 'Interpreter', 'latex', 'FontSize', 18);
grid on;

[h, pValue] = lbqtest(residuals);
if h == 0
    fprintf('Residuals ARE white noise (p=%.3f)\n', pValue);
else
    fprintf('Residuals are NOT white noise (p=%.3f)\n', pValue);
end

%% =========================================================================
%  SECTION 4: Money Supply Rule Estimation (US)
%  =========================================================================
fprintf('\n=== Money Supply Rule (US) ===\n');

% Basic regression
[B, ~, residuals, ~, ~] = regress(dlogM_us(2:end), ...
    [ones(length(sigma_us_lag_lt(1:end-1)), 1), sigma_us_lag_lt(2:end)]);
fprintf('OLS: const=%.4f, beta_sigma=%.4f\n', B(1), B(2));

% ARIMAX model
Mdl = arima('Constant', NaN, 'ARLags', 1:2, 'MALags', 1, 'Beta', NaN, 'D', 0);
EstMdl = estimate(Mdl, dlogM_us(4:end), 'X', sigma_us_lag_lt(2:end));
residuals = infer(EstMdl, dlogM_us(4:end), 'X', sigma_us_lag_lt(2:end));
fittedValues = dlogM_us(4:end) - residuals;

% Plot: Fitted vs Actual
figure('Name', 'Money Supply US ARIMAX', 'NumberTitle', 'off');
plot(1:length(fittedValues), dlogM_us(4:end), '-o', 'MarkerSize', 3); hold on;
plot(1:length(fittedValues), fittedValues, '-r', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Value');
title('ARIMAX - US Money Supply Growth', 'FontSize', 14);
legend('Actual', 'Fitted', 'Location', 'best');
grid on;

% Plot: Residuals
figure('Name', 'Money Supply US Residuals', 'NumberTitle', 'off');
plot(dates(datesperiod(4:end)), residuals, '-o', 'MarkerSize', 3);
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
title('$M^{us}$ Residuals', 'Interpreter', 'latex', 'FontSize', 18);
grid on;

[h, pValue] = lbqtest(residuals);
if h == 0
    fprintf('Residuals ARE white noise (p=%.3f)\n', pValue);
else
    fprintf('Residuals are NOT white noise (p=%.3f)\n', pValue);
end

%% =========================================================================
%  SECTION 5: Money Supply Rule Estimation (EU)
%  =========================================================================
fprintf('\n=== Money Supply Rule (EU) ===\n');

% Basic regression
[B, ~, residuals, ~, ~] = regress(dlogM_eu(2:end), ...
    [ones(length(sigma_eu_lag_lt(1:end-1)), 1), sigma_eu_lag_lt(2:end)]);
fprintf('OLS: const=%.4f, beta_sigma=%.4f\n', B(1), B(2));

% ARIMAX model
Mdl = arima('Constant', NaN, 'ARLags', 1:2, 'MALags', 1:12, 'Beta', NaN, 'D', 0);
EstMdl = estimate(Mdl, dlogM_eu(4:end), 'X', sigma_eu_lag_lt(2:end));
residuals = infer(EstMdl, dlogM_eu(4:end), 'X', sigma_eu_lag_lt(2:end));
fittedValues = dlogM_eu(4:end) - residuals;

% Plot: Fitted vs Actual
figure('Name', 'Money Supply EU ARIMAX', 'NumberTitle', 'off');
plot(1:length(fittedValues), dlogM_eu(4:end), '-o', 'MarkerSize', 3); hold on;
plot(1:length(fittedValues), fittedValues, '-r', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Value');
title('ARIMAX - EU Money Supply Growth', 'FontSize', 14);
legend('Actual', 'Fitted', 'Location', 'best');
grid on; axis tight;

% Plot: Residuals
figure('Name', 'Money Supply EU Residuals', 'NumberTitle', 'off');
plot(dates(datesperiod(4:end)), residuals, '-o', 'MarkerSize', 3);
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
title('$M^{eu}$ Residuals', 'Interpreter', 'latex', 'FontSize', 18);
grid on; axis tight;

[h, pValue] = lbqtest(residuals);
if h == 0
    fprintf('Residuals ARE white noise (p=%.3f)\n', pValue);
else
    fprintf('Residuals are NOT white noise (p=%.3f)\n', pValue);
end

%% =========================================================================
%  SECTION 6: Risk Premium VAR
%  =========================================================================
fprintf('\n=== Risk Premium ARIMAX ===\n');

X = [Delta_im, [NaN; dsigma_us_lt], [NaN; dsigma_eu_lt], ...
     Theta_d_us_lt, Theta_d_eu_lt, dlogM_us, dlogM_eu];

Mdl = arima('Constant', NaN, 'ARLags', 1:2, 'MALags', 1:12, 'Beta', NaN(7,1), 'D', 0);
EstMdl = estimate(Mdl, riskprm_t(4:end), 'X', X(2:end, :));
residuals = infer(EstMdl, riskprm_t(4:end), 'X', X(2:end, :));
fittedValues = riskprm_t(4:end) - residuals;

% Plot: Fitted vs Actual
figure('Name', 'Risk Premium ARIMAX', 'NumberTitle', 'off');
plot(1:length(fittedValues), riskprm_t(4:end), '-o', 'MarkerSize', 3); hold on;
plot(1:length(fittedValues), fittedValues, '-r', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Value');
title('ARIMAX - Risk Premium', 'FontSize', 14);
legend('Actual', 'Fitted', 'Location', 'best');
grid on;

% Plot: Residuals
figure('Name', 'Risk Premium Residuals', 'NumberTitle', 'off');
plot(dates(datesperiod(4:end)), residuals, '-o', 'MarkerSize', 3);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
title('$\xi$ Residuals', 'Interpreter', 'latex', 'FontSize', 18);

[h, pValue] = lbqtest(residuals);
if h == 0
    fprintf('Residuals ARE white noise (p=%.3f)\n', pValue);
else
    fprintf('Residuals are NOT white noise (p=%.3f)\n', pValue);
end

catch ME
    warning('ARIMA estimation failed: %s', ME.message);
    fprintf('Skipping remaining ARIMA sections.\n');
end
end  % has_econometrics

%% =========================================================================
%  SECTION 7: Risk Premium Regression (no Econometrics Toolbox required)
%  =========================================================================
fprintf('\n=== Risk Premium Regression ===\n');

X = [ones(length(riskprm_t), 1), Delta_im, sigma_us_lt, sigma_eu_lt, ...
     Theta_d_us_lt, Theta_d_eu_lt, dlogM_us, dlogM_eu];
[B, ~, R, ~, ~] = regress(riskprm_t, X);

riskprm_fit = B(1) + B(2)*Delta_im + B(3)*sigma_us_lt + B(4)*sigma_eu_lt + ...
              B(5)*Theta_d_us_lt + B(6)*Theta_d_eu_lt + B(7)*dlogM_us + B(8)*dlogM_eu;
riskprm_res = riskprm_t - riskprm_fit;

fprintf('Coefficients:\n');
fprintf('  const=%.4f, Delta_im=%.4f, sigma_us=%.4f, sigma_eu=%.4f\n', B(1), B(2), B(3), B(4));
fprintf('  Theta_us=%.4f, Theta_eu=%.4f, dM_us=%.4f, dM_eu=%.4f\n', B(5), B(6), B(7), B(8));

% Plot: Fitted vs Actual
figure('Name', 'Risk Premium Regression', 'NumberTitle', 'off');
plot(dates(datesperiod), riskprm_fit * abs_scale, 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), riskprm_t * abs_scale, 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Fitted', 'Actual', 'Location', 'best', 'Box', 'off');
title('Risk Premium Regression', 'Interpreter', 'latex', 'FontSize', 18);

% Plot: Residuals
figure('Name', 'Risk Premium Regression Residuals', 'NumberTitle', 'off');
plot(dates(datesperiod), riskprm_res * abs_scale, 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Residual', 'Location', 'best', 'Box', 'off');
title('Risk Premium Residual ($\xi$)', 'Interpreter', 'latex', 'FontSize', 18);

%% =========================================================================
%  SECTION 8: FX Predictability
%  =========================================================================
fprintf('\n=== FX Return Predictability ===\n');

% AR(1) test
X = [ones(length(fx_ret_t), 1), fx_ret_lt];
[B, ~, R, ~, ~] = regress(fx_ret_t, X);

if has_econometrics
    [h, pValue] = lbqtest(R);
    if h == 0
        fprintf('AR(1) residuals ARE white noise (p=%.3f)\n', pValue);
    else
        fprintf('AR(1) residuals are NOT white noise (p=%.3f)\n', pValue);
    end
end

% Full predictability regression
X = [ones(length(fx_ret_t), 1), fx_ret_lt, Delta_im(1:end-1), ...
     sigma_us_lt(1:end-1), sigma_eu_lt(1:end-1), ...
     Theta_d_us_lt(1:end-1), Theta_d_eu_lt(1:end-1), ...
     dlogM_us(1:end-1), dlogM_eu(1:end-1)];
[B, ~, R, ~, ~] = regress(fx_ret_t, X);
ret_fit = fx_ret_t - R;

if has_econometrics
    [h, pValue] = lbqtest(R);
    if h == 0
        fprintf('Full model residuals ARE white noise (p=%.3f)\n', pValue);
    else
        fprintf('Full model residuals are NOT white noise (p=%.3f)\n', pValue);
    end
end

% Plot: FX Returns
figure('Name', 'FX Returns', 'NumberTitle', 'off');
plot(dates(datesperiod(1:end-1)), fx_ret_t, 'LineWidth', 2); hold on;
plot(dates(datesperiod(1:end-1)), fx_ret_lt, 'LineWidth', 1, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Actual', 'Lag', 'Location', 'best', 'Box', 'off');
title('FX Returns', 'Interpreter', 'latex', 'FontSize', 18);

% Plot: Scatter
figure('Name', 'FX Returns Scatter', 'NumberTitle', 'off');
scatter(fx_ret_lt, fx_ret_t, 'filled');
grid on; axis tight;
label_x('Lagged Return');
label_y('Realized Return');
formataxis(gca);
title('FX Returns: Lag vs Realized', 'Interpreter', 'latex', 'FontSize', 18);

% Plot: Fitted vs Actual
figure('Name', 'FX Returns vs Fitted', 'NumberTitle', 'off');
plot(dates(datesperiod(1:end-1)), fx_ret_t, 'LineWidth', 1); hold on;
plot(dates(datesperiod(1:end-1)), ret_fit, 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
legend('Actual', 'Predictable', 'Location', 'best', 'Box', 'off');
title('FX Returns: Actual vs Predictable', 'Interpreter', 'latex', 'FontSize', 18);

%% =========================================================================
%  SECTION 9: Forward Market Diagnostics
%  =========================================================================
fprintf('\n=== Forward Market Diagnostics ===\n');

figure('Name', 'Forward Market', 'NumberTitle', 'off');
plot(dates(datesperiod), zeros(1, length(datesperiod)), 'k--', 'LineWidth', 1); hold on;
plot(dates(datesperiod), (f_t(datesperiod) ./ exp(-inv_e(datesperiod)) - 1) * abs_scale, ...
     'b', 'LineWidth', 2);
plot(dates(datesperiod), riskprm_t * abs_scale, 'r:', 'LineWidth', 1.5);
plot(dates(datesperiod), (exp(im_us(datesperiod)) - exp(im_eu(datesperiod))) * abs_scale, ...
     '--', 'LineWidth', 2);
grid on; axis tight;
datetick('x', 'yyyy-mm', 'keeplimits');
label_x('Time (Year-Month)');
label_y('BPS (annual)');
formataxis(gca);
legend('', 'Forward Return', 'Risk Premium', 'Rate Diff', 'Location', 'best', 'Box', 'off');
title('FX Predictable Returns', 'Interpreter', 'latex', 'FontSize', 18);

%% Restore Dynare paths if removed
if restore_dynare
    path(original_path);
    fprintf('Dynare paths restored.\n');
end

fprintf('\n=== Diagnostics Complete ===\n');
