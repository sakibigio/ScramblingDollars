%% plot_baseline.m - Plot Filter Results
% =========================================================================
% Plots the baseline results from the filtering exercise.
% Called by main_filter.m after the filter completes.
%
% Required variables from workspace:
%   - sigma_us_t, sigma_eu_t, sigma_us_TED_t, sigma_eu_TED_t
%   - dates, datesperiod, abs_scale
%   - Various model outputs (BP_us_t, TED_us_t, theta_us_t, etc.)
%
% (c) Saki Bigio
% =========================================================================

%% Settings
% Set output folder based on machine
foldername = overleaf_dir();

% Define formatting if not already defined
if ~exist('formataxis', 'var')
    formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', 14, 'Box', 'On');
end

if ~exist('desiredDecimalPlaces', 'var')
    desiredDecimalPlaces = 2;
end

if ~exist('printit', 'var')
    printit = 0;
end

if ~exist('matching_type', 'var')
    matching_type = 1;  % default to Cobb-Douglas
end
mt_suffix = '_cd';
if matching_type == 0
    mt_suffix = '_l';
end

if ~exist('label_y', 'var')
    label_y = @(x) ylabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', 14, 'Interpreter', 'latex');
end

%% ========================================================================
%  SIGMA PLOTS
%  ========================================================================

%% Figure 1: Sigma US (TED target)
figure('Name', 'Sigma US (TED)', 'NumberTitle', 'off')
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(sigma_us_t), ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), sigma_us_TED_t(datesperiod), 'LineWidth', 3); 
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
title('Estimated $\sigma^*$ (US TED)', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportfig(gcf, [foldername 'F_sigmaus' mt_suffix], 'color', 'cmyk', 'resolution', 1600);
end

%% Figure 2: Sigma EU (TED target)
figure('Name', 'Sigma EU (TED)', 'NumberTitle', 'off')
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(sigma_eu_t), ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), sigma_eu_TED_t(datesperiod), 'LineWidth', 3); 
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
title('Estimated $\sigma$ (EU TED)', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportfig(gcf, [foldername 'F_sigmaeu' mt_suffix], 'color', 'cmyk', 'resolution', 1600);
end

%% Figure 3: Log Sigma Comparison (US vs EU)
figure('Name', 'Log Sigma Comparison', 'NumberTitle', 'off')
plot(dates(datesperiod), ones(1,length(datesperiod))*log(mean(sigma_us_t)), ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), log(sigma_us_TED_t(datesperiod)), 'LineStyle', ':', 'LineWidth', 2); 
plot(dates(datesperiod), log(sigma_eu_TED_t(datesperiod)), 'LineWidth', 3); 
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('average', 'US', 'EU', 'Location', 'best');
title('Log $\sigma$ Comparison', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportfig(gcf, [foldername 'F_useu_logcomp' mt_suffix], 'color', 'cmyk', 'resolution', 1600);
end

%% Figure 4: Sigma Comparison (US vs EU)
figure('Name', 'Sigma Comparison', 'NumberTitle', 'off')
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(sigma_us_t), ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), sigma_us_TED_t(datesperiod), 'LineStyle', ':', 'LineWidth', 2); 
plot(dates(datesperiod), sigma_eu_TED_t(datesperiod), 'LineWidth', 3); 
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('average', 'US', 'EU', 'Location', 'best');
title('$\sigma$ Comparison', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportfig(gcf, [foldername 'F_useu_comp' mt_suffix], 'color', 'cmyk', 'resolution', 1600);
end

%% Figure 5: Solver Diagnostics (Flags)
figure('Name', 'Solver Diagnostics', 'NumberTitle', 'off')
plot(dates(datesperiod), sigma_us_TED_flag(datesperiod), 'LineWidth', 1, 'LineStyle', '-', 'Color', 'b'); hold on;
plot(dates(datesperiod), sigma_eu_TED_flag(datesperiod), 'LineWidth', 1);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
legend('sigma US flag', 'sigma EU flag', 'Location', 'best');
formataxis(gca);
title('Solver Diagnostics', 'interpreter', 'latex', 'fontsize', 16);

%% ========================================================================
%  SIGMA vs MU CO-MOVEMENT DIAGNOSTICS
%  ========================================================================
% Checks whether the filtered sigma series tracks the input liquidity ratio
% mu too closely (in levels and in first differences).

zscore_local = @(x) (x - mean(x,'omitnan')) ./ std(x,'omitnan');

%% Figure 5b: Sigma vs Mu (US)
ls_us = log(sigma_us_t(datesperiod));
mu_us_d = mu_us(datesperiod);
dls_us = diff(ls_us);
dmu_us = diff(mu_us_d);
rho_lvl_us  = corr_complete(ls_us, mu_us_d);
rho_diff_us = corr_complete(dls_us, dmu_us);

figure('Name', 'Sigma vs Mu (US)', 'NumberTitle', 'off')
subplot(1,3,1)
plot(dates(datesperiod), zscore_local(ls_us), 'LineWidth', 2); hold on;
plot(dates(datesperiod), zscore_local(mu_us_d), 'LineWidth', 2, 'LineStyle', '--');
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
legend('$\log\sigma_{us}$ (z)', '$\mu_{us}$ (z)', 'Location', 'best', 'interpreter', 'latex');
formataxis(gca);
title('Z-scored Time Series', 'interpreter', 'latex', 'fontsize', 14);

subplot(1,3,2)
scatter(mu_us_d, ls_us, 25, 'filled'); hold on;
p = polyfit(mu_us_d, ls_us, 1);
xfit = linspace(min(mu_us_d), max(mu_us_d), 50);
plot(xfit, polyval(p, xfit), 'r-', 'LineWidth', 2);
grid on; axis tight;
label_x('$\mu_{us}$'); label_y('$\log\sigma_{us}$');
formataxis(gca);
title(sprintf('Levels: $\\rho=%.3f$', rho_lvl_us), 'interpreter', 'latex', 'fontsize', 14);

subplot(1,3,3)
scatter(dmu_us, dls_us, 25, 'filled'); hold on;
p = polyfit(dmu_us, dls_us, 1);
xfit = linspace(min(dmu_us), max(dmu_us), 50);
plot(xfit, polyval(p, xfit), 'r-', 'LineWidth', 2);
grid on; axis tight;
label_x('$\Delta\mu_{us}$'); label_y('$\Delta\log\sigma_{us}$');
formataxis(gca);
title(sprintf('First diffs: $\\rho=%.3f$', rho_diff_us), 'interpreter', 'latex', 'fontsize', 14);

%% Figure 5c: Sigma vs Mu (EU)
ls_eu = log(sigma_eu_t(datesperiod));
mu_eu_d = mu_eu(datesperiod);
dls_eu = diff(ls_eu);
dmu_eu = diff(mu_eu_d);
rho_lvl_eu  = corr_complete(ls_eu, mu_eu_d);
rho_diff_eu = corr_complete(dls_eu, dmu_eu);

figure('Name', 'Sigma vs Mu (EU)', 'NumberTitle', 'off')
subplot(1,3,1)
plot(dates(datesperiod), zscore_local(ls_eu), 'LineWidth', 2); hold on;
plot(dates(datesperiod), zscore_local(mu_eu_d), 'LineWidth', 2, 'LineStyle', '--');
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
legend('$\log\sigma_{eu}$ (z)', '$\mu_{eu}$ (z)', 'Location', 'best', 'interpreter', 'latex');
formataxis(gca);
title('Z-scored Time Series', 'interpreter', 'latex', 'fontsize', 14);

subplot(1,3,2)
scatter(mu_eu_d, ls_eu, 25, 'filled'); hold on;
p = polyfit(mu_eu_d, ls_eu, 1);
xfit = linspace(min(mu_eu_d), max(mu_eu_d), 50);
plot(xfit, polyval(p, xfit), 'r-', 'LineWidth', 2);
grid on; axis tight;
label_x('$\mu_{eu}$'); label_y('$\log\sigma_{eu}$');
formataxis(gca);
title(sprintf('Levels: $\\rho=%.3f$', rho_lvl_eu), 'interpreter', 'latex', 'fontsize', 14);

subplot(1,3,3)
scatter(dmu_eu, dls_eu, 25, 'filled'); hold on;
p = polyfit(dmu_eu, dls_eu, 1);
xfit = linspace(min(dmu_eu), max(dmu_eu), 50);
plot(xfit, polyval(p, xfit), 'r-', 'LineWidth', 2);
grid on; axis tight;
label_x('$\Delta\mu_{eu}$'); label_y('$\Delta\log\sigma_{eu}$');
formataxis(gca);
title(sprintf('First diffs: $\\rho=%.3f$', rho_diff_eu), 'interpreter', 'latex', 'fontsize', 14);

fprintf('\n--- Sigma vs Mu co-movement ---\n');
fprintf('US: corr(log sigma, mu) = %.3f (levels), %.3f (diffs)\n', rho_lvl_us, rho_diff_us);
fprintf('EU: corr(log sigma, mu) = %.3f (levels), %.3f (diffs)\n', rho_lvl_eu, rho_diff_eu);

%% ========================================================================
%  RISK PREMIA PLOTS
%  ========================================================================

%% Figure 6: Risk Premium
figure('Name', 'Risk Premia', 'NumberTitle', 'off') 
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(riskprm_t)*abs_scale, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), riskprm_t(datesperiod)*abs_scale, 'LineWidth', 3);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
title('Risk Premium $\xi$ (bps)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 7: Risk Premium vs Policy Rate Differential
figure('Name', 'Risk Premia vs Rate Diff', 'NumberTitle', 'off') 
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(riskprm_t)*abs_scale, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), riskprm_t(datesperiod)*abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), (exp(im_us(datesperiod))-exp(im_eu(datesperiod)))*abs_scale, ...
    'LineWidth', 2, 'LineStyle', ':');
grid on; axis tight;
legend('', 'model', 'rate diff (data)', 'Location', 'best');
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
title('$\xi$ vs Rate Differentials', 'interpreter', 'latex', 'fontsize', 16);

%% ========================================================================
%  MODEL VS DATA PLOTS
%  ========================================================================

%% Figure 8: Bond Premium US (Model vs Data)
figure('Name', 'BP US: Model vs Data', 'NumberTitle', 'off') 
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(BP_us_t(datesperiod))*abs_scale, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), BP_us_t(datesperiod)*abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), Rb_Rm(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('', 'model', 'data', 'Location', 'best');
title('Bond Premium US (bps)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 9: Bond Premium EU (Model vs Data)
figure('Name', 'BP EU: Model vs Data', 'NumberTitle', 'off') 
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(BP_eu_t(datesperiod))*abs_scale, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), BP_eu_t(datesperiod)*abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), Rb_Rm_eu(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('', 'model', 'data', 'Location', 'best');
title('Bond Premium EU (bps)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 10a: TED Spread US (Model vs Data)
figure('Name', 'TED US: Model vs Data', 'NumberTitle', 'off')
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(TED_us_t)*abs_scale, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), TED_us_t(datesperiod)*abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), TED_s_us_t(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('', 'model', 'data', 'Location', 'best');
title('TED Spread US (bps)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 10b: TED Spread EU (Model vs Data)
figure('Name', 'TED EU: Model vs Data', 'NumberTitle', 'off')
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(TED_eu_t)*abs_scale, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), TED_eu_t(datesperiod)*abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), TED_s_eu_t(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('', 'model', 'data', 'Location', 'best');
title('TED Spread EU (bps)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 11: CIP Deviation (Model vs Data)
figure('Name', 'CIP: Model vs Data', 'NumberTitle', 'off') 
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*abs_scale, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), CIP_t(datesperiod)*abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), cip(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('', 'model', 'data', 'Location', 'best');
title('CIP Deviation (bps)', 'interpreter', 'latex', 'fontsize', 16);

%% ========================================================================
%  MARKET STRUCTURE PLOTS
%  ========================================================================

%% Figure 12: Market Tightness
figure('Name', 'Market Tightness', 'NumberTitle', 'off')  
plot(dates(datesperiod), log(theta_us_t(datesperiod)), 'LineWidth', 3); hold on;
plot(dates(datesperiod), log(theta_eu_t(datesperiod)), 'LineWidth', 3); 
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('US', 'EU', 'Location', 'best');
title('Log Market Tightness $\theta$', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 13: Funding Composition US
figure('Name', 'Funding US', 'NumberTitle', 'off')  
plot(dates(datesperiod), Smin_us_t(datesperiod), 'LineWidth', 3); hold on;
plot(dates(datesperiod), DW_us_t(datesperiod), 'LineWidth', 3); 
plot(dates(datesperiod), FF_us_t(datesperiod), 'LineWidth', 1); 
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('Funding Need', 'DW Funded', 'FF Funded', 'Location', 'best');
title('Funding Composition (US)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 14: Funding Composition EU
figure('Name', 'Funding EU', 'NumberTitle', 'off')  
plot(dates(datesperiod), Smin_eu_t(datesperiod), 'LineWidth', 3); hold on;
plot(dates(datesperiod), DW_eu_t(datesperiod), 'LineWidth', 3); 
plot(dates(datesperiod), FF_eu_t(datesperiod), 'LineWidth', 1);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('Funding Need', 'DW Funded', 'FF Funded', 'Location', 'best');
title('Funding Composition (EU)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 15: DW/FF Ratio
figure('Name', 'DW/FF Ratio', 'NumberTitle', 'off')  
plot(dates(datesperiod), DW_us_t(datesperiod)./FF_us_t(datesperiod)*100, 'LineWidth', 3); hold on;
plot(dates(datesperiod), DW_eu_t(datesperiod)./FF_eu_t(datesperiod)*100, 'LineWidth', 3); 
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('US', 'EU', 'Location', 'best');
title('DW/FF Ratio (\%)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 16: DW Borrowing (DISABLED -- flow vs stock mismatch)
% Replaced by Figure 17 (stock vs stock).
% figure('Name', 'DW Borrowing: Model vs Data', 'NumberTitle', 'off')
% plot(dates(datesperiod), DW_us_t(datesperiod)*100, 'LineWidth', 3); hold on;
% plot(dates(datesperiod), DW_n(datesperiod)*100, 'LineWidth', 2, 'LineStyle', ':');
% ylabel('% of checkable deposits');
% datetick('x', 'mmm-yy', 'keeplimits');
% formataxis(gca);
% legend('model', 'data', 'Location', 'best');
% title('DW Borrowing/Deposits (flow)', 'interpreter', 'latex', 'fontsize', 16);

%% Figure 17: DW Loan STOCK (Model vs Data)
% DWS_us_t built in main_filter.m as DWS_t = DW_t + (1-delta)*DWS_{t-1}.
% Data DW_n is already an outstanding-stock series, so this is the like-with-like
% comparison.  Plot in percent of checkable deposits.
if exist('DWS_us_t', 'var')
    figure('Name', 'DW Loan Stock: Model vs Data', 'NumberTitle', 'off')
    plot(dates(datesperiod), DWS_us_t(datesperiod)*100, 'LineWidth', 3); hold on;
    plot(dates(datesperiod), DW_n(datesperiod)*100, 'LineWidth', 2, 'LineStyle', ':');
    ylabel('% of checkable deposits');
    datetick('x', 'mmm-yy', 'keeplimits');
    formataxis(gca);
    legend(sprintf('model stock (\\delta=%.2f)', delta_DWS), 'data', 'Location', 'best');
    title('DW Loan Stock', 'interpreter', 'latex', 'fontsize', 16);
end

%% Figure 17b: EU DW STOCK (Model vs Data, dual y-axes)
% Data DW_eu_n incomplete & on a different scale -> yyaxis decouples levels.
% Stock-vs-stock only (flow plot removed).
if exist('DWS_eu_t', 'var')
    figure('Name', 'EU DW Stock: Model vs Data', 'NumberTitle', 'off')
    if exist('DW_eu_n', 'var')
        yyaxis left
        h_dat = plot(dates(datesperiod), DW_eu_n(datesperiod)*100, ...
            'LineWidth', 2, 'LineStyle', ':');
        ylabel('Data: % of deposits');

        yyaxis right
        h_stock = plot(dates(datesperiod), DWS_eu_t(datesperiod)*100, ...
            'LineWidth', 3);
        ylabel('Model stock: % of deposits');
        legend([h_stock h_dat], ...
            sprintf('model stock (\\delta=%.2f)', delta_DWS), 'data', ...
            'Location', 'best');
    else
        h_stock = plot(dates(datesperiod), DWS_eu_t(datesperiod)*100, ...
            'LineWidth', 3);
        ylabel('% of deposits');
        legend(h_stock, sprintf('model stock (\\delta=%.2f)', delta_DWS), ...
            'Location', 'best');
    end
    datetick('x', 'mmm-yy', 'keeplimits');
    formataxis(gca);
    title('EU DW Stock/Deposits', 'interpreter', 'latex', 'fontsize', 16);
end

%% Figures 18-19: Normalize by reserves (mu = M/D) instead of deposits
% Since DW_n = DW/D and mu = M/D, DW_n/mu = DW/M (per unit reserves).
% Same for FF_n.  On the model side, DW_us_t/mu_us_t and FF_us_t/mu_us_t.
mu_us_t_norm = exp(mu_us);     % reserves/deposits ratio

% Figure 18: DW / Reserves  (flow + stock, model vs data)
figure('Name', 'DW / Reserves: Model vs Data', 'NumberTitle', 'off')
h_flow = plot(dates(datesperiod), DW_us_t(datesperiod) ./ mu_us_t_norm(datesperiod) * 100, ...
    'LineWidth', 3); hold on;
if exist('DWS_us_t', 'var')
    h_stock = plot(dates(datesperiod), DWS_us_t(datesperiod) ./ mu_us_t_norm(datesperiod) * 100, ...
        'LineWidth', 3, 'LineStyle', '--');
end
h_data = plot(dates(datesperiod), DW_n(datesperiod) ./ mu_us_t_norm(datesperiod) * 100, ...
    'LineWidth', 2, 'LineStyle', ':');
ylabel('% of reserves');
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
if exist('DWS_us_t', 'var')
    legend([h_flow h_stock h_data], 'model flow', ...
        sprintf('model stock (\\delta=%.2f)', delta_DWS), 'data', 'Location', 'best');
else
    legend([h_flow h_data], 'model', 'data', 'Location', 'best');
end
title('DW / Reserves (M)', 'interpreter', 'latex', 'fontsize', 16);

% Figure 19: FF / Reserves  (model vs data)
figure('Name', 'FF / Reserves: Model vs Data', 'NumberTitle', 'off')
plot(dates(datesperiod), FF_us_t(datesperiod) ./ mu_us_t_norm(datesperiod) * 100, ...
    'LineWidth', 3); hold on;
plot(dates(datesperiod), FF_n(datesperiod) ./ mu_us_t_norm(datesperiod) * 100, ...
    'LineWidth', 2, 'LineStyle', ':');
ylabel('% of reserves');
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('model', 'data', 'Location', 'best');
title('FF / Reserves (M)', 'interpreter', 'latex', 'fontsize', 16);

%% Summary
fprintf('plot_baseline complete. 20 figures generated.\n');
