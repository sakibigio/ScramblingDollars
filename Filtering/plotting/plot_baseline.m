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
[~, username] = system('whoami');
username = strtrim(username);
if strcmp(username, 'sakibigio')
    foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollars_Revision_Restud/quantfigs/';
elseif strcmp(username, 'sakiclaudia')
    foldername = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/quantfigs/';
else
    warning('Unknown user: %s. Figures will not be saved.', username);
    foldername = './';
    if ~exist('printit', 'var')
        printit = 0;
    end
end

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
    exportfig(gcf, [foldername 'F_sigmaus'], 'color', 'cmyk', 'resolution', 1600);
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
    exportfig(gcf, [foldername 'F_sigmaeu'], 'color', 'cmyk', 'resolution', 1600);
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
    exportfig(gcf, [foldername 'F_useu_logcomp'], 'color', 'cmyk', 'resolution', 1600);
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
    exportfig(gcf, [foldername 'F_useu_comp'], 'color', 'cmyk', 'resolution', 1600);
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

%% Figure 10: TED Spread US (Model vs Data)
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

%% Figure 16: DW Borrowing (Model vs Data)
figure('Name', 'DW Borrowing: Model vs Data', 'NumberTitle', 'off') 
mod_scale  = mean(DW_us_t(datesperiod)./M_us(datesperiod)');
data_scale = mean(DW_t(datesperiod)./M_us(datesperiod));
plot(dates(datesperiod), (DW_us_t(datesperiod)./M_us(datesperiod)')/mod_scale, 'LineWidth', 3); hold on;
plot(dates(datesperiod), (DW_t(datesperiod)./M_us(datesperiod))/data_scale, 'LineWidth', 2, 'LineStyle', ':');
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
legend('model', 'data', 'Location', 'best');
title('DW Borrowing/Reserves (normalized)', 'interpreter', 'latex', 'fontsize', 16);

%% Summary
fprintf('plot_baseline complete. 16 figures generated.\n');
