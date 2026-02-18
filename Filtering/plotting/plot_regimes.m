%% plot_regimes.m - Plot Filter Results with Regime Shading
% =========================================================================
% Plots filter results with Markov regime shading.
% Run AFTER markov_estimation.jl generates MS_sigma_us_prob.csv.
%
% Pipeline: main_filter.m → markov_estimation.jl → plot_regimes.m
%
% Required files:
%   - MS_sigma_us_prob.csv (from Julia Markov estimation)
%
% Required variables from workspace:
%   - sigma_us_t, sigma_eu_t, sigma_us_TED_t, sigma_eu_TED_t
%   - BP_us_t, BP_eu_t, TED_us_t, TED_eu_t, CIP_t
%   - theta_us_t, Smin_us_t, DW_us_t, FF_us_t
%   - dates, datesperiod, abs_scale
%
% (c) Saki Bigio
% =========================================================================

%% Settings
printit = 0;  % Set to 1 to save figures
threshold = 0.5;  % Probability threshold for regime classification

% Colors
eu_color = [0.4 0.1 1.0];
data_color = [0.7 0.1 0.4];

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
    printit = 0;
end

% Define formatting if not already defined
if ~exist('formataxis', 'var')
    formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', 14, 'Box', 'On');
end

%% Load Markov Regime Probabilities
prob_file = 'data/MS_sigma_us_prob.csv';
if ~exist(prob_file, 'file')
    warning('MS_sigma_us_prob.csv not found in data/. Run markov_estimation.jl first.');
    fprintf('Skipping regime plots.\n');
    return;
end

% Read probabilities - use readmatrix for compatibility
try
    sigma_us_stateprob = readmatrix(prob_file, 'NumHeaderLines', 1);
catch
    sigma_us_stateprob = csvread(prob_file, 1, 0);
end
sigma_us_stateprob = [sigma_us_stateprob(:,1); 0];  % Append zero for indexing

% Classify regimes
sigma_us_high_t = sigma_us_stateprob < threshold;

% Regime transitions
sigma_us_high2low_t = [(sigma_us_high_t(1:end-1)==1) .* (sigma_us_high_t(2:end)==0); 0];
sigma_us_low2high_t = [(sigma_us_high_t(1:end-1)==0) .* (sigma_us_high_t(2:end)==1); 0];

% Convert dates to datetime if needed
if isnumeric(dates)
    dates = datetime(dates, 'ConvertFrom', 'datenum');
elseif ~isdatetime(dates)
    dates = datetime(dates);
end

%% ========================================================================
%  REGIME PROBABILITY
%  ========================================================================

figure('Name', 'High Liquidity Regime', 'NumberTitle', 'off')
plot(dates(datesperiod), zeros(1,length(datesperiod)), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), sigma_us_stateprob(datesperiod), 'LineWidth', 3, 'Color', 'r'); 
grid on; axis tight;
xtickformat('MMM-yy');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), sigma_us_stateprob(datesperiod), 'LineWidth', 2, 'Color', 'r'); 
title('Prob. Low Funding Risk State', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_sigmaus_states.pdf'), 'ContentType', 'vector');
end

%% ========================================================================
%  SIGMA ESTIMATES WITH REGIMES
%  ========================================================================

figure('Name', 'Sigma US/EU', 'NumberTitle', 'off')
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(sigma_us_t), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), sigma_us_TED_t(datesperiod), 'LineWidth', 3, 'Color', 'r'); 
plot(dates(datesperiod), sigma_eu_TED_t(datesperiod), 'LineWidth', 3, 'LineStyle', '-.', 'Color', eu_color); 
grid on; axis tight;
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), sigma_us_TED_t(datesperiod), 'LineWidth', 3, 'Color', 'r'); 
plot(dates(datesperiod), sigma_eu_TED_t(datesperiod), 'LineWidth', 3, 'LineStyle', '-.', 'Color', eu_color); 
legend('', 'US', 'EU', 'box', 'off', 'color', 'none');
xtickformat('MMM-yy');
formataxis(gca);
title('$\sigma$ Estimates', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_sigmauseu.pdf'), 'ContentType', 'vector');
end

%% ========================================================================
%  TED SPREADS
%  ========================================================================

figure('Name', 'TED Comparison', 'NumberTitle', 'off') 
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(TED_us_t)*abs_scale, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), TED_us_t(datesperiod)*abs_scale, 'LineWidth', 3, 'Color', 'r');
plot(dates(datesperiod), TED_eu_t(datesperiod)*abs_scale, 'LineWidth', 3, 'LineStyle', '-.', 'Color', eu_color);
grid on; axis tight;
xtickformat('MMM-yy');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), TED_us_t(datesperiod)*abs_scale, 'LineWidth', 3, 'Color', 'r');
plot(dates(datesperiod), TED_eu_t(datesperiod)*abs_scale, 'LineWidth', 3, 'LineStyle', '-.', 'Color', eu_color);
legend('', 'US', 'EU', 'box', 'off', 'color', 'none');
title('TED Spreads (bps)', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_Tedtargets.pdf'), 'ContentType', 'vector');
end

%% ========================================================================
%  BOND PREMIA (MODEL VS DATA)
%  ========================================================================

% Handle NaN values
NaNindex = (Rb_Rm == 0);
Rb_Rm_plot = Rb_Rm;
Rb_Rm_plot(NaNindex) = NaN;

figure('Name', 'US Bond Premium', 'NumberTitle', 'off') 
plot(dates(datesperiod), zeros(1,length(datesperiod)), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), (BP_us_t(datesperiod)-mean(BP_us_t))*abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), (Rb_Rm_plot(datesperiod)-mean(Rb_Rm_plot(~isnan(Rb_Rm_plot))))*abs_scale, 'LineWidth', 3, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), (BP_us_t(datesperiod)-mean(BP_us_t))*abs_scale, 'LineWidth', 3, 'Color', 'r');
plot(dates(datesperiod), (Rb_Rm_plot(datesperiod)-mean(Rb_Rm_plot(~isnan(Rb_Rm_plot))))*abs_scale, 'LineWidth', 3, 'LineStyle', '-.', 'Color', data_color);
legend('', 'Model', 'Data', 'box', 'off', 'color', 'none');
xtickformat('MMM-yy');
formataxis(gca);
title('US Bond Premium (bps, demeaned)', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_BPus_fit.pdf'), 'ContentType', 'vector');
end

% EU Bond Premium
NaNindex = (Rb_Rm_eu == 0);
Rb_Rm_eu_plot = Rb_Rm_eu;
Rb_Rm_eu_plot(NaNindex) = NaN;

figure('Name', 'EU Bond Premium', 'NumberTitle', 'off') 
plot(dates(datesperiod), zeros(1,length(datesperiod)), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), (BP_eu_t(datesperiod)-mean(BP_eu_t))*abs_scale, 'LineWidth', 3);
plot(dates(datesperiod), (Rb_Rm_eu_plot(datesperiod)-mean(Rb_Rm_eu_plot(~isnan(Rb_Rm_eu_plot))))*abs_scale, 'LineWidth', 3, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), (BP_eu_t(datesperiod)-mean(BP_eu_t))*abs_scale, 'LineWidth', 3, 'Color', 'r');
plot(dates(datesperiod), (Rb_Rm_eu_plot(datesperiod)-mean(Rb_Rm_eu_plot(~isnan(Rb_Rm_eu_plot))))*abs_scale, 'LineWidth', 3, 'LineStyle', '-.', 'Color', data_color);
legend('', 'Model', 'Data', 'box', 'off', 'color', 'none');
xtickformat('MMM-yy');
formataxis(gca);
title('EU Bond Premium (bps, demeaned)', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_BPeu_fit.pdf'), 'ContentType', 'vector');
end

%% ========================================================================
%  CIP DEVIATION
%  ========================================================================

figure('Name', 'CIP Deviation', 'NumberTitle', 'off') 
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(CIP_t(datesperiod))*abs_scale, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), cip(datesperiod)*abs_scale, 'LineWidth', 3, 'LineStyle', '-');
grid on; axis tight;
xtickformat('MMM-yy');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), cip(datesperiod)*abs_scale, 'LineWidth', 3, 'LineStyle', '-', 'Color', data_color);
title('CIP Deviation (bps)', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportfig(gcf, [foldername 'F_CIP'], 'color', 'cmyk', 'resolution', 1600);
end

%% ========================================================================
%  INTERBANK DISPERSION
%  ========================================================================

figure('Name', 'Interbank Dispersion', 'NumberTitle', 'off') 
plot(dates(datesperiod), ones(1,length(datesperiod))*mean(Chi_D_US(datesperiod))*abs_scale, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), Chi_D_US(datesperiod)*abs_scale, 'LineWidth', 3, 'LineStyle', '-');
grid on; axis tight;
xtickformat('MMM-yy');
formataxis(gca); 
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), Chi_D_US(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', '-', 'Color', data_color);
title('US Interbank Dispersion (bps)', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_IBdispersion.pdf'), 'ContentType', 'vector');
end

%% ========================================================================
%  FUNDING COMPOSITION
%  ========================================================================

figure('Name', 'Market Tightness', 'NumberTitle', 'off')  
plot(dates(datesperiod), log(theta_us_t(datesperiod)), 'LineWidth', 3); hold on;
plot(dates(datesperiod), log(theta_eu_t(datesperiod)), 'LineWidth', 3); hold off;
grid on; axis tight;
xtickformat('MMM-yy');
formataxis(gca);
legend('US', 'EU', 'Location', 'best');
title('Log Market Tightness $\theta$', 'interpreter', 'latex', 'fontsize', 16);

figure('Name', 'Funding Composition', 'NumberTitle', 'off')  
plot(dates(datesperiod), Smin_us_t(datesperiod), 'LineWidth', 3); hold on;
plot(dates(datesperiod), DW_us_t(datesperiod), 'LineWidth', 3); 
plot(dates(datesperiod), FF_us_t(datesperiod), 'LineWidth', 1); hold off;
grid on; axis tight;
xtickformat('MMM-yy');
formataxis(gca);
legend('Funding Need', 'DW Funded', 'FF Funded', 'Location', 'best');
title('US Funding Composition', 'interpreter', 'latex', 'fontsize', 16);

%% ========================================================================
%  DW BORROWING (MODEL VS DATA)
%  ========================================================================

figure('Name', 'DW Borrowing', 'NumberTitle', 'off') 
mod_scale  = mean(DW_us_t(datesperiod)./M_us(datesperiod)');
data_scale = mean(DW_t(datesperiod)./M_us(datesperiod));
plot(dates(datesperiod), (DW_us_t(datesperiod)./M_us(datesperiod)')/mod_scale, 'LineWidth', 3, 'Color', 'r'); hold on;
plot(dates(datesperiod), (DW_t(datesperiod)./M_us(datesperiod))/data_scale, 'LineWidth', 3, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), (DW_us_t(datesperiod)./M_us(datesperiod)')/mod_scale, 'LineWidth', 3, 'Color', 'r');
plot(dates(datesperiod), (DW_t(datesperiod)./M_us(datesperiod))/data_scale, 'LineWidth', 3, 'LineStyle', '-.', 'Color', data_color);
xtickformat('MMM-yy');
formataxis(gca); 
legend('Model', 'Data', 'box', 'off', 'color', 'none');
title('DW Borrowing / Reserves (normalized)', 'interpreter', 'latex', 'fontsize', 16);
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_DWus_fit.pdf'), 'ContentType', 'vector');
end

%% ========================================================================
%  MULTI-CURRENCY CIP WITH REGIMES
%  ========================================================================

desiredNumXTicks = 3;

figure('Name', 'CIP All Currencies', 'NumberTitle', 'off') 
subplot(3, 3, 1); 
plot(dates(datesperiod), CIP_s_eu_t(datesperiod)*abs_scale, 'LineWidth', 2); hold on;
xtickformat('MM-yy');
grid on; axis tight;
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), CIP_s_eu_t(datesperiod)*abs_scale, 'LineWidth', 2, 'Color', data_color);
ax = gca;
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
title('EU', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    CIP_data = eval(['CIP_s_' curlist{j} '_t(datesperiod)']);
    subplot(3, 3, j+1);
    plot(dates(datesperiod), CIP_data*abs_scale, 'LineWidth', 2); hold on;
    regime_patches(dates, datesperiod, sigma_us_high_t);
    plot(dates(datesperiod), CIP_data*abs_scale, 'LineWidth', 2, 'Color', data_color);
    xtickformat('MM-yy');
    grid on; axis tight;
    ax = gca;
    xLimits = get(ax, 'XLim');
    customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
    set(ax, 'XTick', customXTicks);
    title(conlist{j}, 'interpreter', 'latex', 'Fontsize', 15);
end
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_CIP_all.pdf'), 'ContentType', 'vector');
end

%% ========================================================================
%  FX DEVALUATION WITH REGIMES
%  ========================================================================

figure('Name', 'FX Devaluation All', 'NumberTitle', 'off')
datesperiod_f = datesperiod(2:end);
datesperiod_l = datesperiod(1:end-1);

subplot(3, 3, 1);
temp2 = (ln_eu_us_t(datesperiod_f) - ln_eu_us_t(datesperiod_l)) * abs_scale/100;
plot(dates(datesperiod_f), temp2, 'LineWidth', 1, 'LineStyle', '-', 'Color', data_color); hold on;
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod_f), temp2, 'LineWidth', 1, 'LineStyle', '-', 'Color', data_color);
xtickformat('MM-yy');
ax = gca;
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
grid on; axis tight;
title('EU/USD', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    ln_data = eval(['ln_' curlist{j} '_us_t']);
    temp2 = (ln_data(datesperiod_f) - ln_data(datesperiod_l)) * abs_scale/100;
    plot(dates(datesperiod_f), temp2, 'LineWidth', 1, 'LineStyle', '-', 'Color', data_color); hold on;
    regime_patches(dates, datesperiod, sigma_us_high_t);
    plot(dates(datesperiod_f), temp2, 'LineWidth', 1, 'LineStyle', '-', 'Color', data_color);
    xtickformat('MM-yy');
    grid on; axis tight;
    ax = gca;
    xLimits = get(ax, 'XLim');
    customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
    set(ax, 'XTick', customXTicks);
    title([conlist{j} '/USD'], 'interpreter', 'latex', 'Fontsize', 15);
end
if printit == 1
    exportgraphics(gcf, fullfile(foldername, 'F_devaluation_all.pdf'), 'ContentType', 'vector');
end

%% ========================================================================
%  REGIME MOMENTS (EXPORT TO LATEX)
%  ========================================================================

% Regime indices
index_r2 = sigma_us_high_t;
index_r1 = ~index_r2;
dev_t = ((inv_e_t(2:end)./inv_e_t(1:end-1)).^(-1) - 1) * abs_scale;

% Average regime change
dev_low2high = mean(dev_t(sigma_us_low2high_t == 1));
dev_high2low = mean(dev_t(sigma_us_high2low_t == 1));

% Unconditional moments
E_bp_data = mean(BP_us_t) * abs_scale;
E_cip_data = mean(CIP_s_eu_t) * abs_scale;
E_dev_data = mean(dev_t);

std_bp_data = std(BP_us_t) * abs_scale;
std_cip_data = std(CIP_s_eu_t) * abs_scale;
std_dev_data = std(dev_t);

aux = autocorr(BP_us_t); rho_bp_data = aux(2);
aux = autocorr(CIP_s_eu_t); rho_cip_data = aux(2);
aux = autocorr(dev_t); rho_dev_data = aux(2);

% Conditional moments (Regime 1 = Low risk)
E_bp_data_r1 = mean(BP_us_t(index_r1)) * abs_scale;
E_cip_data_r1 = mean(CIP_s_eu_t(index_r1)) * abs_scale;
E_dev_data_r1 = mean(dev_t(index_r1(1:end-1)));

std_bp_data_r1 = std(BP_us_t(index_r1)) * abs_scale;
std_cip_data_r1 = std(CIP_s_eu_t(index_r1)) * abs_scale;
std_dev_data_r1 = std(dev_t(index_r1(1:end-1)));

aux = autocorr(BP_us_t(index_r1)); rho_bp_data_r1 = aux(2);
aux = autocorr(CIP_s_eu_t(index_r1)); rho_cip_data_r1 = aux(2);
aux = autocorr(dev_t(index_r1(1:end-1))); rho_dev_data_r1 = aux(2);

% Conditional moments (Regime 2 = High risk)
E_bp_data_r2 = mean(BP_us_t(index_r2)) * abs_scale;
E_cip_data_r2 = mean(CIP_s_eu_t(index_r2)) * abs_scale;
E_dev_data_r2 = mean(dev_t(index_r2(1:end-1)));

std_bp_data_r2 = std(BP_us_t(index_r2)) * abs_scale;
std_cip_data_r2 = std(CIP_s_eu_t(index_r2)) * abs_scale;
std_dev_data_r2 = std(dev_t(index_r2(1:end-1)));

aux = autocorr(BP_us_t(index_r2)); rho_bp_data_r2 = aux(2);
aux = autocorr(CIP_s_eu_t(index_r2)); rho_cip_data_r2 = aux(2);
aux = autocorr(dev_t(index_r2(1:end-1))); rho_dev_data_r2 = aux(2);

% Print summary
fprintf('\n=== Regime Moments ===\n');
fprintf('                    Uncond    Low-Risk   High-Risk\n');
fprintf('BP mean (bps):      %6.1f    %6.1f     %6.1f\n', E_bp_data, E_bp_data_r1, E_bp_data_r2);
fprintf('CIP mean (bps):     %6.1f    %6.1f     %6.1f\n', E_cip_data, E_cip_data_r1, E_cip_data_r2);
fprintf('Deval mean (bps):   %6.1f    %6.1f     %6.1f\n', E_dev_data, E_dev_data_r1, E_dev_data_r2);

% Save to LaTeX (if printing)
if printit == 1
    filename = fullfile(foldername, 'Data_CIP_Moments.tex');
    fid = fopen(filename, 'wt');
    fprintf(fid, 'CIP (data) & %.1f & %.2f & %.1f & %.1f & \\{ %.2f, %.2f \\} & %.1f \\\\ \n', ...
        E_cip_data, rho_cip_data, std_cip_data, E_cip_data_r2-E_cip_data_r1, ...
        rho_cip_data_r2, rho_cip_data_r1, std_cip_data_r2/std_cip_data_r1);
    fclose(fid);
    fprintf('Saved: %s\n', filename);
end

%% Summary
fprintf('plot_regimes complete. 13 figures generated.\n');

%% ========================================================================
%  HELPER FUNCTION
%  ========================================================================

function out = regime_patches(dates, datesperiod, sigma_us_high_t)
    % Add shaded patches for high-risk regime periods
    y_limits = ylim;
    for i = 1:length(datesperiod) - 1
        if sigma_us_high_t(i) == 1
            x_patch = [dates(datesperiod(i)), dates(datesperiod(i + 1)), ...
                       dates(datesperiod(i + 1)), dates(datesperiod(i))];
            y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
            patch(x_patch, y_patch, [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        end
    end
    out = [];
end
