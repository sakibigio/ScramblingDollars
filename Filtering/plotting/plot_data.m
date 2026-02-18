%% plot_data.m - Plot Raw Data Series
% =========================================================================
% Generates diagnostic plots of input data series for the Scrambling for
% Dollars model.
%
% INPUTS:
%   - LFX_data.mat (from LFX_getdata.m)
%   - exchange_rate_data.mat (from LFX_getdata.m)
%
% OUTPUTS:
%   - Figures saved to Overleaf folder (if printit=1)
%
% (c) Saki Bigio
% =========================================================================

%% Settings
printit = 0;    % Set to 1 to save figures to PDF
printver = 1;   % Version number for saved files

% Set output folder based on machine
[~, username] = system('whoami');
username = strtrim(username);
if strcmp(username, 'sakibigio')
    foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollars_Revision_Restud/FigsTabs/';
elseif strcmp(username, 'sakiclaudia')
    foldername = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/FigsTabs/';
else
    warning('Unknown user: %s. Figures will not be saved.', username);
    foldername = './';
    printit = 0;
end

%% Load data
load LFX_data;
load exchange_rate_data;

%% Plot formatting
formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
    'Fontsize', 28, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
label_x = @(x) xlabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
    'Fontsize', 28, 'Interpreter', 'latex');
label_y = @(x) ylabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
    'Fontsize', 28, 'Interpreter', 'latex');

%% Date and scale settings
datesperiod = 1:234;
dates = datenum(2001, 1:234, 1);
abs_scale = 12e4;

% Currency lists
curlist = {'au', 'ca', 'jp', 'nz', 'no', 'sw', 'ch', 'uk'};
conlist = {'AUD', 'CAD', 'JPY', 'NZD', 'NOK', 'SWK', 'CHF', 'GBP'};

%% ========================================================================
%  INDIVIDUAL FIGURES
%  ========================================================================

%% Figure 1: Dollar-Euro Exchange Rate
figure('Name', 'Dollar Euro FX', 'NumberTitle', 'off')
plot(dates(datesperiod), exp(-inv_e(datesperiod)), 'LineWidth', 2);
datetick('x', 'yyyy-mm', 'keeplimits');
xlabel('year', 'interpreter', 'latex');
ylabel('level', 'interpreter', 'latex');
formataxis(gca);
grid on; axis tight;
if printit == 1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'raw_data_FX']);
end

%% Figure 2: Liquidity Ratios
figure('Name', 'Liquidity Ratios', 'NumberTitle', 'off')
plot(dates(datesperiod), exp(mu_us(datesperiod)), 'LineWidth', 2); hold on;
plot(dates(datesperiod), exp(mu_eu(datesperiod)), 'LineWidth', 2); hold off;
h = legend('US', 'EU', 'color', 'none', 'box', 'off', 'location', 'northwest');
set(h, 'interpreter', 'latex');
datetick('x', 'yyyy-mm', 'keeplimits');
xlabel('year', 'interpreter', 'latex');
ylabel('level', 'interpreter', 'latex');
formataxis(gca);
grid on; axis tight;
if printit == 1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'raw_data_mu']);
end

%% Figure 3: Policy Rate Differentials
figure('Name', 'Rate Differentials', 'NumberTitle', 'off')
plot(dates(datesperiod), (exp(im_eu(datesperiod))-1)*abs_scale, 'LineWidth', 2); hold on;
plot(dates(datesperiod), (exp(im_us(datesperiod))-1)*abs_scale, 'LineWidth', 2); hold off;
h = legend('EU', 'US', 'color', 'none', 'box', 'off', 'location', 'northeast');
set(h, 'interpreter', 'latex');
datetick('x', 'yyyy-mm', 'keeplimits');
xlabel('year', 'interpreter', 'latex');
ylabel('annual bps', 'interpreter', 'latex');
formataxis(gca);
grid on; axis tight;
if printit == 1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'raw_data_im']);
end

%% Figure 4: Money Supplies
figure('Name', 'Money Supplies', 'NumberTitle', 'off')
plot(dates(datesperiod), exp(M_us(datesperiod)), 'LineWidth', 2); hold on;
plot(dates(datesperiod), exp(M_eu(datesperiod)), 'LineWidth', 2); hold off;
h = legend('US', 'EU', 'color', 'none', 'box', 'off', 'location', 'northwest');
set(h, 'interpreter', 'latex');
datetick('x', 'yyyy-mm', 'keeplimits');
xlabel('year', 'interpreter', 'latex');
ylabel('level', 'interpreter', 'latex');
formataxis(gca);
grid on; axis tight;
if printit == 1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'raw_data_M']);
end

%% Figure 5: Bond Premium and CIP
figure('Name', 'Bond and CIP premia (US)', 'NumberTitle', 'off')
plot(dates(datesperiod), Rb_Rm(datesperiod)*abs_scale, 'LineWidth', 2); hold on;
plot(dates(datesperiod), cip(datesperiod)*abs_scale, 'LineWidth', 2); hold off;
h = legend('$\mathcal{BP}$', '$\mathcal{CIP}$', 'color', 'none', 'box', 'off', 'location', 'northwest');
set(h, 'interpreter', 'latex');
datetick('x', 'yyyy-mm', 'keeplimits');
xlabel('year', 'interpreter', 'latex');
ylabel('annual bps', 'interpreter', 'latex');
formataxis(gca);
grid on; axis tight;
if printit == 1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'raw_data_spreads']);
end

%% Figure 6: Interbank Spreads
figure('Name', 'Interbank Spreads (US)', 'NumberTitle', 'off')
plot(dates(datesperiod), Chi_D_US(datesperiod)*abs_scale, 'LineWidth', 2);
datetick('x', 'yyyy-mm', 'keeplimits');
xlabel('year', 'interpreter', 'latex');
ylabel('annual bps', 'interpreter', 'latex');
formataxis(gca);
grid on; axis tight;
if printit == 1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'raw_data_CIP']);
end

%% ========================================================================
%  MULTI-CURRENCY FIGURES
%  ========================================================================

%% Figure 7: Exchange Rates (All Currencies)
figure('Name', 'Exchange Rates', 'NumberTitle', 'off')
subplot(3, 3, 1);
plot(dates(datesperiod), exp(-inv_e(datesperiod)), 'LineWidth', 2);
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
title('EUR/USD', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    ln_data = eval(['ln_' curlist{j} '_us_t']);
    plot(dates(datesperiod), exp(-ln_data(datesperiod)), 'LineWidth', 2);
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD'], 'interpreter', 'latex', 'Fontsize', 15);
end
if printit == 1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'fig_exchange_rate_path']);
end

%% Figure 8: TED Spreads (All Currencies)
figure('Name', 'TED Spreads', 'NumberTitle', 'off')
abp_scale = 12 * 1e4;

subplot(3, 3, 1);
plot(dates(datesperiod), TED_s_eu_t(datesperiod)*abp_scale, 'LineWidth', 2); hold on;
plot(dates(datesperiod), TED_s_us_t(datesperiod)*abp_scale, 'LineWidth', 2); hold off;
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h = legend('EUR', 'US');
set(h, 'interpreter', 'latex', 'location', 'Northeast', 'Fontsize', 10);
title('TED Spreads', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    TED_data = eval(['TED_s_' curlist{j} '_t']);
    plot(dates(datesperiod), TED_data(datesperiod)*abp_scale, 'LineWidth', 2); hold on;
    plot(dates(datesperiod), TED_s_us_t(datesperiod)*abp_scale, 'LineWidth', 2); hold off;
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} ' TED'], 'interpreter', 'latex', 'Fontsize', 15);
end
if printit == 1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'fig_TED_spreads']);
end

%% Figure 9: Non-Negativity Test
figure('Name', 'Non-Negativity Test', 'NumberTitle', 'off')
abp_scale = 12 * 1e4;

% EUR/USD
aux = (-(im_eu(datesperiod) - im_us(datesperiod)) + Rb_Rm(datesperiod)) * abp_scale;
rsp_aux = -(mean(aux) - mean(cip(datesperiod)) * abp_scale);
m_aux = mean(aux) + rsp_aux;

subplot(3, 3, 1);
plot(dates(datesperiod), aux, 'LineWidth', 2); hold on;
plot(dates(datesperiod), m_aux + 0*dates(datesperiod), 'LineWidth', 1, 'Color', 'r'); hold off;
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
title(['EUR/USD (' num2str(round(m_aux)) ')'], 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    im_data = eval(['im_' curlist{j}]);
    aux = (-(im_data(datesperiod) - im_us(datesperiod)) + Rb_Rm(datesperiod)) * abp_scale;
    m_aux = mean(aux) + rsp_aux;
    
    subplot(3, 3, j+1);
    plot(dates(datesperiod), aux, 'LineWidth', 2); hold on;
    plot(dates(datesperiod), m_aux + 0*dates(datesperiod), 'LineWidth', 1, 'Color', 'r'); hold off;
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD (' num2str(round(m_aux)) ')'], 'interpreter', 'latex', 'Fontsize', 15);
end
if printit == 1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'fig_non_negativity']);
end

%% ========================================================================
%  SUMMARY FIGURES
%  ========================================================================

%% Figure 10: 9-Panel Summary
figure('Name', 'Data Summary (9 panels)', 'NumberTitle', 'off')

subplot(3, 3, 1);
plot(dates(datesperiod), exp(mu_us(datesperiod)), 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$\mu^*$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 3, 2);
plot(dates(datesperiod), exp(-inv_e(datesperiod)), 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$EUR/USD$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 3, 3);
plot(dates(datesperiod), Rb_Rm(datesperiod)*abs_scale, 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$R_b^*-R_m^*$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 3, 4);
plot(dates(datesperiod), exp(mu_eu(datesperiod)), 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$\mu$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 3, 5);
plot(dates(datesperiod), exp(im_eu(datesperiod)), 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$i_m$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 3, 6);
plot(dates(datesperiod), exp(im_us(datesperiod)), 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$i_m^*$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 3, 7);
plot(dates(datesperiod), cip(datesperiod)*abs_scale, 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('CIP', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 3, 8);
plot(dates(datesperiod), ois(datesperiod)*abs_scale, 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('OIS', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 3, 9);
plot(dates(datesperiod), Chi_D_US(datesperiod)*abs_scale, 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$\chi^-_{us}-\chi^+_{us}$', 'interpreter', 'latex', 'Fontsize', 15);

sgtitle('Data', 'interpreter', 'latex', 'Fontsize', 20);

if printit == 1
    orient landscape;
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'raw_data']);
end

%% Figure 11: 6-Panel Summary
figure('Name', 'Data Summary (6 panels)', 'NumberTitle', 'off')

subplot(3, 2, 1);
plot(dates(datesperiod), exp(-inv_e(datesperiod)), 'LineWidth', 2);
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$e$ (EUR/USD)', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 2, 2);
plot(dates(datesperiod), exp(mu_us(datesperiod)), 'LineWidth', 2); hold on;
plot(dates(datesperiod), exp(mu_eu(datesperiod)), 'LineWidth', 2); hold off;
h = legend('c=US', 'c=EU', 'color', 'none', 'box', 'off', 'location', 'northwest');
set(h, 'interpreter', 'latex');
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$\mu^{c}$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 2, 3);
plot(dates(datesperiod), exp(im_eu(datesperiod)), 'LineWidth', 2); hold on;
plot(dates(datesperiod), exp(im_us(datesperiod)), 'LineWidth', 2); hold off;
h = legend('c=EU', 'c=US', 'color', 'none', 'box', 'off', 'location', 'northeast');
set(h, 'interpreter', 'latex');
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$i_m^{c}$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 2, 4);
plot(dates(datesperiod), exp(M_us(datesperiod)), 'LineWidth', 2); hold on;
plot(dates(datesperiod), exp(M_eu(datesperiod)), 'LineWidth', 2); hold off;
h = legend('c=US', 'c=EU', 'color', 'none', 'box', 'off', 'location', 'northwest');
set(h, 'interpreter', 'latex');
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$M^{c}$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 2, 5);
plot(dates(datesperiod), Rb_Rm(datesperiod)*abs_scale, 'LineWidth', 2); hold on;
plot(dates(datesperiod), cip(datesperiod)*abs_scale, 'LineWidth', 2);
plot(dates(datesperiod), Rb_Rm_eu(datesperiod)*abs_scale, 'LineWidth', 2); hold off;
h = legend('$\mathcal{BP}$', '$\mathcal{CIP}$', '$\mathcal{BP}_{eu}$', ...
    'color', 'none', 'box', 'off', 'location', 'northwest');
set(h, 'interpreter', 'latex');
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('$\mathcal{BP}, \mathcal{CIP}$', 'interpreter', 'latex', 'Fontsize', 15);

subplot(3, 2, 6);
plot(dates(datesperiod), ois(datesperiod)*abs_scale, 'LineWidth', 2); hold on;
plot(dates(datesperiod), Chi_D_US(datesperiod)*abs_scale, 'LineWidth', 2); hold off;
h = legend('$\mathcal{OIS}$', '$R^f$ Spread', 'color', 'none', 'box', 'off', 'location', 'northwest');
set(h, 'interpreter', 'latex');
datetick('x', 'yy', 'keeplimits');
grid on; axis tight;
title('Other Spreads', 'interpreter', 'latex', 'Fontsize', 15);

sgtitle('Data', 'interpreter', 'latex', 'Fontsize', 20);

if printit == 1
    set(gcf, 'PaperUnits', 'normalized');
    set(gcf, 'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [foldername 'raw_data_6panel']);
end

fprintf('plot_data complete. %d figures generated.\n', 11);
