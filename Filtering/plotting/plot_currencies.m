%% plot_currencies.m - Plot Multi-Currency Filter Results
% =========================================================================
% Plots filter results for all currencies (AU, CA, JP, NZ, NO, SW, CH, UK).
% Called after main_filter.m completes.
%
% Required variables from workspace:
%   - sigma_*_t, sigma_*_TED_t, sigma_*_TED_flag for each currency
%   - BP_*_t, CIP_*_t, inv_e_*_t, ln_*_us_t for each currency
%   - dates, datesperiod, abs_scale, curlist, conlist
%
% (c) Saki Bigio
% =========================================================================

%% Settings
printit = 0;  % Set to 1 to save figures

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

% Define formatting functions if not already defined
if ~exist('formataxis', 'var')
    formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', 14, 'Box', 'On');
end

if ~exist('label_x', 'var')
    label_x = @(x) xlabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', 14, 'Interpreter', 'latex');
end

%% ========================================================================
%  DIAGNOSTICS
%  ========================================================================

%% Figure 1: Solver Flags (All Currencies)
figure('Name', 'Solver Diagnostics', 'NumberTitle', 'off') 

subplot(3, 3, 1); 
plot(dates(datesperiod), sigma_eu_TED_flag, 'LineWidth', 2);
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h = legend('EU', 'US');
set(h, 'interpreter', 'latex', 'location', 'Northeast', 'Fontsize', 10, 'box', 'off');
title('EU', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    flag_data = eval(['sigma_' curlist{j} '_TED_flag(datesperiod)']);
    plot(dates(datesperiod), flag_data, 'LineWidth', 2);
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title(conlist{j}, 'interpreter', 'latex', 'Fontsize', 15);
end

%% ========================================================================
%  SIGMA ESTIMATES
%  ========================================================================

%% Figure 2: Log Sigma (All Currencies)
figure('Name', 'Currency Sigmas', 'NumberTitle', 'off') 

subplot(3, 3, 1); 
plot(dates(datesperiod), log(sigma_eu_TED_t(datesperiod)), 'LineWidth', 2); hold on;
plot(dates(datesperiod), log(sigma_us_TED_t(datesperiod)), 'LineWidth', 2); hold off;
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h = legend('EU', 'US');
set(h, 'interpreter', 'latex', 'location', 'Northeast', 'Fontsize', 10, 'box', 'off');
title('EU', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    sigma_data = eval(['sigma_' curlist{j} '_t(datesperiod)']);
    plot(dates(datesperiod), log(sigma_data), 'LineWidth', 2); hold on;
    plot(dates(datesperiod), log(sigma_us_TED_t(datesperiod)), 'LineWidth', 1, 'LineStyle', '-'); hold off;
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title(conlist{j}, 'interpreter', 'latex', 'Fontsize', 15);
end

%% ========================================================================
%  BOND PREMIA
%  ========================================================================

%% Figure 3: Bond Premia (Model vs Data)
figure('Name', 'Global Bond Premia', 'NumberTitle', 'off') 

subplot(3, 3, 1); 
plot(dates(datesperiod), BP_eu_t(datesperiod)*abs_scale, 'LineWidth', 2); hold on;
% EU data is saved as BP_s_eu_t or Rb_Rm_eu
try
    plot(dates(datesperiod), BP_s_eu_t(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', '--');
catch
    if exist('Rb_Rm_eu', 'var')
        plot(dates(datesperiod), Rb_Rm_eu(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', '--');
    end
end
hold off;
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h = legend('Model', 'Data');
set(h, 'interpreter', 'latex', 'location', 'Northeast', 'Fontsize', 10);
title('EU', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    BP_model = eval(['BP_' curlist{j} '_t(datesperiod)']);
    plot(dates(datesperiod), BP_model*abs_scale, 'LineWidth', 2); hold on;
    % Try to get data variable (saved as BP_s_*_t)
    try
        BP_data = eval(['BP_s_' curlist{j} '_t(datesperiod)']);
        plot(dates(datesperiod), BP_data*abs_scale, 'LineWidth', 2, 'LineStyle', '--');
        h = legend('Model', 'Data');
        set(h, 'interpreter', 'latex', 'location', 'best', 'Fontsize', 8);
    catch
        % No data available - model only
    end
    hold off;
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title(conlist{j}, 'interpreter', 'latex', 'Fontsize', 15);
end

%% ========================================================================
%  CIP DEVIATIONS
%  ========================================================================

%% Figure 4: CIP Deviations (Model vs Data)
figure('Name', 'CIP Deviations', 'NumberTitle', 'off') 

subplot(3, 3, 1); 
plot(dates(datesperiod), CIP_t(datesperiod)*abs_scale, 'LineWidth', 2); hold on;
% EU CIP data might be saved as 'cip' or 'CIP_s_eu_t'
try
    plot(dates(datesperiod), CIP_s_eu_t(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', '--');
catch
    if exist('cip', 'var')
        plot(dates(datesperiod), cip(datesperiod)*abs_scale, 'LineWidth', 2, 'LineStyle', '--');
    end
end
hold off;
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h = legend('Model', 'Data');
set(h, 'interpreter', 'latex', 'location', 'Northeast', 'Fontsize', 10);
title('EU', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    CIP_model = eval(['CIP_' curlist{j} '_t(datesperiod)']);
    plot(dates(datesperiod), CIP_model*abs_scale, 'LineWidth', 2); hold on;
    % Try to get data variable (saved as CIP_s_*_t)
    try
        CIP_data = eval(['CIP_s_' curlist{j} '_t(datesperiod)']);
        plot(dates(datesperiod), CIP_data*abs_scale, 'LineWidth', 2, 'LineStyle', '--');
        h = legend('Model', 'Data');
        set(h, 'interpreter', 'latex', 'location', 'best', 'Fontsize', 8);
    catch
        % No data available - model only
    end
    hold off;
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title(conlist{j}, 'interpreter', 'latex', 'Fontsize', 15);
end

%% ========================================================================
%  EXCHANGE RATES
%  ========================================================================

%% Figure 5: FX Rates (Model vs Data - Demeaned Dynamics)
figure('Name', 'FX Rates (Demeaned Dynamics)', 'NumberTitle', 'off')

% EU/US exchange rate
subplot(3, 3, 1);
model_eu = -log(inv_e_t(datesperiod));  % Model: -log(1/e) = log(e)
data_eu = -inv_e(datesperiod);          % Data is stored as -log(e), so negate
plot(dates(datesperiod), model_eu - mean(model_eu), 'LineWidth', 2); hold on;
plot(dates(datesperiod), data_eu - mean(data_eu), 'LineWidth', 2, 'LineStyle', '--'); hold off;
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h = legend('Model', 'Data');
set(h, 'interpreter', 'latex', 'location', 'Northeast', 'Fontsize', 10);
title('EUR/USD', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    % Get model
    inv_e_model = eval(['inv_e_' curlist{j} '_t(datesperiod)']);
    
    % Try to get data - use try/catch since exist() doesn't work with string var names
    try
        inv_e_data = eval(['inv_e_' curlist{j} '(datesperiod)']);
        % Demean both for comparison
        plot(dates(datesperiod), inv_e_model - mean(inv_e_model), 'LineWidth', 2); hold on;
        plot(dates(datesperiod), inv_e_data - mean(inv_e_data), 'LineWidth', 2, 'LineStyle', '--'); hold off;
        h = legend('Model', 'Data');
        set(h, 'interpreter', 'latex', 'location', 'best', 'Fontsize', 8);
    catch
        % Only plot model if data doesn't exist
        plot(dates(datesperiod), inv_e_model - mean(inv_e_model), 'LineWidth', 2);
    end
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD'], 'interpreter', 'latex', 'Fontsize', 15);
end

%% Figure 6: FX Changes (Model vs Data with Correlations)
figure('Name', 'FX Changes', 'NumberTitle', 'off')

datesperiod_f = datesperiod(2:end);
datesperiod_l = datesperiod(1:end-1);

% EUR/USD
model_chg = diff(-log(inv_e_t(datesperiod))) * abs_scale/100;
data_chg = diff(-inv_e(datesperiod)) * abs_scale/100;
subplot(3, 3, 1);
plot(dates(datesperiod_f), model_chg, 'LineWidth', 2); hold on;
plot(dates(datesperiod_f), data_chg, 'LineWidth', 1, 'LineStyle', '--'); hold off;
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
rho = corr(model_chg, data_chg);
title(['EUR/USD (corr=' num2str(rho, '%.2f') ')'], 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    inv_e_model = eval(['inv_e_' curlist{j} '_t']);
    model_chg = diff(inv_e_model(datesperiod)) * abs_scale/100;
    
    % Try to get data
    try
        inv_e_data = eval(['inv_e_' curlist{j}]);
        data_chg = diff(inv_e_data(datesperiod)) * abs_scale/100;
        plot(dates(datesperiod_f), model_chg, 'LineWidth', 2); hold on;
        plot(dates(datesperiod_f), data_chg, 'LineWidth', 1, 'LineStyle', '--'); hold off;
        rho = corr(model_chg, data_chg);
        title([conlist{j} '/USD (corr=' num2str(rho, '%.2f') ')'], 'interpreter', 'latex', 'Fontsize', 15);
    catch
        plot(dates(datesperiod_f), model_chg, 'LineWidth', 2);
        title([conlist{j} '/USD (model only)'], 'interpreter', 'latex', 'Fontsize', 15);
    end
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
end

%% Figure 7: Model FX Levels (All Currencies)
figure('Name', 'Model FX Levels', 'NumberTitle', 'off')

subplot(3, 3, 1);
plot(dates(datesperiod), -log(inv_e_t(datesperiod)), 'LineWidth', 2);
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
title('EUR/USD (model)', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    try
        inv_e_model = eval(['inv_e_' curlist{j} '_t(datesperiod)']);
        plot(dates(datesperiod), inv_e_model, 'LineWidth', 2);
    catch
    end
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD (model)'], 'interpreter', 'latex', 'Fontsize', 15);
end

%% Summary
fprintf('plot_currencies complete. 7 figures generated.\n');
