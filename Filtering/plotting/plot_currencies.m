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
if exist('Rb_Rm_eu', 'var')
    plot(dates(datesperiod), Rb_Rm_eu(datesperiod)*abs_scale, 'LineWidth', 2);
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
    % Try to get data variable, use model if not available
    try
        BP_data = eval(['Rb_Rm_' curlist{j} '(datesperiod)']);
    catch
        BP_data = BP_model;  % Fallback to model
    end
    plot(dates(datesperiod), BP_model*abs_scale, 'LineWidth', 2); hold on;
    plot(dates(datesperiod), BP_data*abs_scale, 'LineWidth', 2); hold off;
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
if exist('cip', 'var')
    plot(dates(datesperiod), cip(datesperiod)*abs_scale, 'LineWidth', 2);
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
    % Try to get data variable, use model if not available
    try
        CIP_data = eval(['cip_' curlist{j} '(datesperiod)']);
    catch
        CIP_data = CIP_model;  % Fallback to model
    end
    plot(dates(datesperiod), CIP_model*abs_scale, 'LineWidth', 2); hold on;
    plot(dates(datesperiod), CIP_data*abs_scale, 'LineWidth', 2); hold off;
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title(conlist{j}, 'interpreter', 'latex', 'Fontsize', 15);
end

%% ========================================================================
%  EXCHANGE RATES
%  ========================================================================

%% Figure 5: FX Rates (Model vs Data, Levels)
figure('Name', 'FX Rates (Levels)', 'NumberTitle', 'off')

% EU/US exchange rate
subplot(3, 3, 1);
plot(dates(datesperiod), log(e_euus_t(datesperiod)), 'LineWidth', 2); hold on;
plot(dates(datesperiod), -inv_e(datesperiod), 'LineWidth', 2, 'LineStyle', '--'); hold off;
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
h = legend('Model', 'Data');
set(h, 'interpreter', 'latex', 'location', 'Northeast', 'Fontsize', 10);
title('EUR/USD', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    try
        inv_e_model = eval(['inv_e_' curlist{j} '_t(datesperiod)']);
        inv_e_data = eval(['inv_e_' curlist{j} '(datesperiod)']);
        plot(dates(datesperiod), -log(inv_e_model), 'LineWidth', 2); hold on;
        plot(dates(datesperiod), -inv_e_data, 'LineWidth', 2, 'LineStyle', '--'); hold off;
    catch
        % Skip if variable doesn't exist
    end
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD'], 'interpreter', 'latex', 'Fontsize', 15);
end

%% Figure 6: FX Rates (Demeaned)
figure('Name', 'FX Rates (Demeaned)', 'NumberTitle', 'off')

% Check if FX data variables exist
if exist('inv_e', 'var')
    temp2 = mean(log(inv_e_t(datesperiod)));
    subplot(3, 3, 1);
    plot(dates(datesperiod), -log(inv_e_t(datesperiod)) + temp2, 'LineWidth', 2); hold on;
    plot(dates(datesperiod), -inv_e(datesperiod) + mean(inv_e(datesperiod)), 'LineWidth', 2, 'LineStyle', '--'); hold off;
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    h = legend('Model', 'Data');
    set(h, 'interpreter', 'latex', 'location', 'Northeast', 'Fontsize', 10);
    title('EUR/USD', 'interpreter', 'latex', 'Fontsize', 15);
end

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    try
        inv_e_model = eval(['inv_e_' curlist{j} '_t(datesperiod)']);
        temp2 = mean(log(inv_e_model));
        plot(dates(datesperiod), -log(inv_e_model) + temp2, 'LineWidth', 2);
    catch
        % Skip if variable doesn't exist
    end
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD'], 'interpreter', 'latex', 'Fontsize', 15);
end

%% Figure 7: FX Devaluation Rates (Model Only)
figure('Name', 'FX Devaluation Rates', 'NumberTitle', 'off')

datesperiod_f = datesperiod(2:end);
datesperiod_l = datesperiod(1:end-1);

% EUR/USD
temp1 = (-log(inv_e_t(datesperiod_f)) + log(inv_e_t(datesperiod_l))) * abs_scale/100;
subplot(3, 3, 1);
plot(dates(datesperiod_f), temp1, 'LineWidth', 2);
datetick('x', 'yy');
grid on; axis tight;
label_x('Time (Year)');
formataxis(gca);
title('EUR/USD', 'interpreter', 'latex', 'Fontsize', 15);

for j = 1:length(curlist)
    subplot(3, 3, j+1);
    try
        inv_e_model = eval(['inv_e_' curlist{j} '_t']);
        temp1 = (-log(inv_e_model(datesperiod_f)) + log(inv_e_model(datesperiod_l))) * abs_scale/100;
        plot(dates(datesperiod_f), temp1, 'LineWidth', 2);
    catch
        % Skip if variable doesn't exist
    end
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    formataxis(gca);
    title([conlist{j} '/USD'], 'interpreter', 'latex', 'Fontsize', 15);
end

%% Summary
fprintf('plot_currencies complete. 7 figures generated.\n');
