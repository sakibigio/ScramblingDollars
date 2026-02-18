%% LFX_getdata.m - Load and Process Raw Data
% =========================================================================
% This script reads data from LFX_datainputs.xlsx, calculates steady-state
% values, and estimates AR(1) parameters for interest rates and money bases.
%
% INPUTS:
%   - LFX_datainputs.xlsx (Excel file with raw data)
%
% OUTPUTS:
%   - LFX_data.mat              : Main data file for filtering
%   - LFX_targets.mat            : Steady-state calibration targets
%   - LFX_datamoments.mat        : Data moments (mean, std, autocorr)
%   - dynare_calibration_param.mat : Parameters for Dynare estimation
%   - exchange_rate_data.mat     : Exchange rate data
%
% (c) Saki Bigio
% =========================================================================
clear; 
close all;

%% ========================================================================
%  SETTINGS
%  ========================================================================

% Sample period settings
year_target = (25:72);          % Months 25-72 ≈ 2003-2007 (pre-crisis calibration)
dates = datenum(2001, 1:234, 1); % Monthly dates: Jan 2001 - Jun 2020

% Persistence options
exo_persistence_M = 1;          % 1 = use exogenous persistence for M (0.995)
exo_persistence_i = 0;          % 1 = use exogenous persistence for i (0.99)

% Plotting
plotit = 0;                     % 1 = generate diagnostic plots

% Frequency and scaling
freq = 12;                      % Monthly frequency
abps_factor = freq * 1e4;       % Convert monthly decimal to annual bps

% Liquidity ratio scaling factors
liq_ratio_scale_us = 0.18;      % US reserve ratio scaling
liq_ratio_scale_eu = 0.40;      % EU reserve ratio scaling

% Money supply detrending
detrendM = 1;                   % 1 = detrend money supply
demean_range = year_target;     % Range for computing trend

% Bond premium adjustment
Rb_Rm_scale = 0;                % Shift in Rb-Rm scale (set to 0 = no adjustment)

% Corridor width
iota_ss = 0.1;                  % DW-IOR spread (10% annual)

% Plot formatting
formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal');
label_x = @(x) xlabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
label_y = @(x) ylabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');

%% ========================================================================
%  CURRENCY CONFIGURATION
%  ========================================================================

% Currency list and corresponding Excel sheet names
currencies = struct();
currencies.codes = {'au', 'ca', 'eu', 'jp', 'nz', 'no', 'sw', 'ch', 'uk', 'us'};
currencies.sheets = {'Australia', 'Canada', 'Euro', 'Japan', 'NewZealand', ...
                     'Norway', 'Sweden', 'SwissFranc', 'UK', 'US'};
currencies.col_index = [3, 2, 10, 9, 4, 6, 8, 7, 5, 1];  % Column indices in policyrate matrix

% Rate adjustments (bps, annual) - applied to some countries
currencies.rate_adj = containers.Map();
currencies.rate_adj('au') = -100;  % Australia: subtract 100 bps
currencies.rate_adj('nz') = -200;  % New Zealand: subtract 200 bps
currencies.rate_adj('no') = -50;   % Norway: subtract 50 bps

n_currencies = length(currencies.codes);

%% ========================================================================
%  READ US AND EU LIQUIDITY RATIOS
%  ========================================================================

fprintf('Loading data from LFX_datainputs.xlsx...\n');

% US liquidity ratio (μ_us)
raw_mu_us = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'C27:C260');
mu_us = log(liq_ratio_scale_us * raw_mu_us);
fprintf('  Average US liquidity ratio: %.4f\n', mean(exp(mu_us)));

% EU liquidity ratio (μ_eu)
raw_mu_eu = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'D27:D260');
mu_eu = log(liq_ratio_scale_eu * raw_mu_eu);
fprintf('  Average EU liquidity ratio: %.4f\n', mean(exp(mu_eu)));

%% ========================================================================
%  READ SPREAD DATA
%  ========================================================================

% CIP deviation, OIS spread, Excess Bond Premium
cip = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'E27:E260') / freq / 1e4;
ois = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'F27:F260') / freq / 1e4;
ebp = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'G27:G260') / freq / 1e4;
ebp_eu = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'CG27:CG260') / freq / 1e4;
Chi_D_US = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'H27:H260') / freq / 1e4;

%% ========================================================================
%  READ POLICY RATES AND TED SPREADS (ALL CURRENCIES)
%  ========================================================================

% Policy rates (10 currencies)
policyrate = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'S27:AB260');
policyrate(isnan(policyrate)) = 0;

% LIBOR rates
Liborrate = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'AI27:AQ260');
Liborrate(isnan(Liborrate)) = 0;

% TED spreads (all currencies)
TED_all = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'CQ27:CZ260');

% Bond premia (all currencies)
BP_s_all = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'CG27:CP260');

% CIP deviations (all currencies)
CIP_s_all = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'DA27:DJ260');

% Extract US series
TED_s_us_t = TED_all(:, 1) / freq / 100;
i_us_t = policyrate(:, 1) / freq / 100;
BP_s_us_t = BP_s_all(:, 1) / freq / 100;
CIP_s_us_t = CIP_s_all(:, 1) / freq / 100;

%% ========================================================================
%  READ DISCOUNT WINDOW DATA
%  ========================================================================

DW_t = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'BV27:BV260');

%% ========================================================================
%  READ AND PROCESS MONEY SUPPLY
%  ========================================================================

% US Money Supply
M_us_t = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'AC27:AC260');
log_liq_t = log(M_us_t);
d_log_liq_t = diff(log_liq_t);

% Detrend if requested
if detrendM == 1
    gM_us = mean(d_log_liq_t(demean_range));
    d_log_liq_t = d_log_liq_t - gM_us;
end

% Reconstruct log money supply
M_us = log_liq_t(1) + cumsum([0; d_log_liq_t]);
M_us_ss = exp(mean(M_us));

% Estimate AR(1) for US money
X = [ones(length(M_us)-1, 1), M_us(1:end-1)];
[B, ~, R, ~, ~] = regress(M_us(2:end), X);
rho_M_us = B(2);
sigma_M_us = std(R);

% EU Money Supply
M_eu_t = xlsread('LFX_datainputs.xlsx', 'DataCounterpart', 'AD27:AD260');
log_liq_eu_t = log(M_eu_t);
d_log_liq_eu_t = diff(log_liq_eu_t);

% Detrend if requested
if detrendM == 1
    gM_eu = mean(d_log_liq_eu_t(demean_range));
    d_log_liq_eu_t = d_log_liq_eu_t - gM_eu;
end

% Reconstruct log money supply
M_eu = log_liq_eu_t(1) + cumsum([0; d_log_liq_eu_t]);
M_eu_ss = exp(mean(M_eu));

% Estimate AR(1) for EU money
X = [ones(length(M_eu)-1, 1), M_eu(1:end-1)];
[B, ~, R, ~, ~] = regress(M_eu(2:end), X);
rho_M_eu = B(2);
sigma_M_eu = std(R);

% Override with exogenous persistence if requested
if exo_persistence_M == 1
    rho_M_us = 0.995;
    rho_M_eu = 0.995;
end

fprintf('  Money supply AR(1): rho_M_us=%.3f, rho_M_eu=%.3f\n', rho_M_us, rho_M_eu);

%% ========================================================================
%  READ INFLATION DATA
%  ========================================================================

pi_us_t = xlsread('LFX_datainputs.xlsx', 'EuroData', 'AC28:AC261') / freq / 100;
libor_us_t = xlsread('LFX_datainputs.xlsx', 'EuroData', 'K28:K261') / freq / 100;

%% ========================================================================
%  PROCESS MULTI-CURRENCY DATA (using structs instead of eval)
%  ========================================================================

fprintf('Processing multi-currency data...\n');

% Initialize struct arrays
cur = struct();

% Process non-US currencies (9 currencies)
for j = 1:(length(currencies.codes) - 1)  % Exclude US
    code = currencies.codes{j};
    sheet = [currencies.sheets{j} 'Data'];
    col = currencies.col_index(j);
    
    % Policy rates and LIBOR (from pre-loaded matrices)
    cur.(code).i_t = policyrate(:, col) / freq / 100;
    cur.(code).libor_t = policyrate(:, col) / freq / 100;  % Same as policy rate in original
    cur.(code).TED_s_t = TED_all(:, col) / freq / 100;
    cur.(code).BP_s_t = BP_s_all(:, col) / freq / 100;
    cur.(code).CIP_s_t = CIP_s_all(:, col) / freq / 100;
    
    % Inflation and exchange rates (read from Excel)
    cur.(code).pi_t = xlsread('LFX_datainputs.xlsx', sheet, 'AD28:AD261') / freq / 100;
    cur.(code).ln_us_t = xlsread('LFX_datainputs.xlsx', sheet, 'C28:C261');
end

% Handle US separately (already loaded pi_us_t and libor_us_t above)
cur.us.i_t = policyrate(:, 1) / freq / 100;
cur.us.libor_t = policyrate(:, 1) / freq / 100;
cur.us.TED_s_t = TED_all(:, 1) / freq / 100;
cur.us.BP_s_t = BP_s_all(:, 1) / freq / 100;
cur.us.CIP_s_t = CIP_s_all(:, 1) / freq / 100;
cur.us.pi_t = pi_us_t;  % Already loaded from EuroData sheet

% Also store in individual variables for backward compatibility
for j = 1:n_currencies
    code = currencies.codes{j};
    assignin('caller', ['i_' code '_t'], cur.(code).i_t);
    assignin('caller', ['libor_' code '_t'], cur.(code).libor_t);
    assignin('caller', ['TED_s_' code '_t'], cur.(code).TED_s_t);
    assignin('caller', ['BP_s_' code '_t'], cur.(code).BP_s_t);
    assignin('caller', ['CIP_s_' code '_t'], cur.(code).CIP_s_t);
    assignin('caller', ['pi_' code '_t'], cur.(code).pi_t);
    if j < n_currencies
        assignin('caller', ['ln_' code '_us_t'], cur.(code).ln_us_t);
    end
end

%% ========================================================================
%  CALCULATE STEADY-STATE RATES AND AR(1) PARAMETERS
%  ========================================================================

% Initialize savelist for dynare
savelist_dynare = '';

for j = 1:n_currencies
    code = currencies.codes{j};
    
    % Get data
    i_t = cur.(code).i_t;
    libor_t = cur.(code).libor_t;
    pi_t = cur.(code).pi_t;
    
    % Steady-state rates (real)
    cur.(code).imss = mean(1 + i_t) / mean(1 + pi_t);
    cur.(code).iwss = mean(1 + i_t + iota_ss/freq) / mean(1 + pi_t);
    cur.(code).RLiborss = mean(1 + libor_t) / mean(1 + pi_t);
    
    % Convert to real rates for AR estimation
    i_t_real = (1 + i_t) / mean(1 + pi_t) - 1;
    libor_t_real = (1 + libor_t) / mean(1 + pi_t) - 1;
    
    % AR(1) estimation
    X = [ones(length(i_t_real)-1, 1), i_t_real(1:end-1)];
    [B, ~, R, ~, ~] = regress(i_t_real(2:end), X);
    cur.(code).rho_im = B(2);
    cur.(code).sigma_im = std(R);
    
    % Override persistence if requested
    if exo_persistence_i == 1
        cur.(code).rho_im = 0.99;
    end
    
    % Steady-state inflation
    cur.(code).piss = mean(pi_t) + 1;
    
    % Store transformed rates
    cur.(code).i_t = i_t_real + 1;
    cur.(code).libor_t = libor_t_real + 1;
    
    % Build savelist
    savelist_dynare = [savelist_dynare sprintf('imss_%s iwss_%s rho_im_%s sigma_im_%s RLiborss_%s ', ...
        code, code, code, code, code)];
end

% Apply rate adjustments for specific countries
adj_codes = keys(currencies.rate_adj);
for k = 1:length(adj_codes)
    code = adj_codes{k};
    adj = currencies.rate_adj(code) / 12 / 1e4;  % Convert annual bps to monthly decimal
    cur.(code).imss = cur.(code).imss + adj;
    cur.(code).iwss = cur.(code).iwss + adj;
    cur.(code).RLiborss = cur.(code).RLiborss + adj;
end

% Create individual variables for saving (backward compatibility)
for j = 1:n_currencies
    code = currencies.codes{j};
    assignin('caller', ['imss_' code], cur.(code).imss);
    assignin('caller', ['iwss_' code], cur.(code).iwss);
    assignin('caller', ['rho_im_' code], cur.(code).rho_im);
    assignin('caller', ['sigma_im_' code], cur.(code).sigma_im);
    assignin('caller', ['RLiborss_' code], cur.(code).RLiborss);
    assignin('caller', ['piss_' code], cur.(code).piss);
    assignin('caller', ['i_' code '_t'], cur.(code).i_t);
    assignin('caller', ['libor_' code '_t'], cur.(code).libor_t);
end

%% ========================================================================
%  PROCESS EXCHANGE RATES
%  ========================================================================

savelist_fx = '';

for j = 1:(n_currencies - 1)  % Exclude US
    code = currencies.codes{j};
    
    % Steady-state log exchange rate
    cur.(code).ln_us_ss = mean(cur.(code).ln_us_t);
    
    % Inverse exchange rate
    cur.(code).inv_e = -cur.(code).ln_us_t;
    
    % Create individual variables
    assignin('caller', ['ln_' code '_us_t'], cur.(code).ln_us_t);
    assignin('caller', ['ln_' code '_us_ss'], cur.(code).ln_us_ss);
    assignin('caller', ['inv_e_' code], cur.(code).inv_e);
    
    savelist_fx = [savelist_fx sprintf('ln_%s_us_t ln_%s_us_ss ', code, code)];
end

% Main exchange rate (EUR/USD)
inv_e = cur.eu.inv_e;
inv_e_jp = cur.jp.inv_e;
inv_e_ch = cur.ch.inv_e;

%% ========================================================================
%  PROCESS POLICY RATE TIME SERIES
%  ========================================================================

savelist_rates = '';

for j = 1:n_currencies
    code = currencies.codes{j};
    
    % Log rates
    cur.(code).im = log(cur.(code).i_t);
    cur.(code).RLibor = log(cur.(code).libor_t);
    cur.(code).im_obs = diff(cur.(code).im);
    
    % Create individual variables
    assignin('caller', ['im_' code], cur.(code).im);
    assignin('caller', ['RLibor_' code], cur.(code).RLibor);
    assignin('caller', ['im_' code '_obs'], cur.(code).im_obs);
    
    savelist_rates = [savelist_rates sprintf('im_%s_obs im_%s RLibor_%s TED_s_%s_t BP_s_%s_t CIP_s_%s_t ', ...
        code, code, code, code, code, code)];
end

% Apply log-rate adjustments
for k = 1:length(adj_codes)
    code = adj_codes{k};
    adj = currencies.rate_adj(code) / 12 / 1e4;
    cur.(code).im = cur.(code).im + adj;
    cur.(code).RLibor = cur.(code).RLibor + adj;
    assignin('caller', ['im_' code], cur.(code).im);
    assignin('caller', ['RLibor_' code], cur.(code).RLibor);
end

%% ========================================================================
%  EXCESS BOND PREMIUM
%  ========================================================================

Rb_Rm = ebp;
Rb_Rm_eu = ebp_eu;

% Apply scale adjustment if specified
if Rb_Rm_scale > 0
    Rb_Rm = Rb_Rm - mean(Rb_Rm) + Rb_Rm_scale;
    Rb_Rm_eu = Rb_Rm_eu - mean(Rb_Rm_eu) + Rb_Rm_scale;
end

Rb_us = log(Rb_Rm + exp(cur.us.im));

%% ========================================================================
%  COMPUTE PRE-CRISIS CALIBRATION TARGETS
%  ========================================================================

Ted_us_yt = mean(TED_s_us_t(year_target)) * 1e4 * 12;
Ted_eu_yt = mean(cur.eu.TED_s_t(year_target)) * 1e4 * 12;
ois_yt = mean(ois(year_target)) * 1e4 * 12;
cip_yt = mean(cip(year_target)) * 1e4 * 12;
uip_yt = mean(cur.eu.i_t(year_target) - cur.us.i_t(year_target)) * 1e4 * 12;
mu_us_yt = mean(exp(mu_us(year_target)));
mu_eu_yt = mean(exp(mu_eu(year_target)));
Rb_Rm_yt = mean(Rb_Rm(year_target) * freq * 1e4);
inv_e_yt = mean(exp(inv_e(year_target)));

fprintf('\n=== Pre-Crisis Calibration Targets (2003-2007) ===\n');
fprintf('  TED_us: %.2f bps\n', Ted_us_yt);
fprintf('  TED_eu: %.2f bps\n', Ted_eu_yt);
fprintf('  mu_us:  %.4f\n', mu_us_yt);
fprintf('  mu_eu:  %.4f\n', mu_eu_yt);
fprintf('  Rb_Rm:  %.2f bps\n', Rb_Rm_yt);

%% ========================================================================
%  COMPUTE DATA MOMENTS
%  ========================================================================

bps_scale = 12e4;
moments = struct();

% Time series in bps
ois_ts = ois * bps_scale;
cip_ts = cip * bps_scale;
idiff_ts = (cur.eu.i_t - cur.us.i_t) * bps_scale;
Rb_Rm_ts = Rb_Rm * bps_scale;
mu_us_ts = exp(mu_us);
mu_eu_ts = exp(mu_eu);
e_euus_ts = exp(inv_e);

var_names = {'ois', 'cip', 'idiff', 'Rb_Rm', 'mu_us', 'mu_eu', 'e_euus'};
var_data = {ois_ts, cip_ts, idiff_ts, Rb_Rm_ts, mu_us_ts, mu_eu_ts, e_euus_ts};

for j = 1:length(var_names)
    aux = var_data{j};
    rho = autocorr(aux);
    moments.(var_names{j}).mean = mean(aux(year_target));
    moments.(var_names{j}).std = std(aux);
    moments.(var_names{j}).rho = rho(2);
end

%% ========================================================================
%  DIAGNOSTIC PLOTS
%  ========================================================================

if plotit == 1
    figure('Name', 'Money Supply');
    plot(dates, M_us_t); hold on; 
    plot(dates, M_eu_t); 
    legend('US', 'EU');
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Money Supply');
    
    figure('Name', 'Spreads');
    plot(dates, cip * abps_factor); hold on;
    plot(dates, ois * abps_factor);
    plot(dates, ebp * abps_factor);
    plot(dates, (cur.eu.im - cur.us.im) * abps_factor);
    legend('Bond CIP', 'OIS CIP', 'EBP', 'Interest Differential');
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Spreads (bps)');
    
    figure('Name', 'Bond Premium');
    plot(dates, Rb_Rm * abps_factor); hold on;
    plot(dates, Rb_Rm_eu * abps_factor);
    legend('US', 'EU');
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Excess Bond Premium (bps)');
    
    figure('Name', 'Liquidity Ratios');
    plot(dates, exp(mu_us)); hold on;
    plot(dates, exp(mu_eu));
    legend('US', 'EU');
    datetick('x', 'yy');
    grid on; axis tight;
    label_x('Time (Year)');
    title('Liquidity Ratios');
end

%% ========================================================================
%  SAVE OUTPUT FILES
%  ========================================================================

fprintf('\nSaving output files...\n');

% Rename for compatibility
Ted_us = TED_s_us_t;
Ted_eu = cur.eu.TED_s_t;
im_us = cur.us.im;
im_eu = cur.eu.im;
TED_s_eu_t = cur.eu.TED_s_t;
pi_eu_t = cur.eu.pi_t;

% 1. Main data file
save('LFX_data.mat', 'mu_eu', 'mu_us', 'inv_e', 'Ted_us', 'Ted_eu', ...
    'inv_e_jp', 'inv_e_ch', 'M_us', 'Rb_Rm', 'Rb_Rm_eu', 'Rb_us', 'M_eu', ...
    'ois', 'cip', 'Chi_D_US', 'pi_us_t', 'pi_eu_t', 'DW_t', ...
    'im_us', 'im_eu', 'TED_s_us_t', 'TED_s_eu_t', ...
    'im_au', 'im_ca', 'im_jp', 'im_nz', 'im_no', 'im_sw', 'im_ch', 'im_uk', ...
    'TED_s_au_t', 'TED_s_ca_t', 'TED_s_jp_t', 'TED_s_nz_t', ...
    'TED_s_no_t', 'TED_s_sw_t', 'TED_s_ch_t', 'TED_s_uk_t');
fprintf('  Saved: LFX_data.mat\n');

% 2. Calibration targets
save('LFX_targets.mat', 'Ted_us_yt', 'Ted_eu_yt', 'mu_us_yt', 'mu_eu_yt', ...
    'Rb_Rm_yt', 'inv_e_yt', 'ois_yt', 'cip_yt', 'uip_yt');
fprintf('  Saved: LFX_targets.mat\n');

% 3. Data moments
save('LFX_datamoments.mat', 'moments');
fprintf('  Saved: LFX_datamoments.mat\n');

% 4. Calibration parameters (replaces old calibration.mat and dynare_calibration_param.mat)
save('calibration.mat', 'M_us_ss', 'rho_M_us', 'sigma_M_us', ...
    'M_eu_ss', 'rho_M_eu', 'sigma_M_eu', ...
    'imss_us', 'iwss_us', 'rho_im_us', 'sigma_im_us', 'RLiborss_us', ...
    'imss_eu', 'iwss_eu', 'rho_im_eu', 'sigma_im_eu', 'RLiborss_eu', ...
    'imss_au', 'iwss_au', 'rho_im_au', 'sigma_im_au', 'RLiborss_au', ...
    'imss_ca', 'iwss_ca', 'rho_im_ca', 'sigma_im_ca', 'RLiborss_ca', ...
    'imss_jp', 'iwss_jp', 'rho_im_jp', 'sigma_im_jp', 'RLiborss_jp', ...
    'imss_nz', 'iwss_nz', 'rho_im_nz', 'sigma_im_nz', 'RLiborss_nz', ...
    'imss_no', 'iwss_no', 'rho_im_no', 'sigma_im_no', 'RLiborss_no', ...
    'imss_sw', 'iwss_sw', 'rho_im_sw', 'sigma_im_sw', 'RLiborss_sw', ...
    'imss_ch', 'iwss_ch', 'rho_im_ch', 'sigma_im_ch', 'RLiborss_ch', ...
    'imss_uk', 'iwss_uk', 'rho_im_uk', 'sigma_im_uk', 'RLiborss_uk');
fprintf('  Saved: calibration.mat\n');

% 5. Exchange rate data
save('exchange_rate_data.mat', 'inv_e', 'inv_e_jp', 'inv_e_ch', ...
    'ln_au_us_t', 'ln_ca_us_t', 'ln_eu_us_t', 'ln_jp_us_t', ...
    'ln_nz_us_t', 'ln_no_us_t', 'ln_sw_us_t', 'ln_ch_us_t', 'ln_uk_us_t');
fprintf('  Saved: exchange_rate_data.mat\n');

fprintf('\nData loading complete!\n');
