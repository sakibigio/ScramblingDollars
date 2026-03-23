%% Compare Leontief vs Cobb-Douglas Filtering
% Runs main_filter.m with both matching functions and compares results
%
% (c) Saki Bigio
% 
% Targets: US BP, EU BP, CIP, FF volume, DW volume
% Metrics: corr(level), corr(diff), std, max, min
% Columns: Data, Leontief, Cobb-Douglas

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

% Save data counterparts (same for both runs)
Rb_Rm_data = Rb_Rm;
Rb_Rm_eu_data = Rb_Rm_eu;
cip_data = cip;
DW_n_data = DW_n;          % WLCFLPCL / TCDSL (primary credit / checkable deposits)
FF_n_data = FF_n;           % FF volume / TCDSL

save('data/filter_comparison_leontief.mat', ...
    'sigma_us_L', 'sigma_eu_L', 'BP_us_L', 'BP_eu_L', ...
    'CIP_L', 'TED_us_L', 'TED_eu_L', 'theta_us_L', 'theta_eu_L', ...
    'psi_us_L', 'psi_eu_L', 'riskprm_L', 'DW_us_L', 'FF_us_L', ...
    'Echi_m_us_L', 'Echi_m_eu_L', 'sigma_us_TED_flag_L', 'sigma_eu_TED_flag_L', ...
    'Rb_Rm_data', 'Rb_Rm_eu_data', 'cip_data', 'DW_n_data', 'FF_n_data', ...
    'dates_L', 'datesperiod_L');
fprintf('Leontief results saved to data/filter_comparison_leontief.mat\n');

close all;

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

save('data/filter_comparison_cd.mat', ...
    'sigma_us_CD', 'sigma_eu_CD', 'BP_us_CD', 'BP_eu_CD', ...
    'CIP_CD', 'TED_us_CD', 'TED_eu_CD', 'theta_us_CD', 'theta_eu_CD', ...
    'psi_us_CD', 'psi_eu_CD', 'riskprm_CD', 'DW_us_CD', 'FF_us_CD', ...
    'Echi_m_us_CD', 'Echi_m_eu_CD', 'sigma_us_TED_flag_CD', 'sigma_eu_TED_flag_CD', ...
    'dates_CD', 'datesperiod_CD');
fprintf('Cobb-Douglas results saved to data/filter_comparison_cd.mat\n');

close all;

%% ========================================================================
%  LOAD RESULTS FOR COMPARISON
%  ========================================================================
load('data/filter_comparison_leontief.mat');
load('data/filter_comparison_cd.mat');
dates = dates_L;
datesperiod = datesperiod_L;

%% ========================================================================
%  COMPREHENSIVE COMPARISON TABLE
%  ========================================================================
abs_scale = 12e4;  % BPS annualized

fprintf('\n');
fprintf('================================================================================\n');
fprintf('  COMPREHENSIVE COMPARISON: LEONTIEF vs COBB-DOUGLAS vs DATA\n');
fprintf('================================================================================\n');

% --- Helper: compute stats for a model-data pair ---
% Returns [corr_level, corr_diff, std, max, min]
compute_stats = @(model, data, dp) [ ...
    corr(model(dp), data(dp)), ...
    corr(diff(model(dp)), diff(data(dp))), ...
    std(model(dp)), ...
    max(model(dp)), ...
    min(model(dp))];

% --- Scale to BPS annualized ---
BP_us_L_bps  = BP_us_L  * abs_scale;
BP_us_CD_bps = BP_us_CD * abs_scale;
BP_us_D_bps  = Rb_Rm_data * abs_scale;

BP_eu_L_bps  = BP_eu_L  * abs_scale;
BP_eu_CD_bps = BP_eu_CD * abs_scale;
BP_eu_D_bps  = Rb_Rm_eu_data * abs_scale;

CIP_L_bps  = CIP_L  * abs_scale;
CIP_CD_bps = CIP_CD * abs_scale;
CIP_D_bps  = cip_data * abs_scale;

% --- DW: already deposit-normalized (WLCFLPCL / TCDSL) ---
DW_L_pct  = DW_us_L * 100;   % Model: percent of deposits
DW_CD_pct = DW_us_CD * 100;
DW_D_pct  = DW_n_data * 100; % Data: percent of checkable deposits

% --- FF: already deposit-normalized (FF volume / TCDSL) ---
FF_L_pct  = FF_us_L * 100;
FF_CD_pct = FF_us_CD * 100;
FF_D_pct  = FF_n_data * 100;

% --- Valid data periods (handle NaN gaps) ---
% DW (WLCFLPCL): starts ~2003 (first ~24 months NaN)
% FF (NY Fed volume): starts ~2006 (first ~60 months NaN)
dw_valid = datesperiod(~isnan(DW_n_data(datesperiod)));
ff_valid = datesperiod(~isnan(FF_n_data(datesperiod)));
fprintf('DW data available: %d/%d periods (starts period %d)\n', length(dw_valid), length(datesperiod), dw_valid(1));
fprintf('FF data available: %d/%d periods (starts period %d)\n', length(ff_valid), length(datesperiod), ff_valid(1));

% =====================================================================
%  TABLE 1: US Bond Premium
% =====================================================================
fprintf('\n--- US Bond Premium (bps annualized) ---\n');
fprintf('                    Data        Leontief    Cobb-Douglas\n');
fprintf('                    ----        --------    ------------\n');

L_stats  = compute_stats(BP_us_L_bps, BP_us_D_bps, datesperiod);
CD_stats = compute_stats(BP_us_CD_bps, BP_us_D_bps, datesperiod);

fprintf('corr(level)          —          %8.3f    %12.3f\n', L_stats(1), CD_stats(1));
fprintf('corr(diff)           —          %8.3f    %12.3f\n', L_stats(2), CD_stats(2));
fprintf('Std                %6.1f      %8.1f    %12.1f\n', std(BP_us_D_bps(datesperiod)), L_stats(3), CD_stats(3));
fprintf('Max                %6.1f      %8.1f    %12.1f\n', max(BP_us_D_bps(datesperiod)), L_stats(4), CD_stats(4));
fprintf('Min                %6.1f      %8.1f    %12.1f\n', min(BP_us_D_bps(datesperiod)), L_stats(5), CD_stats(5));

% =====================================================================
%  TABLE 2: EU Bond Premium
% =====================================================================
fprintf('\n--- EU Bond Premium (bps annualized) ---\n');
fprintf('                    Data        Leontief    Cobb-Douglas\n');
fprintf('                    ----        --------    ------------\n');

L_stats  = compute_stats(BP_eu_L_bps, BP_eu_D_bps, datesperiod);
CD_stats = compute_stats(BP_eu_CD_bps, BP_eu_D_bps, datesperiod);

fprintf('corr(level)          —          %8.3f    %12.3f\n', L_stats(1), CD_stats(1));
fprintf('corr(diff)           —          %8.3f    %12.3f\n', L_stats(2), CD_stats(2));
fprintf('Std                %6.1f      %8.1f    %12.1f\n', std(BP_eu_D_bps(datesperiod)), L_stats(3), CD_stats(3));
fprintf('Max                %6.1f      %8.1f    %12.1f\n', max(BP_eu_D_bps(datesperiod)), L_stats(4), CD_stats(4));
fprintf('Min                %6.1f      %8.1f    %12.1f\n', min(BP_eu_D_bps(datesperiod)), L_stats(5), CD_stats(5));

% =====================================================================
%  TABLE 3: CIP Deviation
% =====================================================================
fprintf('\n--- CIP Deviation (bps annualized) ---\n');
fprintf('                    Data        Leontief    Cobb-Douglas\n');
fprintf('                    ----        --------    ------------\n');

L_stats  = compute_stats(CIP_L_bps, CIP_D_bps, datesperiod);
CD_stats = compute_stats(CIP_CD_bps, CIP_D_bps, datesperiod);

fprintf('corr(level)          —          %8.3f    %12.3f\n', L_stats(1), CD_stats(1));
fprintf('corr(diff)           —          %8.3f    %12.3f\n', L_stats(2), CD_stats(2));
fprintf('Std                %6.1f      %8.1f    %12.1f\n', std(CIP_D_bps(datesperiod)), L_stats(3), CD_stats(3));
fprintf('Max                %6.1f      %8.1f    %12.1f\n', max(CIP_D_bps(datesperiod)), L_stats(4), CD_stats(4));
fprintf('Min                %6.1f      %8.1f    %12.1f\n', min(CIP_D_bps(datesperiod)), L_stats(5), CD_stats(5));

% =====================================================================
%  TABLE 4: DW Volume (WLCFLPCL / TCDSL)
% =====================================================================
fprintf('\n--- Discount Window Volume (%% of checkable deposits, %d periods) ---\n', length(dw_valid));
fprintf('                    Data        Leontief    Cobb-Douglas\n');
fprintf('                    ----        --------    ------------\n');

% Correlations with data — restricted to valid DW window
L_stats_dw  = compute_stats(DW_L_pct, DW_D_pct, dw_valid);
CD_stats_dw = compute_stats(DW_CD_pct, DW_D_pct, dw_valid);

fprintf('corr(level)          —          %8.3f    %12.3f\n', L_stats_dw(1), CD_stats_dw(1));
fprintf('corr(diff)           —          %8.3f    %12.3f\n', L_stats_dw(2), CD_stats_dw(2));
fprintf('Std                %6.3f      %8.3f    %12.3f\n', std(DW_D_pct(dw_valid)), std(DW_L_pct(dw_valid)), std(DW_CD_pct(dw_valid)));
fprintf('Max                %6.3f      %8.3f    %12.3f\n', max(DW_D_pct(dw_valid)), max(DW_L_pct(dw_valid)), max(DW_CD_pct(dw_valid)));
fprintf('Min                %6.3f      %8.3f    %12.3f\n', min(DW_D_pct(dw_valid)), min(DW_L_pct(dw_valid)), min(DW_CD_pct(dw_valid)));
fprintf('Mean               %6.3f      %8.3f    %12.3f\n', nanmean(DW_D_pct(dw_valid)), mean(DW_L_pct(dw_valid)), mean(DW_CD_pct(dw_valid)));

% =====================================================================
%  TABLE 5: FF Volume (FF volume / TCDSL)
% =====================================================================
fprintf('\n--- Fed Funds Volume (%% of checkable deposits, %d periods) ---\n', length(ff_valid));
fprintf('                    Data        Leontief    Cobb-Douglas\n');
fprintf('                    ----        --------    ------------\n');

% Correlations with data — restricted to valid FF window
L_stats_ff  = compute_stats(FF_L_pct, FF_D_pct, ff_valid);
CD_stats_ff = compute_stats(FF_CD_pct, FF_D_pct, ff_valid);

fprintf('corr(level)          —          %8.3f    %12.3f\n', L_stats_ff(1), CD_stats_ff(1));
fprintf('corr(diff)           —          %8.3f    %12.3f\n', L_stats_ff(2), CD_stats_ff(2));
fprintf('Std                %6.2f      %8.2f    %12.2f\n', std(FF_D_pct(ff_valid)), std(FF_L_pct(ff_valid)), std(FF_CD_pct(ff_valid)));
fprintf('Max                %6.2f      %8.2f    %12.2f\n', max(FF_D_pct(ff_valid)), max(FF_L_pct(ff_valid)), max(FF_CD_pct(ff_valid)));
fprintf('Min                %6.2f      %8.2f    %12.2f\n', min(FF_D_pct(ff_valid)), min(FF_L_pct(ff_valid)), min(FF_CD_pct(ff_valid)));
fprintf('Mean               %6.2f      %8.2f    %12.2f\n', nanmean(FF_D_pct(ff_valid)), mean(FF_L_pct(ff_valid)), mean(FF_CD_pct(ff_valid)));

% =====================================================================
%  TABLE 6: DW/FF Ratio (restricted to periods where both DW and FF available)
% =====================================================================
fprintf('\n--- DW/FF Ratio (%%) ---\n');
DWFF_L  = DW_us_L ./ FF_us_L * 100;
DWFF_CD = DW_us_CD ./ FF_us_CD * 100;
% Restrict to periods where both data series are available
dwff_valid = intersect(dw_valid, ff_valid);
DWFF_D  = DW_n_data ./ FF_n_data * 100;
fprintf('                    Data        Leontief    Cobb-Douglas\n');
fprintf('                    ----        --------    ------------\n');
fprintf('Mean               %6.2f      %8.2f    %12.2f\n', mean(DWFF_D(dwff_valid)), mean(DWFF_L(dwff_valid)), mean(DWFF_CD(dwff_valid)));
fprintf('Max                %6.2f      %8.2f    %12.2f\n', max(DWFF_D(dwff_valid)), max(DWFF_L(dwff_valid)), max(DWFF_CD(dwff_valid)));
fprintf('Min                %6.2f      %8.2f    %12.2f\n', min(DWFF_D(dwff_valid)), min(DWFF_L(dwff_valid)), min(DWFF_CD(dwff_valid)));

% =====================================================================
%  TABLE 7: Regime Breakdown
% =====================================================================
pre  = 1:72;      % Jan 2001 – Dec 2006
gfc  = 73:120;    % Jan 2007 – Dec 2010
post = 121:234;   % Jan 2011 – Jun 2020

fprintf('\n--- Volume by Regime (%% of checkable deposits) ---\n');
fprintf('                    Pre-GFC      GFC        Post-GFC\n');
fprintf('                    -------    -------      --------\n');
fprintf('DW data:            %.3f%%     %.3f%%     %.3f%%\n', nanmean(DW_n_data(pre))*100, nanmean(DW_n_data(gfc))*100, nanmean(DW_n_data(post))*100);
fprintf('DW Leontief:        %.3f%%     %.3f%%     %.3f%%\n', mean(DW_us_L(pre))*100, mean(DW_us_L(gfc))*100, mean(DW_us_L(post))*100);
fprintf('DW Cobb-Douglas:    %.3f%%     %.3f%%     %.3f%%\n', mean(DW_us_CD(pre))*100, mean(DW_us_CD(gfc))*100, mean(DW_us_CD(post))*100);
fprintf('FF data:            %.2f%%      %.2f%%     %.2f%%\n', nanmean(FF_n_data(pre))*100, nanmean(FF_n_data(gfc))*100, nanmean(FF_n_data(post))*100);
fprintf('FF Leontief:        %.2f%%      %.2f%%     %.2f%%\n', mean(FF_us_L(pre))*100, mean(FF_us_L(gfc))*100, mean(FF_us_L(post))*100);
fprintf('FF Cobb-Douglas:    %.2f%%      %.2f%%     %.2f%%\n', mean(FF_us_CD(pre))*100, mean(FF_us_CD(gfc))*100, mean(FF_us_CD(post))*100);

% =====================================================================
%  SUMMARY SCORECARD
% =====================================================================
fprintf('\n');
fprintf('================================================================================\n');
fprintf('  SCORECARD (closer to data = better)\n');
fprintf('================================================================================\n');
fprintf('                          Leontief    Cobb-Douglas    Winner\n');
fprintf('                          --------    ------------    ------\n');

% Precompute
bp_corr_L  = corr(BP_us_L_bps(datesperiod), BP_us_D_bps(datesperiod));
bp_corr_CD = corr(BP_us_CD_bps(datesperiod), BP_us_D_bps(datesperiod));
bp_std_D   = std(BP_us_D_bps(datesperiod));
bp_std_L   = std(BP_us_L_bps(datesperiod));
bp_std_CD  = std(BP_us_CD_bps(datesperiod));
cip_corr_L  = corr(CIP_L_bps(datesperiod), CIP_D_bps(datesperiod));
cip_corr_CD = corr(CIP_CD_bps(datesperiod), CIP_D_bps(datesperiod));
cip_std_D   = std(CIP_D_bps(datesperiod));
cip_std_L   = std(CIP_L_bps(datesperiod));
cip_std_CD  = std(CIP_CD_bps(datesperiod));
ff_mean_D   = nanmean(FF_D_pct(ff_valid));
ff_mean_L   = mean(FF_L_pct(ff_valid));
ff_mean_CD  = mean(FF_CD_pct(ff_valid));
dw_mean_D   = nanmean(DW_D_pct(dw_valid));
dw_mean_L   = mean(DW_L_pct(dw_valid));
dw_mean_CD  = mean(DW_CD_pct(dw_valid));

w1 = 'Leontief'; if bp_corr_CD > bp_corr_L; w1 = 'CD'; end
w2 = 'Leontief'; if abs(bp_std_CD-bp_std_D) < abs(bp_std_L-bp_std_D); w2 = 'CD'; end
w3 = 'Leontief'; if cip_corr_CD > cip_corr_L; w3 = 'CD'; end
w4 = 'Leontief'; if abs(cip_std_CD-cip_std_D) < abs(cip_std_L-cip_std_D); w4 = 'CD'; end
w5 = 'Leontief'; if abs(ff_mean_CD-ff_mean_D) < abs(ff_mean_L-ff_mean_D); w5 = 'CD'; end
w6 = 'Leontief'; if abs(dw_mean_CD-dw_mean_D) < abs(dw_mean_L-dw_mean_D); w6 = 'CD'; end

fprintf('BP_us corr(level)        %8.3f    %12.3f     %s\n', bp_corr_L, bp_corr_CD, w1);
fprintf('BP_us std (data=%.1f)    %8.1f    %12.1f     %s\n', bp_std_D, bp_std_L, bp_std_CD, w2);
fprintf('CIP corr(level)          %8.3f    %12.3f     %s\n', cip_corr_L, cip_corr_CD, w3);
fprintf('CIP std (data=%.1f)      %8.1f    %12.1f     %s\n', cip_std_D, cip_std_L, cip_std_CD, w4);
fprintf('FF mean (data=%.2f)     %8.2f    %12.2f     %s\n', ff_mean_D, ff_mean_L, ff_mean_CD, w5);
fprintf('DW mean (data=%.3f)    %8.3f    %12.3f     %s\n', dw_mean_D, dw_mean_L, dw_mean_CD, w6);

fprintf('\n--- Solver Convergence ---\n');
fprintf('Leontief US:     %d/%d (100%% fsolve)\n', sum(sigma_us_TED_flag_L > 0), length(sigma_us_TED_flag_L));
fprintf('Cobb-Douglas US: %d/%d (%d via fminbnd)\n', sum(sigma_us_TED_flag_CD > 0), length(sigma_us_TED_flag_CD), sum(sigma_us_TED_flag_CD == 10));

fprintf('\n--- Sigma Statistics ---\n');
fprintf('                        Leontief    Cobb-Douglas    Ratio\n');
fprintf('sigma_us Mean:          %8.4f    %12.4f    %5.2f\n', mean(sigma_us_L), mean(sigma_us_CD), mean(sigma_us_L)/mean(sigma_us_CD));
fprintf('sigma_us Std:           %8.4f    %12.4f    %5.2f\n', std(sigma_us_L), std(sigma_us_CD), std(sigma_us_L)/std(sigma_us_CD));
fprintf('sigma_us Max:           %8.4f    %12.4f    %5.2f\n', max(sigma_us_L), max(sigma_us_CD), max(sigma_us_L)/max(sigma_us_CD));
fprintf('sigma_eu Mean:          %8.4f    %12.4f    %5.2f\n', mean(sigma_eu_L), mean(sigma_eu_CD), mean(sigma_eu_L)/mean(sigma_eu_CD));

%% ========================================================================
%  COMPARISON PLOTS
%  ========================================================================
FSize = 14;
formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', ...
    'Fontsize', FSize, 'Box', 'Off');

%% Figure 1: Main targets — Model vs Data
figure('Name', 'Targets: Model vs Data', 'Position', [50 50 1400 800]);

% US Bond Premium
subplot(2,3,1);
plot(dates(datesperiod), BP_us_D_bps(datesperiod), 'k-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), BP_us_L_bps(datesperiod), 'b--', 'LineWidth', 1.2);
plot(dates(datesperiod), BP_us_CD_bps(datesperiod), 'r-.', 'LineWidth', 1.2);
datetick('x', 'yy'); grid on;
title('US Bond Premium (bps)', 'Interpreter', 'latex', 'FontSize', FSize+1);
legend('Data', 'Leontief', 'CD', 'Location', 'best', 'FontSize', 8);
formataxis(gca);

% EU Bond Premium
subplot(2,3,2);
plot(dates(datesperiod), BP_eu_D_bps(datesperiod), 'k-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), BP_eu_L_bps(datesperiod), 'b--', 'LineWidth', 1.2);
plot(dates(datesperiod), BP_eu_CD_bps(datesperiod), 'r-.', 'LineWidth', 1.2);
datetick('x', 'yy'); grid on;
title('EU Bond Premium (bps)', 'Interpreter', 'latex', 'FontSize', FSize+1);
legend('Data', 'Leontief', 'CD', 'Location', 'best', 'FontSize', 8);
formataxis(gca);

% CIP Deviation
subplot(2,3,3);
plot(dates(datesperiod), CIP_D_bps(datesperiod), 'k-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), CIP_L_bps(datesperiod), 'b--', 'LineWidth', 1.2);
plot(dates(datesperiod), CIP_CD_bps(datesperiod), 'r-.', 'LineWidth', 1.2);
datetick('x', 'yy'); grid on;
title('CIP Deviation (bps)', 'Interpreter', 'latex', 'FontSize', FSize+1);
legend('Data', 'Leontief', 'CD', 'Location', 'best', 'FontSize', 8);
formataxis(gca);

% DW Volume
subplot(2,3,4);
plot(dates(dw_valid), DW_D_pct(dw_valid), 'k-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), DW_L_pct(datesperiod), 'b--', 'LineWidth', 1.2);
plot(dates(datesperiod), DW_CD_pct(datesperiod), 'r-.', 'LineWidth', 1.2);
datetick('x', 'yy'); grid on;
title('DW Volume (\% checkable dep.)', 'Interpreter', 'latex', 'FontSize', FSize+1);
legend('Data', 'Leontief', 'CD', 'Location', 'best', 'FontSize', 8);
formataxis(gca);

% FF Volume
subplot(2,3,5);
plot(dates(ff_valid), FF_D_pct(ff_valid), 'k-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), FF_L_pct(datesperiod), 'b--', 'LineWidth', 1.2);
plot(dates(datesperiod), FF_CD_pct(datesperiod), 'r-.', 'LineWidth', 1.2);
datetick('x', 'yy'); grid on;
title('FF Volume (\% checkable dep.)', 'Interpreter', 'latex', 'FontSize', FSize+1);
legend('Data', 'Leontief', 'CD', 'Location', 'best', 'FontSize', 8);
formataxis(gca);

% DW/FF Ratio
subplot(2,3,6);
plot(dates(dwff_valid), DWFF_D(dwff_valid), 'k-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), DWFF_L(datesperiod), 'b--', 'LineWidth', 1.2);
plot(dates(datesperiod), DWFF_CD(datesperiod), 'r-.', 'LineWidth', 1.2);
datetick('x', 'yy'); grid on;
title('DW/FF Ratio (\%)', 'Interpreter', 'latex', 'FontSize', FSize+1);
legend('Data', 'Leontief', 'CD', 'Location', 'best', 'FontSize', 8);
formataxis(gca);

sgtitle('Model vs Data: Leontief vs Cobb-Douglas', 'Interpreter', 'latex', 'FontSize', 16);

%% Figure 2: Sigma and Theta comparison
figure('Name', 'Sigma and Theta', 'Position', [100 100 1200 500]);

subplot(1,3,1);
plot(dates(datesperiod), sigma_us_L(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_us_CD(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title('$\sigma_{US}$', 'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'CD', 'Location', 'best');
formataxis(gca);

subplot(1,3,2);
plot(dates(datesperiod), sigma_eu_L(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_eu_CD(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title('$\sigma_{EU}$', 'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'CD', 'Location', 'best');
formataxis(gca);

subplot(1,3,3);
plot(dates(datesperiod), theta_us_L(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), theta_us_CD(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title('$\theta_{US}$ (Market Tightness)', 'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'CD', 'Location', 'best');
formataxis(gca);

sgtitle('Latent Variables: Leontief vs Cobb-Douglas', 'Interpreter', 'latex', 'FontSize', 16);

%% Figure 3: Normalized sigma overlay
figure('Name', 'Normalized Dynamics', 'Position', [100 100 1000 400]);

subplot(1,2,1);
sigma_us_L_norm = (sigma_us_L - mean(sigma_us_L)) / std(sigma_us_L);
sigma_us_CD_norm = (sigma_us_CD - mean(sigma_us_CD)) / std(sigma_us_CD);
plot(dates(datesperiod), sigma_us_L_norm(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_us_CD_norm(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title(sprintf('$\\sigma_{US}$ Normalized (corr = %.3f)', corr(sigma_us_L, sigma_us_CD)), ...
    'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'CD', 'Location', 'best');
ylabel('Std. deviations');
formataxis(gca);

subplot(1,2,2);
sigma_eu_L_norm = (sigma_eu_L - mean(sigma_eu_L)) / std(sigma_eu_L);
sigma_eu_CD_norm = (sigma_eu_CD - mean(sigma_eu_CD)) / std(sigma_eu_CD);
plot(dates(datesperiod), sigma_eu_L_norm(datesperiod), 'b-', 'LineWidth', 1.5); hold on;
plot(dates(datesperiod), sigma_eu_CD_norm(datesperiod), 'r--', 'LineWidth', 1.5);
datetick('x', 'yy'); grid on;
title(sprintf('$\\sigma_{EU}$ Normalized (corr = %.3f)', corr(sigma_eu_L, sigma_eu_CD)), ...
    'Interpreter', 'latex', 'FontSize', FSize+2);
legend('Leontief', 'CD', 'Location', 'best');
ylabel('Std. deviations');
formataxis(gca);

%% ========================================================================
%  SAVE ALL RESULTS
%  ========================================================================
save('data/filter_comparison.mat', ...
    'sigma_us_L', 'sigma_eu_L', 'sigma_us_CD', 'sigma_eu_CD', ...
    'BP_us_L', 'BP_eu_L', 'BP_us_CD', 'BP_eu_CD', ...
    'CIP_L', 'CIP_CD', 'TED_us_L', 'TED_eu_L', 'TED_us_CD', 'TED_eu_CD', ...
    'theta_us_L', 'theta_eu_L', 'theta_us_CD', 'theta_eu_CD', ...
    'riskprm_L', 'riskprm_CD', ...
    'DW_us_L', 'FF_us_L', 'DW_us_CD', 'FF_us_CD', ...
    'Rb_Rm_data', 'Rb_Rm_eu_data', 'cip_data', 'DW_n_data', 'FF_n_data', ...
    'dates', 'datesperiod', 'dw_valid', 'ff_valid');

fprintf('\n========================================================\n');
fprintf('All results saved to data/filter_comparison.mat\n');
fprintf('========================================================\n');
