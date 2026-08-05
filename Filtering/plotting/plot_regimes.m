%% plot_regimes.m  –  Filter Results with Markov Regime Shading
% =========================================================================
% Generates the paper figures from the filtering exercise:
%
%   Fig. "Data Variables" (6 panels):
%     F_sigmauseu.pdf        – (a) Implied sigma_t
%     F_Tedtargets.pdf       – (b) TED Spread
%     F_BPus_fit.pdf         – (c) US Bond Premia
%     F_BPeu_fit.pdf         – (d) EU Bond Premia
%     F_DWus_fit.pdf         – (e) US Window / Interbank (% of deposits)
%     F_FFus_fit.pdf         – (f) US Fed Funds (% of deposits)
%     F_DWFFscatter_data.pdf       – DW vs FF scatter, data, regime-colored (paper)
%     F_DWFFscatter_model.pdf      – DW vs FF scatter, model, regime-colored (paper)
%     F_DWFFscatter_combined.pdf   – DW vs FF scatter, data | model combined (internal)
%     F_DWFFscatter_data_detrend.pdf      – DW vs FF scatter, data, FF detrended (paper)
%     F_DWFFscatter_model_detrend.pdf     – DW vs FF scatter, model, FF detrended (paper)
%     F_DWFFscatter_combined_detrend.pdf  – DW vs FF scatter, data | model, FF detrended (internal)
%     F_dDWdFFscatter_data.pdf     – Δlog DW vs Δlog FF, data, regime-colored
%     F_dDWdFFscatter_model.pdf    – Δlog DW vs Δlog FF, model, regime-colored
%     F_dDWdFFscatter_combined.pdf – Δlog DW vs Δlog FF, data | model (internal)
%     F_DWdFFscatter_data.pdf      – log DW vs Δlog FF, data, regime-colored
%     F_DWdFFscatter_model.pdf     – log DW vs Δlog FF, model, regime-colored
%     F_DWdFFscatter_combined.pdf  – log DW vs Δlog FF, data | model (internal)
%     F_IBdispersion.pdf     – Dispersion in US Interbank Rates
%
%   Fig. "Estimation and Filtered Probabilities":
%     F_sigmaus_states.pdf   – Filtered regime probabilities
%
% PIPELINE
%   main_filter.m  →  markov_estimation.jl  →  plot_regimes.m
%
% REQUIRED FILE
%   data/MS_sigma_us_prob.csv   – produced by markov_estimation.jl
%
% REQUIRED WORKSPACE VARIABLES (all produced by main_filter.m)
%   sigma_us_t, sigma_eu_t, sigma_us_TED_t, sigma_eu_TED_t
%   BP_us_t, BP_eu_t, TED_us_t, TED_eu_t, CIP_t
%   DW_us_t, FF_us_t          – model DW/FF as fraction of checkable deposits
%   DW_n, FF_n                – data DW/FF as fraction of checkable deposits
%                               (from LFX_data.mat; NaN before series start)
%   Rb_Rm, Rb_Rm_eu, cip, Chi_D_US
%   dates, datesperiod, abs_scale
%
% (c) Saki Bigio
% =========================================================================

%% ── 0. SETTINGS ─────────────────────────────────────────────────────────

% Inherit printit from main_filter.m if available; default to 0
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

threshold = 0.5;    % Prob. of low-risk state below which = scrambling regime

% Output folder: match convention from other plotting scripts
foldername = overleaf_dir();
if ~exist(foldername, 'dir'), mkdir(foldername); end

% Formatting – guard allows standalone use outside main_filter.m
if ~exist('formataxis', 'var') || ~isa(formataxis, 'function_handle')
    formataxis = @(ax) set(ax, 'Fontname', 'Times', 'FontWeight', 'normal', ...
        'Fontsize', 14, 'Box', 'On');
end

% Colours
eu_color   = [0.40 0.10 1.00];    % blue-violet  – EU model series
data_color = [0.70 0.10 0.40];    % dark magenta  – data series
ff_color   = [0.20 0.50 0.20];    % dark green    – FF interbank

%% ── 1. LOAD MARKOV REGIME PROBABILITIES ─────────────────────────────────

prob_file = 'data/MS_sigma_us_prob.csv';
if ~exist(prob_file, 'file')
    error(['plot_regimes: data/MS_sigma_us_prob.csv not found.\n' ...
           'Run markov_estimation.jl first, then re-run plot_regimes.m.']);
end

% Read by column name so the convention is robust to MarSwitching's
% arbitrary state ordering. CSV has columns: prob_nor, prob_scr.
prob_tbl = readtable(prob_file);
vn = prob_tbl.Properties.VariableNames;
if ismember('prob_nor', vn) && ismember('prob_scr', vn)
    prob_nor_vec = prob_tbl.prob_nor;
elseif ismember('prob_state1', vn) && ismember('prob_state2', vn)
    error(['plot_regimes: %s uses legacy column names (prob_state1/2).\n' ...
           'Re-run markov_estimation.jl to regenerate with prob_nor/prob_scr columns.'], ...
           prob_file);
else
    error('plot_regimes: %s has unexpected columns: %s', prob_file, strjoin(vn, ', '));
end

% sigma_us_stateprob(t) = Prob(normal state at t)  [name kept for back-compat]
% scrambling  <=>  high_t == 1  <=>  normal prob < threshold
sigma_us_stateprob = [prob_nor_vec; 0];   % append 0 for safe end-indexing
sigma_us_high_t    = (sigma_us_stateprob < threshold);

% Load estimated parameters for LaTeX table
params_file = 'data/MS_sigma_us_params.csv';
if ~exist(params_file, 'file')
    warning('plot_regimes: %s not found. LaTeX table will not be generated.', params_file);
    write_table = false;
else
    try
        params_tbl = readtable(params_file);
        pnames = params_tbl.param;
        pvals  = params_tbl.value;
    catch
        % Fallback: lower-level parsing for older MATLAB / quirky envs
        fid = fopen(params_file, 'r');
        raw = textscan(fid, '%s%f', 'Delimiter', ',', 'HeaderLines', 1);
        fclose(fid);
        pnames = raw{1};
        pvals  = raw{2};
    end
    % Build a struct for easy access by param name
    for kk = 1:numel(pnames)
        pname = strrep(pnames{kk}, '.', '_');
        pstruct.(pname) = pvals(kk);
    end
    write_table = true;
end

%% ── 2. HELPER: REGIME SHADING ───────────────────────────────────────────
% Call AFTER data lines are drawn so ylim is stable.
% Patches are pushed to background with uistack.

    function regime_patches(dates_vec, per, high_t)
        % dates_vec : full datenum dates vector
        % per       : index vector for sample period (datesperiod)
        % high_t    : logical(T_full), 1 = scrambling period
        yl  = ylim;
        buf = diff(yl) * 0.02;
        yl(1) = yl(1) - buf;  yl(2) = yl(2) + buf;
        shade = [0.82 0.82 0.82];
        for ii = 1:numel(per) - 1
            if high_t(per(ii)) == 1
                xp = [dates_vec(per(ii))   dates_vec(per(ii+1)) ...
                      dates_vec(per(ii+1)) dates_vec(per(ii))];
                yp = [yl(1) yl(1) yl(2) yl(2)];
                % FaceAlpha=0.3 keeps the bands visible whether they end
                % up under or over the lines (uistack fails on yyaxis figs
                % where patch and lines live in different hidden axes).
                h = patch(xp, yp, shade, 'EdgeColor', 'none', 'FaceAlpha', 0.30);
                try
                    uistack(h, 'bottom');
                catch
                    % yyaxis: uistack can't permute children across the
                    % paired axes objects. Skip silently; FaceAlpha=0.3
                    % keeps the underlying lines legible.
                end
            end
        end
        ylim(yl);
    end

    function plot_regline(x, y, color)
        % OLS y = a + b*x over finite (x,y) pairs. Draws the line across the
        % current x-axis with HandleVisibility off so it doesn't enter legends.
        mask = isfinite(x) & isfinite(y);
        if sum(mask) < 2, return; end
        b  = polyfit(x(mask), y(mask), 1);
        xx = xlim;
        plot(xx, polyval(b, xx), '-', 'Color', color, 'LineWidth', 2, ...
            'HandleVisibility', 'off');
    end

%% ── 3. FIGURE (a)  –  SIGMA ESTIMATES ───────────────────────────────────

figure('Name', 'Sigma US/EU', 'NumberTitle', 'off');
plot(dates(datesperiod), mean(sigma_us_t) * ones(1, numel(datesperiod)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_us = plot(dates(datesperiod), sigma_us_TED_t(datesperiod), 'LineWidth', 2.5, 'Color', 'r');
h_eu = plot(dates(datesperiod), sigma_eu_TED_t(datesperiod), 'LineWidth', 2.5, ...
    'LineStyle', '-.', 'Color', eu_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
legend([h_us h_eu], 'US', 'EU', 'Box', 'off', 'Color', 'none', 'Location', 'northwest');
plot(dates(datesperiod), sigma_us_TED_t(datesperiod), 'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(datesperiod), sigma_eu_TED_t(datesperiod), 'LineWidth', 2.5, ...
    'LineStyle', '-.', 'Color', eu_color, 'HandleVisibility', 'off');
% title('Implied $\sigma^x_t$', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_sigmauseu' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 4. FIGURE (b)  –  TED SPREAD ────────────────────────────────────────

figure('Name', 'TED Spread', 'NumberTitle', 'off');
plot(dates(datesperiod), mean(TED_us_t) * abs_scale * ones(1, numel(datesperiod)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_us = plot(dates(datesperiod), TED_us_t(datesperiod) * abs_scale, 'LineWidth', 2.5, 'Color', 'r');
h_eu = plot(dates(datesperiod), TED_eu_t(datesperiod) * abs_scale, 'LineWidth', 2.5, ...
    'LineStyle', '-.', 'Color', eu_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
legend([h_us h_eu], 'US', 'EU', 'Box', 'off', 'Color', 'none');
plot(dates(datesperiod), TED_us_t(datesperiod) * abs_scale, 'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(datesperiod), TED_eu_t(datesperiod) * abs_scale, 'LineWidth', 2.5, ...
    'LineStyle', '-.', 'Color', eu_color, 'HandleVisibility', 'off');
% title('$\mathcal{TED}$ Spread (bps)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_Tedtargets' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 5. FIGURE (c)  –  US BOND PREMIA ────────────────────────────────────
% Truncated to the window where Rb_Rm data exists (zeros are treated as
% missing). Means are computed on the valid subset to align demeaning.

NaN_us     = (Rb_Rm == 0) | isnan(Rb_Rm);
Rb_Rm_plot = Rb_Rm;  Rb_Rm_plot(NaN_us) = NaN;
Rb_Rm_mean = mean(Rb_Rm_plot(~NaN_us));

bp_first = find(~NaN_us, 1, 'first');
bp_last  = find(~NaN_us, 1, 'last');
bp_per   = max(bp_first, datesperiod(1)) : min(bp_last, datesperiod(end));

figure('Name', 'US Bond Premium', 'NumberTitle', 'off');
plot(dates(bp_per), zeros(1, numel(bp_per)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_mod = plot(dates(bp_per), (BP_us_t(bp_per) - mean(BP_us_t(bp_per))) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r');
h_dat = plot(dates(bp_per), (Rb_Rm_plot(bp_per) - Rb_Rm_mean) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, bp_per, sigma_us_high_t);
legend([h_mod h_dat], 'Model', 'Data', 'Box', 'off', 'Color', 'none');
plot(dates(bp_per), (BP_us_t(bp_per) - mean(BP_us_t(bp_per))) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(bp_per), (Rb_Rm_plot(bp_per) - Rb_Rm_mean) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color, 'HandleVisibility', 'off');
% title('US Bond Premia (bps, demeaned)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_BPus_fit' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 6. FIGURE (d)  –  EU BOND PREMIA ────────────────────────────────────
% Truncated to the window where Rb_Rm_eu data exists.

NaN_eu        = (Rb_Rm_eu == 0) | isnan(Rb_Rm_eu);
Rb_Rm_eu_plot = Rb_Rm_eu;  Rb_Rm_eu_plot(NaN_eu) = NaN;
Rb_Rm_eu_mean = mean(Rb_Rm_eu_plot(~NaN_eu));

bp_eu_first = find(~NaN_eu, 1, 'first');
bp_eu_last  = find(~NaN_eu, 1, 'last');
bp_eu_per   = max(bp_eu_first, datesperiod(1)) : min(bp_eu_last, datesperiod(end));

figure('Name', 'EU Bond Premium', 'NumberTitle', 'off');
plot(dates(bp_eu_per), zeros(1, numel(bp_eu_per)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_mod = plot(dates(bp_eu_per), (BP_eu_t(bp_eu_per) - mean(BP_eu_t(bp_eu_per))) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r');
h_dat = plot(dates(bp_eu_per), (Rb_Rm_eu_plot(bp_eu_per) - Rb_Rm_eu_mean) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, bp_eu_per, sigma_us_high_t);
legend([h_mod h_dat], 'Model', 'Data', 'Box', 'off', 'Color', 'none');
plot(dates(bp_eu_per), (BP_eu_t(bp_eu_per) - mean(BP_eu_t(bp_eu_per))) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(bp_eu_per), (Rb_Rm_eu_plot(bp_eu_per) - Rb_Rm_eu_mean) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color, 'HandleVisibility', 'off');
% title('EU Bond Premia (bps, demeaned)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_BPeu_fit' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 6b. FIGURE  –  CIP DEVIATION ────────────────────────────────────────
% Model: CIP_t = UIP_t + Rm_eu .* riskprm_t  (from main_filter.m).
% Data : cip (column E of LFX_datainputs.xlsx, /freq/1e4 in load_data.m).
% Truncated to periods where the data series is non-zero / non-NaN. Zero
% line is a meaningful reference (no-arbitrage benchmark), so no demeaning.

NaN_cip   = (cip == 0) | isnan(cip);
cip_plot  = cip;  cip_plot(NaN_cip) = NaN;

cip_first = find(~NaN_cip, 1, 'first');
cip_last  = find(~NaN_cip, 1, 'last');
cip_per   = max(cip_first, datesperiod(1)) : min(cip_last, datesperiod(end));

figure('Name', 'CIP Deviation', 'NumberTitle', 'off');
plot(dates(cip_per), zeros(1, numel(cip_per)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_mod = plot(dates(cip_per), CIP_t(cip_per) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r');
h_dat = plot(dates(cip_per), cip_plot(cip_per) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, cip_per, sigma_us_high_t);
legend([h_mod h_dat], 'Model', 'Data', 'Box', 'off', 'Color', 'none');
plot(dates(cip_per), CIP_t(cip_per) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(cip_per), cip_plot(cip_per) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color, 'HandleVisibility', 'off');
% title('CIP Deviation (bps)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_CIP_fit' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7. FIGURE (e)  –  US WINDOW / INTERBANK ─────────────────────────────
%
% DW_us_t, FF_us_t : model DW and FF volumes as fraction of checkable
%                    deposits (units consistent with DW_n, FF_n).
% DW_n, FF_n       : data series from load_data.m (LFX_data.mat):
%                    DW_n = WLCFLPCL / TCDSL  (primary credit / deposits)
%                    FF_n = FF volume / TCDSL  (fed funds / deposits)
%                    Leading NaNs where series start after Jan 2001.
% Multiply by 100 to display as percent of deposits.

if ~exist('DW_n', 'var') || ~exist('FF_n', 'var')
    error('plot_regimes: DW_n and FF_n not found in workspace. Check LFX_data.mat was loaded.');
end

%% ── 7a. FIGURE (e)  –  DW / DEPOSITS ────────────────────────────────────
% Truncated to periods where DW data exists (DW_n starts ~2003). Internal
% gaps are NaN-padded so they render as blanks. Single y-axis (model and
% data are in comparable units after log-percent conversion).

dw_valid = ~isnan(DW_n) & DW_n > 0;
dw_first = find(dw_valid, 1, 'first');
dw_last  = find(dw_valid, 1, 'last');
dw_per   = max(dw_first, datesperiod(1)) : min(dw_last, datesperiod(end));

dw_n_plot = DW_n;  dw_n_plot(~dw_valid) = NaN;

% NOTE: Flow plot disabled -- model DW_us_t is a flow per period while data
% DW_n is an outstanding stock, so the comparison is conceptually mismatched.
% See Section 7a' for the stock-vs-stock figure (the right comparison).
% figure('Name', 'US Discount Window', 'NumberTitle', 'off');
% yyaxis left
% h_dat = plot(dates(dw_per), log(dw_n_plot(dw_per) * 100), ...
%     'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
% ylabel('Data: $\log$(\% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);
% yyaxis right
% h_mod = plot(dates(dw_per), log(DW_us_t(dw_per) * 100), ...
%     'LineWidth', 2.5, 'Color', 'r');
% ylabel('Model: $\log$(\% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);
% grid on; axis tight;
% datetick('x', 'mmm-yy', 'keeplimits');
% formataxis(gca);
% regime_patches(dates, dw_per, sigma_us_high_t);
% legend([h_mod h_dat], 'Model', 'Data', 'Box', 'off', 'Color', 'none');
% if printit
%     exportgraphics(gcf, fullfile(foldername, ['F_DWus_fit' mt_suffix '.pdf']), 'ContentType', 'vector');
% end

%% ── 7a'. FIGURE  –  DW LOAN STOCK / DEPOSITS  (log) ─────────────────────
% Model DWS_us_t built in main_filter.m via DWS_t = DW_t + (1-delta)*DWS_{t-1}.
% Data DW_n is already an outstanding-stock series (WLCFLPCL/TCDSL), so this
% is the apples-to-apples comparison.  Same dw_per / regime shading as 7a.
if exist('DWS_us_t', 'var')
    dws_us_plot = DWS_us_t;  dws_us_plot(~dw_valid | DWS_us_t <= 0) = NaN;
    figure('Name', 'US Discount Window: Loan Stock', 'NumberTitle', 'off');
    yyaxis left
    h_dat = plot(dates(dw_per), log(dw_n_plot(dw_per) * 100), ...
        'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
    ylabel('Data: $\log$(\% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);
    yyaxis right
    h_mod = plot(dates(dw_per), log(dws_us_plot(dw_per) * 100), ...
        'LineWidth', 2.5, 'Color', 'r');
    ylabel('Model: $\log$(\% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);
    grid on; axis tight;
    datetick('x', 'mmm-yy', 'keeplimits');
    formataxis(gca);
    regime_patches(dates, dw_per, sigma_us_high_t);
    legend_lbl = 'Model (stock)';
    if exist('delta_DWS', 'var')
        legend_lbl = sprintf('Model stock (\\delta=%.2f)', delta_DWS);
    end
    legend([h_mod h_dat], legend_lbl, 'Data', 'Box', 'off', 'Color', 'none');
    % title('US Discount Window -- Loan Stock', 'Interpreter', 'latex', 'Fontsize', 16);
    if printit
        exportgraphics(gcf, fullfile(foldername, ['F_DWSus_fit' mt_suffix '.pdf']), 'ContentType', 'vector');
    end
end

%% ── 7a''. FIGURE  –  EU DW (flow) / DEPOSITS  (DISABLED) ────────────────
% Disabled: flow-vs-stock comparison is mismatched (see 7a').  Stock plot
% in Section 7a''' is the right comparison.
%
% dw_eu_valid = DW_eu_t > 0 & isfinite(DW_eu_t);
% if any(dw_eu_valid)
%     ... (yyaxis flow plot)
% end

%% ── 7a'''. FIGURE  –  EU DW LOAN STOCK / DEPOSITS  (log, dual axes) ─────
% Stock built in main_filter.m the same way as the US side.  yyaxis used
% when EU data is available (incomplete series, levels off).
if exist('DWS_eu_t', 'var')
    dws_eu_plot = DWS_eu_t;  dws_eu_plot(DWS_eu_t <= 0 | ~isfinite(DWS_eu_t)) = NaN;
    figure('Name', 'EU Discount Window: Loan Stock', 'NumberTitle', 'off');
    if exist('DW_eu_n', 'var')
        dw_eu_n_valid = ~isnan(DW_eu_n) & DW_eu_n > 0;
        dw_eu_n_plot  = DW_eu_n;  dw_eu_n_plot(~dw_eu_n_valid) = NaN;
        yyaxis left
        h_dat = plot(dates(dw_per), log(dw_eu_n_plot(dw_per) * 100), ...
            'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
        ylabel('Data: $\log$(\% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);
        yyaxis right
        h_mod = plot(dates(dw_per), log(dws_eu_plot(dw_per) * 100), ...
            'LineWidth', 2.5, 'Color', eu_color);
        ylabel('Model: $\log$(\% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);
        legend_lbl_eu = 'Model stock (EU)';
        if exist('delta_DWS', 'var')
            legend_lbl_eu = sprintf('Model stock EU (\\delta=%.2f)', delta_DWS);
        end
        legend([h_mod h_dat], legend_lbl_eu, 'Data (EU)', 'Box', 'off', 'Color', 'none', 'AutoUpdate', 'off');
    else
        h_mod = plot(dates(dw_per), log(dws_eu_plot(dw_per) * 100), ...
            'LineWidth', 2.5, 'Color', eu_color); hold on;
        ylabel('$\log$(\% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);
        legend_lbl_eu = 'Model stock (EU)';
        if exist('delta_DWS', 'var')
            legend_lbl_eu = sprintf('Model stock EU (\\delta=%.2f)', delta_DWS);
        end
        legend(h_mod, legend_lbl_eu, 'Box', 'off', 'Color', 'none', 'AutoUpdate', 'off');
    end
    grid on; axis tight;
    datetick('x', 'mmm-yy', 'keeplimits');
    formataxis(gca);
    regime_patches(dates, dw_per, sigma_us_high_t);
    % title('EU Discount Window -- Loan Stock', 'Interpreter', 'latex', 'Fontsize', 16);
    if printit
        exportgraphics(gcf, fullfile(foldername, ['F_DWSeu_fit' mt_suffix '.pdf']), 'ContentType', 'vector');
    end
end

%% ── 7b. FIGURE (f)  –  FF / DEPOSITS ────────────────────────────────────
% Truncated to periods where FF data exists (FF_n starts ~2006), with
% separate y-axes for data and model so their levels are independently
% legible (model and data differ by an order of magnitude in this fit).

ff_valid = ~isnan(FF_n) & FF_n > 0;
ff_first = find(ff_valid, 1, 'first');
ff_last  = find(ff_valid, 1, 'last');
ff_per   = max(ff_first, datesperiod(1)) : min(ff_last, datesperiod(end));

% NaN-out any internal gaps so the plot leaves them blank
ff_n_plot   = FF_n;        ff_n_plot(~ff_valid)  = NaN;
ff_mod_plot = FF_us_t;     ff_mod_plot(~ff_valid) = NaN;

figure('Name', 'US Fed Funds Volume', 'NumberTitle', 'off');

yyaxis left
h_dat = plot(dates(ff_per), log(ff_n_plot(ff_per) * 100), ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
ylabel('Data: $\log$(FF, \% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);

yyaxis right
h_mod = plot(dates(ff_per), log(ff_mod_plot(ff_per) * 100), ...
    'LineWidth', 2.5, 'Color', 'r');
ylabel('Model: $\log$(FF, \% of Deposits)', 'Interpreter', 'latex', 'Fontsize', 16);

grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, ff_per, sigma_us_high_t);
legend([h_dat h_mod], 'Data', 'Model', 'Box', 'off', 'Color', 'none');
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_FFus_fit' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%{
%  ════════════════════════════════════════════════════════════════════════
%  Sections 7c through 7k DISABLED: only the detrended scatter set
%  (sections 7l/7m/7n) is kept.  Setup variables that 7l-7n depend on are
%  replicated right before 7l (search for "Setup for detrended scatter").
%  ════════════════════════════════════════════════════════════════════════
%% ── 7c. FIGURE  –  DW vs FF SCATTER, DATA (regime-colored) ──────────────
%
% Scatter plot of log(DW) vs log(FF) using the same series as panels (e)
% and (f), colored by Markov regime: scrambling vs normal.

scr_idx = logical(sigma_us_high_t(datesperiod));
nor_idx = ~scr_idx;

scr_color = [0.85 0.20 0.20];   % red    – scrambling
nor_color = [0.20 0.35 0.70];   % blue   – normal

% Data series (handle leading NaNs)
DW_n_per = DW_n(datesperiod) * 100;
FF_n_per = FF_n(datesperiod) * 100;
ok_dat   = ~isnan(DW_n_per) & ~isnan(FF_n_per) & DW_n_per > 0 & FF_n_per > 0;

figure('Name', 'DW vs FF Scatter – Data', 'NumberTitle', 'off');
h_nor = scatter(log(DW_n_per(ok_dat & nor_idx)), log(FF_n_per(ok_dat & nor_idx)), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(log(DW_n_per(ok_dat & scr_idx)), log(FF_n_per(ok_dat & scr_idx)), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on; axis tight;
plot_regline(log(DW_n_per(ok_dat)), log(FF_n_per(ok_dat)), [0 0 0]);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\log(\mathrm{FF}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('Data', 'Interpreter', 'latex', 'Fontsize', 16);
xlim_data = xlim; ylim_data = ylim;
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWFFscatter_data' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7d. FIGURE  –  DW vs FF SCATTER, MODEL (regime-colored) ─────────────

DW_m_per = DW_us_t(datesperiod) * 100;
FF_m_per = FF_us_t(datesperiod) * 100;
ok_mod   = DW_m_per >= 0 & FF_m_per >= 0;

figure('Name', 'DW vs FF Scatter – Model', 'NumberTitle', 'off');
h_nor = scatter(log(DW_m_per(ok_mod & nor_idx)), log(FF_m_per(ok_mod & nor_idx)), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(log(DW_m_per(ok_mod & scr_idx)), log(FF_m_per(ok_mod & scr_idx)), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on; axis tight;
plot_regline(log(DW_m_per(ok_mod)), log(FF_m_per(ok_mod)), [0 0 0]);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\log(\mathrm{FF}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('Model', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWFFscatter_model' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7e. FIGURE  –  DW vs FF SCATTER, COMBINED (internal use) ────────────
%
% Side-by-side panels (data | model) on shared axes for quick internal
% comparison. Paper uses the two separate figures above.

xl_lo = min([log(DW_n_per(ok_dat)); log(DW_m_per(ok_mod))]);
xl_hi = max([log(DW_n_per(ok_dat)); log(DW_m_per(ok_mod))]);
yl_lo = min([log(FF_n_per(ok_dat)); log(FF_m_per(ok_mod))]);
yl_hi = max([log(FF_n_per(ok_dat)); log(FF_m_per(ok_mod))]);
xpad  = 0.05 * (xl_hi - xl_lo);
ypad  = 0.05 * (yl_hi - yl_lo);

figure('Name', 'DW vs FF Scatter – Combined', 'NumberTitle', 'off', ...
    'Position', [100 100 1100 480]);

subplot(1, 2, 1);
h_nor = scatter(log(DW_n_per(ok_dat & nor_idx)), log(FF_n_per(ok_dat & nor_idx)), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(log(DW_n_per(ok_dat & scr_idx)), log(FF_n_per(ok_dat & scr_idx)), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim([xl_lo - xpad, xl_hi + xpad]); ylim([yl_lo - ypad, yl_hi + ypad]);
plot_regline(log(DW_n_per(ok_dat)), log(FF_n_per(ok_dat)), [0 0 0]);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\log(\mathrm{FF}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('(a) Data', 'Interpreter', 'latex', 'Fontsize', 16);

subplot(1, 2, 2);
h_nor = scatter(log(DW_m_per(ok_mod & nor_idx)), log(FF_m_per(ok_mod & nor_idx)), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(log(DW_m_per(ok_mod & scr_idx)), log(FF_m_per(ok_mod & scr_idx)), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim([xl_lo - xpad, xl_hi + xpad]); ylim([yl_lo - ypad, yl_hi + ypad]);
plot_regline(log(DW_m_per(ok_mod)), log(FF_m_per(ok_mod)), [0 0 0]);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\log(\mathrm{FF}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('(b) Model', 'Interpreter', 'latex', 'Fontsize', 16);

if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWFFscatter_combined' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7f. FIGURE  –  ΔlogDW vs ΔlogFF SCATTER, DATA (regime-colored) ──────
%
% First-difference (log-change) scatter: ΔlogX_t = log X_t − log X_{t-1}.
% Point at date t is labeled by the regime at t. Force column orientation
% throughout so logical masks compose without implicit expansion (DW_n,
% dates, sigma_us_high_t are stored with inconsistent shapes in the
% main_filter workspace).

scr_idx = logical(sigma_us_high_t(datesperiod));  scr_idx = scr_idx(:);
nor_idx = ~scr_idx;

scr_color = [0.85 0.20 0.20];   % red    – scrambling
nor_color = [0.20 0.35 0.70];   % blue   – normal

% Levels (% of deposits), forced to columns
DW_n_per = DW_n(datesperiod) * 100;  DW_n_per = DW_n_per(:);
FF_n_per = FF_n(datesperiod) * 100;  FF_n_per = FF_n_per(:);
DW_m_per = DW_us_t(datesperiod) * 100;  DW_m_per = DW_m_per(:);
FF_m_per = FF_us_t(datesperiod) * 100;  FF_m_per = FF_m_per(:);

% Log-differences. log of NaN/non-positive → NaN, propagated by diff.
dlogDW_n = diff(log(DW_n_per));
dlogFF_n = diff(log(FF_n_per));
dlogDW_m = diff(log(DW_m_per));
dlogFF_m = diff(log(FF_m_per));

% Regime indices aligned to differences (drop first date)
scr_idx_d = scr_idx(2:end);
nor_idx_d = nor_idx(2:end);

ok_dat_d = isfinite(dlogDW_n) & isfinite(dlogFF_n);
ok_mod_d = isfinite(dlogDW_m) & isfinite(dlogFF_m);

% Robust axis limits — percentile-based so a few outlier months don't blow
% up the scale (the data Δlog FF has tail jumps that swamp the bulk pattern).
all_dDW = [dlogDW_n(ok_dat_d); dlogDW_m(ok_mod_d)];
all_dFF = [dlogFF_n(ok_dat_d); dlogFF_m(ok_mod_d)];
xdw_lo  = quantile(all_dDW, 0.02);  xdw_hi = quantile(all_dDW, 0.98);
yff_lo  = quantile(all_dFF, 0.02);  yff_hi = quantile(all_dFF, 0.98);
xdw_pad = max(0.05 * (xdw_hi - xdw_lo), 1e-6);
yff_pad = max(0.05 * (yff_hi - yff_lo), 1e-6);
xlim_dDW = [xdw_lo - xdw_pad, xdw_hi + xdw_pad];
ylim_dFF = [yff_lo - yff_pad, yff_hi + yff_pad];

figure('Name', 'dlogDW vs dlogFF Scatter – Data', 'NumberTitle', 'off');
h_nor = scatter(dlogDW_n(ok_dat_d & nor_idx_d), dlogFF_n(ok_dat_d & nor_idx_d), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(dlogDW_n(ok_dat_d & scr_idx_d), dlogFF_n(ok_dat_d & scr_idx_d), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim(xlim_dDW); ylim(ylim_dFF);
plot_regline(dlogDW_n(ok_dat_d & nor_idx_d), dlogFF_n(ok_dat_d & nor_idx_d), nor_color);
plot_regline(dlogDW_n(ok_dat_d & scr_idx_d), dlogFF_n(ok_dat_d & scr_idx_d), scr_color);
formataxis(gca);
xlabel('$\Delta \log(\mathrm{DW}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\Delta \log(\mathrm{FF}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('Data ($\Delta\log$)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_dDWdFFscatter_data' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7g. FIGURE  –  ΔlogDW vs ΔlogFF SCATTER, MODEL (regime-colored) ─────

figure('Name', 'dlogDW vs dlogFF Scatter – Model', 'NumberTitle', 'off');
h_nor = scatter(dlogDW_m(ok_mod_d & nor_idx_d), dlogFF_m(ok_mod_d & nor_idx_d), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(dlogDW_m(ok_mod_d & scr_idx_d), dlogFF_m(ok_mod_d & scr_idx_d), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim(xlim_dDW); ylim(ylim_dFF);
plot_regline(dlogDW_m(ok_mod_d & nor_idx_d), dlogFF_m(ok_mod_d & nor_idx_d), nor_color);
plot_regline(dlogDW_m(ok_mod_d & scr_idx_d), dlogFF_m(ok_mod_d & scr_idx_d), scr_color);
formataxis(gca);
xlabel('$\Delta \log(\mathrm{DW}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\Delta \log(\mathrm{FF}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('Model ($\Delta\log$)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_dDWdFFscatter_model' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7h. FIGURE  –  ΔlogDW vs ΔlogFF SCATTER, COMBINED (internal) ───────

figure('Name', 'dlogDW vs dlogFF Scatter – Combined', 'NumberTitle', 'off', ...
    'Position', [100 100 1100 480]);

subplot(1, 2, 1);
h_nor = scatter(dlogDW_n(ok_dat_d & nor_idx_d), dlogFF_n(ok_dat_d & nor_idx_d), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(dlogDW_n(ok_dat_d & scr_idx_d), dlogFF_n(ok_dat_d & scr_idx_d), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim(xlim_dDW); ylim(ylim_dFF);
plot_regline(dlogDW_n(ok_dat_d & nor_idx_d), dlogFF_n(ok_dat_d & nor_idx_d), nor_color);
plot_regline(dlogDW_n(ok_dat_d & scr_idx_d), dlogFF_n(ok_dat_d & scr_idx_d), scr_color);
formataxis(gca);
xlabel('$\Delta \log(\mathrm{DW}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\Delta \log(\mathrm{FF}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('(a) Data', 'Interpreter', 'latex', 'Fontsize', 16);

subplot(1, 2, 2);
h_nor = scatter(dlogDW_m(ok_mod_d & nor_idx_d), dlogFF_m(ok_mod_d & nor_idx_d), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(dlogDW_m(ok_mod_d & scr_idx_d), dlogFF_m(ok_mod_d & scr_idx_d), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim(xlim_dDW); ylim(ylim_dFF);
plot_regline(dlogDW_m(ok_mod_d & nor_idx_d), dlogFF_m(ok_mod_d & nor_idx_d), nor_color);
plot_regline(dlogDW_m(ok_mod_d & scr_idx_d), dlogFF_m(ok_mod_d & scr_idx_d), scr_color);
formataxis(gca);
xlabel('$\Delta \log(\mathrm{DW}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\Delta \log(\mathrm{FF}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('(b) Model', 'Interpreter', 'latex', 'Fontsize', 16);

if printit
    exportgraphics(gcf, fullfile(foldername, ['F_dDWdFFscatter_combined' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7i. FIGURE  –  logDW vs ΔlogFF SCATTER, DATA (regime-colored) ──────
%
% Tests whether scrambling periods show high co-movement between the
% LEVEL of discount-window borrowing and the FLOW (change) of fed-funds
% activity: high DW + large positive ΔFF would suggest DW-funded
% interbank lending during stress.
%
% X: log(DW_t / Deposits × 100)
% Y: Δlog(FF_t / Deposits) = log FF_t − log FF_{t-1}
% Each point at date t is labeled by the regime at t.

% Align DW level to the difference dates (drop t=1)
logDW_n_t = log(DW_n_per(2:end));
logDW_m_t = log(DW_m_per(2:end));
ok_lvldat = isfinite(logDW_n_t) & isfinite(dlogFF_n);
ok_lvlmod = isfinite(logDW_m_t) & isfinite(dlogFF_m);

% Robust x-axis (log DW). y-axis reuses ylim_dFF from section 7f.
all_lDW  = [logDW_n_t(ok_lvldat); logDW_m_t(ok_lvlmod)];
xldw_lo  = quantile(all_lDW, 0.02);  xldw_hi = quantile(all_lDW, 0.98);
xldw_pad = max(0.05 * (xldw_hi - xldw_lo), 1e-6);
xlim_lDW = [xldw_lo - xldw_pad, xldw_hi + xldw_pad];

figure('Name', 'logDW vs dlogFF Scatter – Data', 'NumberTitle', 'off');
h_nor = scatter(logDW_n_t(ok_lvldat & nor_idx_d), dlogFF_n(ok_lvldat & nor_idx_d), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(logDW_n_t(ok_lvldat & scr_idx_d), dlogFF_n(ok_lvldat & scr_idx_d), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim(xlim_lDW); ylim(ylim_dFF);
plot_regline(logDW_n_t(ok_lvldat & nor_idx_d), dlogFF_n(ok_lvldat & nor_idx_d), nor_color);
plot_regline(logDW_n_t(ok_lvldat & scr_idx_d), dlogFF_n(ok_lvldat & scr_idx_d), scr_color);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\Delta \log(\mathrm{FF}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('Data (level DW vs $\Delta\log$ FF)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWdFFscatter_data' mt_suffix '.pdf']), 'ContentType', 'vector');
end

% Within-regime correlation (data) for quick console diagnostic
if any(ok_lvldat & scr_idx_d)
    rho_scr_dat = corr(logDW_n_t(ok_lvldat & scr_idx_d), dlogFF_n(ok_lvldat & scr_idx_d));
else
    rho_scr_dat = NaN;
end
if any(ok_lvldat & nor_idx_d)
    rho_nor_dat = corr(logDW_n_t(ok_lvldat & nor_idx_d), dlogFF_n(ok_lvldat & nor_idx_d));
else
    rho_nor_dat = NaN;
end
fprintf('\nCorr( log DW_t , dlog FF_t )  [DATA]:\n');
fprintf('  Scrambling: % .3f   Normal: % .3f\n', rho_scr_dat, rho_nor_dat);

%% ── 7j. FIGURE  –  logDW vs ΔlogFF SCATTER, MODEL (regime-colored) ─────

figure('Name', 'logDW vs dlogFF Scatter – Model', 'NumberTitle', 'off');
h_nor = scatter(logDW_m_t(ok_lvlmod & nor_idx_d), dlogFF_m(ok_lvlmod & nor_idx_d), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(logDW_m_t(ok_lvlmod & scr_idx_d), dlogFF_m(ok_lvlmod & scr_idx_d), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim(xlim_lDW); ylim(ylim_dFF);
plot_regline(logDW_m_t(ok_lvlmod & nor_idx_d), dlogFF_m(ok_lvlmod & nor_idx_d), nor_color);
plot_regline(logDW_m_t(ok_lvlmod & scr_idx_d), dlogFF_m(ok_lvlmod & scr_idx_d), scr_color);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\Delta \log(\mathrm{FF}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('Model (level DW vs $\Delta\log$ FF)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWdFFscatter_model' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7k. FIGURE  –  logDW vs ΔlogFF SCATTER, COMBINED (internal) ───────

figure('Name', 'logDW vs dlogFF Scatter – Combined', 'NumberTitle', 'off', ...
    'Position', [100 100 1100 480]);

subplot(1, 2, 1);
h_nor = scatter(logDW_n_t(ok_lvldat & nor_idx_d), dlogFF_n(ok_lvldat & nor_idx_d), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(logDW_n_t(ok_lvldat & scr_idx_d), dlogFF_n(ok_lvldat & scr_idx_d), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim(xlim_lDW); ylim(ylim_dFF);
plot_regline(logDW_n_t(ok_lvldat & nor_idx_d), dlogFF_n(ok_lvldat & nor_idx_d), nor_color);
plot_regline(logDW_n_t(ok_lvldat & scr_idx_d), dlogFF_n(ok_lvldat & scr_idx_d), scr_color);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\Delta \log(\mathrm{FF}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('(a) Data', 'Interpreter', 'latex', 'Fontsize', 16);

subplot(1, 2, 2);
h_nor = scatter(logDW_m_t(ok_lvlmod & nor_idx_d), dlogFF_m(ok_lvlmod & nor_idx_d), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(logDW_m_t(ok_lvlmod & scr_idx_d), dlogFF_m(ok_lvlmod & scr_idx_d), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim(xlim_lDW); ylim(ylim_dFF);
plot_regline(logDW_m_t(ok_lvlmod & nor_idx_d), dlogFF_m(ok_lvlmod & nor_idx_d), nor_color);
plot_regline(logDW_m_t(ok_lvlmod & scr_idx_d), dlogFF_m(ok_lvlmod & scr_idx_d), scr_color);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\Delta \log(\mathrm{FF}/\mathrm{Deposits})$', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('(b) Model', 'Interpreter', 'latex', 'Fontsize', 16);

if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWdFFscatter_combined' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%}

%% ── Setup for detrended scatter (replicates what 7c-7f used to do) ──────
scr_idx   = logical(sigma_us_high_t(datesperiod));  scr_idx = scr_idx(:);
nor_idx   = ~scr_idx;
scr_color = [0.85 0.20 0.20];
nor_color = [0.20 0.35 0.70];
DW_n_per  = DW_n(datesperiod) * 100;     DW_n_per = DW_n_per(:);
FF_n_per  = FF_n(datesperiod) * 100;     FF_n_per = FF_n_per(:);
DW_m_per  = DW_us_t(datesperiod) * 100;  DW_m_per = DW_m_per(:);
FF_m_per  = FF_us_t(datesperiod) * 100;  FF_m_per = FF_m_per(:);
ok_dat    = isfinite(DW_n_per) & isfinite(FF_n_per) & DW_n_per > 0 & FF_n_per > 0;

% Outlier exclusion: drop the single data observation with the lowest DW
% and highest FF (an extreme out-of-sample point that obscures the fit).
% Identified by the largest log(FF) - log(DW) gap within ok_dat.
ratio_test = log(FF_n_per) - log(DW_n_per);
ratio_test(~ok_dat) = -Inf;
[~, idx_outlier_dat] = max(ratio_test);
ok_dat_no_out = ok_dat;
ok_dat_no_out(idx_outlier_dat) = false;
fprintf('Detrended scatter: dropping data outlier at t=%d (log FF - log DW = %.3f)\n', ...
    idx_outlier_dat, max(ratio_test));

%% ── 7l. FIGURE  –  log DW vs DETRENDED log FF, DATA (regime-colored) ────
%
% Same series as 7c, but log(FF) is detrended with a linear time trend
% (OLS of log FF on t = 1..T over the ok_dat sample). Isolates cyclical
% co-movement between DW and FF from any secular trend in FF.

t_idx       = (1:numel(FF_n_per))';
logFF_n_dt  = log(FF_n_per);
b_trend     = polyfit(t_idx(ok_dat_no_out), logFF_n_dt(ok_dat_no_out), 1);   % [slope intercept]
logFF_n_dt  = logFF_n_dt - polyval(b_trend, t_idx);            % residual

figure('Name', 'DW vs detrended FF Scatter – Data', 'NumberTitle', 'off');
h_nor = scatter(log(DW_n_per(ok_dat_no_out & nor_idx)), logFF_n_dt(ok_dat_no_out & nor_idx), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(log(DW_n_per(ok_dat_no_out & scr_idx)), logFF_n_dt(ok_dat_no_out & scr_idx), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim([-10, 5]); ylim([-1, 1]);
plot_regline(log(DW_n_per(ok_dat_no_out)), logFF_n_dt(ok_dat_no_out), [0 0 0]);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\log(\mathrm{FF}/\mathrm{Deposits} \times 100)$ (detrended)', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('Data (FF linearly detrended)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWFFscatter_data_detrend' mt_suffix '.pdf']), 'ContentType', 'vector');
end

% Pooled slope for reference (outlier excluded)
b_pooled = polyfit(log(DW_n_per(ok_dat_no_out)), logFF_n_dt(ok_dat_no_out), 1);
fprintf('\nLevels (log DW, detrended log FF) — pooled OLS slope [DATA]: %.3f\n', b_pooled(1));

%% ── 7m. FIGURE  –  log DW vs DETRENDED log FF, MODEL (regime-colored) ───
%
% Model counterpart to 7l. log(FF_m) detrended with the same linear
% time-trend procedure applied to the simulated model series.

ok_mod_lvl  = DW_m_per > 0 & FF_m_per > 0;   % strict — log-safe
logFF_m_dt  = log(FF_m_per);
b_trend_mod = polyfit(t_idx(ok_mod_lvl), logFF_m_dt(ok_mod_lvl), 1);
logFF_m_dt  = logFF_m_dt - polyval(b_trend_mod, t_idx);

figure('Name', 'DW vs detrended FF Scatter – Model', 'NumberTitle', 'off');
h_nor = scatter(log(DW_m_per(ok_mod_lvl & nor_idx)), logFF_m_dt(ok_mod_lvl & nor_idx), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(log(DW_m_per(ok_mod_lvl & scr_idx)), logFF_m_dt(ok_mod_lvl & scr_idx), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim([-10, 5]); ylim([-1, 1]);
plot_regline(log(DW_m_per(ok_mod_lvl)), logFF_m_dt(ok_mod_lvl), [0 0 0]);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\log(\mathrm{FF}/\mathrm{Deposits} \times 100)$ (detrended)', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('Model (FF linearly detrended)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWFFscatter_model_detrend' mt_suffix '.pdf']), 'ContentType', 'vector');
end

b_pooled_mod = polyfit(log(DW_m_per(ok_mod_lvl)), logFF_m_dt(ok_mod_lvl), 1);
fprintf('Levels (log DW, detrended log FF) — pooled OLS slope [MODEL]: %.3f\n', b_pooled_mod(1));

%% ── 7n. FIGURE  –  log DW vs DETRENDED log FF, COMBINED (internal) ─────
%
% Side-by-side panels (data | model) on shared axes. Paper uses 7l and
% 7m separately; this combined figure is for internal comparison.

xl_lo_dt = min([log(DW_n_per(ok_dat_no_out)); log(DW_m_per(ok_mod_lvl))]);
xl_hi_dt = max([log(DW_n_per(ok_dat_no_out)); log(DW_m_per(ok_mod_lvl))]);
yl_lo_dt = min([logFF_n_dt(ok_dat_no_out); logFF_m_dt(ok_mod_lvl)]);
yl_hi_dt = max([logFF_n_dt(ok_dat_no_out); logFF_m_dt(ok_mod_lvl)]);
xpad_dt  = 0.05 * (xl_hi_dt - xl_lo_dt);
ypad_dt  = 0.05 * (yl_hi_dt - yl_lo_dt);

figure('Name', 'DW vs detrended FF Scatter – Combined', 'NumberTitle', 'off', ...
    'Position', [100 100 1100 480]);

subplot(1, 2, 1);
h_nor = scatter(log(DW_n_per(ok_dat_no_out & nor_idx)), logFF_n_dt(ok_dat_no_out & nor_idx), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(log(DW_n_per(ok_dat_no_out & scr_idx)), logFF_n_dt(ok_dat_no_out & scr_idx), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim([-10, 5]); ylim([-1, 1]);
plot_regline(log(DW_n_per(ok_dat_no_out)), logFF_n_dt(ok_dat_no_out), [0 0 0]);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\log(\mathrm{FF}/\mathrm{Deposits} \times 100)$ (detrended)', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('(a) Data', 'Interpreter', 'latex', 'Fontsize', 16);

subplot(1, 2, 2);
h_nor = scatter(log(DW_m_per(ok_mod_lvl & nor_idx)), logFF_m_dt(ok_mod_lvl & nor_idx), ...
    40, nor_color, 'o', 'filled', 'MarkerFaceAlpha', 0.7); hold on;
h_scr = scatter(log(DW_m_per(ok_mod_lvl & scr_idx)), logFF_m_dt(ok_mod_lvl & scr_idx), ...
    60, scr_color, 'd', 'filled', 'MarkerFaceAlpha', 0.75);
grid on;
xlim([-10, 5]); ylim([-1, 1]);
plot_regline(log(DW_m_per(ok_mod_lvl)), logFF_m_dt(ok_mod_lvl), [0 0 0]);
formataxis(gca);
xlabel('$\log(\mathrm{DW}/\mathrm{Deposits} \times 100)$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\log(\mathrm{FF}/\mathrm{Deposits} \times 100)$ (detrended)', 'Interpreter', 'latex', 'Fontsize', 16);
legend([h_nor h_scr], 'Normal', 'Scrambling', 'Box', 'off', 'Color', 'none', 'Location', 'best');
% title('(b) Model', 'Interpreter', 'latex', 'Fontsize', 16);

if printit
    exportgraphics(gcf, fullfile(foldername, ['F_DWFFscatter_combined_detrend' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 7o. FIGURE  –  MARKET TIGHTNESS vs STL FINANCIAL STRESS INDEX ───────
% Model market tightness theta_us_t (deficit/surplus ratio in interbank
% market) plotted against the St Louis Fed Financial Stress Index.  Dual
% y-axes since they're in different units; visual co-movement test.
if exist('STL_FSI', 'var') && exist('theta_us_t', 'var')
    stl_valid = isfinite(STL_FSI);
    stl_plot  = STL_FSI;  stl_plot(~stl_valid) = NaN;
    stl_first = find(stl_valid, 1, 'first');
    stl_last  = find(stl_valid, 1, 'last');
    stl_per   = max(stl_first, datesperiod(1)) : min(stl_last, datesperiod(end));

    figure('Name', 'Market Tightness vs STL FSI', 'NumberTitle', 'off');
    yyaxis left
    h_dat = plot(dates(stl_per), stl_plot(stl_per), ...
        'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
    ylabel('STL Financial Stress Index', 'Interpreter', 'latex', 'Fontsize', 16);
    yyaxis right
    h_mod = plot(dates(stl_per), log((theta_us_t(stl_per))), ...
        'LineWidth', 2.5, 'Color', 'r');
    ylabel('Model log $\theta_{us,t}$ (market tightness)', 'Interpreter', 'latex', 'Fontsize', 16);
    grid on; axis tight;
    datetick('x', 'mmm-yy', 'keeplimits');
    formataxis(gca);
    regime_patches(dates, stl_per, sigma_us_high_t);
    legend([h_mod h_dat], 'Model log $\theta_{us}$', 'STL FSI', ...
        'Box', 'off', 'Color', 'none', 'Interpreter', 'latex');
    % title('Market Tightness vs Financial Stress', 'Interpreter', 'latex', 'Fontsize', 16);
    if printit
        exportgraphics(gcf, fullfile(foldername, ['F_theta_vs_STL' mt_suffix '.pdf']), 'ContentType', 'vector');
    end

    % Correlation diagnostic
    both = isfinite(STL_FSI(stl_per)) & isfinite(theta_us_t(stl_per));
    if sum(both) > 2
        c_lvl = corr(STL_FSI(stl_per(both)), theta_us_t(stl_per(both)));
        fprintf('corr(STL FSI, theta_us_t): level = %.3f  (n=%d)\n', c_lvl, sum(both));
    end
end

%% ── 8. FIGURE (f)  –  DISPERSION IN US INTERBANK RATES ──────────────────

figure('Name', 'US Interbank Dispersion', 'NumberTitle', 'off');
plot(dates(datesperiod), mean(Chi_D_US(datesperiod)) * abs_scale * ones(1, numel(datesperiod)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), Chi_D_US(datesperiod) * abs_scale, ...
    'LineWidth', 2.5, 'Color', data_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), Chi_D_US(datesperiod) * abs_scale, ...
    'LineWidth', 2.5, 'Color', data_color, 'HandleVisibility', 'off');
% title('Dispersion in US Interbank Rates (bps)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_IBdispersion' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 9. FIGURE (tab-estimated) – FILTERED PROBABILITIES ──────────────────

% Truncate the plotted window by one observation -- the very last value is
% dropped per request (it can be a placeholder when the Markov estimator
% appends a 0 at t = T+1; see line 113).
prob_per = datesperiod(1:end-1);

figure('Name', 'Filtered Regime Probabilities', 'NumberTitle', 'off');
plot(dates(prob_per), 0.5 * ones(1, numel(prob_per)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(prob_per), 1-sigma_us_stateprob(prob_per), ...
    'LineWidth', 2.5, 'Color', 'r');
grid on; axis tight; ylim([0 1]);
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, prob_per, sigma_us_high_t);
plot(dates(prob_per), 1-sigma_us_stateprob(prob_per), ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
% title('Prob. Low Funding Risk State', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, ['F_sigmaus_states' mt_suffix '.pdf']), 'ContentType', 'vector');
end

%% ── 10. REGIME MOMENTS (console) ────────────────────────────────────────

index_scr = logical(sigma_us_high_t(datesperiod));
index_nor = ~index_scr;

fprintf('\n=== Conditional Moments by Regime ===\n');
fprintf('                      Scrambling      Normal\n');
fprintf('E[sigma_us]           %9.4f   %9.4f\n', ...
    mean(sigma_us_t(sigma_us_high_t(1:end-1))), ...
    mean(sigma_us_t(~sigma_us_high_t(1:end-1))));
fprintf('E[BP_us]   (bps)      %9.2f   %9.2f\n', ...
    mean(BP_us_t(index_scr))  * abs_scale, mean(BP_us_t(index_nor))  * abs_scale);
fprintf('E[TED_us]  (bps)      %9.2f   %9.2f\n', ...
    mean(TED_us_t(index_scr)) * abs_scale, mean(TED_us_t(index_nor)) * abs_scale);
fprintf('E[CIP]     (bps)      %9.2f   %9.2f\n', ...
    mean(CIP_t(index_scr))    * abs_scale, mean(CIP_t(index_nor))    * abs_scale);
fprintf('E[DW_us]   (%% dep.)  %9.4f   %9.4f\n', ...
    mean(DW_us_t(index_scr))  * 100,       mean(DW_us_t(index_nor))  * 100);
fprintf('E[FF_us]   (%% dep.)  %9.4f   %9.4f\n', ...
    mean(FF_us_t(index_scr))  * 100,       mean(FF_us_t(index_nor))  * 100);
fprintf('Frac. scrambling:     %.1f%%\n', 100 * mean(index_scr));

%% ── 11. WRITE LATEX TABLE ────────────────────────────────────────────────
if write_table && printit
    tex_file = fullfile(foldername, ['tab_markov_estimates' mt_suffix '.tex']);
    fid = fopen(tex_file, 'wt');
    fprintf(fid, '\\begin{tabular}{lcc}\n');
    fprintf(fid, '    \\toprule\n');
    fprintf(fid, '    \\multicolumn{3}{c}{\\textbf{Within Regime Processes}} \\\\\n');
    fprintf(fid, '    \\midrule\n');
    fprintf(fid, '    \\textbf{Coefficient} & \\textbf{Scrambling} & \\textbf{Normal} \\\\\n');
    fprintf(fid, '    \\midrule\n');
    fprintf(fid, '    $\\hat{\\sigma}_{ss}$ & %.3f & %.3f \\\\\n', ...
        pstruct.sigma_ss_scr, pstruct.sigma_ss_nor);
    fprintf(fid, '                        & (%.3f) & (%.3f) \\\\\n', ...
        pstruct.sigma_ss_scr_se, pstruct.sigma_ss_nor_se);
    fprintf(fid, '    $\\rho^{\\sigma,us}$  & %.3f & %.3f \\\\\n', ...
        pstruct.rho_scr, pstruct.rho_nor);
    fprintf(fid, '                        & (%.3f) & (%.3f) \\\\\n', ...
        pstruct.rho_scr_se, pstruct.rho_nor_se);
    fprintf(fid, '    $\\Sigma^{\\sigma,us}$ & %.3f & %.3f \\\\\n', ...
        pstruct.Sigma_scr, pstruct.Sigma_nor);
    fprintf(fid, '                         & (%.3f) & (%.3f) \\\\\n', ...
        pstruct.Sigma_scr_se, pstruct.Sigma_nor_se);
    % γ (shared μ-loading) — added if MS_sigma_us_params.csv carries it
    if isfield(pstruct, 'gamma_mu') && isfinite(pstruct.gamma_mu)
        fprintf(fid, '    \\midrule\n');
        fprintf(fid, '    $\\gamma^{\\mu,us}$   & \\multicolumn{2}{c}{%.3f} \\\\\n', pstruct.gamma_mu);
        if isfield(pstruct, 'gamma_mu_se') && isfinite(pstruct.gamma_mu_se)
            fprintf(fid, '                          & \\multicolumn{2}{c}{(%.3f)} \\\\\n', pstruct.gamma_mu_se);
        end
    end
    fprintf(fid, '    \\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\vspace{6pt}\n');
    fprintf(fid, '\\begin{tabular}{lcc}\n');
    fprintf(fid, '    \\toprule\n');
    fprintf(fid, '    \\multicolumn{2}{c}{\\textbf{Transition Matrix}} \\\\\n');
    fprintf(fid, '    \\midrule\n');
    fprintf(fid, '    & \\textbf{Scrambling} & \\textbf{Normal} \\\\\n');
    fprintf(fid, '    \\midrule\n');
    fprintf(fid, '    $\\Pr(Z_t = Z_{t-1})$ & %.1f\\%% & %.1f\\%% \\\\\n', ...
        (1 - pstruct.trans_scr) * 100, (1 - pstruct.trans_nor) * 100);
    fprintf(fid, '    \\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);
    fprintf('LaTeX table saved to: %s\n', tex_file);
end

fprintf('\nplot_regimes.m complete. 7 figures saved to:\n  %s\n', foldername);
