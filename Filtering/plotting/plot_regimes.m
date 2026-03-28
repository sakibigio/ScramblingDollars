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
%     F_IBdispersion.pdf     – (f) Dispersion in US Interbank Rates
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
[~, username] = system('whoami');
username = strtrim(username);
if strcmp(username, 'sakibigio')
    foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollars_Revision_Restud/quantfigs/';
elseif strcmp(username, 'sakiclaudia')
    foldername = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/quantfigs/';
else
    warning('Unknown user ''%s''. Figures will not be saved to Overleaf.', username);
    foldername = './quantfigs/';
    printit = 0;
end
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

try
    raw_prob = readmatrix(prob_file, 'NumHeaderLines', 1);
catch
    raw_prob = csvread(prob_file, 1, 0);  %#ok<CSVRD>  % fallback for older MATLAB
end

% sigma_us_stateprob(t) = Prob(low-risk / normal state at t)
% scrambling  <=>  high_t == 1  <=>  low-risk prob < threshold
sigma_us_stateprob = [raw_prob(:,2); 0];   % append 0 for safe end-indexing
sigma_us_high_t    = (sigma_us_stateprob < threshold);

% Load estimated parameters for LaTeX table
params_file = 'data/MS_sigma_us_params.csv';
if ~exist(params_file, 'file')
    warning('plot_regimes: %s not found. LaTeX table will not be generated.', params_file);
    write_table = false;
else
    params_tbl = readtable(params_file);
    % Build a struct for easy access by param name
    for kk = 1:height(params_tbl)
        pname = strrep(params_tbl.param{kk}, '.', '_');
        pstruct.(pname) = params_tbl.value(kk);
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
                h  = patch(xp, yp, shade, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
                uistack(h, 'bottom');
            end
        end
        ylim(yl);
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
title('Implied $\sigma^x_t$', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, 'F_sigmauseu.pdf'), 'ContentType', 'vector');
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
title('$\mathcal{TED}$ Spread (bps)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, 'F_Tedtargets.pdf'), 'ContentType', 'vector');
end

%% ── 5. FIGURE (c)  –  US BOND PREMIA ────────────────────────────────────

NaN_us     = (Rb_Rm == 0);
Rb_Rm_plot = Rb_Rm;  Rb_Rm_plot(NaN_us) = NaN;
Rb_Rm_mean = mean(Rb_Rm_plot(~NaN_us));

figure('Name', 'US Bond Premium', 'NumberTitle', 'off');
plot(dates(datesperiod), zeros(1, numel(datesperiod)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_mod = plot(dates(datesperiod), (BP_us_t(datesperiod) - mean(BP_us_t)) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r');
h_dat = plot(dates(datesperiod), (Rb_Rm_plot(datesperiod) - Rb_Rm_mean) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
legend([h_mod h_dat], 'Model', 'Data', 'Box', 'off', 'Color', 'none');
plot(dates(datesperiod), (BP_us_t(datesperiod) - mean(BP_us_t)) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(datesperiod), (Rb_Rm_plot(datesperiod) - Rb_Rm_mean) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color, 'HandleVisibility', 'off');
title('US Bond Premia (bps, demeaned)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, 'F_BPus_fit.pdf'), 'ContentType', 'vector');
end

%% ── 6. FIGURE (d)  –  EU BOND PREMIA ────────────────────────────────────

NaN_eu        = (Rb_Rm_eu == 0);
Rb_Rm_eu_plot = Rb_Rm_eu;  Rb_Rm_eu_plot(NaN_eu) = NaN;
Rb_Rm_eu_mean = mean(Rb_Rm_eu_plot(~NaN_eu));

figure('Name', 'EU Bond Premium', 'NumberTitle', 'off');
plot(dates(datesperiod), zeros(1, numel(datesperiod)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_mod = plot(dates(datesperiod), (BP_eu_t(datesperiod) - mean(BP_eu_t)) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r');
h_dat = plot(dates(datesperiod), (Rb_Rm_eu_plot(datesperiod) - Rb_Rm_eu_mean) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
legend([h_mod h_dat], 'Model', 'Data', 'Box', 'off', 'Color', 'none');
plot(dates(datesperiod), (BP_eu_t(datesperiod) - mean(BP_eu_t)) * abs_scale, ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(datesperiod), (Rb_Rm_eu_plot(datesperiod) - Rb_Rm_eu_mean) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color, 'HandleVisibility', 'off');
title('EU Bond Premia (bps, demeaned)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, 'F_BPeu_fit.pdf'), 'ContentType', 'vector');
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

figure('Name', 'US Discount Window', 'NumberTitle', 'off');
plot(dates(datesperiod), zeros(1, numel(datesperiod)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_mod = plot(dates(datesperiod), (DW_us_t(datesperiod) - mean(DW_us_t(datesperiod))) * 100, ...
    'LineWidth', 2.5, 'Color', 'r');
h_dat = plot(dates(datesperiod), (DW_n(datesperiod) - nanmean(DW_n(datesperiod))) * 100, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
legend([h_mod h_dat], 'Model', 'Data', 'Box', 'off', 'Color', 'none');
plot(dates(datesperiod), (DW_us_t(datesperiod) - mean(DW_us_t(datesperiod))) * 100, ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(datesperiod), (DW_n(datesperiod) - nanmean(DW_n(datesperiod))) * 100, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color, 'HandleVisibility', 'off');
title('US Discount Window / Deposits (\%, demeaned)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, 'F_DWus_fit.pdf'), 'ContentType', 'vector');
end

%% ── 7b. FIGURE (e2) –  FF / DEPOSITS ────────────────────────────────────

figure('Name', 'US Fed Funds Volume', 'NumberTitle', 'off');
plot(dates(datesperiod), zeros(1, numel(datesperiod)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
h_mod = plot(dates(datesperiod), (FF_us_t(datesperiod) - mean(FF_us_t(datesperiod))) * 100, ...
    'LineWidth', 2.5, 'Color', 'r');
h_dat = plot(dates(datesperiod), (FF_n(datesperiod) - nanmean(FF_n(datesperiod))) * 100, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color);
grid on; axis tight;
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
legend([h_mod h_dat], 'Model', 'Data', 'Box', 'off', 'Color', 'none');
plot(dates(datesperiod), (FF_us_t(datesperiod) - mean(FF_us_t(datesperiod))) * 100, ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
plot(dates(datesperiod), (FF_n(datesperiod) - nanmean(FF_n(datesperiod))) * 100, ...
    'LineWidth', 2.5, 'LineStyle', '-.', 'Color', data_color, 'HandleVisibility', 'off');
title('US Interbank Volume / Deposits (\%, demeaned)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, 'F_FFus_fit.pdf'), 'ContentType', 'vector');
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
title('Dispersion in US Interbank Rates (bps)', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, 'F_IBdispersion.pdf'), 'ContentType', 'vector');
end

%% ── 9. FIGURE (tab-estimated) – FILTERED PROBABILITIES ──────────────────

figure('Name', 'Filtered Regime Probabilities', 'NumberTitle', 'off');
plot(dates(datesperiod), 0.5 * ones(1, numel(datesperiod)), ...
    'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k'); hold on;
plot(dates(datesperiod), sigma_us_stateprob(datesperiod), ...
    'LineWidth', 2.5, 'Color', 'r');
grid on; axis tight; ylim([0 1]);
datetick('x', 'mmm-yy', 'keeplimits');
formataxis(gca);
regime_patches(dates, datesperiod, sigma_us_high_t);
plot(dates(datesperiod), sigma_us_stateprob(datesperiod), ...
    'LineWidth', 2.5, 'Color', 'r', 'HandleVisibility', 'off');
title('Prob. Low Funding Risk State', 'Interpreter', 'latex', 'Fontsize', 16);
if printit
    exportgraphics(gcf, fullfile(foldername, 'F_sigmaus_states.pdf'), 'ContentType', 'vector');
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
    fprintf(fid, '    $\\hat{\\sigma}_{ss}$ & %.1f & %.1f \\\\\n', ...
        pstruct.sigma_ss_scr, pstruct.sigma_ss_nor);
    fprintf(fid, '                        & (%.1f) & (%.1f) \\\\\n', ...
        pstruct.sigma_ss_scr_se, pstruct.sigma_ss_nor_se);
    fprintf(fid, '    $\\rho^{\\sigma,us}$  & %.1f & %.1f \\\\\n', ...
        pstruct.rho_scr, pstruct.rho_nor);
    fprintf(fid, '                        & (%.1f) & (%.1f) \\\\\n', ...
        pstruct.rho_scr_se, pstruct.rho_nor_se);
    fprintf(fid, '    $\\Sigma^{\\sigma,us}$ & %.1f & %.1f \\\\\n', ...
        pstruct.Sigma_scr, pstruct.Sigma_nor);
    fprintf(fid, '                         & (%.1f) & (%.1f) \\\\\n', ...
        pstruct.Sigma_scr_se, pstruct.Sigma_nor_se);
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
