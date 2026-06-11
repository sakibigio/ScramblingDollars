%% plot_robustness.m
%  Overlay plot for the robustness exercise (referee 3, risk aversion):
%  one line per grid point in the (iota, lambda*) or (ploss, lambda*)
%  family.
%
%  Requires in base workspace (set by main_robustness.m):
%    which_grid : 'iota' or 'ploss'
%    dates, datesperiod, abs_scale
%    DW_n, FF_n, Rb_Rm, cip   (data series; loaded by LFX_data.mat)
%
%  Reads:
%    data/robustness_iota.mat  or  data/robustness_ploss.mat
%    data/MS_sigma_us_prob.csv (baseline regimes for shading)
%
%  Output panels:
%    F_robustness_<grid>_sigma_us.pdf      sigma_us, log scale
%    F_robustness_<grid>_TED.pdf           TED in bps
%    F_robustness_<grid>_BP_fit.pdf        BP_us demeaned vs data
%    F_robustness_<grid>_CIP_fit.pdf       CIP vs data
%    F_robustness_<grid>_DWS_fit.pdf       DW loan stock (log) vs data
%    F_robustness_<grid>_FF_fit.pdf        FF (log) vs data
%    F_robustness_<grid>_lambda.pdf        lambda* across grid
%    F_robustness_<grid>_moments.pdf       log-moment fit across grid

if ~exist('which_grid', 'var') || isempty(which_grid)
    which_grid = 'iota';
end
fprintf('\nplot_robustness: which_grid = %s\n', which_grid);

%% Load grid
switch which_grid
    case 'iota'
        S = load('data/robustness_iota.mat');
        rob = S.rob_iota;
        grid_vec = S.iota_grid(:);
        grid_label_tex = '$\iota$';
        grid_label_cb  = '\iota';
        grid_fmt = '%.3f';
        baseline_val = S.iota_baseline;
        fname_tag = 'iota';
    case 'ploss'
        S = load('data/robustness_ploss.mat');
        rob = S.rob_ploss;
        grid_vec = S.ploss_grid(:);
        grid_label_tex = '$p_{\rm loss}$';
        grid_label_cb  = 'p_{loss}';
        grid_fmt = '%.2f';
        baseline_val = S.ploss_baseline;
        fname_tag = 'ploss';
    otherwise
        error('plot_robustness: unknown which_grid ''%s''.', which_grid);
end
N_grid = numel(rob);
[~, base_idx] = min(abs(grid_vec - baseline_val));

%% Baseline regime shading
prob_file = 'data/MS_sigma_us_prob.csv';
high_t = [];
if exist(prob_file, 'file')
    pt = readtable(prob_file);
    if ismember('prob_nor', pt.Properties.VariableNames)
        prob_nor = pt.prob_nor;
        high_t = [(prob_nor < 0.5); false];
    end
end

%% Output folder
[~, username] = system('whoami');  username = strtrim(username);
if strcmp(username, 'sakibigio')
    foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs/';
elseif strcmp(username, 'sakiclaudia')
    foldername = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs/';
else
    foldername = './quantfigs/';
end
if ~exist(foldername, 'dir'), mkdir(foldername); end
do_print = strcmp(username, 'sakibigio') || strcmp(username, 'sakiclaudia');

%% Formatting + monochromatic BLUE palette (Portfolio paper convention)
% Convention:
%   Data     : solid '-',  darkest navy, thick      (primary empirical reference)
%   Baseline : dashed '--',dark navy,    thick      (primary model)
%   Variations: cycled '-.' / ':', light->mid blue, thin
FSize = 16;
formataxis = @(ax) set(ax, 'Fontname', 'Times', 'FontWeight', 'normal', ...
    'Fontsize', FSize, 'Box', 'Off', 'PlotBoxAspectRatio', [1 0.75 1]);
label_y = @(s) ylabel(s, 'Fontname', 'Times', 'Fontsize', FSize, 'Interpreter', 'latex');

data_color     = [0.00 0.00 0.25];          % data: very dark navy
base_color     = [0.10 0.15 0.50];          % baseline: dark navy
base_linestyle = '--';                       % baseline: dashed

% Variations: cycle through two non-baseline line styles, with blue ramp
linestyles_var = {'-.', ':'};
% Blue ramp light -> mid. Lightness across grid (-1) variations.
n_var = max(N_grid - 1, 1);
blue_ramp = [linspace(0.62, 0.25, n_var)', ...   % R: 0.62 (light) -> 0.25 (mid)
             linspace(0.72, 0.35, n_var)', ...   % G
             linspace(0.95, 0.75, n_var)'];      % B: stays high so it reads as blue

% Pre-compute style assignment for each non-baseline grid index
var_idx = setdiff(1:N_grid, base_idx);
var_style = struct();
for kk = 1:numel(var_idx)
    var_style(kk).idx       = var_idx(kk);
    var_style(kk).color     = blue_ramp(kk, :);
    var_style(kk).linestyle = linestyles_var{mod(kk-1, numel(linestyles_var)) + 1};
end

%% Helpers
draw_shading  = @(per) draw_regime_shading(per, dates, high_t);
function_save = @(fig, tag) save_pdf(fig, foldername, fname_tag, tag, do_print);

% Plot one variation line; uses baked-in var_style array
plot_variation = @(x, y, kk) plot(x, y, ...
    'LineWidth', 1.0, ...
    'Color', var_style(kk).color, ...
    'LineStyle', var_style(kk).linestyle);

% Which panels to draw.  Final layout (both exercises): a 2x3 grid in
% the appendix with sigma, BP, EBP, DWS, FF, lambda.  TED/CIP/moments
% are skipped.
do_panel_TED   = false;
do_panel_CIP   = false;
do_panel_lam   = true;
do_panel_mom   = false;
do_panel_EBP   = true;

%% ── Panel 1: sigma_us (log scale) ───────────────────────────────────────
fig = figure('Name', sprintf('Robustness %s: sigma_us', fname_tag), 'NumberTitle','off');
hold on;
% Baseline first (solid black, thick) — Portfolio paper convention
h_base = plot(dates(datesperiod), log(rob(base_idx).series.sigma_us_t(datesperiod)), ...
    'LineWidth', 2.5, 'LineStyle', base_linestyle, 'Color', base_color);
% Then variations: monochromatic gray, cycling line styles
for kk = 1:numel(var_style)
    ii = var_style(kk).idx;
    plot_variation(dates(datesperiod), log(rob(ii).series.sigma_us_t(datesperiod)), kk);
end
grid on; axis tight; datetick('x', 'mmm-yy', 'keeplimits'); formataxis(gca);
label_y('$\log\sigma_{us}$');
legend_with_range(h_base, [], grid_label_cb, baseline_val, grid_vec, grid_fmt);
draw_shading(datesperiod);
function_save(fig, 'sigma_us');

if do_panel_TED
%% ── Panel 2: TED (bps) ──────────────────────────────────────────────────
fig = figure('Name', sprintf('Robustness %s: TED', fname_tag), 'NumberTitle','off');
hold on;
h_base = plot(dates(datesperiod), rob(base_idx).series.TED_us_t(datesperiod) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', base_linestyle, 'Color', base_color);
for kk = 1:numel(var_style)
    ii = var_style(kk).idx;
    plot_variation(dates(datesperiod), rob(ii).series.TED_us_t(datesperiod) * abs_scale, kk);
end
h_dat = [];
if exist('TED_s_us_t', 'var')
    h_dat = plot(dates(datesperiod), TED_s_us_t(datesperiod) * abs_scale, ...
        'LineWidth', 2, 'LineStyle', '-', 'Color', data_color);
end
grid on; axis tight; datetick('x', 'mmm-yy', 'keeplimits'); formataxis(gca);
label_y('TED (bps)');
legend_with_range(h_base, h_dat, grid_label_cb, baseline_val, grid_vec, grid_fmt);
draw_shading(datesperiod);
function_save(fig, 'TED');
end

%% ── Panel 3: BP fit (demeaned, bps) ─────────────────────────────────────
NaN_us = (Rb_Rm == 0) | isnan(Rb_Rm);
Rb_Rm_plot = Rb_Rm; Rb_Rm_plot(NaN_us) = NaN;
Rb_Rm_mean = mean(Rb_Rm_plot(~NaN_us));
bp_first = find(~NaN_us, 1, 'first');
bp_last  = find(~NaN_us, 1, 'last');
bp_per   = max(bp_first, datesperiod(1)) : min(bp_last, datesperiod(end));

fig = figure('Name', sprintf('Robustness %s: BP fit', fname_tag), 'NumberTitle','off');
hold on;
plot(dates(bp_per), zeros(1, numel(bp_per)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
bp = rob(base_idx).series.BP_us_t;
h_base = plot(dates(bp_per), (bp(bp_per) - mean(bp(bp_per))) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', base_linestyle, 'Color', base_color);
for kk = 1:numel(var_style)
    ii = var_style(kk).idx;
    bp = rob(ii).series.BP_us_t;
    plot_variation(dates(bp_per), (bp(bp_per) - mean(bp(bp_per))) * abs_scale, kk);
end
h_dat = plot(dates(bp_per), (Rb_Rm_plot(bp_per) - Rb_Rm_mean) * abs_scale, ...
    'LineWidth', 2, 'LineStyle', '-', 'Color', data_color);
grid on; axis tight; datetick('x', 'mmm-yy', 'keeplimits'); formataxis(gca);
label_y('BP, demeaned (bps)');
legend_with_range(h_base, h_dat, grid_label_cb, baseline_val, grid_vec, grid_fmt);
draw_shading(bp_per);
function_save(fig, 'BP_fit');

if do_panel_EBP
%% ── Panel 3b: EBP (EU Bond Premium) fit ─────────────────────────────────
NaN_eu = (Rb_Rm_eu == 0) | isnan(Rb_Rm_eu);
Rb_Rm_eu_plot = Rb_Rm_eu; Rb_Rm_eu_plot(NaN_eu) = NaN;
Rb_Rm_eu_mean = mean(Rb_Rm_eu_plot(~NaN_eu));
ebp_first = find(~NaN_eu, 1, 'first');
ebp_last  = find(~NaN_eu, 1, 'last');
ebp_per   = max(ebp_first, datesperiod(1)) : min(ebp_last, datesperiod(end));

fig = figure('Name', sprintf('Robustness %s: EBP fit', fname_tag), 'NumberTitle','off');
hold on;
plot(dates(ebp_per), zeros(1, numel(ebp_per)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
bp_eu = rob(base_idx).series.BP_eu_t;
h_base = plot(dates(ebp_per), (bp_eu(ebp_per) - mean(bp_eu(ebp_per))) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', base_linestyle, 'Color', base_color);
for kk = 1:numel(var_style)
    ii = var_style(kk).idx;
    bp_eu = rob(ii).series.BP_eu_t;
    plot_variation(dates(ebp_per), (bp_eu(ebp_per) - mean(bp_eu(ebp_per))) * abs_scale, kk);
end
h_dat = plot(dates(ebp_per), (Rb_Rm_eu_plot(ebp_per) - Rb_Rm_eu_mean) * abs_scale, ...
    'LineWidth', 2, 'LineStyle', '-', 'Color', data_color);
grid on; axis tight; datetick('x', 'mmm-yy', 'keeplimits'); formataxis(gca);
label_y('EBP, demeaned (bps)');
legend_with_range(h_base, h_dat, grid_label_cb, baseline_val, grid_vec, grid_fmt);
draw_shading(ebp_per);
function_save(fig, 'EBP_fit');
end

if do_panel_CIP
%% ── Panel 4: CIP fit ────────────────────────────────────────────────────
NaN_cip = (cip == 0) | isnan(cip);
cip_plot = cip; cip_plot(NaN_cip) = NaN;
cip_first = find(~NaN_cip, 1, 'first');
cip_last  = find(~NaN_cip, 1, 'last');
cip_per = max(cip_first, datesperiod(1)) : min(cip_last, datesperiod(end));

fig = figure('Name', sprintf('Robustness %s: CIP fit', fname_tag), 'NumberTitle','off');
hold on;
plot(dates(cip_per), zeros(1, numel(cip_per)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
h_base = plot(dates(cip_per), rob(base_idx).series.CIP_t(cip_per) * abs_scale, ...
    'LineWidth', 2.5, 'LineStyle', base_linestyle, 'Color', base_color);
for kk = 1:numel(var_style)
    ii = var_style(kk).idx;
    plot_variation(dates(cip_per), rob(ii).series.CIP_t(cip_per) * abs_scale, kk);
end
h_dat = plot(dates(cip_per), cip_plot(cip_per) * abs_scale, ...
    'LineWidth', 2, 'LineStyle', '-', 'Color', data_color);
grid on; axis tight; datetick('x', 'mmm-yy', 'keeplimits'); formataxis(gca);
label_y('CIP (bps)');
legend_with_range(h_base, h_dat, grid_label_cb, baseline_val, grid_vec, grid_fmt);
draw_shading(cip_per);
function_save(fig, 'CIP_fit');
end

%% ── Panel 5: DW fit (loan stock, log) ───────────────────────────────────
dw_valid = ~isnan(DW_n) & DW_n > 0;
dw_first = find(dw_valid, 1, 'first');
dw_last  = find(dw_valid, 1, 'last');
dw_per = max(dw_first, datesperiod(1)) : min(dw_last, datesperiod(end));
dw_n_plot = DW_n; dw_n_plot(~dw_valid) = NaN;

fig = figure('Name', sprintf('Robustness %s: DW fit', fname_tag), 'NumberTitle','off');
hold on;
dws = rob(base_idx).series.DWS_us_t; dws(dws <= 0) = NaN;
h_base = plot(dates(dw_per), log(dws(dw_per) * 100), ...
    'LineWidth', 2.5, 'LineStyle', base_linestyle, 'Color', base_color);
for kk = 1:numel(var_style)
    ii = var_style(kk).idx;
    dws = rob(ii).series.DWS_us_t; dws(dws <= 0) = NaN;
    plot_variation(dates(dw_per), log(dws(dw_per) * 100), kk);
end
h_dat = plot(dates(dw_per), log(dw_n_plot(dw_per) * 100), ...
    'LineWidth', 2, 'LineStyle', '-', 'Color', data_color);
grid on; axis tight; datetick('x', 'mmm-yy', 'keeplimits'); formataxis(gca);
label_y('$\log$(DW stock, \% of Deposits)');
legend_with_range(h_base, h_dat, grid_label_cb, baseline_val, grid_vec, grid_fmt);
draw_shading(dw_per);
function_save(fig, 'DWS_fit');

%% ── Panel 6: FF fit (log) ───────────────────────────────────────────────
ff_valid = ~isnan(FF_n) & FF_n > 0;
ff_first = find(ff_valid, 1, 'first');
ff_last  = find(ff_valid, 1, 'last');
ff_per = max(ff_first, datesperiod(1)) : min(ff_last, datesperiod(end));
ff_n_plot = FF_n; ff_n_plot(~ff_valid) = NaN;

fig = figure('Name', sprintf('Robustness %s: FF fit', fname_tag), 'NumberTitle','off');
hold on;
ff = rob(base_idx).series.FF_us_t; ff(ff <= 0) = NaN;
h_base = plot(dates(ff_per), log(ff(ff_per) * 100), ...
    'LineWidth', 2.5, 'LineStyle', base_linestyle, 'Color', base_color);
for kk = 1:numel(var_style)
    ii = var_style(kk).idx;
    ff = rob(ii).series.FF_us_t; ff(ff <= 0) = NaN;
    plot_variation(dates(ff_per), log(ff(ff_per) * 100), kk);
end
h_dat = plot(dates(ff_per), log(ff_n_plot(ff_per) * 100), ...
    'LineWidth', 2, 'LineStyle', '-', 'Color', data_color);
grid on; axis tight; datetick('x', 'mmm-yy', 'keeplimits'); formataxis(gca);
label_y('$\log$(FF, \% of Deposits)');
legend_with_range(h_base, h_dat, grid_label_cb, baseline_val, grid_vec, grid_fmt);
draw_shading(ff_per);
function_save(fig, 'FF_fit');

if do_panel_lam
%% ── Panel 7: lambda* across grid ────────────────────────────────────────
fig = figure('Name', sprintf('Robustness %s: lambda', fname_tag), 'NumberTitle','off');
lam_vec = arrayfun(@(r) r.lambda, rob);
plot(grid_vec, lam_vec, 'o-', 'LineWidth', 2.5, ...
    'MarkerFaceColor', base_color, 'Color', base_color); hold on; grid on;
xline(baseline_val, '--k', 'baseline');
xlabel(grid_label_tex, 'Interpreter', 'latex', 'Fontsize', FSize);
ylabel('$\lambda^*$', 'Interpreter', 'latex', 'Fontsize', FSize);
formataxis(gca);
function_save(fig, 'lambda');
end

if do_panel_mom
%% ── Panel 8: moment fit ────────────────────────────────────────────────
fig = figure('Name', sprintf('Robustness %s: moment fit', fname_tag), ...
             'NumberTitle','off', 'Position', [100 100 1200 360]);
subplot(1,3,1); hold on; grid on;
m = arrayfun(@(r) r.logmean_FF_m, rob);
plot(grid_vec, m, 'o-', 'LineWidth', 2, 'Color', base_color, 'MarkerFaceColor', base_color);
yline(S.logmean_FF_d, '--k', 'data');
xlabel(grid_label_tex, 'Interpreter', 'latex'); ylabel('$\log\,\overline{FF}$', 'Interpreter', 'latex');
title('mean FF'); formataxis(gca);

subplot(1,3,2); hold on; grid on;
m = arrayfun(@(r) r.logmean_DW_m, rob);
plot(grid_vec, m, 'o-', 'LineWidth', 2, 'Color', base_color, 'MarkerFaceColor', base_color);
yline(S.logmean_DW_d, '--k', 'data');
xlabel(grid_label_tex, 'Interpreter', 'latex'); ylabel('$\log\,\overline{DW}$', 'Interpreter', 'latex');
title('mean DW'); formataxis(gca);

subplot(1,3,3); hold on; grid on;
m = arrayfun(@(r) r.logmax_DW_m, rob);
plot(grid_vec, m, 'o-', 'LineWidth', 2, 'Color', base_color, 'MarkerFaceColor', base_color);
yline(S.logmax_DW_d, '--k', 'data');
xlabel(grid_label_tex, 'Interpreter', 'latex'); ylabel('$\log\,\max\,DW$', 'Interpreter', 'latex');
title('max DW'); formataxis(gca);

function_save(fig, 'moments');
end

fprintf('plot_robustness (%s) complete.\n', which_grid);


%% ===== local helpers =====
function draw_regime_shading(per, dates_vec, high_t)
    if isempty(high_t), return; end
    yl  = ylim;
    buf = diff(yl) * 0.02;
    yl(1) = yl(1) - buf; yl(2) = yl(2) + buf;
    shade = [0.82 0.82 0.82];
    for jj = 1:numel(per) - 1
        t1 = per(jj);
        if t1 <= length(high_t) && high_t(t1)
            t2 = per(jj+1);
            xp = [dates_vec(t1) dates_vec(t2) dates_vec(t2) dates_vec(t1)];
            yp = [yl(1) yl(1) yl(2) yl(2)];
            h = patch(xp, yp, shade, 'EdgeColor', 'none', ...
                      'FaceAlpha', 0.30, 'HandleVisibility', 'off');
            try, uistack(h, 'bottom'); catch, end %#ok<NOSEMI>
        end
    end
    ylim(yl);
end

function legend_with_range(h_base, h_dat, grid_label_cb, baseline_val, grid_vec, grid_fmt)
    base_lbl = sprintf('baseline (%s=%s)', grid_label_cb, sprintf(grid_fmt, baseline_val));
    range_lbl = sprintf('%s \\in [%s,\\,%s]', grid_label_cb, ...
        sprintf(grid_fmt, min(grid_vec)), sprintf(grid_fmt, max(grid_vec)));
    % Proxy line for the range entry (matches variation style: dashed gray)
    hold on;
    h_range = plot(NaN, NaN, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
    if isempty(h_dat) || ~isgraphics(h_dat)
        legend([h_base h_range], {base_lbl, range_lbl}, ...
            'Location', 'northwest', 'Box', 'off', 'Interpreter', 'tex');
    else
        legend([h_base h_dat h_range], {base_lbl, 'Data', range_lbl}, ...
            'Location', 'northwest', 'Box', 'off', 'Interpreter', 'tex');
    end
end

function c = rgb_safe(name)
    % Fall back if rgb.m is missing.
    try
        c = rgb(name);
    catch
        switch lower(name)
            case 'navy',    c = [0 0 128]/255;
            case 'red',     c = [200 0 0]/255;
            case 'green',   c = [0 128 0]/255;
            case 'grey',    c = [0.5 0.5 0.5];
            case 'gray',    c = [0.5 0.5 0.5];
            otherwise,      c = [0 0 0];
        end
    end
end

function save_pdf(fig, foldername, tag, panel, do_print)
    if ~do_print, return; end
    out_pdf = fullfile(foldername, sprintf('F_robustness_%s_%s.pdf', tag, panel));
    try
        exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    catch
        saveas(fig, out_pdf);
    end
    fprintf('  saved %s\n', out_pdf);
end
