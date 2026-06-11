%% estimate_robustness_iota.m
%  Robustness check for referee 3 (risk aversion):
%  Vary iota_ss over a grid; for each value, re-estimate the single
%  best-fitting lambda (eta fixed at its baseline calibration), run the
%  full filter, then re-estimate Markov regimes via markov_estimation.jl
%  for that pair.
%
%  Requires in workspace (set by main_robustness.m):
%    mu_us, mu_eu, im_us, im_eu, TED_s_us_t, TED_s_eu_t
%    FF_n, DW_n, sigma_us, sigma_eu, matching_type, varrho, freq
%    pi_us_ss, pi_eu_ss, use_mu_baseline, mubar_us_t, mubar_eu_t
%    abs_scale, dates, datesperiod
%    lambda_baseline, eta_baseline, iota_baseline, ploss_baseline
%    julia_path   (set by main_robustness.m; '' if Julia not run)
%
%  Output: data/robustness_iota.mat with struct rob_iota
%  Per-pair Julia outputs copied to data/robustness/iota/g{ii}/.

fprintf('\n========================================\n');
fprintf('Robustness exercise: grid over iota_ss\n');
fprintf('Free param: lambda  |  eta fixed at %.4f\n', eta_baseline);
fprintf('========================================\n');

%% Grid
% Lower bound 0.06: below this the TED constraint TED > (1-eta)*iota binds
% and the inversion has no admissible sigma.  10 uniform points.
iota_grid = linspace(0.06, 0.12, 10);
N_grid    = numel(iota_grid);

%% Sample selection for moment objective
valid_FF = ~isnan(FF_n) & FF_n > 0;
valid_DW = ~isnan(DW_n) & DW_n > 0;
idx_FF   = find(valid_FF);
idx_DW   = find(valid_DW);

logmean_FF_d = log(mean(FF_n(idx_FF)));
logmean_DW_d = log(mean(DW_n(idx_DW)));
logmax_DW_d  = log(max(DW_n(idx_DW)));

%% Inputs struct for run_filter_series
in = struct();
in.mu_us = mu_us; in.mu_eu = mu_eu;
in.im_us = im_us; in.im_eu = im_eu;
in.TED_s_us_t = TED_s_us_t; in.TED_s_eu_t = TED_s_eu_t;
in.matching_type = matching_type; in.varrho = varrho; in.freq = freq;
in.pi_us_ss = pi_us_ss; in.pi_eu_ss = pi_eu_ss;
in.use_mu_baseline = use_mu_baseline;
in.mubar_us_t = mubar_us_t; in.mubar_eu_t = mubar_eu_t;
in.sigma_us_init = sigma_us; in.sigma_eu_init = sigma_eu;
in.abs_scale = abs_scale;

%% Output folder for per-pair Julia products
rob_dir = fullfile('data', 'robustness', 'iota');
if ~exist(rob_dir, 'dir'), mkdir(rob_dir); end

%% Pre-allocate result struct array
rob_iota = repmat(struct('iota', NaN, 'lambda', NaN, 'eta', NaN, ...
    'fval', NaN, 'exitflag', NaN, 'series', [], 'prob_nor_us', [], ...
    'logmean_FF_m', NaN, 'logmean_DW_m', NaN, 'logmax_DW_m', NaN), N_grid, 1);

%% Loop over grid
options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 1000, ...
                   'TolX', 1e-6, 'TolFun', 1e-10);

% Save baseline regime CSV so we can restore after per-pair Julia runs.
baseline_prob_csv  = 'data/MS_sigma_us_prob.csv';
baseline_RWshock   = 'RW_shock.csv';
baseline_prob_back = 'data/MS_sigma_us_prob_BASELINE.csv';
baseline_RW_back   = 'RW_shock_BASELINE.csv';
if exist(baseline_prob_csv, 'file') && ~exist(baseline_prob_back, 'file')
    copyfile(baseline_prob_csv, baseline_prob_back);
end
if exist(baseline_RWshock, 'file') && ~exist(baseline_RW_back, 'file')
    copyfile(baseline_RWshock, baseline_RW_back);
end

for ii = 1:N_grid
    iota_g = iota_grid(ii);
    fprintf('\n--- (%d/%d) iota_ss = %.4f ---\n', ii, N_grid, iota_g);

    % 1. Estimate lambda by fminsearch (eta fixed at baseline)
    obj = @(lam) obj_lambda(lam, eta_baseline, iota_g, ploss_baseline, in, ...
        idx_FF, idx_DW, logmean_FF_d, logmean_DW_d, logmax_DW_d);
    x0 = lambda_baseline;
    [lam_opt, fval, ef] = fminsearch(obj, x0, options);

    fprintf('  lambda* = %.6f   obj = %.4e   ef = %d\n', lam_opt, fval, ef);

    % 2. Re-run filter at the optimum, full output
    series = run_filter_series(lam_opt, eta_baseline, iota_g, ploss_baseline, in);

    % 3. Per-pair Markov estimation via Julia
    prob_nor = [];
    if ~isempty(julia_path)
        % Write RW_shock.csv with this pair's sigma series
        T_use = min(length(series.sigma_us_t), length(series.sigma_eu_t));
        sigma_mat = [series.sigma_us_t(1:T_use) series.sigma_eu_t(1:T_use)];
        rw_tbl = array2table(sigma_mat, 'VariableNames', {'sigma_us', 'sigma_eu'});
        writetable(rw_tbl, baseline_RWshock);

        % Invoke Julia
        cmd = sprintf('"%s" "%s"', julia_path, 'markov_estimation.jl');
        [status, jresult] = system(cmd);
        if status == 0
            try
                prob_tbl = readtable('data/MS_sigma_us_prob.csv');
                prob_nor = prob_tbl.prob_nor;
            catch e
                warning('  Failed to read per-pair regime CSV: %s', e.message);
            end
        else
            warning('  Julia failed for iota=%.4f: %s', iota_g, jresult);
        end

        % Archive per-pair output
        pair_dir = fullfile(rob_dir, sprintf('g%02d', ii));
        if ~exist(pair_dir, 'dir'), mkdir(pair_dir); end
        for f = {'MS_sigma_us_prob.csv', 'MS_sigma_us_params.csv', ...
                 'MS_sigma_eu_prob.csv', 'MS_params.csv'}
            src = fullfile('data', f{1});
            if exist(src, 'file'), copyfile(src, fullfile(pair_dir, f{1})); end
        end
    end

    % 4. Store
    rob_iota(ii).iota     = iota_g;
    rob_iota(ii).lambda   = lam_opt;
    rob_iota(ii).eta      = eta_baseline;
    rob_iota(ii).ploss    = ploss_baseline;
    rob_iota(ii).fval     = fval;
    rob_iota(ii).exitflag = ef;
    rob_iota(ii).series   = series;
    rob_iota(ii).prob_nor_us = prob_nor;
    rob_iota(ii).logmean_FF_m = log(mean(series.FF_us_t(idx_FF)));
    rob_iota(ii).logmean_DW_m = log(mean(series.DW_us_t(idx_DW)));
    rob_iota(ii).logmax_DW_m  = log(max(series.DW_us_t(idx_DW)));

    fprintf('  log mean FF: model=%.4f  data=%.4f\n', ...
        rob_iota(ii).logmean_FF_m, logmean_FF_d);
    fprintf('  log mean DW: model=%.4f  data=%.4f\n', ...
        rob_iota(ii).logmean_DW_m, logmean_DW_d);
end

%% Restore baseline regime CSV (so plot_regimes still works)
if exist(baseline_prob_back, 'file')
    copyfile(baseline_prob_back, baseline_prob_csv);
end
if exist(baseline_RW_back, 'file')
    copyfile(baseline_RW_back, baseline_RWshock);
end

%% Save
save('data/robustness_iota.mat', 'rob_iota', 'iota_grid', ...
    'lambda_baseline', 'eta_baseline', 'iota_baseline', 'ploss_baseline', ...
    'logmean_FF_d', 'logmean_DW_d', 'logmax_DW_d', '-v7.3');

fprintf('\nRobustness over iota saved to data/robustness_iota.mat\n');

%% ===== Summary table =====
fprintf('\n%8s %10s %12s %12s %12s\n', 'iota', 'lambda*', 'logFFm-d', 'logDWm-d', 'logmaxDWm-d');
fprintf('%s\n', repmat('-', 1, 60));
for ii = 1:N_grid
    fprintf('%8.4f %10.4f %12.4f %12.4f %12.4f\n', ...
        rob_iota(ii).iota, rob_iota(ii).lambda, ...
        rob_iota(ii).logmean_FF_m - logmean_FF_d, ...
        rob_iota(ii).logmean_DW_m - logmean_DW_d, ...
        rob_iota(ii).logmax_DW_m  - logmax_DW_d);
end

%% ===== Local objective =====
function obj = obj_lambda(lam, eta_fix, iot, ploss_fix, in, ...
        idx_FF, idx_DW, logmean_FF_d, logmean_DW_d, logmax_DW_d)
    pen = 0;
    if lam <= 0.05, pen = pen + 1e6 * (1 + abs(0.05-lam)); end
    if lam >= 10,   pen = pen + 1e6 * (1 + abs(lam-10));   end
    if pen > 0, obj = pen; return; end

    try
        out = run_filter_series(lam, eta_fix, iot, ploss_fix, in);
    catch
        obj = 1e9; return;
    end

    FF = out.FF_us_t; DW = out.DW_us_t;
    if any(FF(idx_FF) <= 0) || any(DW(idx_DW) <= 0) || ...
       any(~isfinite(FF(idx_FF))) || any(~isfinite(DW(idx_DW)))
        obj = 1e8; return;
    end

    m1 = log(mean(FF(idx_FF))) - logmean_FF_d;
    m2 = log(mean(DW(idx_DW))) - logmean_DW_d;
    m3 = log(max (DW(idx_DW))) - logmax_DW_d;
    obj = m1^2 + m2^2 + m3^2;
end
