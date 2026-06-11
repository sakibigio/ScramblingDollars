%% estimate_robustness_ploss.m
%  Robustness check for referee 3 (risk aversion):
%  Vary ploss over a grid; for each value, re-estimate the single
%  best-fitting lambda (eta and iota_ss fixed at baseline), run the full
%  filter, then re-estimate Markov regimes via markov_estimation.jl for
%  that pair.
%
%  Requires the same workspace inputs as estimate_robustness_iota.m
%  Output: data/robustness_ploss.mat with struct rob_ploss
%  Per-pair Julia outputs copied to data/robustness/ploss/g{ii}/.

fprintf('\n========================================\n');
fprintf('Robustness exercise: grid over ploss\n');
fprintf('Free param: lambda  |  eta=%.4f, iota_ss=%.4f fixed\n', ...
    eta_baseline, iota_baseline);
fprintf('========================================\n');

%% Grid
% 10 uniform points spanning the empirically plausible range.
ploss_grid = linspace(0.30, 0.70, 10);
N_grid     = numel(ploss_grid);

%% Sample selection for moment objective
valid_FF = ~isnan(FF_n) & FF_n > 0;
valid_DW = ~isnan(DW_n) & DW_n > 0;
idx_FF   = find(valid_FF);
idx_DW   = find(valid_DW);

logmean_FF_d = log(mean(FF_n(idx_FF)));
logmean_DW_d = log(mean(DW_n(idx_DW)));
logmax_DW_d  = log(max(DW_n(idx_DW)));

%% Inputs struct
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

%% Output folder
rob_dir = fullfile('data', 'robustness', 'ploss');
if ~exist(rob_dir, 'dir'), mkdir(rob_dir); end

rob_ploss = repmat(struct('ploss', NaN, 'lambda', NaN, 'eta', NaN, ...
    'iota', NaN, 'fval', NaN, 'exitflag', NaN, 'series', [], ...
    'prob_nor_us', [], 'logmean_FF_m', NaN, 'logmean_DW_m', NaN, ...
    'logmax_DW_m', NaN), N_grid, 1);

options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 1000, ...
                   'TolX', 1e-6, 'TolFun', 1e-10);

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
    p_g = ploss_grid(ii);
    fprintf('\n--- (%d/%d) ploss = %.3f ---\n', ii, N_grid, p_g);

    obj = @(lam) obj_lambda(lam, eta_baseline, iota_baseline, p_g, in, ...
        idx_FF, idx_DW, logmean_FF_d, logmean_DW_d, logmax_DW_d);
    x0 = lambda_baseline;
    [lam_opt, fval, ef] = fminsearch(obj, x0, options);

    fprintf('  lambda* = %.6f   obj = %.4e   ef = %d\n', lam_opt, fval, ef);

    series = run_filter_series(lam_opt, eta_baseline, iota_baseline, p_g, in);

    prob_nor = [];
    if ~isempty(julia_path)
        T_use = min(length(series.sigma_us_t), length(series.sigma_eu_t));
        sigma_mat = [series.sigma_us_t(1:T_use) series.sigma_eu_t(1:T_use)];
        rw_tbl = array2table(sigma_mat, 'VariableNames', {'sigma_us', 'sigma_eu'});
        writetable(rw_tbl, baseline_RWshock);

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
            warning('  Julia failed for ploss=%.3f: %s', p_g, jresult);
        end

        pair_dir = fullfile(rob_dir, sprintf('g%02d', ii));
        if ~exist(pair_dir, 'dir'), mkdir(pair_dir); end
        for f = {'MS_sigma_us_prob.csv', 'MS_sigma_us_params.csv', ...
                 'MS_sigma_eu_prob.csv', 'MS_params.csv'}
            src = fullfile('data', f{1});
            if exist(src, 'file'), copyfile(src, fullfile(pair_dir, f{1})); end
        end
    end

    rob_ploss(ii).ploss    = p_g;
    rob_ploss(ii).lambda   = lam_opt;
    rob_ploss(ii).eta      = eta_baseline;
    rob_ploss(ii).iota     = iota_baseline;
    rob_ploss(ii).fval     = fval;
    rob_ploss(ii).exitflag = ef;
    rob_ploss(ii).series   = series;
    rob_ploss(ii).prob_nor_us = prob_nor;
    rob_ploss(ii).logmean_FF_m = log(mean(series.FF_us_t(idx_FF)));
    rob_ploss(ii).logmean_DW_m = log(mean(series.DW_us_t(idx_DW)));
    rob_ploss(ii).logmax_DW_m  = log(max(series.DW_us_t(idx_DW)));

    fprintf('  log mean FF: model=%.4f  data=%.4f\n', ...
        rob_ploss(ii).logmean_FF_m, logmean_FF_d);
    fprintf('  log mean DW: model=%.4f  data=%.4f\n', ...
        rob_ploss(ii).logmean_DW_m, logmean_DW_d);
end

if exist(baseline_prob_back, 'file')
    copyfile(baseline_prob_back, baseline_prob_csv);
end
if exist(baseline_RW_back, 'file')
    copyfile(baseline_RW_back, baseline_RWshock);
end

save('data/robustness_ploss.mat', 'rob_ploss', 'ploss_grid', ...
    'lambda_baseline', 'eta_baseline', 'iota_baseline', 'ploss_baseline', ...
    'logmean_FF_d', 'logmean_DW_d', 'logmax_DW_d', '-v7.3');

fprintf('\nRobustness over ploss saved to data/robustness_ploss.mat\n');

fprintf('\n%8s %10s %12s %12s %12s\n', 'ploss', 'lambda*', 'logFFm-d', 'logDWm-d', 'logmaxDWm-d');
fprintf('%s\n', repmat('-', 1, 60));
for ii = 1:N_grid
    fprintf('%8.3f %10.4f %12.4f %12.4f %12.4f\n', ...
        rob_ploss(ii).ploss, rob_ploss(ii).lambda, ...
        rob_ploss(ii).logmean_FF_m - logmean_FF_d, ...
        rob_ploss(ii).logmean_DW_m - logmean_DW_d, ...
        rob_ploss(ii).logmax_DW_m  - logmax_DW_d);
end

%% ===== Local objective =====
function obj = obj_lambda(lam, eta_fix, iot_fix, ploss_g, in, ...
        idx_FF, idx_DW, logmean_FF_d, logmean_DW_d, logmax_DW_d)
    pen = 0;
    if lam <= 0.05, pen = pen + 1e6 * (1 + abs(0.05-lam)); end
    if lam >= 10,   pen = pen + 1e6 * (1 + abs(lam-10));   end
    if pen > 0, obj = pen; return; end

    try
        out = run_filter_series(lam, eta_fix, iot_fix, ploss_g, in);
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
