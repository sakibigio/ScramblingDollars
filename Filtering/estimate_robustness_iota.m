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
% Feasibility floor: the inversion requires max(TED) <= (1-eta)*iota, i.e.
% iota >= 0.0350/(1-0.632) = 0.0951 at the baseline eta -- which is exactly
% where the estimate sits (the bound binds). The grid therefore starts at the
% baseline and explores upward; points below the floor are infeasible by
% construction. 11 uniform points; the first equals the baseline iota_ss.
iota_grid = linspace(0.0951, 0.1401, 11);
N_grid    = numel(iota_grid);

%% Sample selection for moment objective
valid_FF = ~isnan(FF_n) & FF_n > 0;
valid_DW = ~isnan(DW_n) & DW_n > 0;
idx_FF   = find(valid_FF);
idx_DW   = find(valid_DW);

logmean_FF_d = log(mean(FF_n(idx_FF)));
logmean_DW_d = log(mean(DW_n(idx_DW)));
logmax_DW_d  = log(max(DW_n(idx_DW)));

%% Baseline-criterion targets (mirror estimate_params_4d.m; w_cip = 0 there)
NaN_bpus = (Rb_Rm == 0)    | isnan(Rb_Rm);
NaN_bpeu = (Rb_Rm_eu == 0) | isnan(Rb_Rm_eu);
tg = struct();
tg.idx_bpus  = find(~NaN_bpus);
tg.idx_bpeu  = find(~NaN_bpeu);
tg.bpus_d_dm = Rb_Rm(tg.idx_bpus)    - mean(Rb_Rm(tg.idx_bpus));
tg.bpeu_d_dm = Rb_Rm_eu(tg.idx_bpeu) - mean(Rb_Rm_eu(tg.idx_bpeu));
tg.v_bpus    = max(var(tg.bpus_d_dm * abs_scale), 1e-12);
tg.v_bpeu    = max(var(tg.bpeu_d_dm * abs_scale), 1e-12);
tg.idx_ff    = idx_FF;
tg.idx_dw    = idx_DW;
tg.logmean_FF_d = logmean_FF_d;
tg.logmean_DW_d = logmean_DW_d;   % data DW_n is a stock; compared to model DWS
tg.w_bpus = 1.0; tg.w_bpeu = 1.0; tg.w_meanFF = 1.0; tg.w_meanDWS = 1.0;
tg.w_corr = 0.0;   % baseline criterion is penalty-free (nopen_lcr spec, 2026-07-31)

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

    % 1. Estimate lambda by fminsearch (eta fixed at baseline), using the
    %    SAME criterion as the baseline estimation (estimate_params_4d.m)
    %    restricted to lambda.  Multi-start: baseline lambda and, past the
    %    first grid point, the previous grid point's optimum.
    obj = @(lam) obj_lambda_baseline(lam, eta_baseline, iota_g, ploss_baseline, in, tg);
    x0_list = lambda_baseline;
    if ii > 1 && isfinite(rob_iota(ii-1).lambda)
        x0_list = [lambda_baseline, rob_iota(ii-1).lambda];
    end
    fval = Inf; lam_opt = NaN; ef = NaN;
    for jj = 1:numel(x0_list)
        [lam_try, fval_try, ef_try] = fminsearch(obj, x0_list(jj), options);
        if fval_try < fval
            fval = fval_try; lam_opt = lam_try; ef = ef_try;
        end
    end

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

% Objective: functions/obj_lambda_baseline.m (baseline criterion, lambda only).
