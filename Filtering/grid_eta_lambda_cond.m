%% Conditional log-sigma separation + moment fit over (eta, lambda), iota at floor.
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1; pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat'); load('data/calibration.mat');
run('functions/params.m');
varrho = 0; abs_scale = 12e4; use_mu_baseline = 0;
T = length(mu_us);
in = struct('mu_us',mu_us,'mu_eu',mu_eu,'im_us',im_us,'im_eu',im_eu, ...
    'TED_s_us_t',TED_s_us_t,'TED_s_eu_t',TED_s_eu_t,'matching_type',matching_type, ...
    'varrho',varrho,'freq',freq,'pi_us_ss',pi_us_ss,'pi_eu_ss',pi_eu_ss, ...
    'use_mu_baseline',use_mu_baseline,'mubar_us_t',zeros(T,1),'mubar_eu_t',zeros(T,1), ...
    'sigma_us_init',sigma_us,'sigma_eu_init',sigma_eu,'abs_scale',abs_scale);

prob_tbl = readtable('data/MS_sigma_us_prob_cbase.csv');
scr_mask = false(T,1); scr_mask(1 + find(prob_tbl.prob_scr > 0.5)) = true;

valid_FF = ~isnan(FF_n) & FF_n > 0;  idx_FF = find(valid_FF);
valid_DW = ~isnan(DW_n) & DW_n > 0;  idx_DW = find(valid_DW);
logmean_FF_d = log(mean(FF_n(idx_FF)));
logmean_DW_d = log(mean(DW_n(idx_DW)));
max_ted_bps = max(max(TED_s_us_t), max(TED_s_eu_t)) * abs_scale;

eta_grid = [0.45 0.50 0.55 0.611 0.65 0.70];
lam_grid = [1.00 1.25 1.50 1.85 2.20];
fprintf('\n%6s %6s %7s | %9s %8s %8s\n', 'eta', 'iota', 'lambda', 'conddiff', 'FFlvl', 'DWSgap');
for ie = 1:numel(eta_grid)
    et = eta_grid(ie);
    iot = 1.05 * max_ted_bps / ((1 - et) * 1e4);   % feasibility floor
    for il = 1:numel(lam_grid)
        lam = lam_grid(il);
        try
            out = run_filter_series(lam, et, iot, ploss_us, in);
            ls = log(out.sigma_us_t);
            cd_ = mean(ls(scr_mask)) - mean(ls(~scr_mask));
            ffl = mean(out.FF_us_t(idx_FF));
            dwg = log(mean(out.DWS_us_t(idx_DW))) - logmean_DW_d;
            fprintf('%6.3f %6.4f %7.2f | %+9.3f %8.3f %+8.3f\n', et, iot, lam, cd_, ffl, dwg);
        catch
            fprintf('%6.3f %6.4f %7.2f |   FAILED\n', et, iot, lam);
        end
    end
end
