%% Generate filtered sigma series under eta=0.5, iota=0.070 for candidate lambdas.
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
for lam = [1.4, 1.6, 1.8]
    out = run_filter_series(lam, 0.5, 0.070, ploss_us, in);
    tbl = table(out.sigma_us_t, out.sigma_eu_t, 'VariableNames', {'sigma_us','sigma_eu'});
    fn = sprintf('RW_shock_eta50_lam%03.0f.csv', lam*100);
    writetable(tbl, fn);
    fprintf('%s: mean sig_us=%.4f max=%.2f\n', fn, mean(out.sigma_us_t), max(out.sigma_us_t));
end
