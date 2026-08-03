%% gen_est_series.m
%  Generate the filtered sigma series for each estimation version, using the
%  SAME mu the estimation itself used (coherent estimation -> series pairs):
%    raw mu : corr, trend, kpss        (results _estimation_result_<m>.mat)
%    LCR mu : corr, trend, kpss        (results _estimation_result_<m>_lcr.mat)
%  Output: RW_shock_est_<m>[_lcr].csv  (sigma_us, sigma_eu)
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1; pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat'); load('data/calibration.mat');
run('functions/params.m');
varrho = 0; abs_scale = 12e4; use_mu_baseline = 0;
T = length(mu_us);
mu_us_raw = mu_us;
mu_us_lcr = mu_minus_lcr_level(mu_us, LCR_us);

base_in = struct('mu_eu',mu_eu,'im_us',im_us,'im_eu',im_eu, ...
    'TED_s_us_t',TED_s_us_t,'TED_s_eu_t',TED_s_eu_t,'matching_type',matching_type, ...
    'varrho',varrho,'freq',freq,'pi_us_ss',pi_us_ss,'pi_eu_ss',pi_eu_ss, ...
    'use_mu_baseline',use_mu_baseline,'mubar_us_t',zeros(T,1),'mubar_eu_t',zeros(T,1), ...
    'sigma_us_init',sigma_us,'sigma_eu_init',sigma_eu,'abs_scale',abs_scale);

modes = {'corr','trend','kpss'};
sfxs  = {'', '_lcr'};
for is = 1:2
    for im_ = 1:3
        sfx = sfxs{is}; md = modes{im_};
        rf = sprintf('_estimation_result_%s%s.mat', md, sfx);
        if ~exist(rf, 'file')
            fprintf('SKIP %s%s: %s missing\n', md, sfx, rf);
            continue;
        end
        R = load(rf);
        in = base_in;
        if is == 2, in.mu_us = mu_us_lcr; else, in.mu_us = mu_us_raw; end
        out = run_filter_series(R.lambda_opt, R.eta_opt, R.iota_ss_opt, ploss_us, in);
        tbl = table(out.sigma_us_t, out.sigma_eu_t, 'VariableNames', {'sigma_us','sigma_eu'});
        fn = sprintf('RW_shock_est_%s%s.csv', md, sfx);
        writetable(tbl, fn);
        fprintf('%s: lam=%.3f iota=%.4f eta=%.4f | mean sig=%.4f max=%.2f\n', ...
            fn, R.lambda_opt, R.iota_ss_opt, R.eta_opt, mean(out.sigma_us_t), max(out.sigma_us_t));
    end
end
