%% profile_lambda_eta50.m
%  Flatness diagnostic: profile the moment-fit objective in lambda at
%  eta = 0.5 (fixed), iota = 0.070 (feasibility floor), no penalty.
%  Reports each component so we can see WHICH moment identifies lambda.
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1; pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat'); load('data/calibration.mat');
run('functions/params.m');
varrho = 0; abs_scale = 12e4; use_mu_baseline = 0;
T = length(mu_us); mubar_us_t = zeros(T,1); mubar_eu_t = zeros(T,1);

valid_FF = ~isnan(FF_n) & FF_n > 0;  idx_FF = find(valid_FF);
valid_DW = ~isnan(DW_n) & DW_n > 0;  idx_DW = find(valid_DW);
logmean_FF_d = log(mean(FF_n(idx_FF)));
logmean_DW_d = log(mean(DW_n(idx_DW)));
NaN_bpus = (Rb_Rm == 0)    | isnan(Rb_Rm);
NaN_bpeu = (Rb_Rm_eu == 0) | isnan(Rb_Rm_eu);
idx_bpus = find(~NaN_bpus); idx_bpeu = find(~NaN_bpeu);
bpus_d_dm = Rb_Rm(idx_bpus)    - mean(Rb_Rm(idx_bpus));
bpeu_d_dm = Rb_Rm_eu(idx_bpeu) - mean(Rb_Rm_eu(idx_bpeu));
v_bpus = max(var(bpus_d_dm * abs_scale), 1e-12);
v_bpeu = max(var(bpeu_d_dm * abs_scale), 1e-12);

in = struct();
in.mu_us = mu_us; in.mu_eu = mu_eu; in.im_us = im_us; in.im_eu = im_eu;
in.TED_s_us_t = TED_s_us_t; in.TED_s_eu_t = TED_s_eu_t;
in.matching_type = matching_type; in.varrho = varrho; in.freq = freq;
in.pi_us_ss = pi_us_ss; in.pi_eu_ss = pi_eu_ss;
in.use_mu_baseline = use_mu_baseline;
in.mubar_us_t = mubar_us_t; in.mubar_eu_t = mubar_eu_t;
in.sigma_us_init = sigma_us; in.sigma_eu_init = sigma_eu;
in.abs_scale = abs_scale;

eta_fix = 0.5; iot = 0.070; ploss_fix = ploss_us;
lam_grid = [0.7:0.1:2.2, 2.4:0.2:3.6];

fprintf('\nlambda profile at eta=0.5, iota=0.070 (each term unit-free; total = sum)\n');
fprintf('%7s %9s %9s %9s %9s %9s %9s %9s\n', ...
    'lambda', 't_bpus', 't_bpeu', 't_FF', 't_DWS', 'total', 'meansig', 'FFlvl');
res = NaN(numel(lam_grid), 8);
for ii = 1:numel(lam_grid)
    lam = lam_grid(ii);
    try
        out = run_filter_series(lam, eta_fix, iot, ploss_fix, in);
        r_bpus = (out.BP_us_t(idx_bpus) - mean(out.BP_us_t(idx_bpus)) - bpus_d_dm) * abs_scale;
        r_bpeu = (out.BP_eu_t(idx_bpeu) - mean(out.BP_eu_t(idx_bpeu)) - bpeu_d_dm) * abs_scale;
        t_bpus = mean(r_bpus.^2) / v_bpus;
        t_bpeu = mean(r_bpeu.^2) / v_bpeu;
        t_ff   = (log(mean(out.FF_us_t(idx_FF)))  - logmean_FF_d)^2;
        t_dws  = (log(mean(out.DWS_us_t(idx_DW))) - logmean_DW_d)^2;
        tot    = t_bpus + t_bpeu + t_ff + t_dws;
        res(ii,:) = [lam t_bpus t_bpeu t_ff t_dws tot mean(out.sigma_us_t) mean(out.FF_us_t(idx_FF))];
        fprintf('%7.2f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', res(ii,:));
    catch ME
        fprintf('%7.2f   FAILED: %s\n', lam, ME.message);
    end
end
save('profile_lambda_eta50.mat', 'res', 'lam_grid');
