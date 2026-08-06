%% make_insample_cip_row.m
%  Official in-sample model CIP row for Table tab:moddatacomp.
%  Series: filter-implied CIP along the realized path at the promoted
%  calibration (LCR-excess mu). Conventions IDENTICAL to
%  compute_data_moments.m: annualization ((1+r)^12-1)*1e4, regime masks from
%  live MS_sigma_us_prob.csv (prob row p <-> data row p+1; first obs normal),
%  within-regime autocorrelations on concatenated subsets.
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1; pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat'); load('data/calibration.mat');
run('functions/params.m');
T = length(mu_us);

in = struct('mu_us',mu_minus_lcr_level(mu_us,LCR_us),'mu_eu',mu_minus_lcr_level(mu_eu,LCR_us),'im_us',im_us,'im_eu',im_eu, ...
    'TED_s_us_t',TED_s_us_t,'TED_s_eu_t',TED_s_eu_t,'matching_type',1,'varrho',0,'freq',freq, ...
    'pi_us_ss',1,'pi_eu_ss',1,'use_mu_baseline',0,'mubar_us_t',zeros(T,1),'mubar_eu_t',zeros(T,1), ...
    'sigma_us_init',sigma_us,'sigma_eu_init',sigma_eu,'abs_scale',12e4);
out = run_filter_series(lambda_us, eta, iota_ss, ploss_us, in);
Rm_us_v = exp(im_us); Rm_eu_v = exp(im_eu);
Rb_us_t = Rm_us_v + out.BP_us_t;  Rb_eu_t = Rm_eu_v + out.BP_eu_t;
cip_pp  = (Rm_eu_v - Rm_us_v) + Rm_eu_v .* (Rb_us_t ./ Rb_eu_t - 1);   % per-period decimal
CIP_m   = ((1 + cip_pp).^12 - 1) * 1e4;                                % annualized bps

prob_tbl = readtable('data/MS_sigma_us_prob.csv');
in_scr = false(T,1); in_scr(2:end) = prob_tbl.prob_scr > 0.5;
in_r1 = ~in_scr; in_r1(1) = true;

x_v = CIP_m(isfinite(CIP_m));
ac  = @(x) sum((x(1:end-1)-mean(x)).*(x(2:end)-mean(x))) / sum((x-mean(x)).^2);
x1 = CIP_m(in_r1); x2 = CIP_m(in_scr);
row = [mean(x_v), ac(x_v), std(x_v), mean(x2)-mean(x1), ac(x2), ac(x1), std(x2)/std(x1)];

ovl = overleaf_dir();
fid = fopen(fullfile(ovl, 'Mod_CIP_Moments_insample.tex'), 'wt');
fprintf(fid, '$\\mathcal{CIP}$ (model, in-sample) & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', row);
fclose(fid);
fprintf('in-sample row: mean %.1f | AC %.2f | std %.1f | diff %.1f | AC{scr,nor} {%.2f,%.2f} | rel %.1f\n', row);
fprintf('Wrote %s\n', fullfile(ovl, 'Mod_CIP_Moments_insample.tex'));
