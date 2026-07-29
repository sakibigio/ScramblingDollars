%% run_isolate_tv.m
% Does B's quantity-fit gain come from the TIME-VARIATION in iota, or just from
% re-estimating lambda/eta?  Runs a CONSTANT-iota model at B's iota level with
% B's lambda/eta, everything else identical to B.
%   cfg 'A2' : constant iota = mean(B's iota_t), lambda/eta = B's, original mu
% Compare against A (paper params, constant iota) and B (time-varying iota).

clear; close all;
L = load('../Filtering_iota_tv/_calibration_iota_tv.mat');  % flat-iota control
% B's mean iota over the sample (annualised decimal)
S  = load('data/iota_corridor_monthly.mat');
sp = S.iota_sprd_dec; cv = S.cover_mask;
f0=find(cv,1); l0=find(cv,1,'last'); sp(1:f0-1)=sp(f0); sp(l0+1:end)=sp(l0);
iota_mean = mean(sp)*0 + mean(sp) + L.premium_opt;  % constant
fprintf('B: lambda=%.4f eta=%.4f  premium=%.5f  -> mean iota_t = %.5f (%.0f bps)\n', ...
    L.lambda_opt, L.eta_opt, L.premium_opt, iota_mean, iota_mean*1e4);

if exist('_calibration_override.mat','file') && ~exist('_calibration_override_ISOBK.mat','file')
    copyfile('_calibration_override.mat','_calibration_override_ISOBK.mat');
end
lambda_us_override_val = L.lambda_opt;
lambda_eu_override_val = L.lambda_opt;
eta_override_val       = L.eta_opt;
iota_ss_override_val   = iota_mean;     % CONSTANT iota at B's mean level
save('_calibration_override.mat','lambda_us_override_val','lambda_eu_override_val', ...
     'eta_override_val','iota_ss_override_val');

matching_type   = 1;
use_iota_tv     = 0;      % <-- CONSTANT iota (the whole point of this test)
use_new_liq     = 0;      % original mu, same as B
iota_premium    = NaN;
printit         = 0; run_julia = 0; do_regimes = 0; do_plot_baseline = 0;
cfg_tag         = 'FLAT';

main_filter

ok = @(x) x(~isnan(x));
m = struct('cfg','FLAT');
m.lambda=lambda_us; m.eta=eta; m.iota_mean_ann=mean(iota_us_vec)*freq*pi_us_ss;
m.FF_model=mean(FF_us_t)*100;  m.FF_data=mean(ok(FF_n))*100;
if exist('DWS_us_t','var'), m.DWS_model=mean(DWS_us_t)*100; else, m.DWS_model=NaN; end
m.DWS_data=mean(ok(DW_n))*100;
iv=find(~isnan(FF_n)&FF_n>0); id=find(~isnan(DW_n)&DW_n>0);
m.loggap_FF=log(mean(FF_us_t(iv)))-log(mean(FF_n(iv)));
if exist('DWS_us_t','var'), m.loggap_DWS=log(mean(DWS_us_t(id)))-log(mean(DW_n(id))); else, m.loggap_DWS=NaN; end
m.BP_corr=corr(BP_us_t,Rb_Rm); m.BP_std_m=std(BP_us_t)*abs_scale; m.BP_std_d=std(Rb_Rm)*abs_scale;
m.CIP_corr=corr(CIP_t,cip);
vv=isfinite(sigma_us_t)&sigma_us_t>0&isfinite(mu_us);
m.corr_sig_mu=corr(log(sigma_us_t(vv)),mu_us(vv));
m.sig_mean=mean(sigma_us_t); m.sig_std=std(sigma_us_t);
if ~exist('output','dir'), mkdir('output'); end
save(fullfile('output','metrics_FLAT.mat'),'-struct','m');
fprintf('\n[FLAT] constant iota=%.0f bps, B lambda/eta: FF %.2f (data %.2f) | DWS %.2f (data %.2f) | |loggap| %.3f\n', ...
    m.iota_mean_ann*1e4, m.FF_model, m.FF_data, m.DWS_model, m.DWS_data, abs(m.loggap_FF)+abs(m.loggap_DWS));

if exist('_calibration_override_ISOBK.mat','file')
    copyfile('_calibration_override_ISOBK.mat','_calibration_override.mat');
    delete('_calibration_override_ISOBK.mat');
end
