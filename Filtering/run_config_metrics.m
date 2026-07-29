%% run_config_metrics.m
% Run the filtering pipeline for ONE configuration and save its fit metrics
% (emphasis on QUANTITY fit) to output/metrics_<cfg>.mat.
%
% Invoke as:   matlab -batch "cfg='A'; run_config_metrics"
%   cfg = 'A'  baseline      : constant iota, ORIGINAL mu
%   cfg = 'B'  tv iota       : time-varying iota, ORIGINAL mu
%   cfg = 'C'  tv iota + new : time-varying iota, NEW mu (log LiqRatio H_18)

if ~exist('cfg','var'), error('set cfg = ''A''|''B''|''C'' before running'); end

TVDIR = '../Filtering_iota_tv';

% --- back up the live calibration override once ---
if exist('_calibration_override.mat','file') && ~exist('_calibration_override_CMPBK.mat','file')
    copyfile('_calibration_override.mat','_calibration_override_CMPBK.mat');
end

switch cfg
    case 'A'   % baseline: constant iota, original mu -> use the ORIGINAL override
        if exist('_calibration_override_CMPBK.mat','file')
            copyfile('_calibration_override_CMPBK.mat','_calibration_override.mat');
        end
        use_iota_tv = 0;  use_new_liq = 0;  iota_premium = NaN;
    case 'B'   % tv iota, original mu
        L = load(fullfile(TVDIR,'_calibration_tv_origliq.mat'));
        lambda_us_override_val = L.lambda_opt; lambda_eu_override_val = L.lambda_opt;
        eta_override_val = L.eta_opt; iota_ss_override_val = L.premium_opt;
        save('_calibration_override.mat','lambda_us_override_val', ...
             'lambda_eu_override_val','eta_override_val','iota_ss_override_val');
        use_iota_tv = 1;  use_new_liq = 0;  iota_premium = L.premium_opt;
    case 'C'   % tv iota, new mu
        L = load(fullfile(TVDIR,'_calibration_tv_newliq.mat'));
        lambda_us_override_val = L.lambda_opt; lambda_eu_override_val = L.lambda_opt;
        eta_override_val = L.eta_opt; iota_ss_override_val = L.premium_opt;
        save('_calibration_override.mat','lambda_us_override_val', ...
             'lambda_eu_override_val','eta_override_val','iota_ss_override_val');
        use_iota_tv = 1;  use_new_liq = 1;  iota_premium = L.premium_opt;
    otherwise
        error('unknown cfg %s', cfg);
end

matching_type   = 1;
cfg_tag         = cfg;
printit         = 0;   % never touch the Overleaf figures
run_julia       = 0;   % metrics only
do_regimes      = 0;
do_plot_baseline= 0;   % skip figure generation for speed

main_filter

%% ---- collect metrics (workspace was cleared; cfg_tag survived) ----
ok = @(x) x(~isnan(x));
m = struct();
m.cfg      = cfg_tag;
m.lambda   = lambda_us;   m.eta = eta;
m.iota_mean_ann = mean(iota_us_vec)*freq*pi_us_ss;

% --- QUANTITIES (the priority) ---
m.FF_model  = mean(FF_us_t)*100;              m.FF_data  = mean(ok(FF_n))*100;
m.DW_model  = mean(DW_us_t)*100;              m.DW_data  = mean(ok(DW_n))*100;
if exist('DWS_us_t','var'), m.DWS_model = mean(DWS_us_t)*100; else, m.DWS_model = NaN; end
m.DWS_data  = m.DW_data;                      % data DW_n is already a stock
m.DWFF_model = mean(DW_us_t./FF_us_t)*100;
m.DWpFF_model= mean(DW_us_t+FF_us_t)*100;
% like-with-like overlapping-sample log gaps
iv = find(~isnan(FF_n) & FF_n>0);  id = find(~isnan(DW_n) & DW_n>0);
m.loggap_FF  = log(mean(FF_us_t(iv))) - log(mean(FF_n(iv)));
if exist('DWS_us_t','var')
    m.loggap_DWS = log(mean(DWS_us_t(id))) - log(mean(DW_n(id)));
else
    m.loggap_DWS = NaN;
end

% --- PRICES ---
m.BP_corr  = corr(BP_us_t, Rb_Rm);   m.BP_std_m = std(BP_us_t)*abs_scale;  m.BP_std_d = std(Rb_Rm)*abs_scale;
m.TED_corr = corr(TED_us_t, TED_s_us_t);
m.CIP_corr = corr(CIP_t, cip);       m.CIP_std_m= std(CIP_t)*abs_scale;    m.CIP_std_d= std(cip)*abs_scale;

% --- SHOCK ---
vv = isfinite(sigma_us_t) & sigma_us_t>0 & isfinite(mu_us);
m.corr_sig_mu = corr(log(sigma_us_t(vv)), mu_us(vv));
m.sig_mean = mean(sigma_us_t);  m.sig_std = std(sigma_us_t);  m.sig_max = max(sigma_us_t);

if ~exist('output','dir'), mkdir('output'); end
save(fullfile('output',['metrics_' cfg_tag '.mat']), '-struct', 'm');
fprintf('\n[metrics] cfg %s: FF %.2f vs %.2f | DWS %.2f vs %.2f | corr(s,mu)=%.3f\n', ...
    m.cfg, m.FF_model, m.FF_data, m.DWS_model, m.DWS_data, m.corr_sig_mu);

%% restore the live override
if exist('_calibration_override_CMPBK.mat','file')
    copyfile('_calibration_override_CMPBK.mat','_calibration_override.mat');
end
