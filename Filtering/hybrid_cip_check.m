%% hybrid_cip_check.m
%  Audit the filter's random-walk assumption with the nonlinear model:
%  compare data CIP against (a) the filter's CIP (RW: E[dlog e]=0) and
%  (b) a hybrid that adds the MODEL's state-contingent expected devaluation,
%      CIP_hybrid_t = CIP_filter_t + E[dlog e | sigma_t, r_t]   (annualized bps),
%  where expected dollar appreciation (e = EUR per USD rising) adds to the
%  measured dollar wedge. If the model's FX is near-RW, (a) ~ (b).
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1; pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat'); load('data/calibration.mat');
run('functions/params.m');
T = length(mu_us);

% --- model expected devaluation by state ---
G = load('data/global_solcalibration_dynare_sigma_us.mat', 'Q_mat', 'e_euus_vec', 'sigma_us_vec');
ev = G.e_euus_vec(:); sv = G.sigma_us_vec(:);
dev_s = log((G.Q_mat * ev) ./ ev);            % E-approx of dlog e per month, by state
dev_bps = dev_s * 12e4;                        % annualized bps
i1 = 1:76; i2 = 77:152;

% --- empirical state path (promoted live files) ---
RW = readtable('RW_shock.csv');  sig = RW.sigma_us(:);
PR = readtable('data/MS_sigma_us_prob.csv');
p  = [PR.prob_scr(1); PR.prob_scr(:)];
clamp = @(x,g) min(max(x, min(g)), max(g));
dev_nor = interp1(sv(i1), dev_bps(i1), clamp(sig, sv(i1)));
dev_scr = interp1(sv(i2), dev_bps(i2), clamp(sig, sv(i2)));
dev_t = (1 - p).*dev_nor + p.*dev_scr;

% --- filter CIP at the promoted calibration (LCR mu) ---
in = struct('mu_us',mu_minus_lcr_level(mu_us,LCR_us),'mu_eu',mu_eu,'im_us',im_us,'im_eu',im_eu, ...
    'TED_s_us_t',TED_s_us_t,'TED_s_eu_t',TED_s_eu_t,'matching_type',1,'varrho',0,'freq',freq, ...
    'pi_us_ss',1,'pi_eu_ss',1,'use_mu_baseline',0,'mubar_us_t',zeros(T,1),'mubar_eu_t',zeros(T,1), ...
    'sigma_us_init',sigma_us,'sigma_eu_init',sigma_eu,'abs_scale',12e4);
out = run_filter_series(lambda_us, eta, iota_ss, ploss_us, in);
Rm_us_v = exp(im_us); Rm_eu_v = exp(im_eu);
Rb_us_t = Rm_us_v + out.BP_us_t;  Rb_eu_t = Rm_eu_v + out.BP_eu_t;
CIP_f = ((Rm_eu_v - Rm_us_v) + Rm_eu_v .* (Rb_us_t ./ Rb_eu_t - 1)) * 12e4;   % bps

CIP_h = CIP_f + dev_t;

% --- comparison vs data ---
cip_d = cip * 12e4;  ok = isfinite(cip_d) & cip_d ~= 0;
scr = p > 0.5;
fprintf('--- Model expected devaluation along the path (annualized bps) ---\n');
fprintf('mean = %+.2f | std = %.2f | range [%+.1f, %+.1f]\n', mean(dev_t), std(dev_t), min(dev_t), max(dev_t));
fprintf('normal months: mean %+.2f | scrambling months: mean %+.2f\n', mean(dev_t(~scr)), mean(dev_t(scr)));
fprintf('\n--- Fit to data CIP (%d months) ---\n', sum(ok));
fprintf('%-28s corr = %.3f   RMSE = %.2f bps\n', 'filter CIP (RW):', corr(CIP_f(ok), cip_d(ok)), sqrt(mean((CIP_f(ok)-cip_d(ok)).^2)));
fprintf('%-28s corr = %.3f   RMSE = %.2f bps\n', 'hybrid (+model E[dlog e]):', corr(CIP_h(ok), cip_d(ok)), sqrt(mean((CIP_h(ok)-cip_d(ok)).^2)));
fprintf('max |hybrid - filter| = %.2f bps at t=%d\n', max(abs(CIP_h - CIP_f)), find(abs(CIP_h-CIP_f)==max(abs(CIP_h-CIP_f)),1));

save('hybrid_cip_check.mat', 'CIP_f', 'CIP_h', 'dev_t', 'cip_d', 'p');
