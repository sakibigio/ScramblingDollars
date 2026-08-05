%% data_moments_check.m - data counterparts for the nonlinear-model moments
addpath('functions'); addpath('functions/chi'); addpath('data');
S = load('data/LFX_data.mat');
bps = 12e4;

ok = @(x) x(isfinite(x) & x~=0);
ac = @(x) corr_complete(x(1:end-1), x(2:end));

bp  = S.Rb_Rm;  bpv = ok(bp);
cip = S.cip;    cipv = ok(cip);

fprintf('\n===== DATA MOMENTS =====\n');
fprintf('%-10s %10s %10s %10s\n','series','mean(bps)','std(bps)','autocorr');
fprintf('%s\n', repmat('-',1,44));
fprintf('%-10s %10.2f %10.2f %10.3f\n','BP',  mean(bpv)*bps,  std(bpv)*bps,  ac(bpv));
fprintf('%-10s %10.2f %10.2f %10.3f\n','CIP', mean(cipv)*bps, std(cipv)*bps, ac(cipv));

% Exchange rate (EUR/USD log level in LFX_data)
if isfield(S,'ln_eu_us_t')
    e  = exp(S.ln_eu_us_t);
    de = diff(log(e));
    fprintf('%-10s %10.3f %10.3f %10.3f   (level; std/mean = %.4f)\n','FX', ...
        mean(e), std(e), ac(e), std(e)/mean(e));
    fprintf('%-10s %10.2f %10.2f %10.3f   (monthly log change, bps)\n','dFX', ...
        mean(de)*1e4, std(de)*1e4, ac(de));
end

% liquidity ratio actually seen by the filter (H_18 swap) vs raw
Lq = load('data/liq_ratio_monthly.mat');
mn = Lq.mu_us_new; cv = Lq.cover_mask;
f0=find(cv,1); l0=find(cv,1,'last'); mn(1:f0-1)=mn(f0); mn(l0+1:end)=mn(l0);
fprintf('\nliquidity ratio exp(mu_us): H_18 mean=%.4f [%.4f, %.4f]\n', ...
    mean(exp(mn)), min(exp(mn)), max(exp(mn)));
fprintf('                    raw LFX mean=%.4f\n', mean(exp(S.mu_us)));

% What BP does the STATIC chi mapping give at the data's mu and a given sigma?
% (this is exactly what the filter evaluates, so it isolates mu vs sigma)
load('data/calibration.mat');
matching_type = 1; pi_us_ss = 1; pi_eu_ss = 1;
run('functions/params.m');
lam = 1.129; et = 0.662; iot_ann = 0.0855; iot = iot_ann/12;
fprintf('\n--- static BP = Echi_m(mu, ploss, sigma, iota) in bps ---\n');
fprintf('%-22s %10s %10s %10s\n','mu (level)','sig=0.25','sig=0.31','sig=0.40');
for muv = [0.10 0.20 0.30 0.3627 0.50]
    v = zeros(1,3); k=0;
    for sg = [0.25 0.31 0.40]
        k=k+1;
        v(k) = Echi_m(muv, ploss_us, sg, iot, lam, et, 1, 0)*bps;
    end
    lbl = sprintf('%.4f', muv);
    if abs(muv-0.3627)<1e-6, lbl = [lbl ' (H_18 mean)']; end
    fprintf('%-22s %10.2f %10.2f %10.2f\n', lbl, v(1), v(2), v(3));
end
