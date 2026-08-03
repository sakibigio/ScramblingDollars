addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1; pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat'); load('data/calibration.mat');
run('functions/params.m');
iota_us = iota_ss/freq/pi_us_ss; abs_sc = 12e4;
theta_plus_us = ((exp(lambda_us)-1)/(exp(lambda_us)+1))^2;
for tt = [50 100 150 200 250]
    mu_yt = exp(mu_us(tt));
    tgt = TED_s_us_t(tt) * abs_sc;
    smin = find_sigma_min(mu_yt, ploss_us, theta_plus_us, 0);
    sg = linspace(smin+1e-5, 10, 4000);
    r = arrayfun(@(s) Chi_p_psi(mu_yt, ploss_us, s, iota_us, lambda_us, eta, 1, 0)*abs_sc - tgt, sg);
    ok = isfinite(r);
    sc = find(diff(sign(r(ok))) ~= 0);
    sgok = sg(ok); rok = r(ok);
    fprintf('t=%3d  TED=%6.1f bps  sig_min=%.3f  roots at:', tt, tgt, smin);
    for k = sc, fprintf(' %.3f', sgok(k)); end
    fprintf('   (resid at smin+eps=%.1f, at 10=%.1f)\n', rok(1), rok(end));
end
