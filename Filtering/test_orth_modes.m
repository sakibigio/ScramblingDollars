%% Smoke test: run the filter once at the committed calibration and report
%  all three identification-penalty statistics on the filtered log sigma_us,
%  plus wall-clock time per filter run (to budget the estimation).
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1;
pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat');
load('data/calibration.mat');
run('functions/params.m');

fprintf('Committed calibration: lambda_us=%.4f lambda_eu=%.4f iota_ss=%.4f eta=%.4f\n', ...
    lambda_us, lambda_eu, iota_ss, eta);

T_obs = length(mu_us);
mubar_us_t = zeros(T_obs, 1);
mubar_eu_t = zeros(T_obs, 1);
Rm_us_vec = exp(im_us);
Rm_eu_vec = exp(im_eu);
abs_sc = 12e4;

iota_us_loc = iota_ss / freq / pi_us_ss;
theta_plus_us = ((exp(lambda_us) - 1) / (exp(lambda_us) + 1))^2;
theta_plus_eu = ((exp(lambda_eu) - 1) / (exp(lambda_eu) + 1))^2;

fsolve_opt = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12);
sig_us_guess = sigma_us; sig_eu_guess = sigma_eu;
sigma_us_t = zeros(T_obs, 1);

tic;
for tt = 1:T_obs
    mu_us_yt = exp(mu_us(tt)) - mubar_us_t(tt);
    TED_us_tgt = TED_s_us_t(tt) * abs_sc;
    sig_min_us = find_sigma_min(mu_us_yt, ploss_us, theta_plus_us, 0);
    res_fun = @(s) Chi_p_psi(mu_us_yt, ploss_us, s, iota_ss/freq/pi_us_ss, lambda_us, eta, 1, 0) * abs_sc - TED_us_tgt;
    sig_try = max(sig_us_guess, sig_min_us + 1e-3);
    [sig_us, ~, ef] = fsolve(res_fun, sig_try, fsolve_opt);
    if ef <= 0 || sig_us < sig_min_us
        sig_us = fminbnd(@(s) res_fun(s)^2, sig_min_us + 1e-6, 100);
    end
    sig_us_guess = sig_us;
    sigma_us_t(tt) = sig_us;
end
t_filter = toc;
fprintf('One US filter pass over %d periods: %.2f s (full 2-region obj eval ~2x)\n', T_obs, t_filter);

valid = isfinite(sigma_us_t) & sigma_us_t > 0 & isfinite(mu_us);
lsig = log(sigma_us_t(valid)); lsig = lsig(:);
tvec = (1:numel(lsig))';
e = lsig - mean(lsig);

corr_mu   = corr(lsig, mu_us(valid));
corr_trnd = corr(lsig, tvec);
kpss_raw  = sum(cumsum(e).^2) / (numel(e)^2 * var(e));

% LRV (Bartlett l12) version, matching estimate_params_4d.m
T = numel(e);
S = cumsum(e);
l = floor(12 * (T/100)^(1/4));
s2 = (e' * e) / T;
for j = 1:l
    w  = 1 - j/(l+1);
    s2 = s2 + 2 * w * (e(1:T-j)' * e(1+j:T)) / T;
end
kpss_lrv = sum(S.^2) / (T^2 * s2);

fprintf('\nAt the committed calibration (eta=%.3f):\n', eta);
fprintf('  corr(log sigma, mu)   = %+.4f   -> corr penalty  = %.4f\n', corr_mu, corr_mu^2);
fprintf('  corr(log sigma, t)    = %+.4f   -> trend penalty = %.4f\n', corr_trnd, corr_trnd^2);
fprintf('  KPSS raw-variance     = %.4f\n', kpss_raw);
fprintf('  KPSS Bartlett-LRV l=%d = %.4f   (5%% crit = 0.463)\n', l, kpss_lrv);
