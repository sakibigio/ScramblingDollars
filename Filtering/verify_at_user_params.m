%% Verify: at the user's reported params, does the estimation routine's
%  filter give the same DWS/FF as main_filter.m?  Same data, same setup,
%  no optimization -- single function evaluation.
cd(fileparts(mfilename('fullpath')));   % this file's own folder
addpath('functions'); addpath('functions/chi'); addpath('utils'); addpath('data');

matching_type = 1;
pi_us_ss = 1;
pi_eu_ss = 1;

load('data/LFX_data.mat');
load('data/calibration.mat');
run('functions/params.m');

% Override with user's specific params
lambda_us = 1.5;
lambda_eu = 1.5;
iota_ss   = 0.0625;
eta       = 0.5;
varrho    = 0;
ploss_us  = 0.5;
ploss_eu  = 0.5;
iota_us = iota_ss / 12 / pi_us_ss;
iota_eu = iota_ss / 12 / pi_eu_ss;

sigma_us = 0.20;
sigma_eu = 0.15;

T = length(mu_us);
mubar_us_t = zeros(T,1);     % matches main_filter use_mu_baseline ~= 1
mubar_eu_t = zeros(T,1);

abs_sc = 12e4;
fsolve_opt = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12);

DW_us_t = zeros(T,1);  FF_us_t = zeros(T,1);
DW_eu_t = zeros(T,1);  FF_eu_t = zeros(T,1);
sigma_us_t = zeros(T,1); sigma_eu_t = zeros(T,1);

theta_plus_us = ((exp(lambda_us)-1)/(exp(lambda_us)+1))^2;
theta_plus_eu = ((exp(lambda_eu)-1)/(exp(lambda_eu)+1))^2;

sig_us = sigma_us; sig_eu = sigma_eu;
for tt = 1:T
    mu_us_yt = exp(mu_us(tt));
    mu_eu_yt = exp(mu_eu(tt));
    TED_us_tgt = TED_s_us_t(tt) * abs_sc;
    TED_eu_tgt = TED_s_eu_t(tt) * abs_sc;

    sig_min_us = find_sigma_min(mu_us_yt, ploss_us, theta_plus_us, varrho);
    sig_min_eu = find_sigma_min(mu_eu_yt, ploss_eu, theta_plus_eu, varrho);

    res_us = @(s) Chi_p_psi(mu_us_yt, ploss_us, s, iota_us, lambda_us, eta, 1, varrho)*abs_sc - TED_us_tgt;
    res_eu = @(s) Chi_p_psi(mu_eu_yt, ploss_eu, s, iota_eu, lambda_eu, eta, 1, varrho)*abs_sc - TED_eu_tgt;

    if sig_us < sig_min_us + 1e-4, sig_us = sig_min_us + max(1e-4, 2*sig_min_us); end
    if sig_eu < sig_min_eu + 1e-4, sig_eu = sig_min_eu + max(1e-4, 2*sig_min_eu); end

    try
        [s, ~, ef] = fsolve(res_us, sig_us, fsolve_opt);
        if ef > 0 && s >= sig_min_us, sig_us = s; else
            sig_us = fminbnd(@(x) res_us(x)^2, sig_min_us+1e-6, 100);
        end
    catch
        sig_us = fminbnd(@(x) res_us(x)^2, sig_min_us+1e-6, 100);
    end
    try
        [s, ~, ef] = fsolve(res_eu, sig_eu, fsolve_opt);
        if ef > 0 && s >= sig_min_eu, sig_eu = s; else
            sig_eu = fminbnd(@(x) res_eu(x)^2, sig_min_eu+1e-6, 100);
        end
    catch
        sig_eu = fminbnd(@(x) res_eu(x)^2, sig_min_eu+1e-6, 100);
    end

    sigma_us_t(tt) = sig_us;
    sigma_eu_t(tt) = sig_eu;
    [~,~,~,~,DW_us_t(tt),FF_us_t(tt)] = Chi_sys(mu_us_yt, ploss_us, sig_us, iota_us, lambda_us, eta, 1, varrho);
    [~,~,~,~,DW_eu_t(tt),FF_eu_t(tt)] = Chi_sys(mu_eu_yt, ploss_eu, sig_eu, iota_eu, lambda_eu, eta, 1, varrho);
end

% Build DWS stock with delta = 0.2
delta_DWS = 0.2;
DWS_us_t = zeros(T,1);
for tt = 1:T
    if tt == 1, prev = 0; else, prev = DWS_us_t(tt-1); end
    DWS_us_t(tt) = DW_us_t(tt) + (1 - delta_DWS) * prev;
end

% Report in main_filter.m format
fprintf('\n=== Verification at user-supplied params ===\n');
fprintf('  lambda_us=%.4f, lambda_eu=%.4f, iota=%.4f, eta=%.4f, varrho=%.4f, delta=%.2f\n', ...
    lambda_us, lambda_eu, iota_ss, eta, varrho, delta_DWS);
fprintf('FF:   model=%.2f%%,  data=%.2f%%\n', mean(FF_us_t)*100, mean(FF_n(~isnan(FF_n)))*100);
fprintf('DW:   model=%.2f%%,  data=%.2f%%   (flow; mismatched units)\n', ...
    mean(DW_us_t)*100, mean(DW_n(~isnan(DW_n)))*100);
fprintf('DWS:  model=%.2f%%,  data=%.2f%%   (stock; delta=%.2f)\n', ...
    mean(DWS_us_t)*100, mean(DW_n(~isnan(DW_n)))*100, delta_DWS);
fprintf('sigma_us: mean=%.4f, min=%.4f, max=%.4f\n', mean(sigma_us_t), min(sigma_us_t), max(sigma_us_t));
fprintf('sigma_eu: mean=%.4f, min=%.4f, max=%.4f\n', mean(sigma_eu_t), min(sigma_eu_t), max(sigma_eu_t));
fprintf('=============================================\n');
