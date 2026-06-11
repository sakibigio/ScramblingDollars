%% model_dw_ff_corr.m
% Model-implied correlation between discount-window borrowing (DW_t) and
% Fed Funds volume (FF_t), for the reply to Referee 2, Comment 3.
% Produces [X] = unconditional corr and [Y] = corr in the high-sigma
% ("scrambling-for-dollars") regime, to drop into
%   responses/r2q3_interbank_dw_corr.tex
%
% Mirrors diagnose_dw.m's filter loop (same prior-optimum parameters) so the
% numbers are consistent with the rest of the calibration diagnostics.

addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1;
pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat'); load('data/calibration.mat');
run('functions/params.m');

T = length(mu_us);

% Prior optimum (same as diagnose_dw.m)
lam = 1.6415; et = 0.9000; iot = 0.0838; vrh = 0.10;
iota_us_loc = iot / freq / pi_us_ss;
theta_plus  = ((exp(lam) - 1)/(exp(lam) + 1))^2;

sigma_t = zeros(T,1); DW_t = zeros(T,1); FF_t = zeros(T,1);
sig_guess = sigma_us; abs_sc = 12e4;
for tt = 1:T
    mu_t    = exp(mu_us(tt));
    TED_tgt = TED_s_us_t(tt)*abs_sc;
    res_fun = @(sig) Chi_p_psi(mu_t, ploss_us, sig, iota_us_loc, lam, et, 1, vrh)*abs_sc - TED_tgt;
    sig_min = find_sigma_min(mu_t, ploss_us, theta_plus, vrh);
    if sig_guess < sig_min + 1e-4
        sig_guess = sig_min + max(0.0001, 2*sig_min);
    end
    [sig_out,~,ef] = fsolve(res_fun, sig_guess, optimoptions('fsolve','Display','off','TolFun',1e-12));
    if ef<=0 || sig_out<sig_min
        sq_fun  = @(s)(Chi_p_psi(mu_t,ploss_us,s,iota_us_loc,lam,et,1,vrh)*abs_sc - TED_tgt)^2;
        sig_out = fminbnd(sq_fun, sig_min+1e-6, 15);
    end
    sig_guess  = sig_out;
    sigma_t(tt)= sig_out;
    [~, ~, ~, ~, DW_t(tt), FF_t(tt), ~] = ...
        Chi_sys(mu_t, ploss_us, sig_out, iota_us_loc, lam, et, 1, vrh);
end

%% Define the high-sigma ("scrambling-for-dollars") regime
% (a) top quartile of filtered sigma_t
thr_q   = quantile(sigma_t, 0.75);
hi_q    = sigma_t >= thr_q;
% (b) if the MS smoothed-probability file is present, use prob > 0.5
hi_ms = [];
fp = 'data/MS_sigma_us_prob.csv';
if exist(fp,'file')
    p = readmatrix(fp);
    p = p(:,end);                      % last column = smoothed prob of high state
    if numel(p) >= T
        hi_ms = p(1:T) > 0.5;
    end
end

cc = @(a,b,m) corr(a(m), b(m), 'rows','complete');

fprintf('\n==== Model-implied DW_t vs FF_t correlation ====\n');
fprintf('[X] unconditional (full sample)        = %+.3f\n', corr(DW_t, FF_t, 'rows','complete'));
fprintf('[Y] high-sigma regime, top quartile    = %+.3f  (n=%d, sigma>=%.3f)\n', cc(DW_t,FF_t,hi_q), sum(hi_q), thr_q);
if ~isempty(hi_ms)
    fprintf('[Y] high-sigma regime, MS prob>0.5     = %+.3f  (n=%d)\n', cc(DW_t,FF_t,hi_ms), sum(hi_ms));
end

%% Data counterpart over the same sample (for the table / sanity check)
ok = ~(isnan(DW_n) | isnan(FF_n));
fprintf('\n---- Data DW_n vs FF_n (same monthly sample, normalized by deposits) ----\n');
fprintf('unconditional                          = %+.3f  (n=%d)\n', corr(DW_n(ok), FF_n(ok)), sum(ok));
fprintf('high-sigma (top quartile of sigma_t)   = %+.3f  (n=%d)\n', cc(DW_n,FF_n,hi_q & ok), sum(hi_q & ok));
