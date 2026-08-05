%% Diagnose why model DW correlates negatively with data DW in the crisis window.
addpath('functions'); addpath('functions/chi'); addpath('data'); addpath('utils');
matching_type = 1;
pi_eu_ss = 1; pi_us_ss = 1;
load('data/LFX_data.mat'); load('data/calibration.mat');
run('functions/params.m');

% Window: Aug 2007 - Dec 2009 (periods 80-108)
T = length(mu_us);
wstart = (2007-2001)*12 + 8;  % 80
wend   = (2009-2001)*12 + 12; % 108
dates  = datenum(2001, 1:T, 1);

% Run filter at the prior optimum
lam = 1.6415; et = 0.9000; iot = 0.0838; vrh = 0.10;
iota_us_loc = iot / freq / pi_us_ss;
theta_plus = ((exp(lam) - 1)/(exp(lam) + 1))^2;
sigma_t = zeros(T,1); Smin_t=zeros(T,1); DW_t=zeros(T,1); FF_t=zeros(T,1); theta_t=zeros(T,1); psi_t=zeros(T,1);
sig_guess = sigma_us;
abs_sc = 12e4;
for tt=1:T
    mu_t = exp(mu_us(tt));
    TED_tgt = TED_s_us_t(tt)*abs_sc;
    res_fun = @(sig) Chi_p_psi(mu_t, ploss_us, sig, iota_us_loc, lam, et, 1, vrh)*abs_sc - TED_tgt;
    sig_min = find_sigma_min(mu_t, ploss_us, theta_plus, vrh);
    if sig_guess < sig_min + 1e-4
        sig_guess = sig_min + max(0.0001, 2*sig_min);
    end
    [sig_out,~,ef] = fsolve(res_fun, sig_guess, optimoptions('fsolve','Display','off','TolFun',1e-12));
    if ef<=0 || sig_out<sig_min
        sq_fun = @(s)(Chi_p_psi(mu_t,ploss_us,s,iota_us_loc,lam,et,1,vrh)*abs_sc - TED_tgt)^2;
        sig_out = fminbnd(sq_fun, sig_min+1e-6, 15);
    end
    sig_guess = sig_out;
    sigma_t(tt)=sig_out;
    [~, theta_t(tt), psi_t(tt), Smin_t(tt), DW_t(tt), FF_t(tt), ~] = ...
        Chi_sys(mu_t, ploss_us, sig_out, iota_us_loc, lam, et, 1, vrh);
end

%% Print summary
fprintf('\n=== Crisis window (Aug 2007 - Dec 2009) ===\n');
win = wstart:wend;
fprintf('Period range: %d - %d (%d obs)\n', wstart, wend, length(win));
fprintf('\n%-22s %12s %12s %12s\n','variable','full-sample','crisis','pre-crisis');
fprintf('%s\n', repmat('-',1,62));
vars = {'exp(mu_us)','sigma_t','TED_s_us_t','DW_n','DW_t (model)','FF_n','FF_t (model)','Smin_t','psi_t','theta_t'};
data = {exp(mu_us), sigma_t, TED_s_us_t, DW_n, DW_t, FF_n, FF_t, Smin_t, psi_t, theta_t};
pre = 1:(wstart-1);
for i=1:length(vars)
    x = data{i};
    fprintf('%-22s %12.4f %12.4f %12.4f\n', vars{i}, mean(x,'omitnan'), mean(x(win),'omitnan'), mean(x(pre),'omitnan'));
end

%% Correlations
fprintf('\n--- Correlations within crisis window ---\n');
fprintf('corr_complete(DW_n, DW_t)      = %.3f\n', corr(DW_n(win), DW_t(win)));
fprintf('corr_complete(DW_n, mu_us)     = %.3f\n', corr(DW_n(win), mu_us(win)));
fprintf('corr_complete(DW_n, sigma_t)   = %.3f\n', corr(DW_n(win), sigma_t(win)));
fprintf('corr_complete(DW_n, TED)       = %.3f\n', corr(DW_n(win), TED_s_us_t(win)));
fprintf('corr_complete(DW_t, mu_us)     = %.3f\n', corr(DW_t(win), mu_us(win)));
fprintf('corr_complete(DW_t, sigma_t)   = %.3f\n', corr(DW_t(win), sigma_t(win)));
fprintf('corr_complete(mu_us, sigma_t)  = %.3f\n', corr(mu_us(win), sigma_t(win)));

%% Decompose model DW = Smin - psi*Spl
fprintf('\n--- Model DW decomposition (means in crisis) ---\n');
fprintf('mean Smin_t    = %.4f\n', mean(Smin_t(win)));
fprintf('mean psi*Spl  ~  Smin - DW = %.4f\n', mean(Smin_t(win)) - mean(DW_t(win)));
fprintf('mean DW_t     = %.4f  (data: %.4f)\n', mean(DW_t(win)), mean(DW_n(win)));
fprintf('std Smin_t in crisis = %.4f   std DW_t = %.4f\n', std(Smin_t(win)), std(DW_t(win)));
fprintf('corr(Smin_t, DW_t) in crisis = %.3f\n', corr(Smin_t(win), DW_t(win)));

%% When does the data DW peak? When does mu peak? When does sigma peak?
fprintf('\n--- Peaks in crisis window ---\n');
[~, k] = max(DW_n(win)); fprintf('DW_n peak at  period %d  (%s)\n', wstart-1+k, datestr(dates(wstart-1+k),'mmm-yy'));
[~, k] = max(sigma_t(win)); fprintf('sigma peak at period %d  (%s)\n', wstart-1+k, datestr(dates(wstart-1+k),'mmm-yy'));
[~, k] = max(exp(mu_us(win))); fprintf('mu_us peak at period %d  (%s)\n', wstart-1+k, datestr(dates(wstart-1+k),'mmm-yy'));
[~, k] = max(DW_t(win)); fprintf('DW_t peak at period %d  (%s)\n', wstart-1+k, datestr(dates(wstart-1+k),'mmm-yy'));

save('_diagnose_dw.mat','sigma_t','DW_t','FF_t','Smin_t','theta_t','psi_t','win','wstart','wend','dates');
