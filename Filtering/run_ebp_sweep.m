%% run_ebp_sweep.m
% CD-consistent feasibility sweep of the steady-state bond-premium target.
% Replicates main_LFX's preamble, then calls the FULL solve_steady_state at
% each target via the lfx_target_ebp hook (so initial guesses, calibration,
% SS solve and multicurrency block run exactly as in the live pipeline).
% Replaces the VOID legacy probe_ebp (Leontief eta=0.5). For the ROBUSTNESS
% section: how low can the EBP target go with a real, positive-sigma solution?

clear; close all;
matching_type = 1; mt_suffix = '_cd';
addpath('functions'); addpath('functions/chi'); addpath('utils'); addpath('data');

printit = 0; plotit = 0; testplot = 0;
foldername = fullfile(pwd,'output','ebp_sweep'); if ~exist(foldername,'dir'), mkdir(foldername); end
T = 5; LFX_plotprefs; freq = 12; rate_scale = 1e4;

pi_eu_ss = 1; pi_us_ss = 1;
load 'Z_im_us.mat'; load 'Z_im_eu.mat';
load dynare_calibration_param.mat;
params;                                   % lambda/eta/iota from override (paper = A)
assert(abs(lambda_us-1.8565)<0.01 && abs(eta-0.6480)<0.01, ...
    'Override is not the paper calibration (lambda=%.4f, eta=%.4f).', lambda_us, eta);
sigma_eu = 0.15; sigma_us = 0.2;
LFX_nt_0e_eqs_2;                          % CD-delegating equations

% Back up shared outputs solve_steady_state writes (M_curr.mat)
if isfile('M_curr.mat'), copyfile('M_curr.mat','M_curr_SWEEPBK.mat'); end

targets = [200 150 100 75 50 35 25];
out = struct('target',{},'sigma_us',{},'sigma_eu',{},'Theta_d',{}, ...
             'mu_us_ss',{},'mu_eu_ss',{},'EBP_bps',{},'is_real',{},'feasible',{});
fprintf('\n%8s %10s %10s %10s %10s %12s %8s\n', ...
    'target','sigma_us','sigma_eu','Theta_d','mu_us_ss','EBP@sol bps','FEAS?');
fprintf('%s\n', repmat('-',1,76));
for t = targets
    lfx_target_ebp = t;
    try
        solve_steady_state;
        ebp_sol = (Rb_us_ss.^freq - Rm_us.^freq)*1e4;
        isr  = isreal([sigmass_us sigmass_eu mu_us_ss Rb_us_ss]) && all(isfinite([sigmass_us sigmass_eu]));
        feas = isr && sigmass_us>0 && sigmass_eu>0 && abs(real(ebp_sol)-t)<1;
        out(end+1) = struct('target',t,'sigma_us',real(sigmass_us),'sigma_eu',real(sigmass_eu), ...
            'Theta_d',real(Theta_dss_us),'mu_us_ss',real(mu_us_ss),'mu_eu_ss',real(mu_eu_ss), ...
            'EBP_bps',real(ebp_sol),'is_real',isr,'feasible',feas); %#ok<SAGROW>
        fprintf('%8.0f %10.4f %10.4f %10.4f %10.4f %12.2f %8s%s\n', ...
            t, real(sigmass_us), real(sigmass_eu), real(Theta_dss_us), ...
            real(mu_us_ss), real(ebp_sol), string(feas), ternstr(~isr,'  [COMPLEX]'));
    catch ME
        fprintf('%8.0f  FAILED: %s\n', t, ME.message);
        out(end+1) = struct('target',t,'sigma_us',NaN,'sigma_eu',NaN,'Theta_d',NaN, ...
            'mu_us_ss',NaN,'mu_eu_ss',NaN,'EBP_bps',NaN,'is_real',false,'feasible',false); %#ok<SAGROW>
    end
    close all;
end
if ~exist('output','dir'), mkdir('output'); end
save(fullfile('output','ebp_sweep_cd.mat'), 'out', 'targets');
fprintf('\n[saved] output/ebp_sweep_cd.mat  (CD, live eta=%.4f, lambda=%.4f)\n', eta, lambda_us);

% Restore shared file
if isfile('M_curr_SWEEPBK.mat'), copyfile('M_curr_SWEEPBK.mat','M_curr.mat'); delete('M_curr_SWEEPBK.mat'); end
fprintf('[restore] M_curr.mat restored\n');

function s = ternstr(c,a), if c, s=a; else, s=''; end, end
