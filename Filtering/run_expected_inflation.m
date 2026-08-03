%% run_expected_inflation.m
% Run the nonlinear model through the global solution (stages [1]-[3] of
% main_LFX.m, Cobb-Douglas) and compute expected inflation state by state.
%
% Equilibrium objects are saved by solve_global.m to:
%   data/initguess.mat                          - converged x0 (warm start
%                                                 for future runs; loaded
%                                                 automatically if it exists
%                                                 and the grid size matches)
%   data/global_sol<nameplot>.mat               - full solution workspace
%
% Expected inflation output (this script):
%   data/expected_inflation_cd.mat / .csv
%
% Two objects per currency, distinguished by the Jensen correction that was
% introduced in feqm_vec.m:
%   Epi   = E[p'/p | s]          (expected gross inflation; rows 12-13 of x)
%   Epi_h = 1 / E[p/p' | s]      (harmonic mean; the inverse of the object
%                                 that deflates policy rates in feqm_vec)
%
% Invoke:  matlab -batch "run_expected_inflation"

matching_type = 1;   % Cobb-Douglas (live pipeline)
lfx_printit   = 0;   % never print figures from this driver
plotit        = 0;
printit       = 0;
mt_suffix     = '_cd';
foldername    = fullfile('.','output','expinf',filesep);
if ~exist(foldername,'dir'); mkdir(foldername); end
lfx_target_ebp = 200;

addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

tic;

% Plot Preferences (same preamble as main_LFX.m)
T = 5;
LFX_plotprefs;
freq = 12;
rate_scale = 1e4;

% Load Parameters
load 'Z_im_us.mat';
load 'Z_im_eu.mat';
pi_eu_ss = 1;
pi_us_ss = 1;
load dynare_calibration_param.mat;
params;

% Sigma initial guesses for the calibration solver
sigma_eu = 0.15;
sigma_us = 0.2;

% Name Scenario (must match main_LFX so the same global_sol file is reused)
nameplot  = 'calibration_dynare_sigma_us';
xperiment = '$\epsilon^{\lambda^{*}}$';

% Load Eq Equations
LFX_nt_0e_eqs_2;

%% [1] Steady State
solve_steady_state;

%% [2] Markov Regime Setup
setup_markov;

%% [3] Global Solution (saves initguess.mat + global_sol...mat)
solve_global;

%% [4] Expected inflation
% Conditional expected gross inflation (monthly), state by state
Epi_us_grid = pi_us_vec;                              % E[p'/p | s]
Epi_eu_grid = pi_eu_vec;
pi_inv_us_grid = (Q_mat*(1./p_us_vec(:)))'.*p_us_vec; % E[p/p' | s]
pi_inv_eu_grid = (Q_mat*(1./p_eu_vec(:)))'.*p_eu_vec;
Epi_us_harm = 1./pi_inv_us_grid;                      % harmonic-mean version
Epi_eu_harm = 1./pi_inv_eu_grid;

% Ergodic distribution of the chain (power iteration on sparse Q_mat)
erg = ones(1,N_s)/N_s;
for it = 1:1e6
    ergn = erg*Q_mat; ergn = ergn/sum(ergn);
    if max(abs(ergn-erg)) < 1e-14, break; end
    erg = ergn;
end
erg = full(erg(:));

% Regime blocks: states 1:N_s/2 = normal (r1), N_s/2+1:N_s = volatile (r2)
i_r1 = 1:N_s/2; i_r2 = N_s/2+1:N_s;
w_r1 = erg(i_r1)/sum(erg(i_r1));
w_r2 = erg(i_r2)/sum(erg(i_r2));

ann = @(gross) (gross.^freq - 1)*100;   % annualized net inflation, percent
bps = @(gross) (gross - 1)*rate_scale;  % monthly net, basis points

summary = struct();
summary.erg_mass_normal   = sum(erg(i_r1));
summary.erg_mass_volatile = sum(erg(i_r2));
for cc = {'us','eu'}
    c = cc{1};
    G  = eval(['Epi_' c '_grid']);
    Gh = eval(['Epi_' c '_harm']);
    summary.(['Epi_' c '_ann_ergodic'])  = ann(G(:)'*erg);
    summary.(['Epi_' c '_ann_normal'])   = ann(G(i_r1)*w_r1);
    summary.(['Epi_' c '_ann_volatile']) = ann(G(i_r2)*w_r2);
    summary.(['Epi_' c '_harm_ann_ergodic']) = ann(Gh(:)'*erg);
    summary.(['jensen_gap_' c '_bps_ergodic']) = (G(:)'*erg - Gh(:)'*erg)*rate_scale;
end

fprintf('\n================ Expected inflation (annualized, %%) ================\n');
fprintf('%-28s %10s %10s\n','', 'US','EU');
fprintf('%-28s %10.4f %10.4f\n','Ergodic mean, E[pi]',        summary.Epi_us_ann_ergodic,  summary.Epi_eu_ann_ergodic);
fprintf('%-28s %10.4f %10.4f\n','Normal regime mean',          summary.Epi_us_ann_normal,   summary.Epi_eu_ann_normal);
fprintf('%-28s %10.4f %10.4f\n','Volatile regime mean',        summary.Epi_us_ann_volatile, summary.Epi_eu_ann_volatile);
fprintf('%-28s %10.4f %10.4f\n','Ergodic mean, harmonic',      summary.Epi_us_harm_ann_ergodic, summary.Epi_eu_harm_ann_ergodic);
fprintf('%-28s %10.4f %10.4f\n','Jensen gap (bps, monthly)',   summary.jensen_gap_us_bps_ergodic, summary.jensen_gap_eu_bps_ergodic);
fprintf('Ergodic regime mass: normal %.3f / volatile %.3f\n', summary.erg_mass_normal, summary.erg_mass_volatile);
fprintf('=====================================================================\n');

% Persist: .mat with grids + summary, .csv state-by-state
save(fullfile('data',['expected_inflation' mt_suffix '.mat']), ...
    'Epi_us_grid','Epi_eu_grid','Epi_us_harm','Epi_eu_harm', ...
    'pi_inv_us_grid','pi_inv_eu_grid','erg','sigma_us_vec','freq','rate_scale','summary');

Tcsv = table((1:N_s)', [ones(N_s/2,1); 2*ones(N_s/2,1)], sigma_us_vec(:), erg, ...
    Epi_us_grid(:), Epi_us_harm(:), Epi_eu_grid(:), Epi_eu_harm(:), ...
    ann(Epi_us_grid(:)), ann(Epi_eu_grid(:)), ...
    'VariableNames', {'state','regime','sigma_us','erg_prob', ...
    'Epi_us_gross_m','Epi_us_harm_gross_m','Epi_eu_gross_m','Epi_eu_harm_gross_m', ...
    'Epi_us_ann_pct','Epi_eu_ann_pct'});
writetable(Tcsv, fullfile('data',['expected_inflation' mt_suffix '.csv']));

fprintf('Saved: data/expected_inflation%s.mat and .csv\n', mt_suffix);
toc;
