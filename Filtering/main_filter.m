
%% Main Filter Script
% Code for Scrambling for Dollars
% (c) Saki Bigio
%
% Pipeline: main_filter.m → markov_estimation.jl → plot_regimes.m

%% Setup - Preserve matching_type through clear
if ~exist('matching_type', 'var')
    matching_type = 1;  % 0 = Leontief, 1 = Cobb-Douglas (default)
end
save('temp_matching_type.mat', 'matching_type');

%% Setup
clear; close all;

%% Load matching_type
load('temp_matching_type.mat');
delete('temp_matching_type.mat');  % Clean up temp file
% matching_type = saved_matching_type;  % Restore after clear

% Add paths
addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

% Print/plot options
printit =  1;
plotdata = 0;
printver = 0;

% Plotting flags (set to 1 to enable)
do_plot_baseline    = 1;  % Main filter results
do_plot_diagnostics = 0;  % Diagnostic plots
do_counterfactual   = 0;  % Counterfactual analysis (not implemented)
do_sensitivity      = 0;  % Sensitivity analysis (not implemented)
do_regimes          = 1;  % Plot paper results
plot_baseline_curr  = 0;  % Baseline plots for other currencies


% Run markov_estimation.jl automatically after filtering
run_julia = 1;  % Set to 1 to run Julia automatically

% Liquidity-ratio baseline correction:
%   When use_mu_baseline = 1, the filter uses
%       mu_*_yt = exp(mu_*(tt)) - mubar_*(tt)
%   where mubar_* is the annual regulatory baseline from
%   functions/liquidity_baseline.m. m_eff is then "excess reserves above
%   the regulatory baseline" and can be positive or negative -- the chi
%   helpers handle both signs via the piecewise closed form.
use_mu_baseline = 0.1;     % 0 = raw exp(mu); 1 = subtract annual baseline
                          % (must be 1 to match estimate_params_4d.m)

% Stationary-mu variant: HP-detrend the liquidity ratio BEFORE filtering.
%   When mu_hp_detrend = 1, the level ratio exp(mu) is replaced by its
%   (very-low-frequency) HP cyclical component recentered at mu_hp_center,
%   via functions/mu_hp_cyc_level.m. This runs the whole filter against a
%   stationary mu series. Set to 0 to recover the raw (trending) mu used in
%   the locked-in/published run.
mu_hp_detrend = 0;        % 0 = raw mu (published); 1 = HP-cyclical mu
%mu_hp_lambda  = 1e6;      % very-low-frequency HP smoothing (monthly)
% mu_hp_center  = 0.2;      % level at which the cyclical liquidity ratio sits

% Empirical-LCR variant: net out the FEDS-Note LCR/HQLA series (column EM,
% loaded as LCR_us in load_data.m) from the US liquidity ratio IN LEVELS, so the
% filter sees the EXCESS liquidity ratio above the requirement:
%   exp(mu_us) - LCR_us/100  (additive gap, not a ratio).
% Done via functions/mu_minus_lcr_level.m. US only -- there is no EU LCR series,
% so mu_eu is left untouched. Mutually exclusive with mu_hp_detrend for mu_us.
mu_minus_lcr = 1;         % 0 = off; 1 = subtract LCR_us from mu_us in levels

%% Load data
load('data/LFX_data.mat');
if mu_hp_detrend == 1 && mu_minus_lcr == 1
    error('main_filter:muOptionConflict', ...
          'Set only one of mu_hp_detrend / mu_minus_lcr (both detrend mu_us).');
end
if mu_hp_detrend == 1
    mu_us = mu_hp_cyc_level(mu_us, mu_hp_lambda, mu_hp_center);
    mu_eu = mu_hp_cyc_level(mu_eu, mu_hp_lambda, mu_hp_center);
    fprintf(['mu HP-detrend ON: lambda=%.3g, center=%.2f | ' ...
             'exp(mu_us) range=[%.3f,%.3f], exp(mu_eu) range=[%.3f,%.3f]\n'], ...
        mu_hp_lambda, mu_hp_center, min(exp(mu_us)), max(exp(mu_us)), ...
        min(exp(mu_eu)), max(exp(mu_eu)));
end
if mu_minus_lcr == 1
    mu_us = mu_minus_lcr_level(mu_us, LCR_us);
    fprintf('mu-minus-LCR ON: exp(mu_us) excess range=[%.3f,%.3f]\n', ...
        min(exp(mu_us)), max(exp(mu_us)));
end

if plotdata == 1
    run('plotting/plot_data.m');
end

%% Parameters
load('data/calibration.mat');

% Scale: Report Variables in BPS
abs_scale = 12e4;

% List of Currencies
curlist = {'au','ca','jp','nz','no','sw','ch','uk'};
conlist = {'AUD','CAD','JPY','NZD','NOK','SWK','CHF','GBP'};
CURRlist = {'EUR','AUD','CAD','JPY','NZD','NOK','SWK','CHF','GBP'};

% OMO share in MBS 
share_us = 1; 
share_eu = 1;

% Steady state inflation rate
pi_eu_ss = 1;
pi_us_ss = 1;

% Load remaining parameters
run('functions/params.m');

% Additional parameters
epsilon_b = -0.5;
zeta_us = -epsilon_b; 
zeta_eu = -epsilon_b;
Theta_b = 3;

% Date variables
T = length(mu_us);
datesperiod = 1:T;
dates = datenum(2001, 1:T, 1);

% Annual liquidity-ratio baseline (broadcast to monthly)
if use_mu_baseline == 1
    [mubar_us_t, mubar_eu_t] = liquidity_baseline(dates);
    fprintf('Baseline correction ON: mubar_us range=[%.3f,%.3f], mubar_eu range=[%.3f,%.3f]\n', ...
        min(mubar_us_t), max(mubar_us_t), min(mubar_eu_t), max(mubar_eu_t));
else
    mubar_us_t = zeros(T,1);
    mubar_eu_t = zeros(T,1);
end

% Plot formats
FSize = 20;
formataxis = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize, 'Box', 'Off', 'PlotBoxAspectRatio', [1 0.75 1]);
label_x = @(x) xlabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize, 'Interpreter', 'latex');
label_y = @(x) ylabel(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize, 'interpreter', 'latex');
desiredDecimalPlaces = 1;

%% Tests for lower bounds of interbank variables
sigma_vec = (0.01:0.01:10);
N_s = length(sigma_vec);
min_ted = min(TED_s_us_t * abs_scale);
max_ted = max(TED_s_us_t * abs_scale);
if matching_type == 0
    min_test_us = Chi_p_psi(min(exp(mu_us)), ploss_us, min(sigma_vec), iota_us, lambda_us, eta, matching_type, varrho) * abs_scale + 0.00012 * abs_scale;
    min_test_eu = Chi_p_psi(min(exp(mu_eu)), ploss_eu, min(sigma_vec), iota_eu, lambda_eu, eta, matching_type, varrho) * abs_scale + 0.00012 * abs_scale;
    max_test = Chi_p_psi(min(exp(mu_us)), ploss_us, max(sigma_vec), iota_us, lambda_us, eta, matching_type, varrho) * abs_scale;
    fprintf('Leontief matching: TED range in data is [%.2f, %.2f] bps\n', min_ted, max_ted);
    fprintf('Leontief matching: TED range in model is [%.2f, %.2f] bps\n', min_test_us, max_test);
else
    min_test_us = 0 ;
    min_test_eu = 0 ;
    max_test = Chi_p_psi(min(exp(mu_us)), ploss_us, max(sigma_vec), iota_us, lambda_us, eta, matching_type, varrho) * abs_scale;
    fprintf('Cobb-Douglas matching: TED range in data is [%.2f, %.2f] bps\n', min_ted, max_ted);
    fprintf('Cobb-Douglas matching: TED range in model is [%.2f, %.2f] bps\n', min_test_us, max_test);
end


% Test lower bound
Ted_test = ones(N_s, 1);
BP_test = ones(N_s, 1);
Psi_test = ones(N_s, 1);
for ss = 1:N_s
    [Ted_test(ss), ~, Psi_test(ss)] = Chi_p_psi(max(exp(mu_us)), ploss_us, sigma_vec(ss), iota_us, lambda_us, eta, matching_type, varrho);
    Ted_test(ss) = Ted_test(ss) * 1e4 * 12;
    BP_test(ss) = Echi_m(max(exp(mu_us)), ploss_us, sigma_vec(ss), iota_us, lambda_us, eta, matching_type, varrho) * 1e4 * 12;
end

figure('Name', 'Ted test', 'NumberTitle', 'off'); 
plot(sigma_vec, Ted_test); hold on;
plot(sigma_vec, BP_test);
scatter(0 * TED_s_us_t, TED_s_us_t * abs_scale, 'filled'); hold on;
scatter(0 * TED_s_eu_t, TED_s_eu_t * abs_scale, 'filled'); hold on;
title(sprintf('TED Test (matching\\_type = %d)', matching_type));
legend('TED model', 'BP model', 'TED data US', 'TED data EU');

figure('Name', 'Psi test', 'NumberTitle', 'off'); 
plot(sigma_vec, Psi_test);
title(sprintf('Psi Test (matching\\_type = %d)', matching_type));

%% Initialize data
load('data/LFX_data.mat');
if mu_hp_detrend == 1
    mu_us = mu_hp_cyc_level(mu_us, mu_hp_lambda, mu_hp_center);
    mu_eu = mu_hp_cyc_level(mu_eu, mu_hp_lambda, mu_hp_center);
end
if mu_minus_lcr == 1
    mu_us = mu_minus_lcr_level(mu_us, LCR_us);
end
endopath = [mu_us mu_eu Rb_Rm cip];
exopath = [im_eu im_us M_us M_eu];

%% Pre-allocate vectors
T = length(endopath);

% Price variables
sigma_us_t = zeros(T, 1);
sigma_us_bp_t = zeros(T, 1);
sigma_us_TED_t = zeros(T, 1);
sigma_eu_TED_t = zeros(T, 1);
Echi_m_us_t = zeros(T, 1);
Echi_d_us_t = zeros(T, 1);
Chi_p_psi_us_t = zeros(T, 1);
sigma_eu_t = zeros(T, 1);
sigma_eu_bp_t = zeros(T, 1);
Echi_m_eu_t = zeros(T, 1);
Echi_d_eu_t = zeros(T, 1);
Chi_p_psi_eu_t = zeros(T, 1);

% Interbank variables
theta_us_t = zeros(T, 1);
psi_us_t = zeros(T, 1);
Smin_us_t = zeros(T, 1);
DW_us_t = zeros(T, 1);
FF_us_t = zeros(T, 1);
Q_us_t = zeros(T, 1);
theta_eu_t = zeros(T, 1);
psi_eu_t = zeros(T, 1);
Smin_eu_t = zeros(T, 1);
DW_eu_t = zeros(T, 1);
FF_eu_t = zeros(T, 1);
Q_eu_t = zeros(T, 1);

% Diagnostics
sigma_us_flag = zeros(T, 1);
sigma_eu_flag = zeros(T, 1);
sigma_us_TED_flag = zeros(T, 1);
sigma_eu_TED_flag = zeros(T, 1);
rest_flag = zeros(T, 1);

% Initial guesses
sigma_us_BP_guess = sigma_us;
sigma_eu_BP_guess = sigma_eu;
sigma_us_TED_guess = sigma_us;
sigma_eu_TED_guess = sigma_eu;

%% Main Filtering Loop
TED_target_i = 1;
TED_eu_target_i = 1; 
BP_eu_target_i = 0;
BP_us_target_i = 0;

if matching_type == 0
    matching_name = 'Leontief';
else
    matching_name = 'Cobb-Douglas';
end
fprintf('Starting filter with matching_type = %d (%s)\n', matching_type, matching_name);

% For Cobb-Douglas: compute Walrasian threshold theta_plus
if matching_type == 1
    theta_plus_us = ((exp(lambda_us) - 1) / (exp(lambda_us) + 1))^2;
    theta_plus_eu = ((exp(lambda_eu) - 1) / (exp(lambda_eu) + 1))^2;
    fprintf('Cobb-Douglas thresholds: θ⁺_us = %.6f, θ⁺_eu = %.6f\n', theta_plus_us, theta_plus_eu);
end

for tt = 1:T
    % Setup targets
    if matching_type == 0
        BP_us_taget = min_test_us + (Rb_Rm(tt) - min(Rb_Rm)) * abs_scale;
        BP_eu_taget = min_test_eu + (Rb_Rm_eu(tt) - min(Rb_Rm_eu)) * abs_scale;
        TED_us_target = min_test_us + (TED_s_us_t(tt) - min(TED_s_us_t)) * abs_scale; 
        TED_eu_target = min_test_eu + (TED_s_eu_t(tt) - min(TED_s_eu_t)) * abs_scale;
    else
        BP_us_taget   = Rb_Rm(tt)       * abs_scale;
        BP_eu_taget   = Rb_Rm_eu(tt)    * abs_scale;
        TED_us_target = TED_s_us_t(tt)  * abs_scale; 
        TED_eu_target = TED_s_eu_t(tt)  * abs_scale;
    end

    % Data points
    if use_mu_baseline == 1
        mu_us_yt = exp(mu_us(tt)) - mubar_us_t(tt);
        mu_eu_yt = exp(mu_eu(tt)) - mubar_eu_t(tt);
    else
        mu_us_yt = exp(mu_us(tt));
        mu_eu_yt = exp(mu_eu(tt));
    end
    
    % [1] US Sigma - Bond Premia
    if BP_us_target_i == 1
        sigma_us_res = @(sigma) Echi_m(mu_us_yt, ploss_us, sigma, iota_us, lambda_us, eta, matching_type, varrho) * 1e4 * 12 - BP_us_taget;
        [sigma_out, ~, exitflag, ~] = fsolve(@(sigma) sigma_us_res(sigma), 1, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
        sigma_us_bp_t(tt) = sigma_out;
        sigma_us_BP_guess = sigma_out;
        sigma_us_flag(tt) = exitflag;
    end

    % [2] US Sigma - TED Target
    sigma_us_res = @(sigma) Chi_p_psi(mu_us_yt, ploss_us, sigma, iota_us, lambda_us, eta, matching_type, varrho) * 1e4 * 12 - TED_us_target;

    if matching_type == 1
        % Cobb-Douglas: compute sigma_min where theta = theta_plus
        sigma_min_us = find_sigma_min(mu_us_yt, ploss_us, theta_plus_us, varrho);
        
        if tt > 1 && sigma_us_t(tt-1) > 0
                sigma_us_TED_guess = sigma_us_t(tt-1);
            else
                sigma_us_TED_guess = sigma_min_us + max(0.0001, 2*sigma_min_us);  % Start above the cliff
        end
        % Walk initial guess up if Chi_p_psi is NaN there
        sigma_us_TED_guess = max(sigma_us_TED_guess, sigma_min_us + max(1e-4, 2*sigma_min_us));
        for kk = 1:6
            if isfinite(sigma_us_res(sigma_us_TED_guess)), break; end
            sigma_us_TED_guess = sigma_us_TED_guess * 1.5 + 0.1;
        end

        % Try fsolve first, but tolerate it throwing on bad residuals
        try
            [sigma_out, ~, exitflag, ~] = fsolve(sigma_us_res, sigma_us_TED_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
        catch
            sigma_out = NaN; exitflag = -1;
        end

        % If fsolve fails, use bounded solver as backup.
        if exitflag <= 0 || ~isfinite(sigma_out) || sigma_out < sigma_min_us
            sigma_res_sq = @(sig) (Chi_p_psi(mu_us_yt, ploss_us, sig, iota_us, lambda_us, eta, 1, varrho) * 1e4 * 12 - TED_us_target)^2;
            sig_hi = max(sigma_min_us + 1e-3, 200);
            sigma_out = fminbnd(sigma_res_sq, sigma_min_us + 1e-6, sig_hi);
            if isempty(sigma_out) || ~isfinite(sigma_out)
                if tt > 1, sigma_out = sigma_us_t(tt-1); else, sigma_out = sigma_min_us + 1e-3; end
            end
            exitflag = 10;  % Flag 10 = solved via fminbnd
        end
    else
        % Leontief: use standard fsolve
        [sigma_out, ~, exitflag, ~] = fsolve(sigma_us_res, sigma_us_TED_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
    end
    sigma_us_TED_t(tt) = sigma_out;
    sigma_us_TED_guess = sigma_out;
    sigma_us_TED_flag(tt) = exitflag;
    
    % [3] EU Sigma - Bond Premium
    if BP_eu_target_i == 1
        sigma_eu_res = @(sigma) Echi_m(mu_eu_yt, ploss_us, sigma, iota_eu, lambda_eu, eta, matching_type, varrho) * 1e4 * 12 - BP_eu_taget;
        [sigma_out, ~, exitflag, ~] = fsolve(@(sigma) sigma_eu_res(sigma), sigma_eu, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
        sigma_eu_bp_t(tt) = sigma_out;
        sigma_eu_flag(tt) = exitflag;
    end

    % [4] EU Sigma - TED target
    sigma_eu_res = @(sigma) Chi_p_psi(mu_eu_yt, ploss_eu, sigma, iota_eu, lambda_eu, eta, matching_type, varrho) * 1e4 * 12 - TED_eu_target;

    if matching_type == 1
        % Cobb-Douglas: compute sigma_min where theta = theta_plus
        sigma_min_eu = find_sigma_min(mu_eu_yt, ploss_eu, theta_plus_eu, varrho);
        if tt > 1 && sigma_eu_t(tt-1) > 0
                sigma_eu_TED_guess = sigma_eu_t(tt-1);
            else
                sigma_eu_TED_guess = sigma_min_eu + max(0.0001, 2*sigma_min_eu);  % Start above the cliff
        end
        
        % Walk initial guess up if residual is NaN there
        sigma_eu_TED_guess = max(sigma_eu_TED_guess, sigma_min_eu + max(1e-4, 2*sigma_min_eu));
        for kk = 1:6
            if isfinite(sigma_eu_res(sigma_eu_TED_guess)), break; end
            sigma_eu_TED_guess = sigma_eu_TED_guess * 1.5 + 0.1;
        end

        % Try fsolve first, tolerate throw on bad residual
        try
            [sigma_out, ~, exitflag, ~] = fsolve(sigma_eu_res, sigma_eu_TED_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
        catch
            sigma_out = NaN; exitflag = -1;
        end

        % If fsolve fails, use bounded solver as backup.
        if exitflag <= 0 || ~isfinite(sigma_out) || sigma_out < sigma_min_eu
            sigma_res_sq = @(sig) (Chi_p_psi(mu_eu_yt, ploss_eu, sig, iota_eu, lambda_eu, eta, 1, varrho) * 1e4 * 12 - TED_eu_target)^2;
            sig_hi = max(sigma_min_eu + 1e-3, 100);
            sigma_out = fminbnd(sigma_res_sq, sigma_min_eu + 1e-6, sig_hi);
            if isempty(sigma_out) || ~isfinite(sigma_out)
                if tt > 1, sigma_out = sigma_eu_t(tt-1); else, sigma_out = sigma_min_eu + 1e-3; end
            end
            exitflag = 10;  % Flag 10 = solved via fminbnd
        end
    else
        % Leontief: use standard fsolve
        [sigma_out, ~, exitflag, ~] = fsolve(sigma_eu_res, sigma_eu_TED_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
    end
    sigma_eu_TED_t(tt) = sigma_out;
    sigma_eu_TED_guess = sigma_out;
    sigma_eu_TED_flag(tt) = exitflag; 

    % Update based on target choice
    if TED_target_i == 1
        sigma_us_t(tt) = sigma_us_TED_t(tt); 
        sigma_us_TED_guess = sigma_us_TED_t(tt); 
        sigma_eu_t(tt) = sigma_eu_TED_t(tt);
        sigma_eu_TED_guess = sigma_eu_TED_t(tt);
    elseif BP_eu_target_i == 1
        sigma_us_t(tt) = sigma_us_bp_t(tt);
        sigma_eu_t(tt) = sigma_eu_bp_t(tt);
    end

    % [5] Update interbank variables
    Echi_m_us_t(tt) = Echi_m(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type, varrho);
    Echi_d_us_t(tt) = Echi_d(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type, varrho);
    Chi_p_psi_us_t(tt) = Chi_p_psi(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type, varrho);
    Echi_m_eu_t(tt) = Echi_m(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);
    Echi_d_eu_t(tt) = Echi_d(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);
    Chi_p_psi_eu_t(tt) = Chi_p_psi(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);

    % Rest of system
    [~, theta_us_t(tt), psi_us_t(tt), Smin_us_t(tt), DW_us_t(tt), FF_us_t(tt), Q_us_t(tt)] = Chi_sys(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type, varrho);
    [~, theta_eu_t(tt), psi_eu_t(tt), Smin_eu_t(tt), DW_eu_t(tt), FF_eu_t(tt), Q_eu_t(tt)] = Chi_sys(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);

    % Sanity check: with extreme parameters Chi_sys can return imaginary or
    % non-positive DW/FF in pathological periods.  Refuse to proceed: hiding
    % these via NaN would (a) silently corrupt the estimation objective if
    % this code path is reused, (b) corrupt downstream plots and moments.
    bad_us = ~isreal(DW_us_t(tt)) || ~isreal(FF_us_t(tt)) || ...
             ~isfinite(DW_us_t(tt)) || ~isfinite(FF_us_t(tt)) || ...
             DW_us_t(tt) <= 0 || FF_us_t(tt) <= 0;
    bad_eu = ~isreal(DW_eu_t(tt)) || ~isreal(FF_eu_t(tt)) || ...
             ~isfinite(DW_eu_t(tt)) || ~isfinite(FF_eu_t(tt)) || ...
             DW_eu_t(tt) <= 0 || FF_eu_t(tt) <= 0;
    if bad_us || bad_eu
        error(['main_filter: Chi_sys produced invalid output at t=%d ' ...
               '(imag/non-positive/non-finite DW or FF).  ' ...
               'sigma_us=%.4g sigma_eu=%.4g DW_us=%s FF_us=%s DW_eu=%s FF_eu=%s.  ' ...
               'Parameters are outside the filter''s feasible region; ' ...
               're-estimate or relax bounds.'], tt, ...
               sigma_us_t(tt), sigma_eu_t(tt), ...
               num2str(DW_us_t(tt)), num2str(FF_us_t(tt)), ...
               num2str(DW_eu_t(tt)), num2str(FF_eu_t(tt)));
    end

    % Progress
    if mod(tt, 50) == 0
        fprintf('  Period %d/%d complete\n', tt, T);
    end
    
    clear mu_us_yt mu_eu_yt TED_us_target TED_eu_target;
end

% Construct rates
Rm_us = exp(im_us);
Rm_eu = exp(im_eu);
inv_e_yt = exp(inv_e);   
Rb_Rm_yt = Rb_Rm;
    
%% Derive other variables
BP_us_t = Echi_m_us_t;
TED_us_t = Chi_p_psi_us_t;
Rd_us_t = Rm_us + Echi_m_us_t + Echi_d_us_t;
Rb_us_t = Rm_us + Echi_m_us_t;
BP_eu_t = Echi_m_eu_t;
TED_eu_t = Chi_p_psi_eu_t;
Rd_eu_t = Rm_eu + Echi_m_eu_t + Echi_d_eu_t;
Rb_eu_t = Rm_eu + Echi_m_eu_t;

% Risk premium 
riskprm_t = Rb_us_t ./ Rb_eu_t - 1;

% UIP deviation
UIP_t = Rm_eu - Rm_us;

% CIP 
CIP_t = UIP_t + Rm_eu .* riskprm_t;

% Forward FX
f_t = (1 + riskprm_t) ./ exp(inv_e);
CIP_check_t = exp(im_eu) .* f_t .* exp(inv_e) - exp(im_us) - CIP_t;

%% Discount-window LOAN STOCK (model side)
% Data DW_n is an outstanding stock (loans accumulated, depreciating).
% Model DW_*_t is a flow (new loans per period).  Build a stock from the
% flow with depreciation rate delta:  DWS_t = DW_t + (1-delta) * DWS_{t-1}.
delta_DWS = 0.2;     % monthly depreciation (loan repayment) rate
DWS_us_t = zeros(T,1);
DWS_eu_t = zeros(T,1);
for tt = 1:T
    if tt == 1
        prev_us = 0;  prev_eu = 0;
    else
        prev_us = DWS_us_t(tt-1);  prev_eu = DWS_eu_t(tt-1);
    end
    DWS_us_t(tt) = DW_us_t(tt) + (1 - delta_DWS) * prev_us;
    DWS_eu_t(tt) = DW_eu_t(tt) + (1 - delta_DWS) * prev_eu;
end
fprintf('DW stock series built (delta=%.2f). Mean DWS_us=%.4g, DWS_eu=%.4g\n', ...
    delta_DWS, mean(DWS_us_t), mean(DWS_eu_t));

%% Quantity variables
inv_e_t = exp(inv_e);
mu_us_t = exp(mu_us);
mu_eu_t = exp(mu_eu);
M_us_t = exp(M_us);
M_eu_t = exp(M_eu);

nu_t = (M_eu_t ./ M_us_t) .* inv_e_t .* mu_us_t ./ mu_eu_t;
b_t = Theta_b * (Rm_us + Echi_m_us_t).^(1/epsilon_b);
d_us_t = b_t ./ ((1 - mu_us_t) + nu_t .* (1 - mu_eu_t));
d_eu_t = nu_t .* d_us_t;

% Checks
budget_check = b_t + (mu_us_t - 1) .* d_us_t + (mu_eu_t - 1) .* d_eu_t;
Rb_check = Rm_us + Rb_Rm - ((nu_t .* (1 - mu_eu_t) + 1 - mu_us_t) .* d_us_t / Theta_b).^epsilon_b;

% Price system  
p_us_t = M_us_t ./ (d_us_t .* mu_us_t);
p_eu_t = M_eu_t ./ (d_eu_t .* mu_eu_t);
e_euus_t = 1 ./ inv_e_t;
M_eu_t_check = p_eu_t .* (d_eu_t .* mu_eu_t) - M_eu_t;
M_us_t_check = p_us_t .* (d_us_t .* mu_us_t) - M_us_t;
e_check = p_us_t ./ p_eu_t - inv_e_t;

% Deposit funding shocks
N_pre = 72;
MBS_us_ss = mean(d_us_t(1:N_pre));
MBS_eu_ss = mean(d_eu_t(1:N_pre));
Theta_d_us_t = (d_us_t - share_us * (mu_us_t .* d_us_t) + MBS_us_ss) ./ (Rd_us_t).^(1/zeta_us);
Theta_d_eu_t = (d_eu_t - share_eu * (mu_eu_t .* d_eu_t) + MBS_eu_ss) ./ (Rd_eu_t).^(1/zeta_eu);

% Risk premium differential
riskprm_shock_t = riskprm_t(datesperiod) - (exp(im_us(datesperiod)) - exp(im_eu(datesperiod)));

%% Save output for Markov estimation
% Include mu_us and mu_eu so the Julia MS estimator can use them as exogenous
% switching regressors (decomposition of σ into μ-driven + regime-driven parts).
sigma_mat = [sigma_us_t(:) sigma_eu_t(:) mu_us(:) mu_eu(:) LCR_us(:)];
RW_shock = array2table(sigma_mat, ...
    'VariableNames', {'sigma_us', 'sigma_eu', 'mu_us', 'mu_eu', 'lcr_us'});
writetable(RW_shock, 'RW_shock.csv');

fprintf('Filter complete. Output saved to RW_shock.csv\n');
fprintf('Matching type: %d (%s)\n', matching_type, matching_name);

%% Run Markov Estimation in Julia (optional)
if run_julia == 1
    fprintf('Running Markov estimation in Julia...\n');
    
    % Try to find Julia
    [status, julia_path] = system('which julia');
    if status ~= 0
        % Try common locations
        if exist('/usr/local/bin/julia', 'file')
            julia_path = '/usr/local/bin/julia';
        elseif exist('/opt/homebrew/bin/julia', 'file')  % M1 Mac
            julia_path = '/opt/homebrew/bin/julia';
        elseif exist('/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia', 'file')
            julia_path = '/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia';
        elseif exist('/Users/sakibigio/.juliaup/bin/julia', 'file')
            julia_path = '/Users/sakibigio/.juliaup/bin/julia';
        else
            warning('Julia not found. Run markov_estimation.jl manually.');
            julia_path = '';
        end
    else
        julia_path = strtrim(julia_path);
    end
    
    if ~isempty(julia_path)
        julia_script = 'markov_estimation.jl';
        if exist(julia_script, 'file')
            cmd = sprintf('"%s" "%s"', julia_path, julia_script);
            [status, result] = system(cmd);
            if status == 0
                fprintf('Julia estimation complete.\n');
            else
                warning('Julia estimation failed: %s', result);
            end
        else
            warning('Julia script not found: %s', julia_script);
        end
    end
end

%% Plotting
if do_plot_baseline == 1
    run('plotting/plot_baseline.m');
end

if do_plot_diagnostics == 1
    run('plotting/plot_diagnostics.m');
end

if do_counterfactual == 1
    run('plotting/plot_counterfactual.m');
end

if do_sensitivity == 1
    run('plotting/plot_sensitivity.m');
end

if do_regimes == 1
    run('plotting/plot_regimes.m');
end

%% Other currencies
% Pre-allocate
sigma_c_TED_t = zeros(T, 1);
sigma_c_TED_flag = zeros(T, 1);
Echi_m_c_t = zeros(T, 1);
Echi_d_c_t = zeros(T, 1);
Chi_p_psi_c_t = zeros(T, 1);
Echi_m_eu_bp_t = zeros(T, 1);
Echi_d_eu_bp_t = zeros(T, 1);
BP_c_t = zeros(T, 1);
Rd_c_t = zeros(T, 1);
TED_c_t = zeros(T, 1);
theta_c_t = zeros(T, 1);
psi_c_t = zeros(T, 1);
Smin_c_t = zeros(T, 1);
DW_c_t = zeros(T, 1);
FF_c_t = zeros(T, 1);
Q_c_t = zeros(T, 1);

for cc = 1:numel(curlist)
    eval(['Rm_' curlist{cc} '=exp(im_' curlist{cc} ');']);
    eval(['Rm_c=Rm_' curlist{cc} ';']);
    
    for tt = 1:T
        if use_mu_baseline == 1
            mu_eu_yt = exp(mu_eu(tt)) - mubar_eu_t(tt);
        else
            mu_eu_yt = exp(mu_eu(tt));
        end
        eval(['target=min_test_eu+TED_s_' curlist{cc} '_t(tt)*abs_scale;']);        
        
        % For Cobb-Douglas, use mu-dependent initial guess
        if matching_type == 1
            sigma_c_guess = max(sigma_eu_t(tt), mu_eu_yt + 0.15);
        else
            sigma_c_guess = sigma_eu_t(tt);
        end
    
        sigma_res = @(sigma) Chi_p_psi(mu_eu_yt, ploss_eu, sigma, iota_eu, lambda_eu, eta, matching_type, varrho) * 1e4 * 12 - target;
        [sigma_out, ~, exitflag, ~] = fsolve(@(sigma) sigma_res(sigma), sigma_c_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
        
        if exitflag > 0
            sigma_c_TED_t(tt) = sigma_out;
        else
            if tt > 1
                sigma_c_TED_t(tt) = sigma_c_TED_t(tt-1);
            else
                sigma_c_TED_t(tt) = sigma_eu_t(tt);  % Use EU value as fallback
            end
        end
        sigma_c_TED_flag(tt) = exitflag;
        TED_c_t(tt) = Chi_p_psi(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);

        Echi_m_c_t(tt) = Echi_m(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);
        Echi_d_c_t(tt) = Echi_d(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);
        Chi_p_psi_c_t(tt) = Chi_p_psi(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);

        [~, theta_c_t(tt), psi_c_t(tt), Smin_c_t(tt), DW_c_t(tt), FF_c_t(tt), Q_c_t(tt)] = Chi_sys(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type, varrho);
    end
    clear mu_eu_yt target;

    BP_c_t = Echi_m_c_t;
    Rb_c_t = Rm_c + Echi_m_c_t;
    Rd_c_t = Rm_c + Echi_m_c_t + Echi_d_c_t;
    UIP_c_t = Rm_c - Rm_us;
    riskprm_c_t = (Rb_us_t) ./ (Rb_c_t) - 1;
    CIP_c_t = UIP_c_t + Rm_c .* riskprm_c_t;
    
    % Price and FX calculations (match old code formula)
    p_c_t = (M_eu ./ (Rd_c_t.^(1/zeta_eu))) ./ mu_eu_t;
    inv_e_c_t = p_us_t ./ p_c_t;  % Level ratio (same as old code)
    f_c_t = (1 + riskprm_c_t) ./ inv_e_c_t;

    eval(['sigma_' curlist{cc} '_TED_flag=sigma_c_TED_flag;']);
    eval(['sigma_' curlist{cc} '_t=sigma_c_TED_t;']);
    eval(['TED_' curlist{cc} '_t=TED_c_t;']);
    eval(['BP_' curlist{cc} '_t=BP_c_t;']);
    eval(['riskprm_' curlist{cc} '_t=riskprm_c_t;']);
    eval(['UIP_' curlist{cc} '_t=UIP_c_t;']);
    eval(['CIP_' curlist{cc} '_t=CIP_c_t;']);
    eval(['inv_e_' curlist{cc} '_t=inv_e_c_t;']);
    eval(['p_' curlist{cc} '_t=p_c_t;']);
    eval(['f_' curlist{cc} '_t=f_c_t;']);
end

%% Currency plots
if plot_baseline_curr == 1
    run('plotting/plot_currencies.m');
end

%% Steady state averages
sigma_us_av = mean(sigma_us_TED_t);
sigma_eu_av = mean(sigma_eu_TED_t);
Rm_us_av = mean(Rm_us);
Rm_eu_av = mean(Rm_eu);
riskprm_av = mean(riskprm_t);
Theta_d_us_av = mean(Theta_d_us_t);
Theta_d_eu_av = mean(Theta_d_us_t);

%% ========================================================================
%  DIAGNOSTICS
%  ========================================================================
fprintf('\n=== Filter Diagnostics ===\n');
fprintf('lambda=%.2f, eta=%.2f, ploss=%.2f, iota_us=%.6f\n', lambda_us, eta, ploss_us, iota_us);

% --- m_eff effect (reflects baseline correction + varrho) ---
if use_mu_baseline == 1
    m_eff_us = exp(mu_us) - mubar_us_t - varrho;
    m_eff_eu = exp(mu_eu) - mubar_eu_t - varrho;
    fprintf('\n--- m_eff effect (baseline ON, varrho = %.4f) ---\n', varrho);
else
    m_eff_us = exp(mu_us) - varrho;
    m_eff_eu = exp(mu_eu) - varrho;
    fprintf('\n--- m_eff effect (baseline OFF, varrho = %.4f) ---\n', varrho);
end
fprintf('m_eff_us: min=%.4f, mean=%.4f, max=%.4f   (negative in %d/%d periods)\n', ...
    min(m_eff_us), mean(m_eff_us), max(m_eff_us), sum(m_eff_us<0), T);
fprintf('m_eff_eu: min=%.4f, mean=%.4f, max=%.4f   (negative in %d/%d periods)\n', ...
    min(m_eff_eu), mean(m_eff_eu), max(m_eff_eu), sum(m_eff_eu<0), T);
fprintf('theta_us:  mean=%.4f, max=%.4f   theta_eu:  mean=%.4f, max=%.4f\n', ...
    mean(theta_us_t), max(theta_us_t), mean(theta_eu_t), max(theta_eu_t));
fprintf('sigma_us:  mean=%.4f, std=%.4f   sigma_eu:  mean=%.4f, std=%.4f\n', ...
    mean(sigma_us_t), std(sigma_us_t), mean(sigma_eu_t), std(sigma_eu_t));
fprintf('corr(log sigma, mu): US lvl=%.3f diff=%.3f   EU lvl=%.3f diff=%.3f\n', ...
    corr(log(sigma_us_t), mu_us), corr(diff(log(sigma_us_t)), diff(mu_us)), ...
    corr(log(sigma_eu_t), mu_eu), corr(diff(log(sigma_eu_t)), diff(mu_eu)));
fprintf('DW/FF:     US mean=%.4f   EU mean=%.4f\n', ...
    mean(DW_us_t./max(FF_us_t,eps)), mean(DW_eu_t./max(FF_eu_t,eps)));
clear m_eff_us m_eff_eu;

fprintf('sigma_us: mean=%.4f, min=%.4f, max=%.4f, p05=%.4f, p95=%.4f\n', mean(sigma_us_t), min(sigma_us_t), max(sigma_us_t), prctile(sigma_us_t,5), prctile(sigma_us_t,95));
fprintf('sigma_eu: mean=%.4f, min=%.4f, max=%.4f\n', mean(sigma_eu_t), min(sigma_eu_t), max(sigma_eu_t));
fprintf('Post/Pre(48): %.2f\n', mean(sigma_us_t(end-47:end))/mean(sigma_us_t(1:48)));
fprintf('Implied Std(omega): mean=%.1f%%, p95=%.1f%%\n', mean(sigma_us_t)*sqrt(8)*100, prctile(sigma_us_t,95)*sqrt(8)*100);

% Solver
n_fsolve = sum(sigma_us_TED_flag > 0 & sigma_us_TED_flag ~= 10);
n_fminbnd = sum(sigma_us_TED_flag == 10);
n_fail = sum(sigma_us_TED_flag <= 0);
fprintf('Solver: fsolve=%d (%.1f%%), fminbnd=%d (%.1f%%), fail=%d\n', n_fsolve, n_fsolve/T*100, n_fminbnd, n_fminbnd/T*100, n_fail);

% Residuals
model_ted = zeros(T,1); data_ted = zeros(T,1);
for tt2 = 1:T
    if use_mu_baseline == 1
        mu_arg = exp(mu_us(tt2)) - mubar_us_t(tt2);
    else
        mu_arg = exp(mu_us(tt2));
    end
    model_ted(tt2) = Chi_p_psi(mu_arg, ploss_us, sigma_us_t(tt2), iota_us, lambda_us, eta, matching_type, varrho) * abs_scale;
    if matching_type == 0
        data_ted(tt2) = min_test_us + (TED_s_us_t(tt2) - min(TED_s_us_t)) * abs_scale;
    else
        data_ted(tt2) = TED_s_us_t(tt2) * abs_scale;
    end
end
resid = abs(model_ted - data_ted);
fprintf('Residuals(bps): max=%.2f, mean=%.2f, median=%.2f, pct<0.1=%.1f%%\n', max(resid), mean(resid), median(resid), sum(resid<0.1)/T*100);

% Volume ratios
fprintf('\n--- Volume Ratios (Model vs Data) ---\n');
fprintf('FF:   model=%.2f%%,  data=%.2f%%\n', mean(FF_us_t)*100, mean(FF_n(~isnan(FF_n)))*100);
fprintf('DW:   model=%.2f%%,  data=%.2f%%   (flow; mismatched units)\n', ...
    mean(DW_us_t)*100, mean(DW_n(~isnan(DW_n)))*100);
if exist('DWS_us_t', 'var')
    fprintf('DWS:  model=%.2f%%,  data=%.2f%%   (stock; delta=%.2f)\n', ...
        mean(DWS_us_t)*100, mean(DW_n(~isnan(DW_n)))*100, delta_DWS);
end
fprintf('DW+FF=%.2f%%, DW/FF=%.2f%%  (data: 11.9%%)\n', mean(DW_us_t+FF_us_t)*100, mean(DW_us_t./FF_us_t)*100);
if exist('DWS_us_t', 'var')
    fprintf('DWS/FF=%.2f%%  (model stock / model FF, ratio)\n', mean(DWS_us_t./FF_us_t)*100);
end

% Correlations and std devs
fprintf('\n--- Model vs Data ---\n');
dBP_m = diff(BP_us_t(datesperiod)); dBP_d = diff(Rb_Rm(datesperiod));
dTED_m = diff(TED_us_t(datesperiod)); dTED_d = diff(TED_s_us_t(datesperiod));
dCIP_m = diff(CIP_t(datesperiod)); dCIP_d = diff(cip(datesperiod));
fprintf('         corr(level)  corr(diff)  std_model  std_data  dstd_model  dstd_data\n');
fprintf('BP:      %.3f         %.3f        %.1f       %.1f      %.1f        %.1f\n', ...
    corr(BP_us_t(datesperiod), Rb_Rm(datesperiod)), corr(dBP_m, dBP_d), ...
    std(BP_us_t(datesperiod))*abs_scale, std(Rb_Rm(datesperiod))*abs_scale, ...
    std(dBP_m)*abs_scale, std(dBP_d)*abs_scale);
fprintf('TED:     %.3f         %.3f        %.1f       %.1f      %.1f        %.1f\n', ...
    corr(TED_us_t(datesperiod), TED_s_us_t(datesperiod)), corr(dTED_m, dTED_d), ...
    std(TED_us_t(datesperiod))*abs_scale, std(TED_s_us_t(datesperiod))*abs_scale, ...
    std(dTED_m)*abs_scale, std(dTED_d)*abs_scale);
fprintf('CIP:     %.3f         %.3f        %.1f       %.1f      %.1f        %.1f\n', ...
    corr(CIP_t(datesperiod), cip(datesperiod)), corr(dCIP_m, dCIP_d), ...
    std(CIP_t(datesperiod))*abs_scale, std(cip(datesperiod))*abs_scale, ...
    std(dCIP_m)*abs_scale, std(dCIP_d)*abs_scale);

clear n_fsolve n_fminbnd n_fail model_ted data_ted resid tt2 dBP_m dBP_d dTED_m dTED_d dCIP_m dCIP_d;