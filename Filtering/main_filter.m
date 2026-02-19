%% Main Filter Script
% Code for Scrambling for Dollars
% (c) Saki Bigio
%
% Pipeline: main_filter.m → markov_estimation.jl → plot_regimes.m

%% Setup
clear; close all;

% Add paths
addpath('functions');
addpath('functions/chi');
addpath('utils');
addpath('data');
addpath('plotting');

%% Model Settings
% Matching function: 0 = Leontief, 1 = Cobb-Douglas
matching_type = 0;

% Print/plot options
printit = 0;
plotdata = 0;
printver = 0;

% Plotting flags (set to 1 to enable)
do_plot_baseline = 1;      % Main filter results
do_plot_diagnostics = 0;   % Diagnostic plots
do_counterfactual = 0;     % Counterfactual analysis (not implemented)
do_sensitivity = 0;        % Sensitivity analysis (not implemented)
plot_baseline_curr = 1;  % Baseline plots for other currencies

% Run markov_estimation.jl automatically after filtering
run_julia = 0;  % Set to 1 to run Julia automatically

%% Load data
load('data/LFX_data.mat');

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
datesperiod = 1:234;
dates = datenum(2001, 1:234, 1);

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
min_test_us = Chi_p_psi(min(exp(mu_us)), ploss_us, min(sigma_vec), iota_us, lambda_us, eta, matching_type) * abs_scale + 0.00012 * abs_scale;
min_test_eu = Chi_p_psi(min(exp(mu_eu)), ploss_eu, min(sigma_vec), iota_eu, lambda_eu, eta, matching_type) * abs_scale + 0.00012 * abs_scale;
max_test = Chi_p_psi(min(exp(mu_us)), ploss_us, max(sigma_vec), iota_us, lambda_us, eta, matching_type) * abs_scale;

% Test lower bound
Ted_test = ones(N_s, 1);
BP_test = ones(N_s, 1);
Psi_test = ones(N_s, 1);
for ss = 1:N_s
    [Ted_test(ss), ~, Psi_test(ss)] = Chi_p_psi(max(exp(mu_us)), ploss_us, sigma_vec(ss), iota_us, lambda_us, eta, matching_type);
    Ted_test(ss) = Ted_test(ss) * 1e4 * 12;
    BP_test(ss) = Echi_m(max(exp(mu_us)), ploss_us, sigma_vec(ss), iota_us, lambda_us, eta, matching_type) * 1e4 * 12;
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
BP_us_target_i = 1;

if matching_type == 0
    matching_name = 'Leontief';
else
    matching_name = 'Cobb-Douglas';
end
fprintf('Starting filter with matching_type = %d (%s)\n', matching_type, matching_name);

for tt = 1:T
    % Setup targets
    BP_us_taget = min_test_us + (Rb_Rm(tt) - min(Rb_Rm)) * abs_scale;
    BP_eu_taget = min_test_eu + (Rb_Rm_eu(tt) - min(Rb_Rm_eu)) * abs_scale;
    TED_us_target = min_test_us + (TED_s_us_t(tt) - min(TED_s_us_t)) * abs_scale; 
    TED_eu_target = min_test_eu + (TED_s_eu_t(tt) - min(TED_s_eu_t)) * abs_scale;

    % Data points
    mu_us_yt = exp(mu_us(tt));
    mu_eu_yt = exp(mu_eu(tt));
    
    % [1] US Sigma - Bond Premia
    if BP_us_target_i == 1
        sigma_us_res = @(sigma) Echi_m(mu_us_yt, ploss_us, sigma, iota_us, lambda_us, eta, matching_type) * 1e4 * 12 - BP_us_taget;
        [sigma_out, ~, exitflag, ~] = fsolve(@(sigma) sigma_us_res(sigma), 1, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
        sigma_us_bp_t(tt) = sigma_out;
        sigma_us_BP_guess = sigma_out;
        sigma_us_flag(tt) = exitflag;
    end

    % [2] US Sigma - TED Target  
    % For Cobb-Douglas, use mu-dependent initial guess to avoid cliff
    if matching_type == 1
        sigma_us_TED_guess = max(sigma_us_TED_guess, mu_us_yt + 0.15);
    end
    sigma_us_res = @(sigma) Chi_p_psi(mu_us_yt, ploss_us, sigma, iota_us, lambda_us, eta, matching_type) * 1e4 * 12 - TED_us_target;
    [sigma_out, ~, exitflag, ~] = fsolve(@(sigma) sigma_us_res(sigma), sigma_us_TED_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
    sigma_us_TED_t(tt) = sigma_out;
    sigma_us_TED_guess = sigma_out;
    sigma_us_TED_flag(tt) = exitflag;
    
    % [3] EU Sigma - Bond Premium
    if BP_eu_target_i == 1
        sigma_eu_res = @(sigma) Echi_m(mu_eu_yt, ploss_us, sigma, iota_eu, lambda_eu, eta, matching_type) * 1e4 * 12 - BP_eu_taget;
        [sigma_out, ~, exitflag, ~] = fsolve(@(sigma) sigma_eu_res(sigma), sigma_eu, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
        sigma_eu_bp_t(tt) = sigma_out;
        sigma_eu_flag(tt) = exitflag;
    end

    % [4] EU Sigma - TED target
    % For Cobb-Douglas, use mu-dependent initial guess to avoid cliff
    if matching_type == 1
        sigma_eu_TED_guess = max(sigma_eu_TED_guess, mu_eu_yt + 0.15);
    end
    sigma_eu_res = @(sigma) Chi_p_psi(mu_eu_yt, ploss_eu, sigma, iota_eu, lambda_eu, eta, matching_type) * 1e4 * 12 - TED_eu_target;
    [sigma_out, ~, exitflag, ~] = fsolve(@(sigma) sigma_eu_res(sigma), sigma_eu_TED_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
    sigma_eu_TED_t(tt) = sigma_out;
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
    Echi_m_us_t(tt) = Echi_m(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type);
    Echi_d_us_t(tt) = Echi_d(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type);
    Chi_p_psi_us_t(tt) = Chi_p_psi(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type);
    Echi_m_eu_t(tt) = Echi_m(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type);
    Echi_d_eu_t(tt) = Echi_d(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type);
    Chi_p_psi_eu_t(tt) = Chi_p_psi(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type);

    % Rest of system
    [~, theta_us_t(tt), psi_us_t(tt), Smin_us_t(tt), DW_us_t(tt), FF_us_t(tt), Q_us_t(tt)] = Chi_sys(mu_us_yt, ploss_us, sigma_us_t(tt), iota_us, lambda_us, eta, matching_type);
    [~, theta_eu_t(tt), psi_eu_t(tt), Smin_eu_t(tt), DW_eu_t(tt), FF_eu_t(tt), Q_eu_t(tt)] = Chi_sys(mu_eu_yt, ploss_eu, sigma_eu_t(tt), iota_eu, lambda_eu, eta, matching_type);

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
sigma_mat = [sigma_us_t sigma_eu_t];
RW_shock = array2table(sigma_mat, 'VariableNames', {'sigma_us', 'sigma_eu'});
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
        mu_eu_yt = exp(mu_eu(tt));
        eval(['target=min_test_eu+TED_s_' curlist{cc} '_t(tt)*abs_scale;']);        
        
        % For Cobb-Douglas, use mu-dependent initial guess
        if matching_type == 1
            sigma_c_guess = max(sigma_eu_t(tt), mu_eu_yt + 0.15);
        else
            sigma_c_guess = sigma_eu_t(tt);
        end
    
        sigma_res = @(sigma) Chi_p_psi(mu_eu_yt, ploss_eu, sigma, iota_eu, lambda_eu, eta, matching_type) * 1e4 * 12 - target;
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
        TED_c_t(tt) = Chi_p_psi(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);

        Echi_m_c_t(tt) = Echi_m(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);
        Echi_d_c_t(tt) = Echi_d(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);
        Chi_p_psi_c_t(tt) = Chi_p_psi(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);
        
        [~, theta_c_t(tt), psi_c_t(tt), Smin_c_t(tt), DW_c_t(tt), FF_c_t(tt), Q_c_t(tt)] = Chi_sys(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);
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

fprintf('\n=== Summary Statistics ===\n');
fprintf('Mean sigma_us: %.4f\n', sigma_us_av);
fprintf('Mean sigma_eu: %.4f\n', sigma_eu_av);
fprintf('Mean risk premium: %.4f (%.2f bps)\n', riskprm_av, riskprm_av * 1e4);
fprintf('\nNext step: Run markov_estimation.jl, then plot_regimes.m\n');
