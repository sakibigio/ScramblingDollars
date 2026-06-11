%% Update Parameter Values
% continuous_distribution
freq        = 12     ;  % frequency (e.g. 12 monthly) 

%----------- Shock Distribution --------------

sigma_us    = 0.02;    % Was 0.20 (initial guess)
sigma_eu    = 0.015;   % Was 0.15 (initial guess)

%------------- Trading Coefficients ------------
% Lambda depends on matching technology
% matching_type should be defined before calling this script
if ~exist('matching_type', 'var')
    matching_type = 0;  % default to Leontief
    warning('matching_type not defined, defaulting to Leontief (0)');
end

if matching_type == 0
    % Leontief matching
    lambda_us = 3.5;
    lambda_eu = 3.5;
    eta       = 0.5;
    ploss_eu    = 0.75   ;  % probability of an outflow of euro deposits
    ploss_us    = 0.75   ;  % probability of an outflow of dollar deposits
    iota_ss     = 0.1;    % annual, nominal, decimal
elseif matching_type == 1
    %  Cobb-Douglas matching (Run 21 estimates -- common lambda, eta=0.5)
    lambda_us = 1.532;
    lambda_eu = 1.532;
    iota_ss   = 0.055;
    eta       = 0.4341;
    varrho    = 0;
    ploss_eu  = 0.5   ;   % probability of an outflow of euro deposits
    ploss_us  = 0.5   ;   % probability of an outflow of dollar deposits
    % lambda_us   = 0.84 ; % 4.5 nice. eta = 0.625 nice values...
    % lambda_eu   = lambda_us;
    % eta         = 0.7 ;
    % ploss_eu    = 0.5 ;  % probability of an outflow of euro deposits
    % ploss_us    = 0.5 ;  % probability of an outflow of dollar deposits
    % iota_ss     = 0.075;    % annual, nominal, decimal
else
    error('Unknown matching_type: %d. Use 0 (Leontief) or 1 (Cobb-Douglas)', matching_type);
end

varrho      = 0;   % Threshold shift: effective reserve threshold is mu - varrho.
gamma       = 1;

%-----------Policy------------------
% Discount window-IOR spread (single source of truth)


im_eu = imss_eu;
im_us = imss_us;
iota_eu = iota_ss / freq / pi_eu_ss;
iota_us = iota_ss / freq / pi_us_ss;
iw_eu = im_eu + iota_eu * pi_eu_ss;
iw_us = im_us + iota_us * pi_us_ss;

M_eu    = M_eu_ss;
M_us    = M_us_ss;

%----------Supply Demand System------
barb_us        = 0.9;
barb_eu        = 0.9;
barB           = barb_eu+barb_us; 
bard_us        = 1;
bard_eu        = 1;
bard_tot       = bard_us+bard_eu;
bard_tot       = bard_tot;

% ------- Some important ratios ----
nu_us_d        = bard_us/bard_tot;
nu_eu_d        = bard_eu/bard_tot;
nu_b           = barB/bard_tot;

% Steady State Values
Rm_eu  = im_eu/pi_eu_ss;
Rm_us  = im_us/pi_us_ss;

% Transitions One Period Dynamics
M_euus_ratio=M_eu/M_us;

%% Non-financial sector parameters (matching-type independent)
Theta_b = 1;
Theta_d_eu = 1;
Theta_d_us = 1;
epsilon_b = -0.001;
zeta_us = 1000;
zeta_eu = 1000;

% Adjustment to Euro rate (needed for CIP target)
% Targets steady-state CIP = +20 bps/yr to match data (Data_CIP_Moments.tex: 20.6).
% Negative value RAISES Rm_eu by |im_eu_adj| per month. With today's imss_eu=0.999336,
% imss_us=0.999329, solving Rm_eu_ss^12 - Rm_us_ss^12 = 20e-4 gives the value below.
% Earlier value was +0.0006 (lowering Rm_eu by 6 bps/mo), which produced CIP ≈ -70 bps.
im_eu_adj = -0.000161;

%% Session override hook (used by compare_eta_calibrations.sh and estimate_params_4d.m)
% Two override files supported (loaded in this order; later overrides win):
%   _eta_session_override.mat   - just eta (used by eta-sweep driver)
%   _calibration_override.mat   - full calibration snapshot from estimation
% Path is relative to params.m so it works whether called via `params;` or
% `run('functions/params.m')` (which cd's into functions/).
%
% Schema (any subset of these fields):
%   eta_override_val, lambda_us_override_val, lambda_eu_override_val,
%   iota_ss_override_val
tmp_params_dir = fileparts(mfilename('fullpath'));
for tmp_override_file = {'_eta_session_override.mat', '_calibration_override.mat'}
    tmp_override_path = fullfile(tmp_params_dir, '..', tmp_override_file{1});
    if exist(tmp_override_path, 'file') == 2
        tmp_over = load(tmp_override_path);
        if isfield(tmp_over, 'eta_override_val')
            eta = tmp_over.eta_override_val;
            fprintf('  [override %s: eta = %.3f]\n', tmp_override_file{1}, eta);
        end
        if isfield(tmp_over, 'lambda_us_override_val')
            lambda_us = tmp_over.lambda_us_override_val;
            fprintf('  [override %s: lambda_us = %.3f]\n', tmp_override_file{1}, lambda_us);
        end
        if isfield(tmp_over, 'lambda_eu_override_val')
            lambda_eu = tmp_over.lambda_eu_override_val;
            fprintf('  [override %s: lambda_eu = %.3f]\n', tmp_override_file{1}, lambda_eu);
        end
        if isfield(tmp_over, 'iota_ss_override_val')
            iota_ss = tmp_over.iota_ss_override_val;
            % Recompute derived iotas (params.m sets these earlier, but they
            % depend on iota_ss which we just changed).
            iota_eu = iota_ss / freq / pi_eu_ss;
            iota_us = iota_ss / freq / pi_us_ss;
            iw_eu   = im_eu + iota_eu * pi_eu_ss;
            iw_us   = im_us + iota_us * pi_us_ss;
            fprintf('  [override %s: iota_ss = %.4f]\n', tmp_override_file{1}, iota_ss);
        end
        clear tmp_over;
    end
end
clear tmp_override_path tmp_override_file tmp_params_dir;

% Print confirmation
if matching_type == 0
    fprintf('Parameters loaded: Leontief (λ = %.1f)\n', lambda_us);
else
    fprintf('Parameters loaded: Cobb-Douglas (λ = %.1f)\n', lambda_us);
end
