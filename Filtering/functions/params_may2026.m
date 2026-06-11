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
    %  Cobb-Douglas matching
    lambda_us = 1.0 ; % 4.5 nice. eta = 0.625 nice values...
    lambda_eu = 1.0 ;
    eta       = 0.70 ;
    ploss_eu  = 0.5   ;  % probability of an outflow of euro deposits
    ploss_us  = 0.5   ;  % probability of an outflow of dollar deposits
    iota_ss   = 0.065 ;    % annual, nominal, decimal
    % lambda_us   = 0.84 ; % 4.5 nice. eta = 0.625 nice values...
    % lambda_eu   = lambda_us;
    % eta         = 0.7 ;
    % ploss_eu    = 0.5 ;  % probability of an outflow of euro deposits
    % ploss_us    = 0.5 ;  % probability of an outflow of dollar deposits
    % iota_ss     = 0.075;    % annual, nominal, decimal
else
    error('Unknown matching_type: %d. Use 0 (Leontief) or 1 (Cobb-Douglas)', matching_type);
end

varrho      = 0.0   ;   % Threshold shift: effective reserve threshold is mu - varrho.
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
im_eu_adj = 0.0006;

% Print confirmation
if matching_type == 0
    fprintf('Parameters loaded: Leontief (λ = %.1f)\n', lambda_us);
else
    fprintf('Parameters loaded: Cobb-Douglas (λ = %.1f)\n', lambda_us);
end
