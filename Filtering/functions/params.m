%% Update Parameter Values
% continuous_distribution
freq        = 12     ;  % frequency (e.g. 12 monthly) 

%----------- Shock Distribution --------------
ploss_eu    = 0.75   ;  % probability of an outflow of euro deposits
ploss_us    = 0.75   ;  % probability of an outflow of dollar deposits
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
elseif matching_type == 1
    % Cobb-Douglas matching
    lambda_us = 1.0 ; % 4.5 nice. eta = 0.625 nice values...
    lambda_eu = 1.0 ;
    eta       = 0.65;
else
    error('Unknown matching_type: %d. Use 0 (Leontief) or 1 (Cobb-Douglas)', matching_type);
end

varrho      = 0.0   ;
gamma       = 1;

%-----------Policy------------------
im_eu = imss_eu;
im_us = imss_us;
iw_eu = iwss_eu;
iw_us = iwss_us;
iota_eu = (iw_eu-im_eu)/pi_eu_ss;
iota_us = (iw_us-im_us)/pi_us_ss;

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
