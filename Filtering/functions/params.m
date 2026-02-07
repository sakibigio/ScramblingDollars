%% Update Parameter Values
% continuous_distribution
freq        = 12     ;  % fequency (e.g. 25 biweekly) 

%----------- Shock Distribution --------------
ploss_eu    = 0.25   ;  % frac loss Euro deposits
ploss_us    = 0.25   ;  % frac loss $ deposits
% sigma_eu    = 0.168   ;  % scale coefficient 0.168
% sigma_us    = 0.168   ;  % scale coefficient 0.168
% sigma_eu    = 0.1326   ;  % scale coefficient 0.168
% sigma_us    = 0.245   ;  % scale coefficient 0.168
% sigma_eu    = 0.133479   ;  % scale coefficient 0.168
% sigma_us    = 0.24496   ;  % scale coefficient 0.168
sigma_eu   = 0.015   ;  % scale coefficient 0.168
sigma_us   = 0.02   ;  % scale coefficient 0.168


%------------- Trading Coefficients ------------
lambda_us   = 3.5   ;  % efficiency Euro market
lambda_eu   = 3.5   ;  % efficiency $ market
varrho      = 0.0   ;
gamma       = 1;
eta         = 0.5;

%-----------Policy------------------
% im_eu  = 1.02^(1/freq);
% im_us  = 1.02^(1/freq);
% iw_eu  = 1.12^(1/freq);
% iw_us  = 1.12^(1/freq);

% im_eu = (1+mean_im_eu/100)^(1/freq);
% im_us = (1+mean_im_us/100)^(1/freq);
% iw_eu  = (1+mean_im_eu/100+0.1)^(1/freq);
% iw_us  = (1+mean_im_us/100+0.1)^(1/freq);

im_eu = imss_eu;
im_us = imss_us;
iw_eu = imss_eu+(iwss_eu-imss_eu)*1;
iw_us = imss_us+(iwss_us-imss_us)*1;
iota_eu = (iw_eu-im_eu)/pi_eu_ss;
iota_us = (iw_us-im_us)/pi_us_ss;
%M_eu    = 0.5;
%M_us    = 0.5;
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
% Rw_eu  = iw_eu;
% Rw_us  = iw_us;

% Transitions One Period Dynamics
M_euus_ratio=M_eu/M_us;

% Initial guess of loan and deposit demand size (will calibrate at line 63)
Theta_b = 1;
Theta_d_eu = 1^(1/35); 
Theta_d_us = 1^(1/35);

% Loan demand elasticity parameter (elasticity=1/epsilon_b)
epsilon_b = -1/35;

% Deposit demand elasticity parameter (elasticity=1/zeta)
zeta_us = 1;
zeta_eu = 1;

