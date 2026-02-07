%% Load data
% Code for Scrambbling for Dolars
% (c) Saki Bigio
%%%%%%%%%%%%
% Acronyms
% _t variables: deduced from model
% _yt variables within loop

%% Updates
clear; close all;
load('LFX_data3.mat'); % Load Latest Data
printit=0;
plotdata=1;
printver=0;
if plotdata==1
    LFX_plotdata;
end

%% Update paraemters [fix block] 
% Set other predetermined parameters 
load dynare_calibration_param.mat;

% Scale: Report Variables in BPS
abs_scale=12e4;

% List of Currencies
curlist = {'au','ca','jp','nz','no','sw','ch','uk'};
conlist = {'AUD','CAD','JPY','NZD','NOK','SWK','CHF','GBP'};
CURRlist = {'EUR','AUD','CAD','JPY','NZD','NOK','SWK','CHF','GBP'};

% OMO share in MBS 
share_us=1; 
share_eu=1;

% Steady state inflation rate is 1 
pi_eu_ss = 1;
pi_us_ss = 1;

% Update rest of parameters
LFX_params_cd_v3;
% zeta_us = 1/35;
% zeta_eu = 1/35;
% lambda_us=2.5;
% lambda_eu=2.5;
% lambda_us=lambda_us;
% lambda_us=lambda_us*1.5;
% lambda_eu=lambda_eu*1.5;
% ploss_us=0.75; ploss_eu=0.75;
epsilon_b=-0.5;
zeta_us=-epsilon_b; zeta_eu=-epsilon_b;
Theta_b=3;

% Plot Properties - Date Variables
datesperiod = 1:234;
dates=datenum(2001,1:234,1);

% Plot Formats
FSize=20;
formataxis   = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize, 'Box','Off','PlotBoxAspectRatio', [1 0.75 1]);
label_x      = @(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize,'Interpreter','latex');
label_y      = @(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', FSize,'interpreter','latex');
yLimits = get(gca, 'YLim');
desiredDecimalPlaces=1;
%customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
%set(ax, 'YTick', customYTicks);
%desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
%xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
%grid on;

%% Load Eq Equations & Functions
LFX_nt_0e_eqs_3;

%% Tests for Lower Bounds of Interbank Variables
sigma_vec=(0.01:0.01:10);
N_s=length(sigma_vec);
min_ted=min(TED_s_us_t*abs_scale);
max_ted=max(TED_s_us_t*abs_scale);
min_test_us=Chi_p_psi(min(exp(mu_us)),ploss_us,min(sigma_vec),iota_us,lambda_us,eta)*abs_scale+0.00012*abs_scale;
min_test_eu=Chi_p_psi(min(exp(mu_eu)),ploss_eu,min(sigma_vec),iota_eu,lambda_eu,eta)*abs_scale+0.00012*abs_scale;
max_test=Chi_p_psi(min(exp(mu_us)),ploss_us,max(sigma_vec),iota_us,lambda_us,eta)*abs_scale;

% Test Lower Bound
Ted_test=ones(N_s,1);
BP_test=ones(N_s,1);
Psi_test=ones(N_s,1);
for ss=1:N_s
    [Ted_test(ss),~,Psi_test(ss)]=Chi_p_psi(max(exp(mu_us)),ploss_us,sigma_vec(ss),iota_us,lambda_us,eta);
    Ted_test(ss)=Ted_test(ss)*1e4*12;
    BP_test(ss)=Echi_m(max(exp(mu_us)),ploss_us,sigma_vec(ss),iota_us,lambda_us,eta)*1e4*12;
end
figure('Name','Ted test','NumberTitle','off'); 
plot(sigma_vec,Ted_test); hold on;
plot(sigma_vec,BP_test);
scatter(0*TED_s_us_t,TED_s_us_t*abs_scale,'filled'); hold on;
scatter(0*TED_s_eu_t,TED_s_eu_t*abs_scale,'filled'); hold on;

figure('Name','Psi test','NumberTitle','off'); 
plot(sigma_vec,Psi_test);


%% Data Analogues
load LFX_data3.mat;

% initialize variables
endopath=[mu_us mu_eu Rb_Rm cip];
exopath =[im_eu im_us M_us M_eu ];

%% Initializing Vectors - Pre-allocation
T=length(endopath);

% Price Variables
sigma_us_t = zeros(T,1);
sigma_us_bp_t= zeros(T,1);
sigma_us_TED_t= zeros(T,1);
sigma_eu_TED_t= zeros(T,1);
Echi_m_us_t= zeros(T,1);
Echi_d_us_t= zeros(T,1);
Chi_p_psi_us_t= zeros(T,1);
sigma_eu_t = zeros(T,1);
sigma_eu_bp_t = zeros(T,1);
Echi_m_eu_t= zeros(T,1);
Echi_d_eu_t= zeros(T,1);
Chi_p_psi_eu_t= zeros(T,1);
%Echi_m_eu_bp_t= zeros(T,1);
% Echi_d_eu_bp_t= zeros(T,1);
% BP_us_t    = zeros(T,1);
% BP_eu_t    = zeros(T,1);
% Rd_us_t    = zeros(T,1);
% Rd_eu_t    = zeros(T,1);
% Rb_us_t    = zeros(T,1);
% Rb_eu_t    = zeros(T,1);
% TED_us_t   = zeros(T,1);
% TED_eu_t   = zeros(T,1);
% riskprm_t  = zeros(T,1);

% Quantity Variables
% nu_t       = zeros(T,1);
% d_us_t     = zeros(T,1);
% Theta_d_us_t=zeros(T,1);
% Theta_d_eu_t=zeros(T,1);
% d_eu_t     = zeros(T,1);
% b_t        = zeros(T,1);
% mu_us_t    = zeros(T,1);
% mu_eu_t    = zeros(T,1);
% p_eu_t     = zeros(T,1);
% inv_e_t    = zeros(T,1);
% p_us_t     = zeros(T,1);
% e_euus_t   = zeros(T,1);
% M_eu_t     = zeros(T,1);
% nu_t  = zeros(T,1);
% M_us_t= zeros(T,1);
% f_t   = zeros(T,1);

% UIP and CIP deviations Variables
% UIP_t     = zeros(T,1);
% CIP_t     = zeros(T,1);

% Interbank variables
theta_us_t= zeros(T,1);
psi_us_t  = zeros(T,1);
Smin_us_t = zeros(T,1);
DW_us_t   = zeros(T,1);
FF_us_t   = zeros(T,1);
Q_us_t    = zeros(T,1);
theta_eu_t= zeros(T,1);
psi_eu_t  = zeros(T,1);
Smin_eu_t = zeros(T,1);
DW_eu_t   = zeros(T,1);
FF_eu_t   = zeros(T,1);
Q_eu_t    = zeros(T,1);

% Code Diagnostics
sigma_us_flag     = zeros(T,1);
sigma_eu_flag     = zeros(T,1);
sigma_us_TED_flag = zeros(T,1);
sigma_eu_TED_flag = zeros(T,1);
rest_flag     = zeros(T,1);

% Start Guess
sigma_us_BP_guess=sigma_us;
sigma_eu_BP_guess=sigma_eu;
sigma_us_TED_guess=sigma_us;
sigma_eu_TED_guess=sigma_eu;

%% Main Filtering Step
TED_target_i=1; % turn of for other targets
TED_eu_target_i=1; 
BP_eu_target_i=0;
BP_us_target_i=1;

for tt=1:T
   % Other Targets 
   % target = [Rb_Rm(tt)*abs_scale; exp(mu_us(tt)); exp(mu_eu(tt)); cip(tt)*abs_scale; Rb_Rm_eu(tt)*abs_scale; min_test_us+(TED_s_us_t(tt)-min(TED_s_us_t))*abs_scale; min_test_eu+(TED_s_eu_t(tt)-min(TED_s_eu_t))*abs_scale;];
   % Distance=@(Rb_us,mu_us,mu_eu,riskprm) [(Rb_us-Rm_us)*abs_scale-target(1); (mu_us)-target(2); (mu_eu)-target(3); (Rm_eu-Rm_us+riskprm)*abs_scale-target(4)];

    % Setup baseline targets
    BP_us_taget=min_test_us+(Rb_Rm(tt)-min(Rb_Rm))*abs_scale;
    BP_eu_taget=min_test_eu+(Rb_Rm_eu(tt)-min(Rb_Rm_eu))*abs_scale;
    TED_us_target=min_test_us+(TED_s_us_t(tt)-min(TED_s_us_t))*abs_scale; 
    TED_eu_target=min_test_eu+(TED_s_eu_t(tt)-min(TED_s_eu_t))*abs_scale;

    % Update Data Points
    mu_us_yt=exp(mu_us(tt));
    mu_eu_yt=exp(mu_eu(tt));
    
    % [1] Backing Out US Sigma
    % Target #2: US Sigma - Bond Premia
    if BP_us_target_i==1
        sigma_us_res=@(sigma) Echi_m(mu_us_yt,ploss_us,sigma,iota_us,lambda_us,eta)*1e4*12-BP_us_taget;
        [sigma_out,~,exitflag,~]=fsolve(@(sigma) sigma_us_res(sigma),1,optimoptions('fsolve','Display','iter','TolFun',1e-12,'MaxFunEval',1e9,'MaxIter',1e6));
        sigma_us_bp_t(tt)=sigma_out; sigma_us_BP_guess= sigma_out; clear sigma_out;
        sigma_us_flag(tt)     = exitflag;
    end

    % Target #1: US Sigma - TED Target  
    sigma_us_res=@(sigma) Chi_p_psi(mu_us_yt,ploss_us,sigma,iota_us,lambda_us,eta)*1e4*12-TED_us_target;
    [sigma_out,~,exitflag,~]=fsolve(@(sigma) sigma_us_res(sigma),sigma_us_TED_guess,optimoptions('fsolve','Display','iter','TolFun',1e-12,'MaxFunEval',1e9,'MaxIter',1e6));
    sigma_us_TED_t(tt)=sigma_out; sigma_us_TED_guess= sigma_out; clear sigma_out;
    sigma_us_TED_flag(tt)     = exitflag;
    
    % [2] Backing Out EU Sigma
    % Target #2: EU Sigma - Bond Premium Target
    if BP_eu_target_i==1
        sigma_eu_res=@(sigma) Echi_m(mu_eu_yt,ploss_us,sigma,iota_eu,lambda_eu,eta)*1e4*12-BP_eu_taget;
        [sigma_out,~,exitflag,~]=fsolve(@(sigma) sigma_eu_res(sigma),sigma_eu,optimoptions('fsolve','Display','iter','TolFun',1e-12,'MaxFunEval',1e9,'MaxIter',1e6));
        sigma_eu_bp_t(tt)=sigma_out; clear sigma_out;
        sigma_eu_flag(tt)     = exitflag;
    end

    % EU Sigma - CIP Target
    % sigma_eu_res=@(sigma_eu) Echi_m(mu_eu_yt,ploss_eu,sigma_eu,iota_eu,lambda_eu)*1e4*12-(target(1)-target(4));
    % [sigma_out,~,exitflag,~]=fsolve(@(sigma) sigma_eu_res(sigma),sigma_us_t(tt),optimoptions('fsolve','Display','iter','TolFun',1e-12,'MaxFunEval',1e9,'MaxIter',1e6));
    % sigma_eu_t(tt)=sigma_out; clear sigma_out;
    % sigma_eu_flag(tt)     = exitflag;

    % Target #1: EU Sigma - TED target
    sigma_eu_res=@(sigma) Chi_p_psi(mu_eu_yt,ploss_eu,sigma,iota_eu,lambda_eu,eta)*1e4*12-TED_eu_target;
    [sigma_out,~,exitflag,~]=fsolve(@(sigma) sigma_eu_res(sigma),sigma_eu_TED_guess,optimoptions('fsolve','Display','iter','TolFun',1e-12,'MaxFunEval',1e9,'MaxIter',1e6));
    sigma_eu_TED_t(tt)=sigma_out; % sigma_eu_TED_guess= sigma_out; clear sigma_out;
    sigma_eu_TED_flag(tt)     = exitflag; 

    % Update Targets    
    if  TED_target_i==1
        sigma_us_t(tt)=sigma_us_TED_t(tt); 
        sigma_us_TED_guess=sigma_us_TED_t(tt); 
        sigma_eu_t(tt)=sigma_eu_TED_t(tt);
        sigma_eu_TED_guess=sigma_eu_TED_t(tt);
    elseif BP_eu_target==1
        sigma_us_t(tt)=sigma_us_bp_t(tt);
        sigma_eu_t(tt)=sigma_eu_bp_t(tt);
      %  TED_eu_t(tt)=Chi_p_psi(mu_eu_yt,ploss_eu,sigma_eu_TED_t(tt),iota_eu,lambda_eu,eta);
    end

    % [3] Update Interbank Variables
    Echi_m_us_t(tt)=Echi_m(mu_us_yt,ploss_us,sigma_us_t(tt),iota_us,lambda_us,eta);
    Echi_d_us_t(tt)=Echi_d(mu_us_yt,ploss_us,sigma_us_t(tt),iota_us,lambda_us,eta);
    Chi_p_psi_us_t(tt)=Chi_p_psi(mu_us_yt,ploss_us,sigma_us_t(tt),iota_us,lambda_us,eta);
    Echi_m_eu_t(tt)=Echi_m(mu_eu_yt,ploss_eu,sigma_eu_t(tt),iota_eu,lambda_eu,eta);
    Echi_d_eu_t(tt)=Echi_d(mu_eu_yt,ploss_eu,sigma_eu_t(tt),iota_eu,lambda_eu,eta);
    Chi_p_psi_eu_t(tt)=Chi_p_psi(mu_eu_yt,ploss_eu,sigma_eu_t(tt),iota_eu,lambda_eu,eta);

    % Rest of the System
    [~,theta_us_t(tt),psi_us_t(tt),Smin_us_t(tt),DW_us_t(tt),FF_us_t(tt),Q_us_t(tt)]=Chi_sys(mu_us_yt,ploss_us,sigma_us_t(tt),iota_us,lambda_us,eta);
    [~,theta_eu_t(tt),psi_eu_t(tt),Smin_eu_t(tt),DW_eu_t(tt),FF_eu_t(tt),Q_eu_t(tt)]=Chi_sys(mu_eu_yt,ploss_eu,sigma_eu_t(tt),iota_eu,lambda_eu,eta);

    % Destroy stuff
    clear mu_us_yt mu_eu_yt TED_us_target TED_eu_target;
end
% Construction of variables
Rm_us=exp(im_us);
Rm_eu=exp(im_eu);
inv_e_yt=exp(inv_e);   
Rb_Rm_yt=Rb_Rm;
    
%% Deducing Other Variables
BP_us_t=Echi_m_us_t;
TED_us_t=Chi_p_psi_us_t;
Rd_us_t=Rm_us+Echi_m_us_t+Echi_d_us_t;
Rb_us_t=Rm_us+Echi_m_us_t;
BP_eu_t=Echi_m_eu_t;
TED_eu_t=Chi_p_psi_eu_t;
Rd_eu_t=Rm_eu+Echi_m_eu_t+Echi_d_eu_t;
Rb_eu_t=Rm_eu+Echi_m_eu_t;

% Risk Premium 
riskprm_t=Rb_us_t./Rb_eu_t-1;

% UIP Deviation Model
UIP_t=Rm_eu-Rm_us;

% CIP 
CIP_t=UIP_t+Rm_eu.*riskprm_t;

% Forward FX
f_t=(1+riskprm_t)./exp(inv_e);
CIP_check_t=exp(im_eu).*f_t.*exp(inv_e)-exp(im_us)-CIP_t;

%% Computing Quantity Variables
inv_e_t=exp(inv_e);

%     Rb_Rm_yt=Rb_Rm;
mu_us_t=exp(mu_us);
mu_eu_t=exp(mu_eu);
M_us_t=exp(M_us);
M_eu_t=exp(M_eu);

% Independent block, should have its own code
nu_t=(M_eu_t./M_us_t).*inv_e_t.*mu_us_t./mu_eu_t;
b_t=Theta_b*(Rm_us+Echi_m_us_t).^(1/epsilon_b);
d_us_t=b_t./((1-mu_us_t)+nu_t.*(1-mu_eu_t));
d_eu_t=nu_t.*d_us_t;

% Checks
budget_check=b_t+(mu_us_t-1).*d_us_t+(mu_eu_t-1).*d_eu_t;
Rb_check=Rm_us+Rb_Rm-((nu_t.*(1-mu_eu_t)+1-mu_us_t).*d_us_t/Theta_b).^epsilon_b;

% Price System  
p_us_t  = M_us_t./(d_us_t.*mu_us_t)    ;
p_eu_t  = M_eu_t./(d_eu_t.*mu_eu_t)    ;
e_euus_t= 1./inv_e_t;
M_eu_t_check  = p_eu_t.*(d_eu_t.*mu_eu_t)-M_eu_t         ;
M_us_t_check  = p_us_t.*(d_us_t.*mu_us_t)-M_us_t         ;
e_check=p_us_t./p_eu_t-inv_e_t;
%     p_eu_t(tt)  = M_us(tt)/(d_us_t(tt)*mu_us_t(tt))/inv_e_yt    ;
%     inv_e_t(tt) = M_us(tt)/(d_us_t(tt)*mu_us_t(tt))/p_eu_t(tt)  ;

% Backout Deposit Funding Shocks
N_pre=72;
MBS_us_ss  =mean(d_us_t(1:N_pre)); % ---> check number from SS
MBS_eu_ss  =mean(d_eu_t(1:N_pre));
share_us=1;
share_eu=1;
Theta_d_us_t=(d_us_t-share_us*(mu_us_t.*d_us_t)+MBS_us_ss)./(Rd_us_t).^(1/zeta_us);
Theta_d_eu_t=(d_eu_t-share_eu*(mu_eu_t.*d_eu_t)+MBS_eu_ss)./(Rd_eu_t).^(1/zeta_eu);

% Differential of Risk-Premium vis-a-vis rates
riskprm_shock_t=riskprm_t(datesperiod)-(exp(im_us(datesperiod))-exp(im_eu(datesperiod)));

%% Saving Outcome for External Analysis
sigma_mat=[sigma_us_t sigma_eu_t];
RW_shock = array2table(sigma_mat, 'VariableNames', {'sigma_us', 'sigma_eu'});
writetable(RW_shock, 'RW_shock.csv');

%% Plotting Results
close all;

plot_baseline=0;
if plot_baseline==1
    RW_filter_plot_baseline;
end

run_counterfactual=0;
if run_counterfactual==1
    LFX_rw_filter_counterfactual;
end

run_sensitivity=0;
if run_sensitivity==1
    LFX_sensitivity;
end

%% Diagnostics - Shocks
% Plot Data
plot_diagnostics=0;
if plot_diagnostics==1
    LFX_rw_shock_diagnostics;
end

% 
% 
% 
% 
% %% Baseline Regressoins
% % Do log-regressions
% log_sigma_us_t  = log(sigma_us_rw);
% % d_log_sigma_us_t= diff(log_sigma_us_t);
% Gsigma_us       = mean(log_sigma_us_t(estimate_range));
% log_sigma_us_t = log_sigma_us_t-Gsigma_us;
% lag_aux=log_sigma_us_t(estimate_range(1:end-1));
% X = [ones(length(lag_aux),1) lag_aux];
% [B,~,R,~,~] = regress(log_sigma_us_t(estimate_range(2:end)),X);
% rho_sigma_us   = B(2);
% Sigma_sigma_us = std(R);
% 
% fprintf('******\n');
% fprintf('sigma us\n');
% fprintf('AR(1) coeff: %f\n', B(2));
% 
% log_sigma_eu_t  = log(sigma_eu_rw);
% % d_log_sigma_us_t= diff(log_sigma_us_t);
% Gsigma_eu       = mean(log_sigma_eu_t(estimate_range));
% log_sigma_eu_t = log_sigma_eu_t-Gsigma_eu;
% lag_aux=log_sigma_eu_t(estimate_range(1:end-1));
% X = [ones(length(lag_aux),1) lag_aux];
% [B,~,R,~,~] = regress(log_sigma_eu_t(estimate_range(2:end)),X);
% rho_sigma_eu   = B(2);
% Sigma_sigma_eu = std(R);   
% fprintf('******\n');
% fprintf('sigma eu\n');
% fprintf('AR(1) coeff: %f\n', B(2));
% 
% % Do log-regressions
% log_Theta_d_us_t  = log(Theta_d_us_rw);
% % d_log_sigma_us_t= diff(log_sigma_us_t);
% GTheta_d_us       = mean(log_Theta_d_us_t(estimate_range));
% log_Theta_d_us_t = log_Theta_d_us_t-GTheta_d_us;
% lag_aux=log_Theta_d_us_t(estimate_range(1:end-1));
% X = [ones(length(lag_aux),1) lag_aux];
% [B,~,R,~,~] = regress(log_Theta_d_us_t(estimate_range(2:end)),X);
% rho_Theta_d_us   = B(2);
% Sigma_Theta_d_us = std(R);
% fprintf('******\n');
% fprintf('Theta d us\n');
% fprintf('AR(1) coeff: %f\n', B(2));
% 
% log_Theta_d_eu_t  = log(Theta_d_eu_rw);
% % d_log_sigma_us_t= diff(log_sigma_us_t);
% GTheta_d_eu       = mean(log_Theta_d_eu_t(estimate_range));
% log_Theta_d_eu_t = log_Theta_d_eu_t-GTheta_d_eu;
% lag_aux=log_Theta_d_eu_t(estimate_range(1:end-1));
% X = [ones(length(lag_aux),1) lag_aux];
% [B,~,R,~,~] = regress(log_Theta_d_eu_t(estimate_range(2:end)),X);
% rho_Theta_d_eu   = B(2);
% Sigma_Theta_d_eu = std(R);
% fprintf('******\n');
% fprintf('Theta d eu\n');
% fprintf('AR(1) coeff: %f\n', B(2));
% 
% Griskprm_rw       = mean(riskprm_rw(estimate_range));
% friskprm_rw_t     = riskprm_rw-Griskprm_rw;
% lag_aux=friskprm_rw_t(estimate_range(1:end-1));
% X = [ones(length(lag_aux),1) lag_aux];
% [B,~,R,~,~] = regress(friskprm_rw_t(estimate_range(2:end)),X);
% rho_riskprm   = B(2);
% Sigma_riskprm = std(R);
% fprintf('******\n');
% fprintf('risk premia\n');
% fprintf('AR(1) coeff: %f\n', B(2));
% 
% % Random Walk Tests
% % Is sigma_us a random walk?
% [h, pValue, stat, cValue] = adftest(sigma_us_rw, 'model', 'TS', 'lags', 2);
% fprintf('******\n');
% 
% fprintf('sigma us\n');
% % Display the results
% if h == 0
%     fprintf('The null hypothesis of a unit root (random walk) cannot be rejected.\n');
% else
%     fprintf('The null hypothesis of a unit root (random walk) is rejected.\n');
% end
% 
% fprintf('p-Value: %f\n', pValue);
% fprintf('Test Statistic: %f\n', stat);
% fprintf('Critical Value: %f\n', cValue);
% 
% % Is sigma_eu a random walk?
% fprintf('******\n');
% 
% fprintf('sigma eu\n');
% [h, pValue, stat, cValue] = adftest(sigma_eu_rw, 'model', 'TS', 'lags', 2);
% 
% % Display the results
% if h == 0
%     fprintf('The null hypothesis of a unit root (random walk) cannot be rejected.\n');
% else
%     fprintf('The null hypothesis of a unit root (random walk) is rejected.\n');
% end
% 
% fprintf('p-Value: %f\n', pValue);
% fprintf('Test Statistic: %f\n', stat);
% fprintf('Critical Value: %f\n', cValue);
% 
% % Is sigma_eu a random walk?
% fprintf('******\n');
% 
% fprintf('Theta d us\n');
% [h, pValue, stat, cValue] = adftest(Theta_d_us_rw, 'model', 'TS', 'lags', 2);
% 
% % Display the results
% if h == 0
%     fprintf('The null hypothesis of a unit root (random walk) cannot be rejected.\n');
% else
%     fprintf('The null hypothesis of a unit root (random walk) is rejected.\n');
% end
% 
% fprintf('p-Value: %f\n', pValue);
% fprintf('Test Statistic: %f\n', stat);
% fprintf('Critical Value: %f\n', cValue);
% 
% % Is sigma_eu a random walk?
% fprintf('******\n');
% 
% fprintf('Theta d eu\n');
% [h, pValue, stat, cValue] = adftest(Theta_d_eu_rw, 'model', 'TS', 'lags', 2);
% 
% % Display the results
% if h == 0
%     fprintf('The null hypothesis of a unit root (random walk) cannot be rejected.\n');
% else
%     fprintf('The null hypothesis of a unit root (random walk) is rejected.\n');
% end
% 
% fprintf('p-Value: %f\n', pValue);
% fprintf('Test Statistic: %f\n', stat);
% fprintf('Critical Value: %f\n', cValue);
% 
% % Is sigma_eu a random walk?
% % fprintf('******\n');
% fprintf('Risk premia');
% [h, pValue, stat, cValue] = adftest(riskprm_rw, 'model', 'TS', 'lags', 2);
% 
% % Display the results
% if h == 0
%     fprintf('The null hypothesis of a unit root (random walk) cannot be rejected.\n');
% else
%     fprintf('The null hypothesis of a unit root (random walk) is rejected.\n');
% end
% 
% fprintf('p-Value: %f\n', pValue);
% fprintf('Test Statistic: %f\n', stat);
% fprintf('Critical Value: %f\n', cValue);
% 
% % Is sigma_eu a random walk?
% 
% fprintf('******\n');
% fprintf('Fx\n');
% [h, pValue, stat, cValue] = adftest(inv_e, 'model', 'TS', 'lags', 2);
% 
% % Display the results
% if h == 0
%     fprintf('The null hypothesis of a unit root (random walk) cannot be rejected.\n');
% else
%     fprintf('The null hypothesis of a unit root (random walk) is rejected.\n');
% end
% 
% fprintf('p-Value: %f\n', pValue);
% fprintf('Test Statistic: %f\n', stat);
% fprintf('Critical Value: %f\n', cValue);



%% Sigma in Other Currencies 
% preallocation
sigma_c_TED_t = zeros(T,1); sigma_c_TED_flag=zeros(T,1);
Echi_m_c_t= zeros(T,1);
Echi_d_c_t= zeros(T,1);
Chi_p_psi_c_t= zeros(T,1);
Echi_m_eu_bp_t= zeros(T,1);
Echi_d_eu_bp_t= zeros(T,1);
BP_c_t    = zeros(T,1);
Rd_c_t    = zeros(T,1);
TED_c_t   = zeros(T,1);
theta_c_t = zeros(T,1);
psi_c_t   = zeros(T,1);
Smin_c_t  = zeros(T,1);
DW_c_t    = zeros(T,1);
FF_c_t    = zeros(T,1);
Q_c_t     = zeros(T,1);

% By currency loops
for cc=1:numel(curlist)
    % Update rates
    eval(['Rm_' curlist{cc} '=exp(im_' curlist{cc} ');']);
    eval(['Rm_c=Rm_' curlist{cc} ';']); % Stopped here...save rates...
    for tt=1:T
        mu_eu_yt=exp(mu_eu(tt));
        eval(['target=min_test_eu+TED_s_' curlist{cc} '_t(tt)*abs_scale;']);        
    
        % Target #1: EU Sigma - TED target
        sigma_res=@(sigma) Chi_p_psi(mu_eu_yt,ploss_eu,sigma,iota_eu,lambda_eu,eta)*1e4*12-target;
        [sigma_out,~,exitflag,~]=fsolve(@(sigma) sigma_res(sigma),sigma_eu_t(tt),optimoptions('fsolve','Display','iter','TolFun',1e-12,'MaxFunEval',1e9,'MaxIter',1e6));
        if exitflag>0
            sigma_c_TED_t(tt)=sigma_out; % sigma_eu_TED_guess= sigma_out; clear sigma_out;
        else
            sigma_c_TED_t(tt)=sigma_c_TED_t(tt-1);
        end
        sigma_c_TED_flag(tt)     = exitflag;

        TED_c_t(tt)=Chi_p_psi(mu_eu_yt,ploss_eu,sigma_c_TED_t(tt),iota_eu,lambda_eu,eta);

        % Collecting rest
        Echi_m_c_t(tt)=Echi_m(mu_eu_yt,ploss_eu,sigma_c_TED_t(tt),iota_eu,lambda_eu,eta);
        Echi_d_c_t(tt)=Echi_d(mu_eu_yt,ploss_eu,sigma_c_TED_t(tt),iota_eu,lambda_eu,eta);
        Chi_p_psi_c_t(tt)=Chi_p_psi(mu_eu_yt,ploss_eu,sigma_c_TED_t(tt),iota_eu,lambda_eu,eta);
        
        % Rest of the System
        [~,theta_c_t(tt),psi_c_t(tt),Smin_c_t(tt),DW_c_t(tt),FF_c_t(tt),Q_c_t(tt)]=Chi_sys(mu_eu_yt,ploss_eu,sigma_c_TED_t(tt),iota_eu,lambda_eu,eta);
    end
    % Destroy stuff
    clear mu_eu_yt target;

    % Construct Other Variables
    BP_c_t=Echi_m_c_t;
    Rb_c_t=Rm_c+Echi_m_c_t;
    Rd_c_t=Rm_c+Echi_m_c_t+Echi_d_c_t;
    UIP_c_t=Rm_c-Rm_us;
    riskprm_c_t=(Rb_us_t)./(Rb_c_t)-1;
    CIP_c_t=UIP_c_t+Rm_c.*riskprm_c_t;
    
    % Construct Counterfactual FX
    p_c_t = (M_eu./(Rd_c_t.^(1/zeta_eu)))./mu_eu_t;
    inv_e_c_t = p_us_t./p_c_t;

    % Forward FX
    f_c_t=(1+riskprm_c_t)./exp(inv_e_c_t);

    % UIP Deviation Model
    eval(['sigma_' curlist{cc} '_TED_flag'  '=sigma_c_TED_flag;']);
    eval(['sigma_' curlist{cc} '_t'  '=sigma_c_TED_t;']);
    eval(['TED_' curlist{cc} '_t'  '=TED_c_t;']);
    eval(['BP_' curlist{cc} '_t'  '=BP_c_t;']);
    eval(['riskprm_' curlist{cc} '_t=riskprm_c_t;']);
    eval(['UIP_' curlist{cc} '_t=UIP_c_t;']);
    eval(['CIP_' curlist{cc} '_t=CIP_c_t;']);
    eval(['inv_e_' curlist{cc} '_t=inv_e_c_t;']);
    eval(['p_' curlist{cc} '_t=p_c_t;']);
    eval(['f_' curlist{cc} '_t=f_c_t;']);
end

%% Plots for Other Currencies
plot_baseline_curr=0;
if plot_baseline_curr==1
    RW_filter_plot_baseline_curr;
end

%% Data Plots
figure
plot(dates(datesperiod),(M_us(datesperiod)./M_eu(datesperiod)-1)*12e2,'LineWidth',3); hold on;
grid on; axis tight;
datetick('x','yyyy-mm','keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title(['Relative Money Supply']);
set(h,'interpreter','latex','fontsize',20);

figure
plot(dates(datesperiod(1:end-1)),(p_eu_t(datesperiod(2:end))./p_eu_t(datesperiod(1:end-1))-1)*12e2,'LineWidth',3); hold on;
plot(dates(datesperiod(1:end-1)),(p_us_t(datesperiod(2:end))./p_us_t(datesperiod(1:end-1))-1)*12e2,'LineWidth',3);
grid on; axis tight;
datetick('x','yyyy-mm','keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]);  
h=title('Inflation');
set(h,'interpreter','latex','fontsize',20);

figure
plot(dates(datesperiod(1:end-1)),(e_euus_t(datesperiod(2:end))./e_euus_t(datesperiod(1:end-1))-1)*12e2,'LineWidth',3); hold on;
plot(dates(datesperiod(1:end-1)),(exp(inv_e(datesperiod(1:end-1)))./exp(inv_e(datesperiod(2:end)))-1)*12e2,'LineWidth',1,'LineStyle','--'); hold on;
grid on; axis tight;
datetick('x','yyyy-mm','keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
legend('Model','Data');
h=title('FX');
set(h,'interpreter','latex','fontsize',20);

figure
plot(dates(datesperiod),im_us(datesperiod)*12e2,'LineWidth',3); hold on;
plot(dates(datesperiod),im_eu(datesperiod)*12e2,'LineWidth',1,'LineStyle','--'); hold on;
grid on; axis tight;
datetick('x','yyyy-mm','keeplimits');
label_x('Time (Year-Month)');
formataxis(gca);
orient landscape;
set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 1 1]); 
legend('im_{us}','im_{eu}');
h=title('Interest Rate');
set(h,'interpreter','latex','fontsize',20);

%% Computing Steady State using estimates
% mean shocks
sigma_us_av= mean(sigma_us_TED_t);
sigma_eu_av= mean(sigma_eu_TED_t);
Rm_us_av   = mean(Rm_us);
Rm_eu_av   = mean(Rm_eu);
riskprm_av = mean(riskprm_t);
Theta_d_us_av =mean(Theta_d_us_t);
Theta_d_eu_av =mean(Theta_d_us_t);

%% Saving results from non-linear filter
% save LFX_rwfilter.mat sigma_us_rw sigma_eu_rw Theta_d_eu_rw Theta_d_us_rw riskprm_rw;