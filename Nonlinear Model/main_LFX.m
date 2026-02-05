%% Liquidity Exchange Rate Model
%
% Version 1: Los Angeles, July, 2019
% Version 9: Los Angeles, February, 2025

% TO DO LIST:
% [1] Solves for Steady State of the model..
% [2] check interior solutions
clear; close all;%

if strcmp(getenv('HOME'),'/Users/sakiclaudia')
    foldername='/Users/sakiclaudia/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/quantfigs/';
else
    foldername='/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion/quantfigs/';
end

printit=1;
plotit= 0;

tic;

% Plot Preferences
T = 5;
LFX_plotprefs;
freq = 12;
rate_scale=1e4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'Z_im_us.mat';
load 'Z_im_eu.mat';

%load 'im_parameter.mat';
%load 'inf_mean.mat';
% LFX_empty_mats;

% Parameters
% pi_eu_ss = 1;
pi_eu_ss = 1;
pi_us_ss = 1;
load dynare_calibration_param.mat;
LFX_params_cd_v3;

% Theta_b = 100^(-1/-35);
Theta_b = 1;
%Theta_b = 1.5;
Theta_d_eu = 1; %300
Theta_d_us = 1;
% Theta_d_eu = 300^(-1/35);
% Theta_d_us = 300^(-1/35);
%epsilon_b = -0.1;
% epsilon_b = -35;
% epsilon_b = -1/35;
epsilon_b = -0.001;
% epsilon_b  = -inf;

zeta_us = 1000;
zeta_eu = 1000;

% Adjustment to Euro rate
im_eu_adj=0.0006; % needed to get cip

%% Name Scenario
nameplot='calibration_dynare_sigma_us';
xperiment='$\epsilon^{\lambda^{*}}$';

%% Load Eq Equations
LFX_nt_0e_eqs_2;

%% Plot Partial Equilibrium Equations
N_mu=50;
max_mu_us=(1-nu_b)/nu_us_d;
max_mu_eu=(1-nu_b)/nu_eu_d;
mu_us_vec=linspace(0.00,max_mu_us,N_mu); % Vector of Liquidity Ratio in Dollars
mu_eu_vec=mu_eu_ame(mu_us_vec);

% Plot Euro balances as function of $
figure(1)
plot(mu_us_vec,mu_eu_vec,'LineWidth',3,'Color',[0.5 0.5 0.5]);  hold on;
xlabel('$\mu_{us}$','interpreter','latex');
ylabel('$\mu_{eu}$','interpreter','latex');
title('Iso-Dollar Premium ($R^{m}-R^{*,m})$ (BPS)','interpreter','latex')
line([0 max_mu_us],[0 max_mu_eu],'Color','k','LineStyle','--'); hold on;

grid on;
% xlim([0 1]); ylim([0 1]);

%% Graphic for Steady State
N_mu=50;
mu_vec=linspace(0,max_mu_us,N_mu);
res_mat=NaN(N_mu,N_mu);
res2_mat=NaN(N_mu,N_mu);
for mm=1:length(mu_vec) % Euro
    for nn=1:length(mu_vec) % dollar
        res2_mat(mm,nn)=(1-F(-mu_vec(nn),ploss_us,sigma_us))*chi_p(theta(mu_vec(nn),ploss_us,sigma_us),iota_us,lambda_us)...
            +F(-mu_vec(nn),ploss_us,sigma_us)*chi_m(theta(mu_vec(nn),ploss_us,sigma_us),iota_us,lambda_us)...
            -(1-F(-mu_vec(mm),ploss_eu,sigma_eu))*chi_p(theta(mu_vec(mm),ploss_eu,sigma_eu),iota_eu,lambda_eu)...
            -F(-mu_vec(mm),ploss_eu,sigma_eu)*chi_m(theta(mu_vec(mm),ploss_eu,sigma_eu),iota_eu,lambda_eu);
    end
end

% Figures
figure(2)
% imagesc(mu_vec,mu_vec,-res_mat*rate_scale); hold on;
surf(mu_vec,mu_vec,res2_mat*rate_scale*freq); hold on;
alpha(0.5); axis tight;

% Figures
figure(1)
% imagesc(mu_vec,mu_vec,-res_mat*rate_scale); hold on;
contour(mu_vec,mu_vec,res2_mat*rate_scale*freq,'ShowText','on','LineWidth',1,'LineStyle',':'); hold on;
alpha(0.5); axis tight;

% Orientation
if printit==2
    orient landscape
    print('F_LFX_inelastic_indifference','-dpdf','-fillpage')
end

%% [MZ- please tranform into MPCC format...] Solve Steady State Equilibrium
% I'm assumming stuff from above produces SS.
% test if solution was found:
mu_us_guess=0.01;
res_init=mind_res_eq(mu_us_guess);
testplot=0;
if testplot
    figure; fplot(mind_res_eq,[-1 2]);
end

% Solving for Steady State
[mu_us_ss,fval,exitflag,~]=mu_us_star_f();
mu_eu_ss=mu_eu_ame(mu_us_ss);
scatter(mu_us_ss,mu_eu_ss,40,'r','filled'); drawnow;

% Test1=Rm_eu-Rm_us-0.5*(chi_p(theta(mu_us,delta_us),iota_us,lambda_us)+chi_m(theta(mu_us,delta_us),iota_us,lambda_us));
Test0=mind_res_eq(mu_us_ss);
Test0b=mind_res(mu_eu_ss,mu_us_ss);

% Interest Rates
Rd_us_ss=Rd_us_f(mu_us_ss);
Rd_eu_ss=Rd_eu_f(mu_us_ss);
Rb_us_ss   =Rb_us_f(mu_us_ss);
Rb_eu_ss   =Rb_eu_f(mu_us_ss);

% Calibration
load exchange_rate_data.mat ln_eu_us_ss;
load LFX_data2.mat mu_us;
share = 1;
target = [200;ln_eu_us_ss;mean((mu_us))];
%target = [56.5403;ln_eu_us_ss];
x0 = [mu_eu_ss;mu_us_ss;Rd_eu_ss;Rd_us_ss;Rb_us_ss;bard_us;bard_eu/bard_us];
[x,fval,exitflag,~]=fsolve(@(x) ...
    feqm_calibrate(x,Echi_d,Echi_m,Rm_eu,Rm_us,ploss_eu,ploss_us,iota_eu,iota_us,lambda_eu,lambda_us,Theta_b,epsilon_b,Theta_d_eu,Theta_d_us,zeta_eu,zeta_us,M_eu,M_us,target,share),...
    [x0;sigma_eu;sigma_us;Theta_d_us],optimoptions('fsolve','Display','off','TolFun',1e-12,'MaxFunEval',1e9,'MaxIter',1e6));
sigma_eu = x(8);
sigma_us = x(9);
Theta_d_eu = x(10);
Theta_d_us = x(10);

% Calculate SS
[x,fval,exitflag,~]=fsolve(@(x) ...
    feqm(x,Echi_d,Echi_m,Rm_eu,Rm_us,ploss_eu,ploss_us,sigma_eu,sigma_us,iota_eu,iota_us,lambda_eu,lambda_us,Theta_b,epsilon_b,Theta_d_eu,Theta_d_us,zeta_eu,zeta_us,share),...
    x0,optimoptions('fsolve','Display','off','TolFun',1e-12));
mu_eu_ss = x(1);
mu_us_ss = x(2);
Rd_eu_ss = x(3);
Rd_us_ss = x(4);
Rb_us_ss = x(5);
Rb_eu_ss = x(5);
d_us_ss = x(6);
nu_ss = x(7);
d_eu_ss = Theta_d_eu*(Rd_eu_ss)^(1/zeta_eu);
p_eu_ss = M_eu/(d_eu_ss*mu_eu_ss);
inv_e_ss = M_us/(d_us_ss*mu_us_ss)/p_eu_ss;
MBS_ss = d_us_ss*share;
sigmass_eu = sigma_eu;
sigmass_us = sigma_us;
Theta_dss_eu = Theta_d_eu;
Theta_dss_us = Theta_d_us;

% Price System
%p_eu_ss=p_eu_f(mu_eu_ss)          ;
%inv_e_ss=inv_e_f(p_eu_ss,mu_us_ss);
p_us_ss=p_us_f(p_eu_ss,inv_e_ss)  ;
e_euus_ss=e_euus_f(inv_e_ss)      ;

% Steady State Liquidity Values
theta_us_ss=theta(mu_us_ss,ploss_us,sigma_us);
theta_eu_ss= theta(mu_eu_ss,ploss_eu,sigma_eu);
chi_p_us_ss                   = chi_p(theta(mu_us_ss,ploss_us,sigma_us),iota_us,lambda_us);
chi_m_us_ss                   = chi_m(theta(mu_us_ss,ploss_us,sigma_us),iota_us,lambda_us);
chi_p_eu_ss                   = chi_p(theta(mu_eu_ss,ploss_eu,sigma_eu),iota_eu,lambda_eu);
chi_m_eu_ss                   = chi_m(theta(mu_eu_ss,ploss_eu,sigma_eu),iota_eu,lambda_eu);
Echi_d_us_ss                  = Echi_d(mu_us_ss,ploss_us,sigma_us,iota_us,lambda_us);
Echi_d_eu_ss                  = Echi_d(mu_eu_ss,ploss_eu,sigma_eu,iota_eu,lambda_eu);
Echi_m_us_ss                  = Echi_m(mu_us_ss,ploss_us,sigma_us,iota_us,lambda_us);
Echi_m_eu_ss                  = Echi_m(mu_eu_ss,ploss_eu,sigma_eu,iota_eu,lambda_eu);
Echi_eu_ss                    = Echi(mu_eu_ss,ploss_eu,sigma_eu,iota_eu,lambda_eu);
Echi_us_ss                    = Echi(mu_us_ss,ploss_us,sigma_us,iota_us,lambda_us);

% Saving Additional Rates
RBond_us_ss  = RBond_us(mu_us_ss) ;
RBond_eu_ss  = RBond_eu(mu_eu_ss) ;
RLibor_us_ss = RLibor_us(mu_us_ss);
RLibor_eu_ss = RLibor_eu(mu_eu_ss);

%load im_multicur_param;
load exchange_rate_data;
curlist = {'au','ca','jp','nz','no','sw','ch','uk'};
x0 = [mu_us_ss;Rd_us_ss];
for j=1:length(curlist)
    eval(['Rm_temp = imss_' curlist{j} ';']);
    eval(['iota_temp = iwss_' curlist{j} '-imss_' curlist{j} ';']);
    eval(['target = ln_' curlist{j} '_us_ss;']);
    [x,fval] = fsolve(@(x) ...
    feqm_multicur(x,Echi_d,Echi_m,Rm_temp,Rm_us,ploss_eu,ploss_us,sigma_eu,sigma_us,iota_temp,iota_us,lambda_eu,lambda_us,Rd_us_ss,mu_us_ss),...
    x0,optimoptions('fsolve','Display','off','TolFun',1e-12,'MaxFunEval',1e9,'MaxIter',1e6));
%     [~,temp1,temp2] = feqm_multicur(x0,Echi_d,Echi_m,Rm_temp,Rm_us,ploss_eu,ploss_us,sigma_eu,sigma_us,iota_temp,iota_us,lambda_eu,lambda_us,Rd_us_ss,mu_us_ss);
%     [x,fval,exitlag] = fsolve(@(x) ...
%     fsigma_multicur(x,target,Echi_d,Echi_m,Rm_temp,Rm_us,ploss_eu,ploss_us,sigma_eu,sigma_us,iota_temp,iota_us,lambda_eu,lambda_us,Rd_us_ss,mu_us_ss,p_us_ss,zeta_eu),...
%     x0,optimoptions('fsolve','Display','off','TolFun',1e-12));
    eval(['Rm_' curlist{j} ' = Rm_temp;']);
    eval(['mu_' curlist{j} '_ss = x(1);']);
    eval(['Rd_' curlist{j} '_ss = x(2);']);
%     eval(['M_' curlist{j} ' = x(3);']);
%     eval(['p_' curlist{j} '_ss = x(4);']);
    %eval(['p_' curlist{j} '_ss = M_eu/mu_' curlist{j} '_ss/(Rd_' curlist{j} '_ss^(1/zeta_eu));']);
    eval(['p_' curlist{j} '_ss = p_us_ss*exp(target);']);
    eval(['M_' curlist{j} ' = p_' curlist{j} '_ss*x(1)*x(2)^(1/zeta_eu);']);
    eval(['inv_e_' curlist{j} '_ss = p_us_ss/p_' curlist{j} '_ss;']);
    eval(['e_' curlist{j} 'us_ss = 1/inv_e_' curlist{j} '_ss;']);
end

%% Test Solutions: all returns are equal:
Rm_us_ss=Rm_us;
Rm_eu_ss=Rm_eu-im_eu_adj;
CIP_approx=(Rm_eu_ss^12-Rm_us_ss^12)*1e4;
Rm_us=Rm_us;
Rm_eu=Rm_eu-0.0006;

% Test for solutions
Test1=Rm_us_ss+(1-F(-mu_us_ss,ploss_us,sigma_us))*chi_p_us_ss+F(-mu_us_ss,ploss_us,sigma_us)*chi_m_us_ss...
    -Rm_eu_ss-(1-F(-mu_eu_ss,ploss_eu,sigma_eu))*chi_p_eu_ss-F(-mu_eu_ss,ploss_eu,sigma_eu)*chi_m_eu_ss;
Test_Echi_m=Echi_m_us_ss-(F(-mu_us_ss,ploss_us,sigma_us)*chi_m_us_ss...
    +(1-F(-mu_us_ss,ploss_us,sigma_us))*chi_p_us_ss);
Test_Echi_d=Echi_d_us_ss...
    -(MassUnd(-mu_us_ss,ploss_us,sigma_us)*chi_m(theta(mu_us_ss,ploss_us,sigma_us),iota_us,lambda_us)...
    +(-MassUnd(-mu_us_ss,ploss_us,sigma_us))*chi_p(theta(mu_us_ss,ploss_us,sigma_us),iota_us,lambda_us));
Test2=Rm_us_ss+...
    +F(-mu_us_ss,ploss_us,sigma_us)*chi_m_us_ss+(1-F(-mu_us_ss,ploss_us,sigma_us))*chi_p_us_ss...
    -Rd_eu_ss...
    +MassUnd(-mu_eu_ss,ploss_eu,sigma_eu)*chi_m(theta(mu_eu_ss,ploss_eu,sigma_eu),iota_eu,lambda_eu)...
    +(-MassUnd(-mu_eu_ss,ploss_eu,sigma_eu))*chi_p(theta(mu_eu_ss,ploss_eu,sigma_eu),iota_eu,lambda_eu);
Test3=Rm_us_ss...
    +(1-F(-mu_us_ss,ploss_us,sigma_us))*chi_p_us_ss+F(-mu_us_ss,ploss_us,sigma_us)*chi_m_us_ss...
    -Rd_us_ss+Echi_d_us_ss;
Test_d     = -Rd_us_ss+Echi_d_us_ss+Rd_eu_ss-Echi_d_eu_ss;
Test4=Echi_us_ss-(Echi_m_us_ss*mu_us_ss+Echi_d_us_ss);

% Construct Bank Profits
Portfolio_ss=Rb_eu_ss*nu_b-Rd_us_ss*nu_us_d-Rd_eu_ss*nu_eu_d+Rm_us_ss*mu_us_ss*nu_us_d+Rm_eu_ss*mu_eu_ss*nu_eu_d;
settle_us_ss=Echi_us_ss*nu_us_d;
settle_eu_ss=Echi_eu_ss*nu_eu_d;
BankProfits_ss=Portfolio_ss+settle_us_ss+settle_eu_ss;

% Recording Steady State
steady.p_us_ss = p_us_ss       ;
steady.p_eu_ss = p_eu_ss       ;
steady.mu_us_ss = mu_us_ss     ;
steady.mu_eu_ss = mu_eu_ss     ;
steady.e_euus_ss = e_euus_ss   ;
steady.Rd_us_ss = Rd_us_ss     ;
steady.Rd_eu_ss = Rd_eu_ss     ;
steady.Rb_us_ss = Rb_us_ss     ;
steady.Rb_eu_ss = Rb_eu_ss     ;

%% Report Variables - Steady State Variables
disp('--------------------------------------------');
disp('-------------- Baseline Solution ---------');
disp('--------------------------------------------');
disp(['Eq. USD deposits: ' num2str(bard_us)]);
disp(['Eq. EU deposits: ' num2str(bard_eu)]);
disp(['Eq. USD loans: ' num2str(barB)]);
disp(['Eq. Euro loans: ' num2str(barB)]);
disp(['Eq. USD reserves: ' num2str(mu_us_ss)]);
disp(['Eq. EU reserves: ' num2str(mu_eu_ss)]);
disp(['Eq. EU/USD deposits: ' num2str(e_euus_ss)]);
disp(['Eq. US bank price: ' num2str(p_us_ss)]);
disp(['Eq. Eud bank price: ' num2str(p_eu_ss)]);
disp(['Eq. USD deposits: ' num2str(Rd_us_ss)]);
disp(['Eq. EU deposit rate: ' num2str(Rd_eu_ss)]);
disp(['Eq. USD loan rate: ' num2str(Rb_us_ss)]);
disp(['Eq. EU  loan rate: ' num2str(Rb_eu_ss)]);
disp(['EBP: ' num2str((Rb_us_ss-Rm_us_ss)*1e4*12) 'bps']);
disp(['Rb_eu-Rd_eu: ' num2str((Rb_eu_ss-Rd_eu_ss)*100*12) '%']);
disp(['Rb_us-Rd_us: ' num2str((Rb_us_ss-Rd_us_ss)*100*12) '%']);
disp(['LP: ' num2str((Rm_eu_ss-Rm_us_ss)*1e4*12) 'bps']);
disp(['Average log(e_euus_ss): ' num2str(log(e_euus_ss))]);

curlist = {'au','ca','jp','nz','no','sw','ch','uk'};
conlist = {'AUS','CAN','JPN','NZL','NOR','SWE','CHE','GBR'};
savelist = '';
for j=1:length(curlist)
    disp(['M_' curlist{j} ': ' num2str(eval(['M_' curlist{j}]))]);
    disp(['Eq. ' conlist{j} ' deposits: ' num2str(eval(['Rd_' curlist{j} '_ss']))]);
    disp(['Eq. ' conlist{j} ' reserves: ' num2str(eval(['mu_' curlist{j} '_ss']))]);
    disp(['Eq. ' conlist{j} ' inv exchg: ' num2str(eval(['inv_e_' curlist{j} '_ss']))]);
    disp(['Eq. ' conlist{j} ' exchg: ' num2str(eval(['e_' curlist{j} 'us_ss']))]);
    disp(['Eq. ' conlist{j} ' price: ' num2str(eval(['p_' curlist{j} '_ss']))]);
    savelist = [savelist 'M_' curlist{j} ' '];
end
eval(['save M_curr.mat ' savelist ' MBS_ss sigmass_us sigmass_eu Theta_dss_eu Theta_dss_us;']);


%% Global Solution
N_sigma_us = 102;

% Rates
Z_im_us = (imss_us);
Zprob_im_us = 1    ;
Z_im_eu = (imss_eu);
Zprob_im_eu = 1    ;

% Index by States
index1=1:N_sigma_us/2;
index2=N_sigma_us/2+1:N_sigma_us;

% % % Normal Regime:
mu_sigma_us_r1 = -1.27 ;
rho_sigma_us_r1 = 0.88 ;
Sigma_sigma_us_r1 = 0.2;
m_sigma_us_r1 = 8.0   ;
% 
% % Volatile Regime
mu_sigma_us_r2 = -0.799;
rho_sigma_us_r2 = 0.61 ;
Sigma_sigma_us_r2 = 0.9;
m_sigma_us_r2 = 2.5    ;

% % Transition Probabilities
P=[0.984 0.016; 0.061 0.939];

% % % Normal Regime:
% mu_sigma_us_r1 = -1.0;
% rho_sigma_us_r1 = 0.975;
% Sigma_sigma_us_r1 = 0.2;
% m_sigma_us_r1 = 3.5;
% 
% % Volatile Regime
% mu_sigma_us_r2 = -1.0;
% rho_sigma_us_r2 = 0.0;
% Sigma_sigma_us_r2 = 0.2;
% m_sigma_us_r2 = 14.5;
% 
% % Transition Probabilities
% P=[0.999 0.001; 0.001 0.999];

% Construction of process
Tauchen_method=1; Rouwenhorst_method=0;
if Tauchen_method==1
    [Z_sigma_us_r1,Zprob_sigma_us_r1] = tauchen(N_sigma_us/2,mu_sigma_us_r1,rho_sigma_us_r1,Sigma_sigma_us_r1,m_sigma_us_r1);
    [Z_sigma_us_r2,Zprob_sigma_us_r2] = tauchen(N_sigma_us/2,mu_sigma_us_r2,rho_sigma_us_r2,Sigma_sigma_us_r2,m_sigma_us_r2);

elseif Rouwenhorst_method==1
    [Z_sigma_us_r1,Zprob_sigma_us_r1] =discretizeAR1_Rouwenhorst(mu_sigma_us_r1*(1-rho_sigma_us_r1),rho_sigma_us_r1,Sigma_sigma_us_r1,N_sigma_us/2);
    [Z_sigma_us_r2,Zprob_sigma_us_r2] =discretizeAR1_Rouwenhorst(mu_sigma_us_r2,rho_sigma_us_r2,Sigma_sigma_us_r2,N_sigma_us/2);
end
[eig_vecs1,eigs1]=eig(Zprob_sigma_us_r1');
invp1=eig_vecs1(:,1)/(sum(eig_vecs1(:,1)));
[eig_vecs2,eigs2]=eig(Zprob_sigma_us_r2');
invp2=eig_vecs2(:,1)/(sum(eig_vecs2(:,1)));
figure('Name','Invariant Distributions')
area(exp(Z_sigma_us_r1),invp1,'FaceAlpha',0.5,'FaceColor','b'); hold on;
area(exp(Z_sigma_us_r2),invp2,'FaceAlpha',0.5,'FaceColor','r'); hold on;

% Test location
figure('Name','Locations')
scatter(Z_sigma_us_r1,Z_sigma_us_r1*0+1,'Filled','Color',[0.1 0.5 0.8]); hold on;
scatter(Z_sigma_us_r2,Z_sigma_us_r2*0+1.1,'Filled','Color',[0.8 0.2 0.8]); grid on;
legend('low','high');

% Index Key - Find nearest location for construction of process:
for ii=1:N_sigma_us/2
    [~,index_r1_to_r2(ii)]=min(abs(Z_sigma_us_r2-Z_sigma_us_r1(ii)));
    [~,index_r2_to_r1(ii)]=min(abs(Z_sigma_us_r1-Z_sigma_us_r2(ii)));
end

% Single Chain
Z_sigma_us=[Z_sigma_us_r1; Z_sigma_us_r2];
Zprob_sigma_us=zeros(N_sigma_us);
for ii=1:N_sigma_us/2
    Zprob_sigma_us(ii,:)=[P(1,1)*Zprob_sigma_us_r1(ii,:) P(1,2)*Zprob_sigma_us_r2(index_r1_to_r2(ii),:)];
    Zprob_sigma_us(N_sigma_us/2+ii,:)=[P(2,1)*Zprob_sigma_us_r1(index_r2_to_r1(ii),:) P(2,2)*Zprob_sigma_us_r2(ii,:)];
end
[eig_vecs,eigs]=eig(Zprob_sigma_us');
invp=eig_vecs(:,1)/(sum(eig_vecs(:,1)));
invpp1=invp(1:N_sigma_us/2); invp1=invpp1/sum(invpp1);
invpp2=invp(N_sigma_us/2+1:end); invp2=invpp2/sum(invpp2);
figure('Name','Conditional Distributions')
area(exp(Z_sigma_us_r1),invp1,'FaceAlpha',0.5,'FaceColor','b'); hold on;
area(exp(Z_sigma_us_r2),invp2,'FaceAlpha',0.5,'FaceColor','r'); hold on;

% Other Parameters
N_Theta_d_us = 1;
mu_Theta_d_us = 0;
%rho_Theta_d_us = 0.98;
rho_Theta_d_us = 0.9828;
%sigma_Theta_d_us = 0.05;
sigma_Theta_d_us = 0.073;
m_Theta_d_us = 3;
[Z_Theta_d_us,Zprob_Theta_d_us] = tauchen(N_Theta_d_us,mu_Theta_d_us,rho_Theta_d_us,sigma_Theta_d_us,m_Theta_d_us);
Z_Theta_d_us = 0;
Zprob_Theta_d_us = 1;

% N_sigma_eu = 11;
% mu_sigma_eu = log(sigma_eu);
% rho_sigma_eu = 0.95;
% sigma_sigma_eu = 0.2;
% m_sigma_eu = 3;
% [Z_sigma_eu,Zprob_sigma_eu] = tauchen(N_sigma_eu,mu_sigma_eu,rho_sigma_eu,sigma_sigma_eu,m_sigma_eu);

sigmashock_us_vec = zeros(1,length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
imshock_us_vec = zeros(1,length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
imshock_eu_vec = zeros(1,length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
Thetadshock_us_vec = zeros(1,length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
Q_mat = zeros(length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us),length(Z_sigma_us)*length(Z_im_us)*length(Z_im_eu)*length(Z_Theta_d_us));
id = 0;
for i1=1:length(Z_im_us)
    for i2=1:length(Z_sigma_us)
        for i3=1:length(Z_im_eu)
            for i4=1:length(Z_Theta_d_us)
                id = id+1;
                sigmashock_us_vec(id) = exp(Z_sigma_us(i2));
                imshock_us_vec(id) = Z_im_us(i1);
                imshock_eu_vec(id) = Z_im_eu(i3);
                Thetadshock_us_vec(id) = exp(Z_Theta_d_us(i4));
                temp = Zprob_sigma_us(i2,:)'*Zprob_im_us(i1,:);
                temp = Zprob_im_eu(i3,:)'*temp(:)';
                temp = Zprob_Theta_d_us(i4,:)'*temp(:)';
                Q_mat(id,:) = temp(:)';
            end
        end
    end
end

Q_mat(Q_mat<1e-9) = 0;
Q_mat = sparse(Q_mat);
N_s=size(Q_mat,1);

if N_s==length(Z_sigma_us)
    Q_mat=Zprob_sigma_us;
end
% Initiate Endogenous Variables
LFX_empty_vecs_state; % initiates times series

%% [III.a] Create Exogenous Paths for shocks
lambda_eu_vec(1:N_s) = lambda_eu;
lambda_us_vec(1:N_s) = lambda_us;
ploss_eu_vec(1:N_s)  = ploss_eu ;
ploss_us_vec(1:N_s)  = ploss_us ;
sigma_eu_vec(1:N_s)  = sigma_eu ;
sigma_us_vec(1:N_s)  = sigmashock_us_vec(:)' ;
iw_eu_vec(1:N_s)     = iw_eu    ;
iw_us_vec(1:N_s)     = iw_us-im_us+imshock_us_vec(:)' ;
im_eu_vec(1:N_s)     = im_eu   - im_eu_adj;
im_us_vec(1:N_s)     = im_us-im_us+imshock_us_vec(:)'    ;
M_eu_vec(1:N_s)      = M_eu     ;
M_us_vec(1:N_s)      = M_us     ;
bard_eu_vec(1:N_s)   = bard_eu  ;
bard_us_vec(1:N_s)   = bard_us  ;
Theta_d_us_vec(1:N_s)= Theta_d_us*Thetadshock_us_vec(:)';
Theta_d_eu_vec(1:N_s)= Theta_d_eu;
Theta_b_vec(1:N_s)   = Theta_b;

% Paths for exogenous
techpath_list={'lambda_eu_vec','lambda_us_vec','ploss_eu_vec','ploss_us_vec','sigma_eu_vec','sigma_us_vec'};
polpath_list={'iw_eu_vec','iw_us_vec','im_eu_vec','im_us_vec','M_eu_vec','M_us_vec'};
scalepath_list={'bard_eu_vec','bard_us_vec'};
paths_list={techpath_list{:} polpath_list{:} scalepath_list{:}};
techparam_list={'lambda_eu','lambda_us','ploss_eu','ploss_us','sigma_eu','sigma_us'};
polparam_list={'iw_eu','iw_us','im_eu','im_us','M_eu','M_us'};
scaleparam_list={'bard_eu','bard_us'};
param_list={techparam_list{:} polparam_list{:} scaleparam_list{:}};

for ii=1:numel(paths_list)
    eval(['paths.' paths_list{ii} '=' paths_list{ii} ';']);
end

% [III.a.a] Run a code that reports the shocks
% [MZ fix]  LFX_report_shocks; % _>Copy from MPPC_report_shocks; [Done]
%     LFX_report_shocks;

% Reconstruction of normalized variables...
iota_eu_vec = iw_eu_vec-im_eu_vec;
iota_us_vec = iw_us_vec-im_us_vec;

%---------- Supply D --------------------
bard_tot_vec       = bard_us_vec+bard_eu_vec;

% ------- Some important ratios ----
nu_us_d_vec        = bard_us_vec./bard_tot_vec;
nu_eu_d_vec        = bard_eu_vec./bard_tot_vec;
nu_b_vec           = barB/bard_tot;

% Transitions One Period Dynamics
M_euus_ratio_vec=M_eu_vec./M_us_vec;

% proceed backwards... (not needed if use fsolve)
close all;
greedout=0.03;
greedin=0.01;
tol=1e-7; condout=2*tol;
iterout=0;
maxiterations=5000;
p_us_vec_in=p_us_vec;
p_eu_vec_in=p_eu_vec;
%x0 = [mu_eu_ss;mu_us_ss;Rd_eu_ss;Rd_us_ss;Rb_us_ss;bard_us;bard_eu/bard_us];
pi_eu_ss = 1;
pi_us_ss = 1;
x0 = [mu_eu_ss;mu_us_ss;Rd_eu_ss;Rd_us_ss;Rb_us_ss;bard_us;bard_eu/bard_us;Rm_eu_ss;Rm_us_ss;iota_eu_vec(1);iota_us_vec(1);pi_eu_ss;pi_us_ss;p_eu_ss;p_us_ss];
x0 = kron(x0,ones(1,N_s));
%load initguess.mat;
% solve transition path
%{
options = optimoptions('fsolve','Display','off');
pout = fsolve(@(p) LFX_notrade_inelastic_path_solve(p(1:N_t),p(N_t+1:end),param_list,paths,steady),[p_us_vec_in;p_eu_vec_in],options);
p_us_vec = pout(1:N_t);
p_eu_vec = pout(N_t+1:end);
%}
clear mind_res;
clear mu_eu_ame;
clear mu_us_star_f;
LFX_nt_0e_eqs_2;

%% Computing Solutions
%[mu_us_vec(ss),fval,exitflag(ss),~]=mu_us_star_f();
%mu_eu_vec(ss)=mu_eu_ame(mu_us_vec(ss));
%load 'initguess.mat';
[x,fval,exitflag,~]=fsolve(@(x) ...
    feqm_vec(x,Echi_d,Echi_m,p_us_f,ploss_eu_vec,ploss_us_vec,sigma_eu_vec,sigma_us_vec,lambda_eu_vec,lambda_us_vec,Theta_b_vec,epsilon_b,Theta_d_eu_vec,Theta_d_us_vec,iw_eu_vec,iw_us_vec,im_eu_vec,im_us_vec,zeta_eu,zeta_us,M_eu,M_us,Q_mat,N_s),...
    x0,optimoptions('fsolve','Display','iter','TolFun',1e-15,'MaxFunctionEvaluations',1e9,'MaxIterations',1e9));
mu_eu_vec = x(1,:);
mu_us_vec = x(2,:);
Rd_eu_vec = x(3,:);
Rd_us_vec = x(4,:);
Rb_us_vec = x(5,:);
Rb_eu_vec = x(5,:);
d_us_vec = x(6,:);
nu_vec = x(7,:);
Rm_eu_vec = x(8,:);
Rm_us_vec = x(9,:);
iota_eu_vec = x(10,:);
iota_us_vec = x(11,:);
pi_eu_vec = x(12,:);
pi_us_vec = x(13,:);
p_eu_vec = x(14,:);
p_us_vec = x(15,:);
d_eu_vec = Theta_d_eu_vec.*(Rd_eu_vec).^(1./zeta_eu);
inv_e_vec = M_us./(d_us_vec.*mu_us_vec)./p_eu_vec;
e_euus_vec=1./(inv_e_vec);

% Liquidity Values
theta_us_vec =theta(mu_us_vec,ploss_us_vec,sigma_us_vec);
theta_eu_vec = theta(mu_eu_vec,ploss_eu_vec,sigma_eu_vec);
chi_p_us_vec = chi_p(theta(mu_us_vec,ploss_us_vec,sigma_us_vec),iota_us_vec,lambda_us_vec);
chi_m_us_vec = chi_m(theta(mu_us_vec,ploss_us_vec,sigma_us_vec),iota_us_vec,lambda_us_vec);
chi_p_eu_vec = chi_p(theta(mu_eu_vec,ploss_eu_vec,sigma_eu_vec),iota_eu_vec,lambda_eu_vec);
chi_m_eu_vec = chi_m(theta(mu_eu_vec,ploss_eu_vec,sigma_eu_vec),iota_eu_vec,lambda_eu_vec);
Echi_d_us_vec = Echi_d(mu_us_vec,ploss_us_vec,sigma_us_vec,iota_us_vec,lambda_us_vec)          ;
Echi_d_eu_vec = Echi_d(mu_eu_vec,ploss_eu_vec,sigma_eu_vec,iota_eu_vec,lambda_eu_vec)          ;
Echi_m_us_vec = Echi_m(mu_us_vec,ploss_us_vec,sigma_us_vec,iota_us_vec,lambda_us_vec)          ;
Echi_m_eu_vec = Echi_m(mu_eu_vec,ploss_eu_vec,sigma_eu_vec,iota_eu_vec,lambda_eu_vec)          ;
Echi_eu_vec = Echi(mu_eu_vec,ploss_eu_vec,sigma_eu_vec,iota_eu_vec,lambda_eu_vec)            ;
Echi_us_vec = Echi(mu_us_vec,ploss_us_vec,sigma_us_vec,iota_us_vec,lambda_us_vec)            ;

% Saving Additional Rates
RBond_us_vec  = RBond_us(mu_us_vec) ;
RBond_eu_vec  = RBond_eu(mu_eu_vec) ;
RLibor_us_vec = RLibor_us(mu_us_vec);
RLibor_eu_vec = RLibor_eu(mu_eu_vec);

% Test iteration
Test1_vec=Rm_us_vec+(1-F(-mu_us_vec,ploss_us_vec,sigma_us_vec)).*chi_p_us_vec+F(-mu_us_vec,ploss_us_vec,sigma_us_vec).*chi_m_us_vec...
    -Rm_eu_vec-(1-F(-mu_eu_vec,ploss_eu_vec,sigma_eu_vec)).*chi_p_eu_vec-F(-mu_eu_vec,ploss_eu_vec,sigma_eu_vec).*chi_m_eu_vec;
% Test of deposit premium
Test_Echi_m_vec=Echi_m_us_vec-(F(-mu_us_vec,ploss_us_vec,sigma_us_vec).*chi_m_us_vec+(1-F(-mu_us_vec,ploss_us_vec,sigma_us_vec)).*chi_p_us_vec);
Test_Echi_d_vec=Echi_d_us_vec...
    -(MassUnd(-mu_us_vec,ploss_us_vec,sigma_us_vec).*chi_m(theta(mu_us_vec,ploss_us_vec,sigma_us_vec),iota_us_vec,lambda_us_vec)...
    +(-MassUnd(-mu_us_vec,ploss_us_vec,sigma_us_vec)).*chi_p(theta(mu_us_vec,ploss_us_vec,sigma_us_vec),iota_us_vec,lambda_us_vec));
Test2_vec=Rm_us_vec+...
    +F(-mu_us_vec,ploss_us_vec,sigma_us_vec).*chi_m_us_vec+(1-F(-mu_us_vec,ploss_us_vec,sigma_us_vec)).*chi_p_us_vec...
    -Rd_eu_vec...
    +MassUnd(-mu_eu_vec,ploss_eu_vec,sigma_eu_vec).*chi_m(theta(mu_eu_vec,ploss_eu_vec,sigma_eu_vec),iota_eu_vec,lambda_eu_vec)...
    +(-MassUnd(-mu_eu_vec,ploss_eu_vec,sigma_eu_vec)).*chi_p(theta(mu_eu_vec,ploss_eu_vec,sigma_eu_vec),iota_eu_vec,lambda_eu_vec);
Test3_vec=Rm_us_vec...
    +(1-F(-mu_us_vec,ploss_us_vec,sigma_us_vec)).*chi_p_us_vec+F(-mu_us_vec,ploss_us_vec,sigma_us_vec).*chi_m_us_vec...
    -Rd_us_vec+Echi_d_us_vec;
Test_d_vec     = -Rd_us_vec+Echi_d_us_vec+Rd_eu_vec-Echi_d_eu_vec;
Test4_vec=Echi_us_vec-(Echi_m_us_vec.*mu_us_vec+Echi_d_us_vec);

% Construct Bank Profits
nu_b_vec = (nu_vec.*(1-mu_eu_vec)+1-mu_us_vec).*d_us_vec./(d_us_vec+d_eu_vec);
nu_us_d_vec = d_us_vec./(d_us_vec+d_eu_vec);
nu_eu_d_vec = d_eu_vec./(d_us_vec+d_eu_vec);
Portfolio_vec=Rb_eu_vec.*nu_b_vec-Rd_us_vec.*nu_us_d_vec-Rd_eu_vec.*nu_eu_d_vec+Rm_us_vec.*mu_us_vec.*nu_us_d_vec+Rm_eu_vec.*mu_eu_vec.*nu_eu_d_vec;
settle_us_vec=Echi_us_vec.*nu_us_d_vec;
settle_eu_vec=Echi_eu_vec.*nu_eu_d_vec;
BankProfits_mat_vec=Portfolio_vec+settle_us_vec+settle_eu_vec;

% Resource Test
TestClear_vec=nu_b_vec-(1-mu_us_vec).*nu_us_d_vec-(1-mu_eu_vec).*nu_eu_d_vec;
TestMus_vec=M_us_vec./p_us_vec-mu_us_vec;
TestMeu_vec=M_eu_vec./p_eu_vec-mu_eu_vec;
            
%{
while condout>tol&&iterout<maxiterations
    iterout=iterout+1;
    p_us_vec_in=greedout.*p_us_vec+(1-greedout).*p_us_vec_in;
    p_eu_vec_in=greedout.*p_eu_vec+(1-greedout).*p_eu_vec_in;
    iterin=0;
    condin=2*tol;
    Rm_eu_vec_in=Rm_eu_vec;
    Rm_us_vec_in=Rm_us_vec;
    iota_eu_vec_in=iota_eu_vec;
    iota_us_vec_in=iota_us_vec;
    while condin>tol&&iterin<maxiterations
        iterin=iterin+1;
        Rm_eu_vec_in=greedin.*Rm_eu_vec+(1-greedin).*Rm_eu_vec_in;
        Rm_us_vec_in=greedin.*Rm_us_vec+(1-greedin).*Rm_us_vec_in;
        iota_eu_vec_in=greedin.*iota_eu_vec+(1-greedin).*iota_eu_vec_in;
        iota_us_vec_in=greedin.*iota_us_vec+(1-greedin).*iota_us_vec_in;
        for ss=1:N_s
            % [MZ fix - make this via a string loop, in case we add more variables]
            
            lambda_eu = lambda_eu_vec(ss) ;
            lambda_us = lambda_us_vec(ss) ;
            ploss_eu= ploss_eu_vec(ss)    ;
            ploss_us= ploss_us_vec(ss)    ;
            sigma_eu= sigma_eu_vec(ss)    ;
            sigma_us= sigma_us_vec(ss)    ;
            iw_eu   = iw_eu_vec(ss)       ;
            iw_us   = iw_us_vec(ss)       ;
            im_eu   = im_eu_vec(ss)       ;
            im_us   = im_us_vec(ss)       ;
            M_eu    = M_eu_vec(ss)        ;
            M_us    = M_us_vec(ss)        ;
            bard_eu = bard_eu_vec(ss)     ;
            bard_us = bard_us_vec(ss)     ;
            Theta_d_us = Theta_d_us_vec(ss);
            Theta_d_eu = Theta_d_eu_vec(ss);
            
            %                 for ii=1:numel(param_list)
            %                     eval([param_list{ii} '=' param_list{ii} '_t(tt);']);
            %                 end
            
            % Reconstruction of normalized variables...
            
            %---------- Supply D --------------------
            bard_tot       = bard_us+bard_eu;
            
            % ------- Some important ratios ----
            nu_us_d        = bard_us/bard_tot;
            nu_eu_d        = bard_eu/bard_tot;
            nu_b           = barB/bard_tot;
            
            % Transitions One Period Dynamics
            M_euus_ratio=M_eu/M_us;
            
            % Update future prices prices...
            %                 p_us_prime=p_us_vec_in(ss+1);
            %                 p_eu_prime=p_us_vec_in(ss+1);
            %                 e_euus_prime=p_us_vec_in(ss+1);
            %                 p_prime_eu=p_us_vec_in(ss+1);
            
            % Guess
            p_us=p_us_ss;
            p_eu=p_eu_ss;
            
            % Inflation Rates
            %             pi_us_vec(ss)=p_us_prime/p_us_vec_in(ss);
            %             pi_eu_vec(ss)=p_eu_prime/p_eu_vec_in(ss);
            
            % Initial Values
            Rm_eu         = Rm_eu_vec_in(ss);
            Rm_us         = Rm_us_vec_in(ss);
            %                 Rm_eu_t(ss,1)   = im_eu/pi_eu_vec(ss);
            %                 Rm_us_t(ss,1)   = im_us/pi_us_vec(ss);
            iota_eu = iota_eu_vec_in(ss);
            iota_us = iota_us_vec_in(ss);
            
            % Updating Equations
            clear mind_res;
            clear mu_eu_ame;
            clear mu_us_star_f;
            LFX_nt_0e_eqs_2;
            
            %% Computing Solutions
            %[mu_us_vec(ss),fval,exitflag(ss),~]=mu_us_star_f();
            %mu_eu_vec(ss)=mu_eu_ame(mu_us_vec(ss));
            [x,fval,exitflag(ss),~]=fsolve(@(x) ...
                feqm(x,Echi_d,Echi_m,Rm_eu,Rm_us,ploss_eu,ploss_us,sigma_eu,sigma_us,iota_eu,iota_us,lambda_eu,lambda_us,Theta_b,epsilon_b,Theta_d_eu,Theta_d_us,zeta_eu,zeta_us),...
                x0,optimoptions('fsolve','Display','off','TolFun',1e-12));
            mu_eu_vec(ss) = x(1);
            mu_us_vec(ss) = x(2);
            Rd_eu_vec(ss) = x(3);
            Rd_us_vec(ss) = x(4);
            Rb_us_vec(ss) = x(5);
            Rb_eu_vec(ss) = x(5);
            d_us_vec(ss) = x(6);
            nu_vec(ss) = x(7);
            d_eu_vec(ss) = Theta_d_eu*(Rd_eu_vec(ss))^(1/zeta_eu);
            p_eu_vec(ss) = M_eu/(d_eu_vec(ss)*mu_eu_vec(ss));
            inv_e_vec(ss) = M_us/(d_us_vec(ss)*mu_us_vec(ss))/p_eu_vec(ss);
            x0 = x;
            % Test1=Rm_eu-Rm_us-0.5*(chi_p(theta(mu_us,delta_us),iota_us,lambda_us)+chi_m(theta(mu_us,delta_us),iota_us,lambda_us));
            %Test0=mind_res_eq(mu_us_vec(ss));
            %Test0b=mind_res(mu_eu_vec(ss),mu_us_vec(ss));
            
            % Interest Rates
            %             Rd_us_vec(ss)=Rd_us_f(mu_us_vec(ss));
            %             Rd_eu_vec(ss)=Rd_eu_f(mu_us_vec(ss));
            %             Rb_us_vec(ss)=Rb_us_f(mu_us_vec(ss));
            %             Rb_eu_vec(ss)=Rb_eu_f(mu_us_vec(ss));
            
            % Price System
            %             p_eu_vec(ss)=p_eu_f(mu_eu_vec(ss))             ;
            %             inv_e_vec(ss)=inv_e_f(p_eu_vec(ss),mu_us_vec(ss));
            p_us_vec(ss)=p_us_f(p_eu_vec(ss),inv_e_vec(ss))  ;
            e_euus_vec(ss)=e_euus_f(inv_e_vec(ss))         ;
            %p_us_t1(ss)=p_us_f2(mu_us_vec(ss))             ;
            pi_us_vec(ss)=Q_mat(ss,:)*p_us_vec_in(:)/p_us_vec(ss);
            pi_eu_vec(ss)=Q_mat(ss,:)*p_eu_vec_in(:)/p_eu_vec(ss);
            
            % Policy Rates
            Rm_eu         = im_eu/pi_eu_vec(ss);
            Rm_us         = im_us/pi_us_vec(ss);
            Rm_eu_vec(ss)=im_eu/pi_eu_vec(ss);
            Rm_us_vec(ss)=im_us/pi_us_vec(ss);
            iota_eu_vec(ss) = (iw_eu-im_eu)/pi_eu_vec(ss);
            iota_us_vec(ss) = (iw_us-im_us)/pi_us_vec(ss);
            
            % Steady State Liquidity Values
            theta_us_vec(ss)=theta(mu_us_vec(ss),ploss_us,sigma_us);
            theta_eu_vec(ss)= theta(mu_eu_vec(ss),ploss_eu,sigma_eu);
            chi_p_us_vec(ss)                   = chi_p(theta(mu_us_vec(ss),ploss_us,sigma_us),iota_us_vec(ss),lambda_us);
            chi_m_us_vec(ss)                   = chi_m(theta(mu_us_vec(ss),ploss_us,sigma_us),iota_us_vec(ss),lambda_us);
            chi_p_eu_vec(ss)                   = chi_p(theta(mu_eu_vec(ss),ploss_eu,sigma_eu),iota_eu_vec(ss),lambda_eu);
            chi_m_eu_vec(ss)                   = chi_m(theta(mu_eu_vec(ss),ploss_eu,sigma_eu),iota_eu_vec(ss),lambda_eu);
            Echi_d_us_vec(ss)                  = Echi_d(mu_us_vec(ss),ploss_us,sigma_us,iota_us_vec(ss),lambda_us)          ;
            Echi_d_eu_vec(ss)                  = Echi_d(mu_eu_vec(ss),ploss_eu,sigma_eu,iota_eu_vec(ss),lambda_eu)          ;
            Echi_m_us_vec(ss)                  = Echi_m(mu_us_vec(ss),ploss_us,sigma_us,iota_us_vec(ss),lambda_us)          ;
            Echi_m_eu_vec(ss)                  = Echi_m(mu_eu_vec(ss),ploss_eu,sigma_eu,iota_eu_vec(ss),lambda_eu)          ;
            Echi_eu_vec(ss)                    = Echi(mu_eu_vec(ss),ploss_eu,sigma_eu,iota_eu_vec(ss),lambda_eu)            ;
            Echi_us_vec(ss)                    = Echi(mu_us_vec(ss),ploss_us,sigma_us,iota_us_vec(ss),lambda_us)            ;
            
            % Saving Additional Rates
            RBond_us_vec(ss)  = RBond_us(mu_us_vec(ss)) ;
            RBond_eu_vec(ss)  = RBond_eu(mu_eu_vec(ss)) ;
            RLibor_us_vec(ss) = RLibor_us(mu_us_vec(ss));
            RLibor_eu_vec(ss) = RLibor_eu(mu_eu_vec(ss));
            
            % Test iteration
            Test1_vec(ss)=Rm_us+(1-F(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss)))*chi_p_us_vec(ss)+F(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss))*chi_m_us_vec(ss)...
                -Rm_eu-(1-F(-mu_eu_vec(ss),ploss_eu_vec(ss),sigma_eu_vec(ss)))*chi_p_eu_vec(ss)-F(-mu_eu_vec(ss),ploss_eu_vec(ss),sigma_eu_vec(ss))*chi_m_eu_vec(ss);
            % Test of deposit premium
            Test_Echi_m_vec(ss)=Echi_m_us_vec(ss)-(F(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss))*chi_m_us_vec(ss)+(1-F(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss)))*chi_p_us_vec(ss));
            Test_Echi_d_vec(ss)=Echi_d_us_vec(ss)...
                -(MassUnd(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss))*chi_m(theta(mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss)),iota_us_vec(ss),lambda_us_vec(ss))...
                +(-MassUnd(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss)))*chi_p(theta(mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss)),iota_us_vec(ss),lambda_us_vec(ss)));
            Test2_vec(ss)=Rm_us+...
                +F(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss))*chi_m_us_vec(ss)+(1-F(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss)))*chi_p_us_vec(ss)...
                -Rd_eu_vec(ss)...
                +MassUnd(-mu_eu_vec(ss),ploss_eu_vec(ss),sigma_eu_vec(ss))*chi_m(theta(mu_eu_vec(ss),ploss_eu_vec(ss),sigma_eu_vec(ss)),iota_eu_vec(ss),lambda_eu_vec(ss))...
                +(-MassUnd(-mu_eu_vec(ss),ploss_eu_vec(ss),sigma_eu_vec(ss)))*chi_p(theta(mu_eu_vec(ss),ploss_eu_vec(ss),sigma_eu_vec(ss)),iota_eu_vec(ss),lambda_eu_vec(ss));
            Test3_vec(ss)=Rm_us...
                +(1-F(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss)))*chi_p_us_vec(ss)+F(-mu_us_vec(ss),ploss_us_vec(ss),sigma_us_vec(ss))*chi_m_us_vec(ss)...
                -Rd_us_vec(ss)+Echi_d_us_vec(ss);
            Test_d_vec(ss)     = -Rd_us_vec(ss)+Echi_d_us_vec(ss)+Rd_eu_vec(ss)-Echi_d_eu_vec(ss);
            Test4_vec(ss)=Echi_us_vec(ss)-(Echi_m_us_vec(ss)*mu_us_vec(ss)+Echi_d_us_vec(ss));
            
            % Construct Bank Profits
            nu_b_vec(ss) = (nu_vec(ss)*(1-mu_eu_vec(ss))+1-mu_us_vec(ss))*d_us_vec(ss)/(d_us_vec(ss)+d_eu_vec(ss));
            nu_us_d_vec(ss) = d_us_vec(ss)/(d_us_vec(ss)+d_eu_vec(ss));
            nu_eu_d_vec(ss) = d_eu_vec(ss)/(d_us_vec(ss)+d_eu_vec(ss));
            Portfolio_vec(ss)=Rb_eu_vec(ss)*nu_b_vec(ss)-Rd_us_vec(ss)*nu_us_d_vec(ss)-Rd_eu_vec(ss)*nu_eu_d_vec(ss)+Rm_us*mu_us_vec(ss)*nu_us_d_vec(ss)+Rm_eu*mu_eu_vec(ss)*nu_eu_d_vec(ss);
            settle_us_vec(ss)=Echi_us_vec(ss)*nu_us_d_vec(ss);
            settle_eu_vec(ss)=Echi_eu_vec(ss)*nu_eu_d_vec(ss);
            BankProfits_mat_vec(ss)=Portfolio_vec(ss)+settle_us_vec(ss)+settle_eu_vec(ss);
            
            % Resource Test
            TestClear_vec(ss)=nu_b_vec(ss)-(1-mu_us_vec(ss))*nu_us_d_vec(ss)-(1-mu_eu_vec(ss))*nu_eu_d_vec(ss);
            TestMus_vec(ss)=M_us_vec(ss)/p_us_vec(ss)-mu_us_vec(ss);
            TestMeu_vec(ss)=M_eu_vec(ss)/p_eu_vec(ss)-mu_eu_vec(ss);
        end
        % Find Conditions
        cond_us_in=max(abs([Rm_us_vec_in(:);iota_us_vec_in(:)]./[Rm_us_vec(:);iota_us_vec(:)]-1));
        cond_eu_in=max(abs([Rm_eu_vec_in(:);iota_eu_vec_in(:)]./[Rm_eu_vec(:);iota_eu_vec(:)]-1));
        condin=max(cond_us_in,cond_eu_in);
        if iterout==1
            disp(['Condition at innerloop iteration: ' num2str(condin)]);
        end
    end
    disp(['Final condition at inner loop iteration: ' num2str(condin)]);
    
    % Find Conditions
    cond_us_out=max(abs(p_us_vec_in./p_us_vec-1));
    cond_eu_out=max(abs(p_eu_vec_in./p_eu_vec-1));
    condout=max(cond_us_out,cond_eu_out);
    disp(['Condition at outerloop iteration: ' num2str(condout)])
    
    %     figure(2)
    %     subplot(2,1,1); plot(p_us_vec); hold on;
    %     subplot(2,1,2); plot(p_eu_vec); hold on; drawnow;
    %     if mod(iterout,100)==0
    %         close all;
    %     end
end
%}
toc;

b_vec = (nu_vec.*(1-mu_eu_vec)+1-mu_us_vec).*d_us_vec;
x0 = x;
save 'initguess.mat' 'x0';
eval(['save global_sol' nameplot '.mat']);


% compute transition path
% LFX_notrade_inelastic_path;

% Some transformations
% TedDelta_us_t= RLibor_us_vec-RBond_us_vec;
% TedDelta_eu_t= RLibor_eu_vec-RBond_eu_vec;
% BondBasis_t  = RBond_eu_vec -RBond_us_vec;
% LiborBasis_t = RLibor_eu_vec-RLibor_us_vec;

% Plot Solutions
% LFX_update_mat;
% LFX_plotresults;

% LFX_plot_scenarios;
% figure
% subplot(2,1,1); plot(s_vec,p_us_vec);
% subplot(2,1,2); plot(s_vec,p_eu_vec);
% Plot 3-D Solutions
close all;
cc=0;

X_vec=sigma_eu_vec;
Y_vec=sigma_us_vec;

% NX=length(s_vec);
% NY=length(s_vec);
% X_mat=reshape(X_vec,NY,NX);
% Y_mat=reshape(Y_vec,NY,NX);
Z_vec=p_us_vec;
% Z_mat=reshape(Z_vec,NY,NX);

% cc=cc+1; plottype='_price_us';
% Z_vec=p_us_vec;
% figure(cc)
% splot1(Y_vec,Z_vec);
% %surf(X_mat,Y_mat,Z_mat);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$p_us$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Price Level U.S.','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end
% 
% cc=cc+1; plottype='_exchangerate';
% Z_vec=e_euus_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% figure(cc)
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('Exchange rate','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Exchange rate','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end
% 
% %cc=cc+1; plottype='_uip_dep';
% cc=cc+1; plottype='_uip';
% figure(cc)
% subplot(2,2,1);
% Z_vec=Rd_eu_vec-Rd_us_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{d}-R^{us}_{d}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{d}-R^{us}_{d}$ (UIP deposit deviation)','interpreter','latex');
% % if printit==1
% %     orient landscape;
% %     printsb(['fig' nameplot plottype]);
% % end
% 
% % cc=cc+1; plottype='_uip_loan';
% subplot(2,2,2);
% Z_vec=Rb_eu_vec-Rb_us_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% % figure(cc)
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{b}-R^{us}_{b}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{b}-R^{us}_{b}$ (UIP loans deviation)','interpreter','latex');
% % if printit==1
% %     orient landscape;
% %     printsb(['fig' nameplot plottype]);
% % end
% 
% % cc=cc+1; plottype='_uip_bond';
% subplot(2,2,3);
% Z_vec=RBond_eu_vec-RBond_us_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% % figure(cc)
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{gov}-R^{us}_{gov}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{gov}-R^{us}_{gov}$ (UIP gov bond deviation)','interpreter','latex');
% % if printit==1
% %     orient landscape;
% %     printsb(['fig' nameplot plottype]);
% % end
% 
% % cc=cc+1; plottype='_uip_libor';
% subplot(2,2,4);
% Z_vec=RLibor_eu_vec-RLibor_us_vec;
% %Z_mat=reshape(Z_vec,NY,NX);
% % figure(cc)
% %surf(X_mat,Y_mat,Z_mat);
% splot1(Y_vec,Z_vec);
% %xlabel('$\sigma_{eu}$','interpreter','latex');
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{libor}-R^{us}_{libor}$ (UIP libor deviation)','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end
% 
% cc=cc+1; plottype='_uip_Rm';
% figure(cc)
% Z_vec=Rm_eu_vec-Rm_us_vec;
% splot1(Y_vec,Z_vec);
% xlabel('$\sigma_{us}$','interpreter','latex');
% ylabel('$R^{eu}_{m}-R^{us}_{m}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{m}-R^{us}_{m}$ (UIP Rm deviation)','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end

% Simulate Markov Chain
rng('default');
start_value = round(size(Q_mat,1)/2);
chain_length = 5000;
chain = zeros(1,chain_length);
draw = zeros(1,chain_length);
chain(1) = start_value;

for t=2:chain_length
    this_step_dist = Q_mat(chain(t-1),:);
    cum_dist = cumsum(this_step_dist);
    draw(t) = rand();
    chain(t) = find(cum_dist>draw(t),1);
end


e_euus_t = e_euus_vec(chain(:));
bp_us_t  = Rb_us_vec(chain(:))-Rm_us_vec(chain(:));
bp_eu_t  = Rb_eu_vec(chain(:))-Rm_eu_vec(chain(:));
uip_Rm_t = Rm_eu_vec(chain(:))-Rm_us_vec(chain(:));
uip_dep_t = Rd_eu_vec(chain(:))-Rd_us_vec(chain(:));
uip_loan_t = Rb_eu_vec(chain(:))-Rb_us_vec(chain(:));
uip_bond_t = RBond_eu_vec(chain(:))-RBond_us_vec(chain(:));
uip_libor_t= RLibor_eu_vec(chain(:))-RLibor_us_vec(chain(:));
d_us_t = d_us_vec(chain(:));
d_eu_t = d_eu_vec(chain(:));
mu_us_t = mu_us_vec(chain(:));
mu_eu_t = mu_eu_vec(chain(:));
pidiff_t = pi_eu_vec(chain(:))-pi_us_vec(chain(:));
im_us_t = im_us_vec(chain(:));
im_eu_t = im_eu_vec(chain(:));
imdiff_t = im_eu_t-im_us_t;
pi_us_t  =pi_us_vec((chain(:)));

disp(['Corr. exchange rate vs UIP Rm: ' num2str(corr(e_euus_t',uip_Rm_t'))]);
disp(['Corr. exchange rate vs UIP Libor: ' num2str(corr(e_euus_t',uip_libor_t'))]);
disp(['Corr. exchange rate vs UIP Deposit: ' num2str(corr(e_euus_t',uip_dep_t'))]);
disp(['Corr. exchange rate vs UIP Loan: ' num2str(corr(e_euus_t',uip_loan_t'))]);
disp(['Corr. exchange rate vs UIP Bond: ' num2str(corr(e_euus_t',uip_bond_t'))]);



% Regression
%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),diff(uip_libor_t(:)),(pidiff_t(2:end)'),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)),ones(999,1)];
%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),diff(uip_libor_t(:)),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)),ones(999,1)];
% X_t = [diff(log(mu_us_t(:))),diff(uip_libor_t(:)),(pidiff_t(2:end)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(mu_us_t(:))),diff(uip_libor_t(:)),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];

%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),diff(imdiff_t(:)),(pidiff_t(2:end)'),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),diff(imdiff_t(:)),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(mu_us_t(:))),diff(imdiff_t(:)),(pidiff_t(2:end)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(mu_us_t(:))),diff(imdiff_t(:)),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(mu_us_t(:))),(pidiff_t(2:end)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];

%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),(pidiff_t(2:end)'),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
%X_t = [diff(log(d_us_t(:))),diff(log(mu_us_t(:))),log(d_us_t(1:end-1)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
X_t = [diff(log(mu_us_t(:))),(pidiff_t(2:end)'),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];
% X_t = [diff(log(mu_us_t(:))),log(mu_us_t(1:end-1)'),ones(chain_length-1,1)];

Y_t = diff(log(e_euus_t)');
sample = chain_length-191+1;
coef = (X_t(sample:end,:)'*X_t(sample:end,:))\(X_t(sample:end,:)'*Y_t(sample:end,:));
[coefb,cib,~,~,stats] = regress(Y_t(sample:end,:),X_t(sample:end,:));
talphaup = 1-0.05/2;
tnu = length(Y_t(sample:end,:))-length(coefb);
tvalue = tinv(talphaup,tnu);
sdb = (coefb-cib(:,1))/tvalue;
disp('         Coef   Std.Err.');
for i=1:length(coef)-1
    %disp(['beta',num2str(i),': ',num2str(coef(i)),',']);
    fprintf('beta%d: %0.4f  %0.4f\n',i,coefb(i),sdb(i));
end
%disp(['const: ',num2str(coef(end))]);
fprintf('const: %0.4f  %0.4f\n',coefb(end),sdb(end));
fprintf('R-sq: %0.4f\n',stats(1));

% moment table
ln_e_sample = log(e_euus_t(sample:end));
mean(ln_e_sample)
std(ln_e_sample)
corrcoef(ln_e_sample(2:end),ln_e_sample(1:end-1))
X=[ones(length(ln_e_sample)-1,1) ln_e_sample(1:end-1)'];
[B,BINT,R,RINT,STATS] = regress(ln_e_sample(2:end)' ,X);

disp('AR 1 coef of ln_e')
B(2)

disp('std residual')
sqrt( std(R)^2/(1-B(2)^2))

mean(log(d_us_t(sample:end)))
std(log(d_us_t(sample:end)))

mean(log(mu_us_t(sample:end)))
std(log(mu_us_t(sample:end)))
corrcoef(log(mu_us_t(sample:end-1))',log(mu_us_t(sample+1:end))')

mean(log(d_eu_t(sample:end)))
std(log(d_eu_t(sample:end)))

mean(log(mu_eu_t(sample:end)))
std(log(mu_eu_t(sample:end)))

mean(pidiff_t(sample:end))*100*12
std(pidiff_t(sample:end))*100*12
corrcoef(pidiff_t(sample:end-1)',pidiff_t(sample+1:end)')

b = zeros(size(Q_mat,1),1);
b(1) = 1;
row = [zeros(1,0),1,zeros(1,size(Q_mat,1)-1)];
A = eye(size(Q_mat))-Q_mat';
A(1,:) = row;
f_ss = A\b;
f_sum = sum(f_ss);
f_ss = f_ss./f_sum;
disp(['E[Rb_us-Rm_us] = ',num2str((Rb_us_vec.^(freq)-Rm_us_vec.^(freq))*f_ss(:)*1e4),' bps']);
disp(['E[Rm_eu-Rm_us] = ',num2str((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*f_ss(:)*1e4),' bps']);
%%%%
% cc=cc+1; plottype='_simexchguipRm';
% figure(cc)
% yyaxis left;plot(e_euus_t(500:end),'LineWidth',2);hold on;
% yyaxis right;plot(uip_Rm_t(500:end),'LineWidth',2);
% xlabel('Period','interpreter','latex');
% h=legend('Exchange Rate','$R^{eu}_{m}-R^{us}_{m}$ UIP Rm deviation','location','best');
% set(h,'interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Exchange Rate and UIP Rm deviation','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end
% 
% cc=cc+1; plottype='_simuiplibor';
% figure(cc)
% plot(uip_dep_t(500:end),'LineWidth',2);xlabel('Period','interpreter','latex');
% ylabel('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex');
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('$R^{eu}_{libor}-R^{us}_{libor}$ (UIP libor deviation)','interpreter','latex');
% if printit==1
%     orient landscape;
%     printsb(['fig' nameplot plottype]);
% end

% %% Plot figures together
% % policy functions
% cc = 0;
% cc = cc+1; plottype='_policyfun';
% x_vec = sigma_us_vec(:);
% x_lab = '$\sigma_{us}$';
% % x_vec = (im_us_vec.^(freq)-1)*100;
% % x_lab = '$i^{us}_m(\%)$';
% figure(cc)
% subplot(4,4,1);
% splot1(x_vec,mu_us_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\mu_{us}$','interpreter','latex');
% scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
% ylim(mean(mu_us_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (US)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,2);
% splot1(x_vec,mu_eu_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\mu_{eu}$','interpreter','latex');
% scale = max(0.001,max(abs(mu_eu_vec-mean(mu_eu_vec))));
% ylim(mean(mu_eu_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (EU)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,3);
% splot1(x_vec,e_euus_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('Exchange rate','interpreter','latex');
% scale = max(0.01,max(abs(e_euus_vec-mean(e_euus_vec))));
% ylim([mean(e_euus_vec)-scale,mean(e_euus_vec)+scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Exchange rate','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,4);
% splot1(x_vec,(Rm_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)]),max([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)])]-1)*100);
% %scale = 10;
% %ylim(mean(([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)]-1)*100)+[-scale,+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,5);
% splot1(x_vec,(Rm_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)]),max([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)])]-1)*100);
% %scale = 0.001;
% %ylim(mean(([Rm_us_vec(:).^(freq);Rm_eu_vec(:).^(freq)]-1)*100)+[-scale,+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,6);
% splot1(x_vec,(Rd_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rd_us_vec(:).^(freq);Rd_eu_vec(:).^(freq)]),max([Rd_us_vec(:).^(freq);Rd_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,7);
% splot1(x_vec,(Rd_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rd_us_vec(:).^(freq);Rd_eu_vec(:).^(freq)]),max([Rd_us_vec(:).^(freq);Rd_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,8);
% splot1(x_vec,(RLibor_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([RLibor_us_vec(:).^(freq);RLibor_eu_vec(:).^(freq)]),max([RLibor_us_vec(:).^(freq);RLibor_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,9);
% splot1(x_vec,(RLibor_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([RLibor_us_vec(:).^(freq);RLibor_eu_vec(:).^(freq)]),max([RLibor_us_vec(:).^(freq);RLibor_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,10);
% splot1(x_vec,(RBond_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([RBond_us_vec(:).^(freq);RBond_eu_vec(:).^(freq)]),max([RBond_us_vec(:).^(freq);RBond_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,11);
% splot1(x_vec,(RBond_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([RBond_us_vec(:).^(freq);RBond_eu_vec(:).^(freq)]),max([RBond_us_vec(:).^(freq);RBond_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,12);
% splot1(x_vec,(Rb_us_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rb_us_vec(:).^(freq);Rb_eu_vec(:).^(freq)]),max([Rb_us_vec(:).^(freq);Rb_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,13);
% splot1(x_vec,(Rb_eu_vec.^(freq)-1)*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% ylim(([min([Rb_us_vec(:).^(freq);Rb_eu_vec(:).^(freq)]),max([Rb_us_vec(:).^(freq);Rb_eu_vec(:).^(freq)])]-1)*100);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,14);
% splot1(x_vec,(Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% scale=max(0.001,max(abs((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100-mean((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100))));
% %ylim(([min([Rm_eu_vec-Rm_us_vec]),max([Rm_eu_vec-Rm_us_vec])])*100);
% ylim([mean((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100)-scale,mean((Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*100)+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{m}-R^{us}_{m}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,15);
% splot1(x_vec,(Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% %ylim(([min([Rd_eu_vec-Rd_us_vec]),max([Rd_eu_vec-Rd_us_vec])])*100);
% scale=max(0.001,max(abs((Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100-mean((Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100))));
% ylim([mean((Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100)-scale,mean((Rd_eu_vec.^(freq)-Rd_us_vec.^(freq))*100)+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{d}-R^{us}_{d}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,16);
% splot1(x_vec,(RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% %ylim(([min([RLibor_eu_vec-RLibor_us_vec]),max([RLibor_eu_vec-RLibor_us_vec])])*100);
% scale=max(0.001,max(abs((RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100-mean((RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100))));
% ylim([mean((RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100)-scale,mean((RLibor_eu_vec.^(freq)-RLibor_us_vec.^(freq))*100)+scale]);
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.30 -0.30 1.60 1.60]);
%     printsb(['fig' nameplot plottype]);
% end
% 
% % Quantities
% cc = cc+1;
% figure(cc); plottype='_policyquant';
% subplot(2,2,1);
% splot1(x_vec,d_us_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('US real deposits $d_{us}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,2);
% splot1(x_vec,d_eu_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('EU real deposits $d_{eu}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,3);
% splot1(x_vec,b_vec);xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real loans $b$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end
% 
% % Simulations
% cc=cc+1; plottype='_simpaths';
% Periods = 1:200;
% figure(cc)
% subplot(4,4,1);
% plot(Periods,mu_us_vec(chain(Periods)),'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\mu_{us}$','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (US)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,2);
% plot(Periods,mu_eu_vec(chain(Periods)),'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\mu_{eu}$','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (EU)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,3);
% plot(Periods,e_euus_vec(chain(Periods)),'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('Exchange rate','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Exchange rate','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,4);
% plot(Periods,(Rm_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,5);
% plot(Periods,(Rm_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,6);
% plot(Periods,(Rd_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,7);
% plot(Periods,(Rd_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,8);
% plot(Periods,(RLibor_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,9);
% plot(Periods,(RLibor_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,10);
% plot(Periods,(RBond_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,11);
% plot(Periods,(RBond_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,12);
% plot(Periods,(Rb_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,13);
% plot(Periods,(Rb_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,14);
% plot(Periods,(Rm_eu_vec(chain(Periods)).^(freq)-Rm_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{m}-R^{us}_{m}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,15);
% plot(Periods,(Rd_eu_vec(chain(Periods)).^(freq)-Rd_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{d}-R^{us}_{d}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,16);
% plot(Periods,(RLibor_eu_vec(chain(Periods)).^(freq)-RLibor_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% ylabel('$\%$','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end

% sim paths of quantities
% cc=cc+1;plottype='_simpathsquant0';
% figure(cc);
% subplot(2,2,1);
% plot(Periods,d_us_vec(chain(Periods)),'Linewidth',2);
% ylabel('$d_{us}$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('US real deposits $d_{us}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,2);
% plot(Periods,d_eu_vec(chain(Periods)),'Linewidth',2);
% ylabel('$d_{eu}$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('EU real deposits $d_{eu}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,3);
% plot(Periods,b_vec(chain(Periods)),'Linewidth',2);
% ylabel('$b$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real loans $b$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end
% 
% % Simulated shock paths
% 
% 
% cc=cc+1; plottype='_shockpath';
% figure(cc)
% subplot(1,2,1);
% plot(Periods,sigma_us_vec(chain(Periods)),'Linewidth',2);
% ylabel('$\sigma_{us}$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Dollar Volatility $\sigma_{us}$','interpreter','latex','Fontsize',15);
% 
% subplot(1,2,2);
% plot(Periods,(im_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$i^{us}_m(\%)$','interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Dollar Libor Rate $i^{us}_m(\%)$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'Units','Inches','PaperUnits','Inches');
%     set(gcf, 'PaperSize', [11,4.5]);
%     set(gcf, 'PaperPosition', [0 0.25 11 4]);
%     printsb(['fig' nameplot plottype]);
%     print(gcf,'-dpdf',[path_g,filesep,['fig' nameplot plottype]])
% end
% 
% % Simulation 2
% % Simulations
% cc=cc+1; plottype='_simpaths2';
% Periods = 1:200;
% figure(cc)
% subplot(4,4,1);
% yyaxis left;plot(Periods,mu_us_vec(chain(Periods)),'Linewidth',2);
% ylabel('$\mu_{us}$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (US)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,2);
% yyaxis left;plot(Periods,mu_eu_vec(chain(Periods)),'Linewidth',2);
% ylabel('$\mu_{eu}$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity ratio (EU)','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,3);
% yyaxis left;plot(Periods,e_euus_vec(chain(Periods)),'Linewidth',2);
% ylabel('Exchange rate','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Exchange rate','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,4);
% yyaxis left;plot(Periods,(Rm_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,5);
% yyaxis left;plot(Periods,(Rm_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_m$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,6);
% yyaxis left;plot(Periods,(Rd_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,7);
% yyaxis left;plot(Periods,(Rd_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_d$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,8);
% yyaxis left;plot(Periods,(RLibor_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,9);
% yyaxis left;plot(Periods,(RLibor_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,10);
% yyaxis left;plot(Periods,(RBond_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,11);
% yyaxis left;plot(Periods,(RBond_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{gov}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,12);
% yyaxis left;plot(Periods,(Rb_us_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{us}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,13);
% yyaxis left;plot(Periods,(Rb_eu_vec(chain(Periods)).^(freq)-1)*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_b$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,14);
% yyaxis left;plot(Periods,(Rm_eu_vec(chain(Periods)).^(freq)-Rm_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{m}-R^{us}_{m}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,15);
% yyaxis left;plot(Periods,(Rd_eu_vec(chain(Periods)).^(freq)-Rd_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{d}-R^{us}_{d}$','interpreter','latex','Fontsize',15);
% 
% subplot(4,4,16);
% yyaxis left;plot(Periods,(RLibor_eu_vec(chain(Periods)).^(freq)-RLibor_us_vec(chain(Periods)).^(freq))*100,'Linewidth',2);
% ylabel('$\%$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('$R^{eu}_{libor}-R^{us}_{libor}$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end
% 
% % sim paths of quantities
% cc=cc+1;plottype='_simpathsquant';
% figure(cc);
% subplot(2,2,1);
% yyaxis left;plot(Periods,d_us_vec(chain(Periods)),'Linewidth',2);
% ylabel('$d_{us}$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('US real deposits $d_{us}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,2);
% yyaxis left;plot(Periods,d_eu_vec(chain(Periods)),'Linewidth',2);
% ylabel('$d_{eu}$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('EU real deposits $d_{eu}$','interpreter','latex','Fontsize',15);
% 
% subplot(2,2,3);
% yyaxis left;plot(Periods,b_vec(chain(Periods)),'Linewidth',2);
% ylabel('$b$','interpreter','latex');
% yyaxis right;plot(Periods,x_vec(chain(Periods)),'Linewidth',2);
% ylabel(x_lab,'interpreter','latex');
% xlim([Periods(1) Periods(end)]);
% xlabel('Periods','interpreter','latex');
% grid on;
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real loans $b$','interpreter','latex','Fontsize',15);
% 
% if printit==1
%     orient landscape;
%     set(gcf, 'PaperUnits', 'normalized');
%     set(gcf, 'PaperPosition', [-0.06 -0.06 1.12 1.12]);    
%     printsb(['fig' nameplot plottype]);
% end

%% Final Plots for paper
% 1) Combine dollar and euro liquidity ratio in one plot [x]
% 2) Combine Rm,Rm* in one plot
% 3) Exchange rate
% 4) BP.
% 5) CIP
% 6) Dollar and Euro deposit rates in one plot

% policy functions
cc = cc+1; 
x_vec = sigma_us_vec(:)/mean(sigma_us_vec);
x_lab = '$\sigma_{us}$';
% x_vec = (im_us_vec.^(freq)-1)*100;
% x_lab = '$i^{us}_m(\%)$';



% Begin Figures
% figure('Name','Liquidity Ratios','NumberTitle','off'); % Liquidity Ratios
% splot1(x_vec,mu_us_vec); hold on;
% splot2(x_vec,mu_eu_vec);
% xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% % ylabel('$liquidity$','interpreter','latex');
% formataxis(gca);
% h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
% set(h,'interpreter','latex');
% %scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
% %ylim(mean(mu_us_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% % title('Liquidity','interpreter','latex','Fontsize',15);
% ax = gca;
% % Specify the number of desired ticks on the x-axis
% desiredNumXTicks = 5; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
% yLimits = get(ax, 'YLim');
% customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
% set(ax, 'YTick', customYTicks);
% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% % Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
% print(gcf,'-dpdf',[foldername 'nlmod_mu']);

% Limit Graph

%limsigmas=[sigma_us_vec(1) sigma_us_vec(1)+0.6*(sigma_us_vec(end)-sigma_us_vec(1))];
limsigmas=[sigma_us_vec(1) sigma_us_vec(end-10)];


figure('Name','US Liquidity Ratio','NumberTitle','off'); % Liquidity Ratios
splot1(sigma_us_vec(index1),mu_us_vec(index1)); hold on;
splot2(sigma_us_vec(index2),mu_us_vec(index2));
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
% ylabel('$liquidity$','interpreter','latex');
formataxis(gca);
% h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_mu_us.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','EU Liquidity Ratio','NumberTitle','off'); % Liquidity Ratios
splot1(sigma_us_vec(index1),mu_eu_vec(index1)); hold on;
splot2(sigma_us_vec(index2),mu_eu_vec(index2));
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
% ylabel('$liquidity$','interpreter','latex');
formataxis(gca);
h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_mu.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

% figure('Name','Real Rates Dollar','NumberTitle','off'); % Endogenous policy rates
% splot1(x_vec(index1),(Rm_us_vec(index1).^(freq)-1)*rate_scale); hold on;
% splot2(x_vec(index1),(Rm_us_vec(index1).^(freq)-1)*rate_scale); 
% xlim([x_vec(1) x_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('bps','interpreter','latex');
% formataxis(gca);
% h=legend('$R^{m,us}$','$R^{m,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
% %scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
% %ylim(mean(mu_us_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% ax = gca;
% % Specify the number of desired ticks on the x-axis
% desiredNumXTicks = 5; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
% yLimits = get(ax, 'YLim');
% customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
% set(ax, 'YTick', customYTicks);
% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% % Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
% print(gcf,'-dpdf',[foldername 'nlmod_Rm']);

figure('Name','Real Dollar Rate','NumberTitle','off'); % Endogenous policy rates
splot1(sigma_us_vec(index1),(Rm_us_vec(index1).^(freq)-1)*rate_scale); hold on;
splot2(sigma_us_vec(index2),(Rm_us_vec(index2).^(freq)-1)*rate_scale); 
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_Rm_us.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Dollar Inflation','NumberTitle','off'); % Endogenous policy rates
splot1(sigma_us_vec(index1),(pi_us_vec(index1).^(freq)-1)*rate_scale); hold on;
splot2(sigma_us_vec(index2),(pi_us_vec(index2).^(freq)-1)*rate_scale); 
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_Rm_us.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Real Euro Rate','NumberTitle','off'); % Endogenous policy rates
splot1(sigma_us_vec(index1),(Rm_eu_vec(index1).^(freq)-1)*rate_scale); hold on;
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
splot2(sigma_us_vec(index2),(Rm_eu_vec(index2).^(freq)-1)*rate_scale); 
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_Rm_eu.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','all rates','NumberTitle','off'); % Endogenous policy rates
splot1(sigma_us_vec(index1),(Rm_us_vec(index1).^(freq)-1)*rate_scale); hold on;
splot2(sigma_us_vec(index2),(Rm_us_vec(index2).^(freq)-1)*rate_scale); 
splot3(sigma_us_vec(index1),(Rm_eu_vec(index1).^(freq)-1)*rate_scale); hold on;
splot4(sigma_us_vec(index2),(Rm_eu_vec(index2).^(freq)-1)*rate_scale); 
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20);
set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 5; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  


figure('Name','FX','NumberTitle','off');
e_mean=mean(e_euus_vec);
splot1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean); hold on;
splot2(sigma_us_vec(index2),e_euus_vec(index2)/e_mean); 
xlabel(x_lab,'interpreter','latex');
formataxis(gca);
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','east','FontSize',20,'AutoUpdate', 'off');
set(h,'interpreter','latex');
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_e.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end
area(sigma_us_vec(index1),yLimits(1)+(yLimits(2)-yLimits(1))*invp1/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base, 'DisplayName', 'Normal Regime'); hold on;
area(sigma_us_vec(index2),yLimits(1)+(yLimits(2)-yLimits(1))*invp2/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base2); hold on;
ylim([yLimits(1) yLimits(2)]);
xlim(limsigmas);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_FX.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Invariants','NumberTitle','off');
area(sigma_us_vec(index1),invp1/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base, 'DisplayName', 'Normal Regime'); hold on;
area(sigma_us_vec(index2),invp2/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base2); hold on;
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
formataxis(gca);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','northeast','FontSize',20,'AutoUpdate', 'off');
set(h,'interpreter','latex');
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
ylim([yLimits(1) yLimits(2)]);
if printit==1
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    exportgraphics(gcf, fullfile(foldername, 'nlmod_invariants.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

% figure('Name','Deposit Rates','NumberTitle','off');
% splot1(sigma_us_vec,(Rd_us_vec.^(freq)-1)*rate_scale); hold on;
% splot2(sigma_us_vec,(Rd_eu_vec.^(freq)-1)*rate_scale); 
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('bps','interpreter','latex');
% formataxis(gca);
% h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','southwest','FontSize',20);
% set(h,'interpreter','latex');
% %scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
% %ylim(mean(mu_us_vec)+[-scale,scale]);
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% ax = gca;
% % Specify the number of desired ticks on the x-axis
% desiredNumXTicks = 5; % Change this to the number of ticks you want
% xLimits = get(ax, 'XLim');
% customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
% set(ax, 'XTick', customXTicks);
% yLimits = get(ax, 'YLim');
% customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
% set(ax, 'YTick', customYTicks);
% desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% % Format the x-axis tick labels
% xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
% grid on;
% %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% % title('Real Policy Rates','interpreter','latex','Fontsize',15);
% set(gcf, 'PaperUnits', 'normalized');
% orient landscape;
% set(gcf, 'PaperPosition', [0 0 1 1]);  
% %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
% print(gcf,'-dpdf',[foldername 'nlmod_Rd']);

figure('Name','Bond Premium','NumberTitle','off')
splot1(sigma_us_vec(index1),(Rb_us_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale); hold on;
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
splot2(sigma_us_vec(index2),(Rb_us_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_BP.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Real Rate Differential','NumberTitle','off')
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
aux1=(Rm_eu_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale;
aux2=(Rm_eu_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale;
splot1(sigma_us_vec(index1),aux1); hold on;
splot2(sigma_us_vec(index2),aux2);
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlim(limsigmas);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
if printit==1
    %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
    % title('Real Policy Rates','interpreter','latex','Fontsize',15);
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    exportgraphics(gcf, fullfile(foldername, 'nlmod_UIP.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Real Rate Differential shade','NumberTitle','off')
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
aux1=(Rm_eu_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale;
aux2=(Rm_eu_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale;
diff=aux2-aux1;
for ii=1:length(index1)-1
    alphaVal1 = invp1(ii)/max(invp1);
    alphaVal2 = invp2(ii)/max(invp2);
    plot(sigma_us_vec(index1(ii:ii+1)),aux1(ii:ii+1),'Color',[0 0 1 alphaVal1],'LineWidth',2); hold on;
    plot(sigma_us_vec(index2(ii:ii+1)),aux2(ii:ii+1),'Color',[1 0 1 alphaVal2],'LineWidth',2);
end
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;

%% Intensity Plots
figure('Name','Real Rate Differential shade','NumberTitle','off')
plot(sigma_us_vec(index1),0*e_euus_vec(index1)+1, 'LineStyle', '--','Color','k','LineWidth',1); hold on;
aux1=e_euus_vec(index1)/mean(e_euus_vec);
aux2=e_euus_vec(index2)/mean(e_euus_vec);
for ii=1:length(index1)-1
    alphaVal1 = sqrt(invp1(ii)/max(invp1));
    alphaVal2 = sqrt(invp2(ii)/max(invp2));
    plot(sigma_us_vec(index1(ii:ii+1)),aux1(ii:ii+1),'Color',[0 0 1 alphaVal1],'LineWidth',2); hold on;
    plot(sigma_us_vec(index2(ii:ii+1)),aux2(ii:ii+1),'Color',[1 0 1 alphaVal2],'LineWidth',2);
end
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;

figure('Name','Real Rate Differential shade','NumberTitle','off')
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
aux1=(Rm_eu_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale;
aux2=(Rm_eu_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale;
diff=aux2-aux1;
for ii=1:length(index1)-1
    alphaVal1 = sqrt(invp1(ii)/max(invp1));
    alphaVal2 = sqrt(invp2(ii)/max(invp2));
    plot(sigma_us_vec(index1(ii:ii+1)),aux1(ii:ii+1),'Color',[0 0 1 alphaVal1],'LineWidth',2); hold on;
    plot(sigma_us_vec(index2(ii:ii+1)),aux2(ii:ii+1),'Color',[1 0 1 alphaVal2],'LineWidth',2);
end
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;


%%  Single Regime Plots
figure('Name','FX','NumberTitle','off');
% e_mean=mean(e_euus_vec);
% splot1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean); hold on;
splot1(sigma_us_vec(index2),e_euus_vec(index2)/e_mean); 
xlabel(x_lab,'interpreter','latex');
formataxis(gca);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
%h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','east','FontSize',20,'AutoUpdate', 'off');
% set(h,'interpreter','latex');
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_e_scramble.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end
area(sigma_us_vec(index1),yLimits(1)+(yLimits(2)-yLimits(1))*invp1/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base, 'DisplayName', 'Normal Regime'); hold on;
area(sigma_us_vec(index2),yLimits(1)+(yLimits(2)-yLimits(1))*invp2/max([invp1; invp2]),'FaceAlpha',0.5,'FaceColor',color_base2); hold on;
ylim([yLimits(1) yLimits(2)]);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_FX_scramble.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','Real Rate Differential','NumberTitle','off')
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
%aux1=(Rm_eu_vec(index1).^(freq)-Rm_us_vec(index1).^(freq))*rate_scale;
aux2=(Rm_eu_vec(index2).^(freq)-Rm_us_vec(index2).^(freq))*rate_scale;
%diff=aux2-aux1;
splot1(sigma_us_vec(index1),aux2); hold on;
% splot2(sigma_us_vec(index2),aux2);
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
ylabel('bps','interpreter','latex');
formataxis(gca);
%h=legend('$R^{d,us}$','$R^{d,eu}$','color','none','box','off','location','northeast','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
% Get the current axes handle
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
if printit==1
    %set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
    % title('Real Policy Rates','interpreter','latex','Fontsize',15);
    set(gcf, 'PaperUnits', 'normalized');
    orient landscape;
    set(gcf, 'PaperPosition', [0 0 1 1]);  
    %  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
    exportgraphics(gcf, fullfile(foldername, 'nlmod_UIP_scramble.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

figure('Name','US Liquidity Ratio','NumberTitle','off'); % Liquidity Ratios
% splot1(sigma_us_vec(index1),mu_us_vec(index1)); hold on;
splot1(sigma_us_vec(index2),mu_us_vec(index2));
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
xlabel(x_lab,'interpreter','latex');
% ylabel('$liquidity$','interpreter','latex');
formataxis(gca);
% h=legend('$\mu^{us}$','$\mu^{eu}$','color','none','box','off','location','southwest','FontSize',20);
% set(h,'interpreter','latex');
%scale = max(0.001,max(abs(mu_us_vec-mean(mu_us_vec))));
%ylim(mean(mu_us_vec)+[-scale,scale]);
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 10);
% title('Liquidity','interpreter','latex','Fontsize',15);
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
set(gcf, 'PaperPosition', [0 0 1 1]);  
%  print(gcf,'-dpdf',['raw_data_' num2str(printver)]);
if printit==1
    exportgraphics(gcf, fullfile(foldername, 'nlmod_mu_us_scramble.pdf'), 'ContentType', 'vector', 'Resolution', 1600);
end

%% Moments - Theoretical
bp_vec=(Rb_us_vec.^(freq)-Rm_us_vec.^(freq))*rate_scale;
cip_vec=(Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*rate_scale;
pi_us_ret_vec=(pi_us_vec.^(freq)-1)*rate_scale;

% CIP check
figure('Name','CIP+Ret comp','NumberTitle','off');
plot(sigma_us_vec(index1),0*Rm_eu_vec(index1), 'LineStyle', '--','Color','k','LineWidth',1); hold on;
splot1(sigma_us_vec(index1),cip_vec(index1)); hold on;
splot2(sigma_us_vec(index1),pi_us_ret_vec(index1));
splot3(sigma_us_vec(index2),cip_vec(index2)); 
splot4(sigma_us_vec(index2),pi_us_ret_vec(index2));
hold off;
%splot3(sigma_us_vec(index2),diff);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);

% Construct the regime specific vectors
e_euus_vec_r1=e_euus_vec(index1);
bp_vec_r1=bp_vec(index1);
cip_vec_r1=cip_vec(index1);
pi_us_ret_vec_r1=pi_us_ret_vec(index1);
e_euus_vec_r2=e_euus_vec(index2);
bp_vec_r2=bp_vec(index2);
cip_vec_r2=cip_vec(index2);
pi_us_ret_vec_r2=pi_us_ret_vec(index2);

% means:
E_e=e_euus_vec*invp;
E_bp=bp_vec*invp;
E_cip=cip_vec*invp;
E_pi_ret_us=pi_us_ret_vec*invp;
E_e_r1=e_euus_vec_r1*invp1;
E_bp_r1=bp_vec_r1*invp1;
E_cip_r1=cip_vec_r1*invp1;
E_pi_ret_us_r1=pi_us_ret_vec_r1*invp1;
E_e_r2=e_euus_vec_r2*invp2;
E_bp_r2=bp_vec_r2*invp2;
E_cip_r2=cip_vec_r2*invp2;
E_pi_ret_us_r2=pi_us_ret_vec_r2*invp2;

% standard deviations:
std_e=sqrt((e_euus_vec-E_e).^2*invp);
std_bp=sqrt((bp_vec-E_bp).^2*invp);
std_cip=sqrt((cip_vec-E_cip).^2*invp);
std_pi_ret_us=sqrt((pi_us_ret_vec-E_pi_ret_us).^2*invp);
std_e_r1=sqrt((e_euus_vec_r1-E_e_r1).^2*invp1);
std_bp_r1=sqrt((bp_vec_r1-E_bp_r1).^2*invp1);
std_cip_r1=sqrt((cip_vec_r1-E_cip_r1).^2*invp1);
std_pi_ret_us_r1=sqrt((pi_us_ret_vec_r1-E_pi_ret_us_r1).^2*invp1);
std_e_r2=sqrt((e_euus_vec_r2-E_e_r2).^2*invp2);
std_bp_r2=sqrt((bp_vec_r2-E_bp_r2).^2*invp2);
std_cip_r2=sqrt((cip_vec_r2-E_cip_r2).^2*invp2);
std_pi_ret_us_r2=sqrt((pi_us_ret_vec_r2-E_pi_ret_us_r2).^2*invp2);

% Moments Conditional on Jump
Dev_Mat=(((e_euus_vec').^(-1)*e_euus_vec)-1)*rate_scale;
Jump_mat=(Zprob_sigma_us.*Dev_Mat);
Jump_test=Zprob_sigma_us;
Dev_r1_to_r2=invp1'*(Jump_mat(1:N_sigma_us/2,N_sigma_us/2+1:end))/P(1,2)*ones(N_sigma_us/2,1);
%Dev_test=invp1'*(Jump_test(1:N_sigma_us/2,N_sigma_us/2+1:end))/P(1,2)*ones(N_sigma_us/2,1);
Dev_r2_to_r1=invp2'*(Jump_mat(N_sigma_us/2+1:end,1:N_sigma_us/2))/P(2,1)*ones(N_sigma_us/2,1);

filename = fullfile(foldername, 'Mod_SwitchDev_Moments.tex');
fid = fopen(filename, 'wt');
% devswith_moments = {...
    %'FX', 1, rho_e_sim, std_e/E_e, (E_e_r2./E_e_r1-1)*10000, rho_e_sim_r2, rho_e_sim_r1, (std_e_r2./std_e_r1); 
    %             '$\Delta$ FX', E_pi_ret_us, rho_pi_us_sim, std_pi_ret_us, E_pi_ret_us_r2-E_pi_ret_us_r1, rho_e_sim_r2, rho_e_sim_r1, (std_pi_ret_us_r2./std_pi_ret_us_r1);
    %             'BP', E_bp, rho_bp_sim, std_bp, E_bp_r2-E_bp_r1, rho_bp_sim_r2, rho_bp_sim_r1, (std_bp_r2./std_bp_r1);                 
               % , };

% fprintf(fid, '\\bottomrule\n');
%for ii = 1:size(model_moments, 1)
    fprintf(fid, '%.1f & %.1f \\\\', Dev_r1_to_r2,Dev_r2_to_r1);
%end
fclose(fid);

%% Moments - Simulation Based
% means:
E_e_sim=mean(e_euus_t);
E_bp_sim=mean(bp_us_t);
E_cip_sim=mean(uip_Rm_t);

% standard deviations:
std_e_sim=std(e_euus_t);
std_bp_sim=std(bp_us_t);
std_cip_sim=std(uip_Rm_t);

% autocorrelation:
aux=autocorr(e_euus_t);
rho_e_sim=aux(2);
aux=autocorr(bp_us_t);
rho_bp_sim=aux(2);
aux=autocorr(uip_Rm_t);
rho_cip_sim=aux(2);
aux=autocorr(pi_us_t);
rho_pi_us_sim=aux(2);

% Autocorrelation by Regime
index_r1=chain<=N_sigma_us/2;
index_r2=chain>N_sigma_us/2;

% e, bp, cip
cip_t=uip_Rm_t;
aux=autocorr(e_euus_t(index_r1));
rho_e_sim_r1=aux(2);
aux=autocorr(e_euus_t(index_r2));
rho_e_sim_r2=aux(2);
aux=autocorr(bp_us_t(index_r1));
rho_bp_sim_r1=aux(2);
aux=autocorr(bp_us_t(index_r2));
rho_bp_sim_r2=aux(2);
aux=autocorr(cip_t(index_r1));
rho_cip_sim_r1=aux(2);
aux=autocorr(cip_t(index_r2));
rho_cip_sim_r2=aux(2);
aux=autocorr(pi_us_t(index_r1));
rho_pi_us_sim_r1=aux(2);
aux=autocorr(pi_us_t(index_r2));
rho_pi_us_sim_r2=aux(2);



%% Print Model Moments to Table
% Define folder and file name
% filename = fullfile(foldername, 'Mod_Moments.tex');
% fid = fopen(filename, 'wt');

% Write the LaTeX table header
% fprintf(fid, '\\begin{table}[h!]\n');
% fprintf(fid, '\\centering\n');
% fprintf(fid, '\\begin{threeparttable}\n');
% fprintf(fid, '\\begin{tabular}{lcccccc}\n');
% fprintf(fid, '\\toprule\n');
% fprintf(fid, '\\textbf{Variable} & \\multicolumn{6}{c}{\\textbf{Moment}} \\\\\n');
% fprintf(fid, '\\cmidrule(lr){2-7}\n');
% fprintf(fid, ' & \\multicolumn{3}{c}{Unconditional} & \\multicolumn{3}{c}{Regime 1 / Regime 2} \\\\\n');
% fprintf(fid, '\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}\n');
% fprintf(fid, ' & Mean & Autocorrelation & Std & Mean & Autocorrelation & Std \\\\\n');
% fprintf(fid, '\\midrule\n');

% % Section for Data Moments
% fprintf(fid, '\\rowcolor{gray!20}\n');
% fprintf(fid, '\\multicolumn{7}{c}{\\textbf{Data Moments}} \\\\\n');

% % Replace the following with your actual data
% data_moments = {'FX', 0.123, 0.456, 0.789, 0.123, 0.456, 0.789; 
%                 'BP', 0.223, 0.556, 0.889, 0.223, 0.556, 0.889;
%                 'CIP', 0.323, 0.656, 0.989, 0.323, 0.656, 0.989};

% for ii = 1:size(data_moments, 1)
%     fprintf(fid, '%s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n', data_moments{ii, :});
% end

% Section for Model Moments
% fprintf(fid, '\\midrule\n');
% fprintf(fid, '\\rowcolor{gray!20}\n');
% fprintf(fid, '\\multicolumn{7}{c}{\\textbf{Model Moments}} \\\\\n');

% Replace the following with your actual model moments
filename = fullfile(foldername, 'Mod_Moments.tex');
fid = fopen(filename, 'wt');
model_moments = {'FX', 1, rho_e_sim, std_e/E_e, (E_e_r2./E_e_r1-1)*10000, rho_e_sim_r2, rho_e_sim_r1, (std_e_r2./std_e_r1); 
                 '$\Delta$ FX', E_pi_ret_us, rho_pi_us_sim, std_pi_ret_us, E_pi_ret_us_r2-E_pi_ret_us_r1, rho_e_sim_r2, rho_e_sim_r1, (std_pi_ret_us_r2./std_pi_ret_us_r1);
                 'BP', E_bp, rho_bp_sim, std_bp, E_bp_r2-E_bp_r1, rho_bp_sim_r2, rho_bp_sim_r1, (std_bp_r2./std_bp_r1);                 
                 'CIP', E_cip, rho_cip_sim, std_cip, E_cip_r2-E_cip_r1, rho_cip_sim_r2, rho_cip_sim_r1, (std_cip_r2./std_cip_r1)};

% fprintf(fid, '\\bottomrule\n');
for ii = 1:size(model_moments, 1)
    fprintf(fid, '%s & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', model_moments{ii, :});
end

% Close the table
% fprintf(fid, '\\bottomrule \n');
fclose(fid);

filename = fullfile(foldername, 'Mod_CIP_Moments.tex');
fid = fopen(filename, 'wt');
model_moments = {...
    %'FX', 1, rho_e_sim, std_e/E_e, (E_e_r2./E_e_r1-1)*10000, rho_e_sim_r2, rho_e_sim_r1, (std_e_r2./std_e_r1); 
    %             '$\Delta$ FX', E_pi_ret_us, rho_pi_us_sim, std_pi_ret_us, E_pi_ret_us_r2-E_pi_ret_us_r1, rho_e_sim_r2, rho_e_sim_r1, (std_pi_ret_us_r2./std_pi_ret_us_r1);
    %             'BP', E_bp, rho_bp_sim, std_bp, E_bp_r2-E_bp_r1, rho_bp_sim_r2, rho_bp_sim_r1, (std_bp_r2./std_bp_r1);                 
                 'CIP (model)', E_cip, rho_cip_sim, std_cip, E_cip_r2-E_cip_r1, rho_cip_sim_r2, rho_cip_sim_r1, (std_cip_r2./std_cip_r1)};

% fprintf(fid, '\\bottomrule\n');
for ii = 1:size(model_moments, 1)
    fprintf(fid, '%s & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \n', model_moments{ii, :});
end

% Close the table
% fprintf(fid, '\\bottomrule \n');
fclose(fid);
% fprintf(fid, '\\end{tabular}\n');

% % Add the units and notes
% fprintf(fid, '\\begin{tablenotes}\n');
% fprintf(fid, '\\small\n');
% fprintf(fid, '\\item \\textbf{Units:} \\newline Mean: Percentage (\\%%), Autocorrelation: Unitless, Std: Percentage (\\%%)\n');
% fprintf(fid, '\\end{tablenotes}\n');
% fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '\\caption{Data vs. Model Moments.}\n');
% fprintf(fid, '\\label{tab:example}\n');
% fprintf(fid, '\\end{table}\n');

% Close the file

%% Phase Diagram-Like Plots
% Interpolation
N_r1_1=1;
N_r2_1=60;
N_r1_2=20;
N_sim =N_r1_1+N_r2_1+N_r1_2;
lsigma_0=mu_sigma_us_r2+2;
lsigma_typath=zeros(N_sim,1);
lsigma_typath(1)=lsigma_0;
for ii=1:N_sim
    if ii<=N_r1_1
        lsigma_typath(ii+1)=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*lsigma_typath(ii,1);
    elseif ii>(N_r1_1+N_r2_1)
        lsigma_typath(ii+1)=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*lsigma_typath(ii,1);
    else
        lsigma_typath(ii+1)=mu_sigma_us_r2*(1-rho_sigma_us_r2)+rho_sigma_us_r2*lsigma_typath(ii,1);
    end
end
sigma_typath=exp(lsigma_typath);
e_euus_typath=0*sigma_typath;
for ii=1:N_sim
    if ii<=N_r1_1
        e_euus_typath(ii)=interp1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean,sigma_typath(ii),'spline');
    elseif ii>N_r1_1+N_r2_1
        e_euus_typath(ii)=interp1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean,sigma_typath(ii),'spline');
    else
        e_euus_typath(ii)=interp1(sigma_us_vec(index2),e_euus_vec(index2)/e_mean,sigma_typath(ii),'spline');
    end
end

figure('Name','FX Path','NumberTitle','off');
e_mean=mean(e_euus_vec);
splot1(sigma_us_vec(index1),e_euus_vec(index1)/e_mean); hold on;
splot2(sigma_us_vec(index2),e_euus_vec(index2)/e_mean); 
for ii=1:N_sim-1
   plot(sigma_typath(ii:ii+1),e_euus_typath(ii:ii+1),'Color',[0.5 0 0.5],'LineWidth',1,'Marker','<','MarkerFaceColor', [0.5 0 0.5], 'MarkerEdgeColor', [0.5 0 0.5], 'MarkerSize', 10); hold on;
   % annotation('arrow',[sigma_typath(ii:ii+1)],[e_euus_typath(ii:ii+1)]); hold on;
end
xlabel(x_lab,'interpreter','latex');
formataxis(gca);
xlim([sigma_us_vec(1) sigma_us_vec(end)]);
h=legend('Normal Regime','Scrambling Regime','color','none','box','off','location','east','FontSize',20,'AutoUpdate', 'off');
set(h,'interpreter','latex');
ax = gca;
% Specify the number of desired ticks on the x-axis
desiredNumXTicks = 3; % Change this to the number of ticks you want
xLimits = get(ax, 'XLim');
customXTicks = linspace(xLimits(1), xLimits(2), desiredNumXTicks);
set(ax, 'XTick', customXTicks);
yLimits = get(ax, 'YLim');
customYTicks = linspace(yLimits(1), yLimits(2), desiredNumXTicks);
set(ax, 'YTick', customYTicks);
desiredDecimalPlaces = 1; % Change this to your desired number of decimal places
% Format the x-axis tick labels
xtickformat(ax, sprintf('%%.%df', desiredDecimalPlaces));
grid on;
%set(gca, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
% title('Real Policy Rates','interpreter','latex','Fontsize',15);
set(gcf, 'PaperUnits', 'normalized');
orient landscape;
hold on;

% figure('Name','FX Typical Path','NumberTitle','off')
% lsigma_0=mu_sigma_us_r1;
% lsigma_typath=zeros(36,1);
% lsigma_typath(1)=lsigma_0;
% %mu_sigma_us_r1 ;
% %rho_sigma_us_r1;
% for ii=1:5
%     aux=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*lsigma_typath(ii,1);
%     [~,index]=min(abs(sigma_us_vec(1:N_sigma_us/2)-exp(aux)));
%     lsigma_typath(ii+1)=log(sigma_us_vec(index));
% end
% for ii=6:22
%     aux=mu_sigma_us_r2*(1-rho_sigma_us_r2)+rho_sigma_us_r2*lsigma_typath(ii,1);
%     [~,index]=min(abs(sigma_us_vec(N_sigma_us/2+1:end)-exp(aux)));
%     lsigma_typath(ii+1)=log(sigma_us_vec(N_sigma_us/2+index));
% end
% for ii=23:30
%     aux=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*lsigma_typath(ii,1);
%     [~,index]=min(abs(sigma_us_vec(1:N_sigma_us/2)-exp(aux)));
%     lsigma_typath(ii+1)=log(sigma_us_vec(index));
% end
% %splot3(sigma_us_vec(index2),diff);
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('bps','interpreter','latex');
% formataxis(gca);
% 
% figure('Name','FX Typical Path','NumberTitle','off')
% sigma_0=mu_sigma_us_r1;
% sigma_typath=zeros(30,1);
% index=zeros(30,1);
% sigma_typath(1)=sigma_0;
% aux=mu_sigma_us_r1*(1-rho_sigma_us_r1)+rho_sigma_us_r1*sigma_typath(1);
% [~,index0]=min(abs(sigma_us_vec(1:N_sigma_us/2)-exp(aux)));
% index(1)=index0;
% sigma_typath(1)=sigma_us_vec(index0);
% %mu_sigma_us_r1 ;
% %rho_sigma_us_r1;
% for ii=1:5
%     [~,index_aux]=max(Zprob_sigma_us(index(ii),1:N_sigma_us/2));
%     index(ii+1)=index_aux;
%     sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% end
% [~,index0]=min(abs(sigma_us_vec(N_sigma_us/2+1:end)-(sigma_typath(ii+1))));
% index(ii+1)=N_sigma_us/2+index0;
% sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% for ii=6:22
%     [~,index_aux]=max(Zprob_sigma_us(index(ii),N_sigma_us/2+1:end));
%     index(ii+1)=N_sigma_us/2+index_aux;
%     sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% end
% [~,index0]=min(abs(sigma_us_vec(1:N_sigma_us/2)-(sigma_typath(ii+1))));
% index(ii+1)=index0;
% sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% for ii=23:30
%     [~,index_aux]=max(Zprob_sigma_us(index(ii),1:N_sigma_us/2));
%     index(ii+1)=index_aux;
%     sigma_typath(ii+1)=sigma_us_vec(index(ii+1));
% end
% figure
% for ii=1:5
%     line(sigma_us_vec(index1(ii:ii+1)),e_euus_vec(index1(ii:ii+1)),'Color',[1 0 0],'LineWidth',2,'Marker','<'); hold on;
% end
% for ii=6
%     line(sigma_us_vec(index1(ii)),e_euus_vec(index2(ii:ii+1)),'Color',[0 0 1],'LineWidth',2,'Marker','>'); hold on;
% end
% for ii=7:22
%     line(sigma_us_vec(index2(ii:ii+1)),e_euus_vec(index2(ii:ii+1)),'Color',[0 0 1],'LineWidth',2,'Marker','>'); hold on;
% end
% for ii=23:30
%     line(sigma_us_vec(index1(ii:ii+1)),e_euus_vec2(index(ii:ii+1)),'Color',[1 0 0],'LineWidth',2,'Marker','<'); hold on;
% end
% %splot3(sigma_us_vec(index2),diff);
% xlim([sigma_us_vec(1) sigma_us_vec(end)]);
% xlabel(x_lab,'interpreter','latex');
% ylabel('bps','interpreter','latex');
% formataxis(gca);