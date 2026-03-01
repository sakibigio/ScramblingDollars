% solve_steady_state.m
% Extracted from main_LFX.m (lines 93–336) on 2026-02-28
% Called by: main_LFX.m
%
% Contents:
%   - Partial equilibrium diagnostics (iso-dollar premium plot)
%   - Steady state graphic (contour of E[chi] difference)
%   - Calibration via feqm_calibrate (200bps target, exchange rate, mu_us)
%   - Steady state solution via feqm
%   - Multi-currency steady states (8 currencies)
%   - Solution tests (return equalization, deposit premium, bank profits)
%   - Report steady state variables
%   - Save M_curr.mat to data/

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

%% Solve Steady State Equilibrium
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
    eval(['Rm_' curlist{j} ' = Rm_temp;']);
    eval(['mu_' curlist{j} '_ss = x(1);']);
    eval(['Rd_' curlist{j} '_ss = x(2);']);
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
eval(['save ' fullfile('data','M_curr.mat') ' ' savelist ' MBS_ss sigmass_us sigmass_eu Theta_dss_eu Theta_dss_us;']);
