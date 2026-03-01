% solve_global.m
% Extracted from main_LFX.m (lines 478–846) on 2026-02-28
% Called by: main_LFX.m
%
% Contents:
%   - Create exogenous shock paths (sigma, rates, deposits)
%   - Re-load equations for global solution
%   - Solve global equilibrium via fsolve (feqm_vec)
%   - Compute liquidity values, rates, bank profits
%   - Solution tests
%   - Save initial guess and global solution to data/

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
save(fullfile('data','initguess.mat'), 'x0');
eval(['save ' fullfile('data',['global_sol' nameplot '.mat'])]);
