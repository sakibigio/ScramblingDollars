%% Equilibrium System
% Equilibrium in Asset Markets
% barB=bard_us*(1-mu_us)+bard_eu*(1-mu_eu);

% Indifference in both currencies:
% Rm_eu+Echi_m_eu-Rm_us-Echi_m_us=;
% Equilibrium in Asset Markets - Expressed in Ratios
% bard_us+mu_us+mu_eu*nu_d=(1+nu_d) -> ((1+nu_d)- (nu_b+mu_us))/nu
mu_eu_ame=@(mu_us) (1-(nu_b+mu_us*nu_us_d))./nu_eu_d;

% -------------------------------------------------------------------------
% Interbank Definitions
% -------------------------------------------------------------------------
% Requires: eta, matching_type (from params.m)
%
% This file defines anonymous function handles that the rest of main_LFX.m
% uses. The matching-dependent functions (chi_p, chi_m, Echi_m, Echi_d,
% psi) now delegate to the standalone .m files (Chi_p.m, Chi_m.m, etc.)
% which handle both Leontief and Cobb-Douglas matching.
%
% The distribution functions (F, Smin, Spl, theta, MassUnd, Econd) do NOT
% depend on matching type and remain as anonymous functions.

% -------------------------------------------------------------------------
% Distribution Functions (matching-type independent)
% -------------------------------------------------------------------------
F=@(omega,p,sigma) (omega<=0).*exp(p./sigma.*omega).*p+...
    (omega>0).*(p+(1-p).*(1-exp(-(1-p)./sigma.*omega)));
Smin=@(mu,p,sigma) sigma.*exp(-p./sigma.*mu)            ;
Spl =@(mu,p,sigma) mu+sigma.*exp(-p./sigma.*mu)         ;
theta=@(mu,p,sigma) Smin(mu,p,sigma)./Spl(mu,p,sigma);
MassUnd=@(omega,p,sigma) (omega<=0).*(omega-sigma./p).*exp(p./sigma.*omega).*p+...
   (omega>0).*(-(1-p).*(omega+sigma./(1-p)).*exp(-(1-p)./sigma.*omega));
Econd=@(omega,p,sigma) MassUnd(omega,p,sigma)./F(omega,p,sigma);

% -------------------------------------------------------------------------
% Matching-Dependent Functions — Wrappers to .m files
% -------------------------------------------------------------------------
% chi_p, chi_m: 3-arg handles (theta, iota, lambda) -> Chi_p.m, Chi_m.m
% These close over eta and matching_type from the workspace.
chi_p=@(th,iota,lambda) Chi_p(th, iota, lambda, eta, matching_type);
chi_m=@(th,iota,lambda) Chi_m(th, iota, lambda, eta, matching_type);

% psi: 2-arg handle (theta, lambda) -> Psi_p.m
psi=@(th,lambda) Psi_p(th, lambda, matching_type);

% Echi_m, Echi_d: 5-arg handles (mu, p, sigma, iota, lambda)
% NOTE: We first capture file function handles to avoid name collision,
% since the anonymous variables will shadow the .m file names.
Echi_d_file = @Echi_d;
Echi_m_file = @Echi_m;

Echi_m=@(mu,p,sigma,iota,lambda) Echi_m_file(mu, p, sigma, iota, lambda, eta, matching_type);
Echi_d=@(mu,p,sigma,iota,lambda) Echi_d_file(mu, p, sigma, iota, lambda, eta, matching_type);

% Echi: Total expected liquidity (uses chi_p, chi_m handles defined above)
Echi =@(mu,p,sigma,iota,lambda)  F(-mu,p,sigma).*chi_m(theta(mu,p,sigma),iota,lambda).*...
    (mu-varrho+(1-varrho).*Econd(-mu,p,sigma))+...
    (1-F(-mu,p,sigma)).*chi_p(theta(mu,p,sigma),iota,lambda).*...
    (mu-varrho+(1-varrho).*(-Econd(-mu,p,sigma).*F(-mu,p,sigma)./(1-F(-mu,p,sigma))))...
    ;

% -------------------------------------------------------------------------
% LFX Equilibrium Equations
% -------------------------------------------------------------------------
% Residual in M indifference
mind_res=@(mu_eu,mu_us) Rm_eu-Rm_us...
    -(1-F(-mu_us,ploss_us,sigma_us)).*chi_p(theta(mu_us,ploss_us,sigma_us),iota_us,lambda_us)...
    -F(-mu_us,ploss_us,sigma_us).*chi_m(theta(mu_us,ploss_us,sigma_us),iota_us,lambda_us)...
    +(1-F(-mu_eu,ploss_eu,sigma_eu)).*chi_p(theta(mu_eu,ploss_eu,sigma_eu),iota_eu,lambda_eu)...
    +F(-mu_eu,ploss_eu,sigma_eu).*chi_m(theta(mu_eu,ploss_eu,sigma_eu),iota_eu,lambda_eu);

mu_eu_ind=@(mu_us) fsolve(@(x) mind_res(x,mu_us),0.1,optimoptions('fsolve','TolFun',1e-16));

% Equilibrium Condition
mind_res_eq=@(mu_us) mind_res(mu_eu_ame(mu_us),mu_us);
mu_us_star_f=@(z) fsolve(@(x) mind_res_eq(x),0.1,optimoptions('fsolve','Display','off','TolFun',1e-16));

% Equilibrium Returns
Rd_us_f=@(mu_us) Rm_us+Echi_m(mu_us,ploss_us,sigma_us,iota_us,lambda_us)...
    +Echi_d(mu_us,ploss_us,sigma_us,iota_us,lambda_us);

Rd_eu_f=@(mu_us) Rm_us+Echi_m(mu_us,ploss_us,sigma_us,iota_us,lambda_us)...
    +Echi_d(mu_eu_ame(mu_us),ploss_eu,sigma_eu,iota_eu,lambda_eu);

% Loan condition
Rb_us_f=@(mu_us) Rm_us+Echi_m(mu_us,ploss_us,sigma_us,iota_us,lambda_us);
Rb_eu_f=@(mu_us) Rm_us+Echi_m(mu_us,ploss_us,sigma_us,iota_us,lambda_us);

% Price System Equations
p_eu_f=@(mu_eu) M_eu./(bard_eu.*mu_eu);
inv_e_f=@(p_eu,mu_us) M_us./(bard_us.*mu_us)./p_eu;
e_euus_f=@(inv_e) inv(inv_e);
p_us_f=@(p_eu,inv_e) p_eu.*inv_e;
p_us_f2=@(mu_us) M_us./(bard_us.*mu_us);

% Side Rates
RBond_us=@(mu_us) Rm_us+chi_p(theta(mu_us,ploss_us,sigma_us),iota_us,lambda_us);
RBond_eu=@(mu_eu) Rm_eu+chi_p(theta(mu_eu,ploss_eu,sigma_eu),iota_eu,lambda_eu);
RLibor_us=@(mu_us) Rm_us+chi_p(theta(mu_us,ploss_us,sigma_us),iota_us,lambda_us)./...
    psi(theta(mu_us,ploss_us,sigma_us),lambda_us);
RLibor_eu=@(mu_eu) Rm_eu+chi_p(theta(mu_eu,ploss_eu,sigma_eu),iota_eu,lambda_eu)./...
    psi(theta(mu_eu,ploss_eu,sigma_eu),lambda_eu);
