function F=feqm_calibrate(x,Echi_d,Echi_m,Rm_eu,Rm_us,ploss_eu,ploss_us,iota_eu,iota_us,lambda_eu,lambda_us,Theta_b,epsilon_b,Theta_d_eu,Theta_d_us,zeta_eu,zeta_us,M_eu,M_us,target,share)
mu_eu = x(1);
mu_us = x(2);
Rd_eu = x(3);
Rd_us = x(4);
Rb_us = x(5);
d_us = x(6);
nu = x(7);
sigma_eu = x(8);
sigma_us = x(9);
Theta_d_eu = x(10);
Theta_d_us = x(10);
F = zeros(10,1);
F(1) = Rd_eu-Echi_d(mu_eu,ploss_eu,sigma_eu,iota_eu,lambda_eu)-Rd_us+Echi_d(mu_us,ploss_us,sigma_us,iota_us,lambda_us);
F(2) = Rm_eu+Echi_m(mu_eu,ploss_eu,sigma_eu,iota_eu,lambda_eu)-Rm_us-Echi_m(mu_us,ploss_us,sigma_us,iota_us,lambda_us);
F(3) = Rb_us-Rd_us+Echi_d(mu_us,ploss_us,sigma_us,iota_us,lambda_us);
F(4) = Rm_eu+Echi_m(mu_eu,ploss_eu,sigma_eu,iota_eu,lambda_eu)-Rd_us+Echi_d(mu_us,ploss_us,sigma_us,iota_us,lambda_us);
F(5) = Rb_us-((nu*(1-mu_eu)+1-mu_us)*d_us/Theta_b)^epsilon_b;
if epsilon_b==-Inf
    F(5) = 1-((nu*(1-mu_eu)+1-mu_us)*d_us/Theta_b);
end
%F(6) = nu-Theta_d_us^(-1)*(Rd_us)^(-1/zeta_us)*Theta_d_eu*(Rd_eu)^(1/zeta_eu);
F(6) = nu-Theta_d_eu*(Rd_eu)^(1/zeta_eu)/d_us;
F(7) = d_us-Theta_d_us*(Rd_us)^(1/zeta_us)/(1+share);
if zeta_us==Inf
    F(6) = nu-Theta_d_us^(-1)*Theta_d_eu*(Rd_eu)^(1/zeta_eu);
    F(7) = d_us-Theta_d_us;
end
if zeta_eu==Inf
    F(6) = nu-Theta_d_us^(-1)*(Rd_us)^(-1/zeta_us)*Theta_d_eu;
end
if zeta_us==Inf && zeta_eu==Inf
    F(6) = nu-Theta_d_us^(-1)*Theta_d_eu;
    F(7) = d_us-Theta_d_us;
end
d_eu = Theta_d_eu*(Rd_eu)^(1/zeta_eu);
p_eu = M_eu/(d_eu*mu_eu);
inv_e = M_us/(d_us*mu_us)/p_eu;
e_euus = inv(inv_e);
F(8) = (Rb_us-Rm_us)*1e4*12-target(1);
%F(8) = (Rb_us-Rd_us)*1e4*12-target(1);
F(9) = log(e_euus)-target(2);
F(10) = log(mu_us)-target(3);
end