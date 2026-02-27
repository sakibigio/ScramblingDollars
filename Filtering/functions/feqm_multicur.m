function [F,temp1,temp2]=feqm_multicur(x,Echi_d,Echi_m,Rm,Rm_us,ploss,ploss_us,sigma,sigma_us,iota,iota_us,lambda,lambda_us,Rd_us,mu_us)
mu = x(1);
Rd = x(2);
F = zeros(2,1);
F(1) = Rd-Echi_d(mu,ploss,sigma,iota,lambda)-Rd_us+Echi_d(mu_us,ploss_us,sigma_us,iota_us,lambda_us);
F(2) = Rm+Echi_m(mu,ploss,sigma,iota,lambda)-Rm_us-Echi_m(mu_us,ploss_us,sigma_us,iota_us,lambda_us);
mugrid = 0:0.01:1;
temp1 = Echi_m(mugrid,ploss,sigma,iota,lambda);
temp2 = Rm_us+Echi_m(mu_us,ploss_us,sigma_us,iota_us,lambda_us)-Rm;
end