function F=feqm_vec(x,Echi_d,Echi_m,p_us_f,ploss_eu_vec,ploss_us_vec,sigma_eu_vec,sigma_us_vec,lambda_eu_vec,lambda_us_vec,Theta_b_vec,epsilon_b,Theta_d_eu_vec,Theta_d_us_vec,iw_eu_vec,iw_us_vec,im_eu_vec,im_us_vec,zeta_eu,zeta_us,M_eu,M_us,Q_mat,N_s)
    mu_eu_vec = x(1,:);
    mu_us_vec = x(2,:);
    Rd_eu_vec = x(3,:);
    Rd_us_vec = x(4,:);
    Rb_us_vec = x(5,:);
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
    
    F = zeros(15,N_s);
    F(1,:) = Rd_eu_vec-Echi_d(mu_eu_vec,ploss_eu_vec,sigma_eu_vec,iota_eu_vec,lambda_eu_vec)-Rd_us_vec+Echi_d(mu_us_vec,ploss_us_vec,sigma_us_vec,iota_us_vec,lambda_us_vec);
    F(2,:) = Rm_eu_vec+Echi_m(mu_eu_vec,ploss_eu_vec,sigma_eu_vec,iota_eu_vec,lambda_eu_vec)-Rm_us_vec-Echi_m(mu_us_vec,ploss_us_vec,sigma_us_vec,iota_us_vec,lambda_us_vec);
    F(3,:) = Rb_us_vec-Rd_us_vec+Echi_d(mu_us_vec,ploss_us_vec,sigma_us_vec,iota_us_vec,lambda_us_vec);
    F(4,:) = Rm_eu_vec+Echi_m(mu_eu_vec,ploss_eu_vec,sigma_eu_vec,iota_eu_vec,lambda_eu_vec)-Rd_us_vec+Echi_d(mu_us_vec,ploss_us_vec,sigma_us_vec,iota_us_vec,lambda_us_vec);
    F(5,:) = Rb_us_vec-((nu_vec.*(1-mu_eu_vec)+1-mu_us_vec).*d_us_vec./Theta_b_vec).^epsilon_b;
    if epsilon_b==-Inf
        F(5,:) = 1-((nu_vec.*(1-mu_eu_vec)+1-mu_us_vec).*d_us_vec./Theta_b_vec);
    end
    F(6,:) = nu_vec-Theta_d_us_vec.^(-1).*(Rd_us_vec).^(-1./zeta_us).*Theta_d_eu_vec.*(Rd_eu_vec).^(1./zeta_eu);
    F(7,:) = d_us_vec-Theta_d_us_vec.*(Rd_us_vec).^(1./zeta_us);
    if zeta_us==Inf
        F(6,:) = nu_vec-Theta_d_us_vec.^(-1).*Theta_d_eu_vec.*(Rd_eu_vec).^(1./zeta_eu);
        F(7,:) = d_us_vec-Theta_d_us_vec;
    end
    if zeta_eu==Inf
        F(6,:) = nu_vec-Theta_d_us_vec.^(-1).*(Rd_us_vec).^(-1./zeta_us).*Theta_d_eu_vec;
    end
    if zeta_us==Inf && zeta_eu==Inf
        F(6,:) = nu_vec-Theta_d_us_vec.^(-1).*Theta_d_eu_vec;
        F(7,:) = d_us_vec-Theta_d_us_vec;
    end
    F(8,:)    = Rm_eu_vec-im_eu_vec./pi_eu_vec;
    F(9,:)    = Rm_us_vec-im_us_vec./pi_us_vec;
    F(10,:)   = iota_eu_vec-(iw_eu_vec-im_eu_vec)./pi_eu_vec;
    F(11,:)   = iota_us_vec-(iw_us_vec-im_us_vec)./pi_us_vec;
    F(12,:)   = pi_eu_vec-(Q_mat*p_eu_vec(:))'./p_eu_vec; % p_eu_vec_in
    F(13,:)   = pi_us_vec-(Q_mat*p_us_vec(:))'./p_us_vec; % p_us_vec_in
    d_eu_vec  = Theta_d_eu_vec.*(Rd_eu_vec).^(1./zeta_eu);
    F(14,:)   = p_eu_vec-M_eu./(d_eu_vec.*mu_eu_vec);
    inv_e_vec = M_us./(d_us_vec.*mu_us_vec)./p_eu_vec;
    F(15,:)   = p_us_vec-p_us_f(p_eu_vec,inv_e_vec);
end