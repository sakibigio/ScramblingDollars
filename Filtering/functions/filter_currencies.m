% filter_currencies.m
% Extracted from main_filter.m (lines 450–529) on 2026-02-28
% Called by: main_filter.m
%
% Contents:
%   - Filter sigma for 8 other currencies (AU,CA,JP,NZ,NO,SW,CH,UK)
%   - Compute interbank variables, rates, FX for each
%   - Store results in currency-specific variables

%% Pre-allocate
sigma_c_TED_t = zeros(T, 1);
sigma_c_TED_flag = zeros(T, 1);
Echi_m_c_t = zeros(T, 1);
Echi_d_c_t = zeros(T, 1);
Chi_p_psi_c_t = zeros(T, 1);
Echi_m_eu_bp_t = zeros(T, 1);
Echi_d_eu_bp_t = zeros(T, 1);
BP_c_t = zeros(T, 1);
Rd_c_t = zeros(T, 1);
TED_c_t = zeros(T, 1);
theta_c_t = zeros(T, 1);
psi_c_t = zeros(T, 1);
Smin_c_t = zeros(T, 1);
DW_c_t = zeros(T, 1);
FF_c_t = zeros(T, 1);
Q_c_t = zeros(T, 1);

for cc = 1:numel(curlist)
    eval(['Rm_' curlist{cc} '=exp(im_' curlist{cc} ');']);
    eval(['Rm_c=Rm_' curlist{cc} ';']);
    
    for tt = 1:T
        mu_eu_yt = exp(mu_eu(tt));
        eval(['target=min_test_eu+TED_s_' curlist{cc} '_t(tt)*abs_scale;']);        
        
        % For Cobb-Douglas, use mu-dependent initial guess
        if matching_type == 1
            sigma_c_guess = max(sigma_eu_t(tt), mu_eu_yt + 0.15);
        else
            sigma_c_guess = sigma_eu_t(tt);
        end
    
        sigma_res = @(sigma) Chi_p_psi(mu_eu_yt, ploss_eu, sigma, iota_eu, lambda_eu, eta, matching_type) * 1e4 * 12 - target;
        [sigma_out, ~, exitflag, ~] = fsolve(@(sigma) sigma_res(sigma), sigma_c_guess, optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'MaxFunEval', 1e9, 'MaxIter', 1e6));
        
        if exitflag > 0
            sigma_c_TED_t(tt) = sigma_out;
        else
            if tt > 1
                sigma_c_TED_t(tt) = sigma_c_TED_t(tt-1);
            else
                sigma_c_TED_t(tt) = sigma_eu_t(tt);  % Use EU value as fallback
            end
        end
        sigma_c_TED_flag(tt) = exitflag;
        TED_c_t(tt) = Chi_p_psi(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);

        Echi_m_c_t(tt) = Echi_m(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);
        Echi_d_c_t(tt) = Echi_d(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);
        Chi_p_psi_c_t(tt) = Chi_p_psi(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);
        
        [~, theta_c_t(tt), psi_c_t(tt), Smin_c_t(tt), DW_c_t(tt), FF_c_t(tt), Q_c_t(tt)] = Chi_sys(mu_eu_yt, ploss_eu, sigma_c_TED_t(tt), iota_eu, lambda_eu, eta, matching_type);
    end
    clear mu_eu_yt target;

    BP_c_t = Echi_m_c_t;
    Rb_c_t = Rm_c + Echi_m_c_t;
    Rd_c_t = Rm_c + Echi_m_c_t + Echi_d_c_t;
    UIP_c_t = Rm_c - Rm_us;
    riskprm_c_t = (Rb_us_t) ./ (Rb_c_t) - 1;
    CIP_c_t = UIP_c_t + Rm_c .* riskprm_c_t;
    
    % Price and FX calculations
    p_c_t = (M_eu ./ (Rd_c_t.^(1/zeta_eu))) ./ mu_eu_t;
    inv_e_c_t = p_us_t ./ p_c_t;
    f_c_t = (1 + riskprm_c_t) ./ inv_e_c_t;

    eval(['sigma_' curlist{cc} '_TED_flag=sigma_c_TED_flag;']);
    eval(['sigma_' curlist{cc} '_t=sigma_c_TED_t;']);
    eval(['TED_' curlist{cc} '_t=TED_c_t;']);
    eval(['BP_' curlist{cc} '_t=BP_c_t;']);
    eval(['riskprm_' curlist{cc} '_t=riskprm_c_t;']);
    eval(['UIP_' curlist{cc} '_t=UIP_c_t;']);
    eval(['CIP_' curlist{cc} '_t=CIP_c_t;']);
    eval(['inv_e_' curlist{cc} '_t=inv_e_c_t;']);
    eval(['p_' curlist{cc} '_t=p_c_t;']);
    eval(['f_' curlist{cc} '_t=f_c_t;']);
end
