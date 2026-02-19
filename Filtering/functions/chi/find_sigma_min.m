function sigma_min = find_sigma_min(mu, ploss, theta_target)
% Find sigma such that theta(mu, ploss, sigma) = theta_target
% Uses bisection since theta is monotonic in sigma

    % theta function
    theta_func = @(sig) (sig * exp(-ploss/sig * mu)) / (mu + sig * exp(-ploss/sig * mu));
    
    % Bisection search
    sig_lo = 0.01;
    sig_hi = 20.0;
    
    for iter = 1:100
        sig_mid = (sig_lo + sig_hi) / 2;
        theta_mid = theta_func(sig_mid);
        
        if theta_mid < theta_target
            sig_lo = sig_mid;
        else
            sig_hi = sig_mid;
        end
        
        if (sig_hi - sig_lo) < 1e-6
            break;
        end
    end
    sigma_min = sig_hi;
end