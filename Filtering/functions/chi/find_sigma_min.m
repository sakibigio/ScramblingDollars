function sigma_min = find_sigma_min(mu, ploss, theta_target, varrho)
% Find sigma such that theta(mu, ploss, sigma; varrho) = theta_target.
% Uses bisection. For the effective threshold m_eff = mu - varrho > 0,
% theta is monotonically increasing in sigma (theta -> 0 as sigma -> 0,
% theta -> 1 as sigma -> inf), so the standard bisection direction holds.

    if nargin < 4 || isempty(varrho)
        varrho = 0;
    end

    % theta from distribution_aggregates (handles both signs of m_eff)
    theta_func = @(sig) get_theta(mu, ploss, sig, varrho);

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

function th = get_theta(mu, ploss, sigma, varrho)
    [~, ~, th] = distribution_aggregates(mu, ploss, sigma, varrho);
end
