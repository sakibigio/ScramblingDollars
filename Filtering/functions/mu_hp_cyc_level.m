function mu_out = mu_hp_cyc_level(mu_in, lam, center)
%MU_HP_CYC_LEVEL Replace the liquidity ratio by its HP cyclical component + center.
%
%   mu_out = MU_HP_CYC_LEVEL(mu_in, lam, center)
%
%   mu_in is the LOG liquidity ratio (as stored in LFX_data.mat). This routine
%   works on the LEVEL ratio L = exp(mu_in): it extracts a very-low-frequency
%   two-sided HP trend (smoothing parameter `lam`), keeps the cyclical deviation
%   cyc = L - trend, recenters it at `center` (e.g. 0.2 -- a typical observed
%   liquidity-ratio level), and returns log(center + cyc) so the downstream
%   filter (which consumes exp(mu)) sees a stationary liquidity-ratio series.
%
%   Used by main_filter.m when mu_hp_detrend == 1 to run the filter against a
%   stationary mu instead of the raw (trending) liquidity ratio.
%
%   Inputs:
%       mu_in  - vector of log liquidity ratios
%       lam    - HP smoothing parameter (e.g. 1e6 for a very-low-freq trend)
%       center - level around which the cyclical ratio is recentered (e.g. 0.2)
%   Output:
%       mu_out - log of the recentered, detrended level, same size as mu_in

    L     = exp(mu_in(:));
    tau   = local_hp_trend(L, lam);     % very-low-frequency trend of the level
    cyc   = L - tau;                    % cyclical component (mean ~ 0)
    L_new = center + cyc;               % stationary level centered at `center`

    nbad = sum(L_new <= 0);
    if nbad > 0
        warning('mu_hp_cyc_level:nonpos', ...
            '%d non-positive level(s) after detrend; flooring at 1e-6.', nbad);
        L_new(L_new <= 0) = 1e-6;
    end

    mu_out = reshape(log(L_new), size(mu_in));
end

function tau = local_hp_trend(y, lam)
%LOCAL_HP_TREND Two-sided Hodrick-Prescott trend via direct linear solve.
%   Returns tau minimizing sum((y - tau).^2) + lam * sum(diff(tau,2).^2).
    y = y(:);
    T = numel(y);
    D = diff(speye(T), 2);              % (T-2) x T second-difference operator
    tau = (speye(T) + lam * (D.' * D)) \ y;
end
