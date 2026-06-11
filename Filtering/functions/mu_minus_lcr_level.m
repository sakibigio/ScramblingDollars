function mu_out = mu_minus_lcr_level(mu_in, lcr_pct)
%MU_MINUS_LCR_LEVEL Liquidity ratio in excess of the LCR/HQLA requirement.
%
%   mu_out = MU_MINUS_LCR_LEVEL(mu_in, lcr_pct)
%
%   mu_in is the LOG liquidity ratio (as stored in LFX_data.mat). Like
%   mu_hp_cyc_level, this routine works on the LEVEL ratio L = exp(mu_in), but
%   instead of an HP trend it nets out the empirical regulatory baseline: the
%   FEDS-Note LCR / HQLA-to-assets series (supplied in PERCENT). The filter then
%   sees the EXCESS liquidity ratio above the requirement,
%
%       L_excess = exp(mu_in) - lcr_pct/100,
%
%   and the function returns log(L_excess) so the downstream filter (which
%   consumes exp(mu)) sees that excess series. This is ADDITIVE in levels (a
%   coverage GAP), not a ratio/multiple.
%
%   Inputs:
%       mu_in    - vector of log liquidity ratios
%       lcr_pct  - LCR / HQLA-to-assets series, PERCENT, same length as mu_in
%   Output:
%       mu_out   - log of the excess level, same size as mu_in

    L     = exp(mu_in(:));
    lcr   = lcr_pct(:) / 100;           % percent -> ratio
    L_new = L - lcr;                    % excess liquidity ratio above requirement

    nbad = sum(L_new <= 0 | isnan(L_new));
    if nbad > 0
        warning('mu_minus_lcr_level:nonpos', ...
            '%d non-positive/NaN excess level(s); flooring at 1e-6.', nbad);
        L_new(L_new <= 0 | isnan(L_new)) = 1e-6;
    end

    mu_out = reshape(log(L_new), size(mu_in));
end
