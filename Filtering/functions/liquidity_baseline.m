function [mubar_us_t, mubar_eu_t] = liquidity_baseline(dates)
% LIQUIDITY_BASELINE - Annual baseline liquidity ratios for US and EU.
% Broadcasts a year-level table to a monthly series matching `dates`.
%
% Captures the secular increase in liquidity ratios driven by regulatory
% regime (Basel III LCR, APP, QE/QT cycles). Used in main_filter.m to
% subtract a year-specific reference from exp(mu) so the filter sees
% deviations from the regulatory baseline rather than raw levels.
%
% Input:
%   dates - vector of MATLAB datenums (any monthly grid)
% Output:
%   mubar_us_t, mubar_eu_t - vectors of length numel(dates), with the
%                            baseline broadcast to each month's calendar year.
%
% Years outside [2001, 2025] are clamped to the nearest endpoint.

    table_year   = (2001:2025)';
    table_mu_us  = [0.16; 0.16; 0.15; 0.15; 0.15; 0.15; 0.15; 0.17; 0.18; ...
                    0.18; 0.18; 0.18; 0.19; 0.20; 0.21; 0.22; 0.23; 0.24; ...
                    0.24; 0.27; 0.27; 0.25; 0.24; 0.24; 0.24];
    table_mu_eu  = [0.14; 0.14; 0.13; 0.13; 0.13; 0.13; 0.13; 0.15; 0.18; ...
                    0.18; 0.17; 0.18; 0.19; 0.20; 0.22; 0.24; 0.25; 0.26; ...
                    0.27; 0.30; 0.32; 0.30; 0.28; 0.26; 0.25];

    dates = dates(:);
    yrs = year(dates);
    mubar_us_t = zeros(size(dates));
    mubar_eu_t = zeros(size(dates));
    for k = 1:numel(dates)
        if yrs(k) < table_year(1)
            idx = 1;
        elseif yrs(k) > table_year(end)
            idx = numel(table_year);
        else
            idx = find(table_year == yrs(k), 1);
        end
        mubar_us_t(k) = table_mu_us(idx);
        mubar_eu_t(k) = table_mu_eu(idx);
    end
end
