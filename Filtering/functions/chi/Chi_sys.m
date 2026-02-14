function [echi_m, theta_val, psi_val, Smin_val, DW_val, FF_val, Q_val] = Chi_sys(mu, ploss, sigma, iota, lambda, eta, matching_type)
% CHI_SYS - Full interbank system calculation
%
% Inputs:
%   mu            : Reserve ratio μ
%   ploss         : Probability of loss (p)
%   sigma         : Volatility of withdrawals
%   iota          : Interest spread (R^w - R^m)
%   lambda        : Matching efficiency
%   eta           : Bargaining power (default = 0.5)
%   matching_type : 0 = Leontief, 1 = Cobb-Douglas (default = 0)
%
% Outputs:
%   echi_m    : Expected marginal liquidity E[χ^m]
%   theta_val : Market tightness θ = S⁻/S⁺
%   psi_val   : Matching probability Ψ⁺
%   Smin_val  : Deficit side measure S⁻
%   DW_val    : Discount window usage (1 - Ψ⁻) * F(-μ)
%   FF_val    : Fed Funds volume Ψ⁺ * (1-F(-μ))
%   Q_val     : Total interbank quantity

if nargin < 6 || isempty(eta)
    eta = 0.5;
end
if nargin < 7 || isempty(matching_type)
    matching_type = 0;
end

% Distribution functions (don't depend on matching type)
F = @(omega, p, s) (omega <= 0) .* exp(p./s .* omega) .* p + ...
    (omega > 0) .* (p + (1-p) .* (1 - exp(-(1-p)./s .* omega)));

Smin = @(m, p, s) s .* exp(-p./s .* m);
Spl  = @(m, p, s) m + s .* exp(-p./s .* m);

% Compute distribution values
Smin_val = Smin(mu, ploss, sigma);
Spl_val  = Spl(mu, ploss, sigma);
F_val    = F(-mu, ploss, sigma);

% Market tightness from distribution
theta_val = Smin_val ./ Spl_val;

% Chi and Psi functions (depend on matching type)
chi_p = Chi_p(theta_val, iota, lambda, eta, matching_type);
chi_m = Chi_m(theta_val, iota, lambda, eta, matching_type);
psi_p = Psi_p(theta_val, lambda, matching_type);
psi_m = Psi_m(theta_val, lambda, matching_type);

% Expected marginal liquidity
echi_m = (1 - F_val) .* chi_p + F_val .* chi_m;

% Matching probability (surplus side)
psi_val = psi_p;

% Discount window usage: banks in deficit who don't find a match
DW_val = (1 - psi_m) .* F_val;

% Fed Funds volume: matched trades
% Volume = Ψ⁺ * S⁺ = Ψ⁻ * S⁻ (by matching identity)
FF_val = psi_p .* Spl_val;

% Total interbank quantity
Q_val = FF_val;

end
