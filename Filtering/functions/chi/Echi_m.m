function echi_m = Echi_m(mu, ploss, sigma, iota, lambda, eta, matching_type)
% ECHI_M - Expected marginal liquidity yield for reserves
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
% Output:
%   echi_m : E[χ^m] = (1-F(-μ)) * χ⁺(θ) + F(-μ) * χ⁻(θ)
%
% This is the bond premium in the model.

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

% Market tightness from distribution
th = Smin(mu, ploss, sigma) ./ Spl(mu, ploss, sigma);

% Chi functions (depend on matching type)
chi_p = Chi_p(th, iota, lambda, eta, matching_type);
chi_m = Chi_m(th, iota, lambda, eta, matching_type);

% Expected marginal liquidity
% Banks in surplus (prob 1-F(-μ)) earn χ⁺
% Banks in deficit (prob F(-μ)) pay χ⁻
echi_m = (1 - F(-mu, ploss, sigma)) .* chi_p + F(-mu, ploss, sigma) .* chi_m;

end
