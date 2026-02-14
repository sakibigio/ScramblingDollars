function echi_d = Echi_d(mu, ploss, sigma, iota, lambda, eta, matching_type)
% ECHI_D - Expected marginal liquidity yield for deposits
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
%   echi_d : E[χ^d] = MassUnd(-μ) * χ⁻(θ) + (-MassUnd(-μ)) * χ⁺(θ)
%
% This is the deposit premium in the model.

if nargin < 6 || isempty(eta)
    eta = 0.5;
end
if nargin < 7 || isempty(matching_type)
    matching_type = 0;
end

% Distribution functions (don't depend on matching type)
MassUnd = @(omega, p, s) (omega <= 0) .* (omega - s./p) .* exp(p./s .* omega) .* p + ...
    (omega > 0) .* (-(1-p) .* (omega + s./(1-p)) .* exp(-(1-p)./s .* omega));

Smin = @(m, p, s) s .* exp(-p./s .* m);
Spl  = @(m, p, s) m + s .* exp(-p./s .* m);

% Market tightness from distribution
th = Smin(mu, ploss, sigma) ./ Spl(mu, ploss, sigma);

% Chi functions (depend on matching type)
chi_p = Chi_p(th, iota, lambda, eta, matching_type);
chi_m = Chi_m(th, iota, lambda, eta, matching_type);

% Expected deposit liquidity
echi_d = MassUnd(-mu, ploss, sigma) .* chi_m + (-MassUnd(-mu, ploss, sigma)) .* chi_p;

end
