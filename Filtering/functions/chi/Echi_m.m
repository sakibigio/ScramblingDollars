function echi_m = Echi_m(mu, ploss, sigma, iota, lambda, eta, matching_type, varrho)
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
%   varrho        : Threshold shift (default = 0). Effective threshold is mu-varrho.
%
% Output:
%   echi_m : E[χ^m] = (1-F(-(μ-ϱ))) * χ⁺(θ) + F(-(μ-ϱ)) * χ⁻(θ)
%
% This is the bond premium in the model.

if nargin < 6 || isempty(eta)
    eta = 0.5;
end
if nargin < 7 || isempty(matching_type)
    matching_type = 0;
end
if nargin < 8 || isempty(varrho)
    varrho = 0;
end

% Distribution aggregates at effective threshold (mu - varrho)
[~, ~, th, F_val, ~] = distribution_aggregates(mu, ploss, sigma, varrho);

% Chi functions (depend on matching type)
chi_p = Chi_p(th, iota, lambda, eta, matching_type);
chi_m = Chi_m(th, iota, lambda, eta, matching_type);

% Expected marginal liquidity
% Banks in surplus (prob 1-F) earn χ⁺
% Banks in deficit (prob F)   pay χ⁻
echi_m = (1 - F_val) .* chi_p + F_val .* chi_m;

end
