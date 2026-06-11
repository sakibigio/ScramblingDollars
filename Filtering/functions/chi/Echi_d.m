function echi_d = Echi_d(mu, ploss, sigma, iota, lambda, eta, matching_type, varrho)
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
%   varrho        : Threshold shift (default = 0). Effective threshold is mu-varrho.
%
% Output:
%   echi_d : E[χ^d] = MassUnd(-(μ-ϱ)) * χ⁻(θ) + (-MassUnd(-(μ-ϱ))) * χ⁺(θ)
%
% This is the deposit premium in the model.

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
[~, ~, th, ~, MassUnd_val] = distribution_aggregates(mu, ploss, sigma, varrho);

% Chi functions (depend on matching type)
chi_p = Chi_p(th, iota, lambda, eta, matching_type);
chi_m = Chi_m(th, iota, lambda, eta, matching_type);

% Expected deposit liquidity
echi_d = MassUnd_val .* chi_m + (-MassUnd_val) .* chi_p;

end
