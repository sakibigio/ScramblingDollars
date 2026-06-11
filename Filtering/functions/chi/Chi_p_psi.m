function [ted, chi_p_val, psi_p_val] = Chi_p_psi(mu, ploss, sigma, iota, lambda, eta, matching_type, varrho)
% CHI_P_PSI - TED spread calculation
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
% Outputs:
%   ted       : TED spread = χ⁺(θ) / Ψ⁺(θ)
%   chi_p_val : χ⁺(θ) value
%   psi_p_val : Ψ⁺(θ) value
%
% The TED spread measures the funding premium in the interbank market.

if nargin < 6 || isempty(eta)
    eta = 0.5;
end
if nargin < 7 || isempty(matching_type)
    matching_type = 0;
end
if nargin < 8 || isempty(varrho)
    varrho = 0;
end

% Market tightness from distribution aggregates at effective threshold (mu - varrho)
[~, ~, th] = distribution_aggregates(mu, ploss, sigma, varrho);

% Chi and Psi functions (depend on matching type)
chi_p_val = Chi_p(th, iota, lambda, eta, matching_type);
psi_p_val = Psi_p(th, lambda, matching_type);

% TED spread = chi_p / psi (NOT multiplication!)
ted = chi_p_val ./ psi_p_val;

end
