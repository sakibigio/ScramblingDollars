function [ted, chi_p_val, psi_p_val] = Chi_p_psi(mu, ploss, sigma, iota, lambda, eta, matching_type)
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
%
% Outputs:
%   ted       : TED spread = χ⁺(θ) * Ψ⁺(θ)
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

% Distribution functions (don't depend on matching type)
Smin = @(m, p, s) s .* exp(-p./s .* m);
Spl  = @(m, p, s) m + s .* exp(-p./s .* m);

% Market tightness from distribution
th = Smin(mu, ploss, sigma) ./ Spl(mu, ploss, sigma);

% Chi and Psi functions (depend on matching type)
chi_p_val = Chi_p(th, iota, lambda, eta, matching_type);
psi_p_val = Psi_p(th, lambda, matching_type);

% TED spread
ted = chi_p_val .* psi_p_val;

end
