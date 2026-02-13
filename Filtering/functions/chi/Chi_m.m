function chi_m = Chi_m(theta, iota, lambda, eta, matching_type)
% CHI_M - Liquidity yield for deficit positions
%
% Inputs:
%   theta         : Market tightness θ (can be vector)
%   iota          : Interest spread (R^w - R^m)
%   lambda        : Matching efficiency
%   eta           : Bargaining power (default = 0.5)
%   matching_type : 0 = Leontief, 1 = Cobb-Douglas (default = 0)
%
% Output:
%   chi_m : Liquidity yield coefficient χ⁻
%
% Formula:
%   χ⁻ = ι * (θ̄/θ)^η * (1 - θ^η * θ̄^{1-η}) / (1 - θ̄)

if nargin < 4 || isempty(eta)
    eta = 0.5;
end
if nargin < 5 || isempty(matching_type)
    matching_type = 0;  % Default to Leontief
end

% Get terminal tightness
barth = bartheta(theta, lambda, matching_type);

% Preallocate
chi_m = zeros(size(theta));

% General case: theta ~= 1
idx = abs(theta - 1) > 1e-10;
if any(idx(:))
    th = theta(idx);
    bt = barth(idx);
    chi_m(idx) = iota .* (bt./th).^eta .* (1 - th.^eta .* bt.^(1-eta)) ./ (1 - bt);
end

% Special case: theta = 1
idx_one = abs(theta - 1) <= 1e-10;
if any(idx_one(:))
    chi_m(idx_one) = iota .* (1 - eta.*(1 - exp(-lambda)));
end

end
