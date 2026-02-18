function chip = Chi_p(theta, iota, lambda, eta, matching_type)
% CHI_P - Liquidity yield for surplus positions
%
% Inputs:
%   theta         : Market tightness θ (can be vector)
%   iota          : Interest spread (R^w - R^m)
%   lambda        : Matching efficiency
%   eta           : Bargaining power (default = 0.5)
%   matching_type : 0 = Leontief, 1 = Cobb-Douglas (default = 0)
%
% Output:
%   chip : Liquidity yield coefficient χ⁺
%
% Formula (Leontief):
%   χ⁺ = ι * (θ̄ - θ̄^η * θ^{1-η}) / (θ̄ - 1)

if nargin < 4 || isempty(eta)
    eta = 0.5;
end
if nargin < 5 || isempty(matching_type)
    matching_type = 0;  % Default to Leontief
end

% Get terminal tightness
barth = bartheta(theta, lambda, matching_type);

% Preallocate
chip = zeros(size(theta));

if matching_type == 0
    %% LEONTIEF - Original formula from old code
    % chi_p = iota * (bartheta - bartheta^eta * theta^(1-eta)) / (bartheta - 1)
    
    % General case: bartheta ~= 1
    idx = abs(barth - 1) > 1e-10;
    if any(idx(:))
        th = theta(idx);
        bt = barth(idx);
        chip(idx) = iota .* (bt - bt.^eta .* th.^(1-eta)) ./ (bt - 1);
    end
    
    % Special case: bartheta = 1 (use L'Hopital or limit)
    idx_one = abs(barth - 1) <= 1e-10;
    if any(idx_one(:))
        chip(idx_one) = iota .* (1 - eta);
    end

elseif matching_type == 1
    %% COBB-DOUGLAS
    % Similar structure but with CD matching
    idx = abs(barth - 1) > 1e-10;
    if any(idx(:))
        th = theta(idx);
        bt = barth(idx);
        chip(idx) = iota .* (bt - bt.^eta .* th.^(1-eta)) ./ (bt - 1);
    end
    
    idx_one = abs(barth - 1) <= 1e-10;
    if any(idx_one(:))
        chip(idx_one) = iota .* (1 - eta);
    end
end

end
