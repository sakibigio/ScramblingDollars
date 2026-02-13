function psi_m = Psi_m(theta, lambda, matching_type)
% PSI_M - Matching probability for deficit side
%
% Inputs:
%   theta         : Market tightness θ (can be vector)
%   lambda        : Matching efficiency
%   matching_type : 0 = Leontief, 1 = Cobb-Douglas (default = 0)
%
% Output:
%   psi_m : Matching probability Ψ⁻

if nargin < 3 || isempty(matching_type)
    matching_type = 0;
end

psi_m = zeros(size(theta));

if matching_type == 0
    %% LEONTIEF
    base_prob = 1 - exp(-lambda);
    
    idx_low = theta <= 1;
    psi_m(idx_low) = base_prob;
    
    idx_high = theta > 1;
    psi_m(idx_high) = base_prob ./ theta(idx_high);

elseif matching_type == 1
    %% COBB-DOUGLAS
    idx_low = theta < 1;
    if any(idx_low(:))
        th = theta(idx_low);
        sqrt_th = sqrt(th);
        alpha = (1 + sqrt_th) ./ (1 - sqrt_th);
        T_max = min(log(abs(alpha)) ./ lambda, 1);
        psi_m(idx_low) = 1 - exp(-lambda .* T_max) .* ...
            ((alpha - exp(lambda .* T_max)) ./ (alpha - 1)).^2;
    end
    
    idx_high = theta > 1;
    if any(idx_high(:))
        th = theta(idx_high);
        sqrt_th = sqrt(th);
        alpha = (1 + sqrt_th) ./ (1 - sqrt_th);
        T_max = min(log(abs(alpha)) ./ lambda, 1);
        psi_m(idx_high) = 1 - exp(-lambda .* T_max) .* ...
            ((alpha - exp(lambda .* T_max)) ./ (alpha - 1)).^2;
    end
    
    idx_one = theta == 1;
    psi_m(idx_one) = 1 - exp(-lambda);
end

end
