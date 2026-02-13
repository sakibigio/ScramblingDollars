function [barth, theta_t] = bartheta(theta, lambda, matching_type, t)
% BARTHETA - Terminal market tightness for OTC matching
%
% Inputs:
%   theta         : Initial market tightness θ₀ (can be vector)
%   lambda        : Matching efficiency parameter
%   matching_type : 0 = Leontief, 1 = Cobb-Douglas
%   t             : Time in [0,1] (optional, default = 0)
%
% Outputs:
%   barth   : Terminal tightness θ̄ (at τ = 1 or T_max for CD)
%   theta_t : Tightness at time t
%
% Reference: Table 1 in Bianchi-Bigio (JET 2025)

if nargin < 4
    t = 0;
end

% Preallocate
barth = zeros(size(theta));
theta_t = zeros(size(theta));

if matching_type == 0
    %% LEONTIEF
    % θ < 1: θ̄ = θ / (θ + (1-θ)e^λ)
    % θ > 1: θ̄ = 1 + (θ-1)e^λ
    
    idx_low = theta < 1;
    if any(idx_low(:))
        th = theta(idx_low);
        barth(idx_low) = th ./ (th + (1-th).*exp(lambda));
        theta_t(idx_low) = th ./ (th + (1-th).*exp(lambda.*t));
    end
    
    idx_high = theta > 1;
    if any(idx_high(:))
        th = theta(idx_high);
        barth(idx_high) = 1 + (th-1).*exp(lambda);
        theta_t(idx_high) = 1 + (th-1).*exp(lambda.*t);
    end
    
    idx_one = theta == 1;
    if any(idx_one(:))
        barth(idx_one) = 1;
        theta_t(idx_one) = 1;
    end

elseif matching_type == 1
    %% COBB-DOUGLAS
    % θ̄ = [(1+√θ)e^{-λT} - (1-√θ)]² / [(1+√θ)e^{-λT} + (1-√θ)]²
    % T = min(T*, 1), T* = (1/λ)log|(1+√θ)/(1-√θ)|
    
    idx_low = theta < 1;
    if any(idx_low(:))
        th = theta(idx_low);
        sqrt_th = sqrt(th);
        alpha = (1 + sqrt_th) ./ (1 - sqrt_th);
        T_star = log(abs(alpha)) ./ lambda;
        T_max = min(T_star, 1);
        
        % Terminal tightness
        exp_T = alpha .* exp(-lambda .* T_max);
        barth(idx_low) = ((exp_T - 1) ./ (exp_T + 1)).^2;
        barth(idx_low & T_star < 1) = 0;  % Early stop
        
        % Tightness at time t
        t_eff = min(t, T_max);
        exp_t = alpha .* exp(-lambda .* t_eff);
        theta_t(idx_low) = ((exp_t - 1) ./ (exp_t + 1)).^2;
    end
    
    idx_high = theta > 1;
    if any(idx_high(:))
        th = theta(idx_high);
        sqrt_th = sqrt(th);
        alpha = (1 + sqrt_th) ./ (1 - sqrt_th);
        T_star = log(abs(alpha)) ./ lambda;
        T_max = min(T_star, 1);
        
        exp_T = alpha .* exp(-lambda .* T_max);
        barth(idx_high) = ((exp_T - 1) ./ (exp_T + 1)).^2;
        barth(idx_high & T_star < 1) = Inf;  % Early stop
        
        t_eff = min(t, T_max);
        exp_t = alpha .* exp(-lambda .* t_eff);
        theta_t(idx_high) = ((exp_t - 1) ./ (exp_t + 1)).^2;
    end
    
    idx_one = theta == 1;
    if any(idx_one(:))
        barth(idx_one) = 1;
        theta_t(idx_one) = 1;
    end
    
else
    error('matching_type must be 0 (Leontief) or 1 (Cobb-Douglas)');
end

end
