function [Smin_val, Spl_val, theta, F_meff, MassUnd_meff] = distribution_aggregates(mu, ploss, sigma, varrho)
% DISTRIBUTION_AGGREGATES - Closed-form aggregates of the asymmetric-Laplace
% liquidity-shock distribution evaluated at the effective threshold (mu - varrho).
%
% The bank's net position is (mu - varrho) + omega, with omega drawn from an
% asymmetric Laplace with left-tail mass ploss and scale sigma:
%   f(w) = (p^2/s) * exp(p/s * w)        for w <= 0
%   f(w) = ((1-p)^2/s) * exp(-(1-p)/s*w) for w >  0
% (mean zero; reduces to symmetric Laplace at p = 1/2).
%
% Inputs:
%   mu     : Reserve ratio (scalar or vector)
%   ploss  : Probability of withdrawal loss (left-tail mass)
%   sigma  : Scale of the shock distribution
%   varrho : Threshold shift (optional, default 0). The effective reserve
%            threshold for deficit/surplus is (mu - varrho).
%
% Outputs:
%   Smin_val     : E[((mu-varrho)+omega)^-]   (aggregate deficit mass)
%   Spl_val      : E[((mu-varrho)+omega)^+]   (aggregate surplus mass)
%   theta        : Market tightness Smin/Spl
%   F_meff       : F(-(mu-varrho), ploss, sigma) -- prob of deficit
%   MassUnd_meff : MassUnd(-(mu-varrho), ploss, sigma) -- deposit-side kernel

if nargin < 4 || isempty(varrho)
    varrho = 0;
end

m_eff = mu - varrho;

% Aggregate masses: piecewise closed form on sign of m_eff.
% (For m_eff>=0 the integration stays in the left-density region;
%  for m_eff<0 the surplus integral stays in the right-density region.
%  Both branches reduce to (Smin, Spl) = (sigma, sigma) at m_eff = 0.)
expo_lo = sigma .* exp(-ploss     ./ sigma .* m_eff);   % valid when m_eff >= 0
expo_hi = sigma .* exp((1-ploss)  ./ sigma .* m_eff);   % valid when m_eff <  0
idx_pos = (m_eff >= 0);
Smin_val = idx_pos .* expo_lo               + (~idx_pos) .* (-m_eff + expo_hi);
Spl_val  = idx_pos .* (m_eff + expo_lo)     + (~idx_pos) .* expo_hi;

% Market tightness
theta = Smin_val ./ Spl_val;

% Distribution functions evaluated at omega = -m_eff
om = -m_eff;
F_meff = (om <= 0) .* exp(ploss./sigma .* om) .* ploss + ...
         (om >  0) .* (ploss + (1-ploss) .* (1 - exp(-(1-ploss)./sigma .* om)));

MassUnd_meff = (om <= 0) .* (om - sigma./ploss) .* exp(ploss./sigma .* om) .* ploss + ...
               (om >  0) .* (-(1-ploss) .* (om + sigma./(1-ploss)) .* exp(-(1-ploss)./sigma .* om));

end
