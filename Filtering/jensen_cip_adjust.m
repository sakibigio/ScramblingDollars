%% jensen_cip_adjust.m
%  Recompute the model CIP with the state-by-state Jensen adjustment, WITHOUT
%  re-solving the model. Uses the committed-run expected-inflation objects
%  (data/expected_inflation_cd.mat, computed by run_expected_inflation.m):
%    Epi_x_grid    = E[(1+pi_x')|s]        (arithmetic)
%    pi_inv_x_grid = E[1/(1+pi_x')|s]      (harmonic side)
%  The current tables' CIP uses the arithmetic convention per leg; the
%  risk-neutral nominal CIP deflates each leg by its own-currency harmonic
%  mean. The additive Jensen adjustment per state is
%    J(s) = [ (E[1/(1+pi_us')|s] - 1/E[(1+pi_us')|s])
%           - (E[1/(1+pi_eu')|s] - 1/E[(1+pi_eu')|s]) ] * 12e4   (annualized bps)
%  (up to a (1+i(s)) multiplier ~ 1). Being additive, its regime-conditional
%  means add directly to the existing table's conditional CIP entries.
addpath('functions'); addpath('data');

S = load('data/expected_inflation_cd.mat');
G = load('data/global_solcalibration_dynare_sigma_us.mat', ...
         'index1', 'index2', 'sigma_us_vec', 'invp', 'invp1', 'invp2');

% Grid consistency: regime indexing is structural (fixed N), but verify sigma grids agree.
if isfield(S, 'sigma_us_vec') && numel(S.sigma_us_vec) == numel(G.sigma_us_vec)
    gd = max(abs(S.sigma_us_vec(:) - G.sigma_us_vec(:)));
    fprintf('sigma grid max |diff| across files: %.2e\n', gd);
end

J_us = S.pi_inv_us_grid(:) - 1 ./ S.Epi_us_grid(:);   % per-period (monthly) net
J_eu = S.pi_inv_eu_grid(:) - 1 ./ S.Epi_eu_grid(:);
J    = (J_us - J_eu) * 12e4;                           % annualized bps
J_us_bps = J_us * 12e4;  J_eu_bps = J_eu * 12e4;

erg = S.erg(:);  erg = erg / sum(erg);
i1 = G.index1; i2 = G.index2;                          % r1 = normal, r2 = scrambling
p1 = erg(i1) / sum(erg(i1));                           % conditional invariant dists
p2 = erg(i2) / sum(erg(i2));

fprintf('\n--- Jensen wedge per leg (annualized bps) ---\n');
fprintf('            ergodic     normal    scrambling\n');
fprintf('US leg:   %8.3f  %9.3f  %11.3f\n', erg'*J_us_bps, p1'*J_us_bps(i1), p2'*J_us_bps(i2));
fprintf('EU leg:   %8.3f  %9.3f  %11.3f\n', erg'*J_eu_bps, p1'*J_eu_bps(i1), p2'*J_eu_bps(i2));
fprintf('CIP  J:   %8.3f  %9.3f  %11.3f\n', erg'*J,        p1'*J(i1),        p2'*J(i2));

dJ = p2'*J(i2) - p1'*J(i1);
fprintf('\nConditional Jensen contribution to CIP (scr - nor): %+.2f bps\n', dJ);
fprintf('Baseline table conditional CIP diff: 14.5 bps -> adjusted: %.1f bps  (data: 45.4)\n', 14.5 + dJ);

% Profile in sigma: where does the wedge live?
sv = S.sigma_us_vec(:);
[~, ord] = sort(sv(i1));
fprintf('\nJ (bps) at selected sigma states (normal block):\n');
pick = round(linspace(1, numel(i1), 8));
for k = pick(:)'
    st = i1(ord(k));
    fprintf('  sigma=%7.3f   J_us=%8.3f  J_eu=%8.3f  J_cip=%8.3f\n', sv(st), J_us_bps(st), J_eu_bps(st), J(st));
end
[~, ord2] = sort(sv(i2));
fprintf('J (bps) at selected sigma states (scrambling block):\n');
for k = pick(:)'
    st = i2(ord2(k));
    fprintf('  sigma=%7.3f   J_us=%8.3f  J_eu=%8.3f  J_cip=%8.3f\n', sv(st), J_us_bps(st), J_eu_bps(st), J(st));
end

save('data/jensen_cip_adjust.mat', 'J', 'J_us_bps', 'J_eu_bps', 'dJ', 'erg', 'i1', 'i2');
