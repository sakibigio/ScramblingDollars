%% CIP = DLP / E[1/(1+pi_eu')]: the proposition's comparability scaling,
%  state-by-state from the promoted (nopen) solution, along the filtered path.
addpath('functions'); addpath('data');
G = load('data/global_solcalibration_dynare_sigma_us.mat', 'Q_mat', 'p_eu_vec', 'sigma_us_vec');
pv = G.p_eu_vec(:); sv = G.sigma_us_vec(:);
pi_inv_eu = (G.Q_mat * (1 ./ pv)) .* pv;      % E[p/p' | s] = E[1/(1+pi_eu')|s]
i1 = 1:76; i2 = 77:152;

H = load('hybrid_cip_check.mat');              % CIP_f (filter DLP mapping, bps), cip_d, p
RW = readtable('RW_shock.csv'); sig = RW.sigma_us(:);
p  = H.p(:);
clamp = @(x,g) min(max(x, min(g)), max(g));
d_nor = interp1(sv(i1), pi_inv_eu(i1), clamp(sig, sv(i1)));
d_scr = interp1(sv(i2), pi_inv_eu(i2), clamp(sig, sv(i2)));
defl_t = (1 - p).*d_nor + p.*d_scr;            % E[1/(1+pi_eu')] along path

CIP_adj = H.CIP_f ./ defl_t;                   % correct comparability scaling
adj_bps = CIP_adj - H.CIP_f;

cip_d = H.cip_d; ok = isfinite(cip_d) & cip_d ~= 0; scr = p > 0.5;
fprintf('Euro deflator along path: mean %.6f | range [%.6f, %.6f]\n', mean(defl_t), min(defl_t), max(defl_t));
fprintf('Adjustment (bps): mean %+.3f | std %.3f | range [%+.2f, %+.2f]\n', mean(adj_bps), std(adj_bps), min(adj_bps), max(adj_bps));
fprintf('  normal months mean %+.3f | scrambling months mean %+.3f\n', mean(adj_bps(~scr)), mean(adj_bps(scr)));
fprintf('\nFit vs data CIP:\n');
fprintf('  unadjusted:  corr %.4f  RMSE %.3f bps\n', corr(H.CIP_f(ok), cip_d(ok)), sqrt(mean((H.CIP_f(ok)-cip_d(ok)).^2)));
fprintf('  eu-deflated: corr %.4f  RMSE %.3f bps\n', corr(CIP_adj(ok), cip_d(ok)), sqrt(mean((CIP_adj(ok)-cip_d(ok)).^2)));
save('euro_deflator_adjust.mat', 'CIP_adj', 'defl_t', 'adj_bps');
