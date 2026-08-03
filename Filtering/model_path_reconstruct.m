%% model_path_reconstruct.m
%  Nonlinear-model-consistent CIP and FX time series along the empirical
%  state path (nopen_lcr spec), by interpolating the solved model's
%  state-contingent pricing functions at (sigma_t, regime_t):
%      X_t = (1 - p_t) * X_nor(sigma_t) + p_t * X_scr(sigma_t)
%  with p_t the Hamilton P(scrambling). Pure post-processing -- no solve.
%  Conditional moments are then computed EXACTLY like compute_data_moments
%  (split by p_t > 0.5), giving the apples-to-apples model counterpart.
addpath('functions'); addpath('data');

G = load('data/global_solcalibration_dynare_sigma_us.mat', ...
    'sigma_us_vec', 'Rm_us_vec', 'Rm_eu_vec', 'e_euus_vec', 'freq', 'rate_scale');
sv   = G.sigma_us_vec(:);
cipv = (G.Rm_eu_vec(:).^G.freq - G.Rm_us_vec(:).^G.freq) * G.rate_scale;  % annualized bps
ev   = G.e_euus_vec(:);
i1 = 1:76; i2 = 77:152;

RW = readtable('RW_shock_est_nopen_lcr.csv');
PR = readtable('data/MS_sigma_us_prob_est_nopen_lcr.csv');
sig = RW.sigma_us(:);                       % 295, from 2001-01
p   = [PR.prob_scr(1); PR.prob_scr(:)];     % probs start 2001-02; hold-first to align
T = numel(sig);

% Interpolate within each regime block (clamped at grid edges)
gnor = sv(i1); gscr = sv(i2);
clamp = @(x, g) min(max(x, min(g)), max(g));
n_clamp_lo = sum(sig < min(gnor)); n_clamp_hi = sum(sig > max(gnor));
cip_nor = interp1(gnor, cipv(i1), clamp(sig, gnor));
cip_scr = interp1(gscr, cipv(i2), clamp(sig, gscr));
e_nor   = interp1(gnor, ev(i1),  clamp(sig, gnor));
e_scr   = interp1(gscr, ev(i2),  clamp(sig, gscr));

cip_t = (1 - p) .* cip_nor + p .* cip_scr;
e_t   = (1 - p) .* e_nor   + p .* e_scr;
FX_t  = e_t / mean(e_t);                    % normalize like the data side

fprintf('Grid: sigma in [%.3f, %.3f]; path in [%.3f, %.3f]; clamped: %d low, %d high (of %d)\n', ...
    min(gnor), max(gnor), min(sig), max(sig), n_clamp_lo, n_clamp_hi, T);

% Conditional moments, data-side construction (p > 0.5)
scr = p > 0.5; nor = ~scr;
fprintf('\nScrambling months: %d | normal: %d\n', sum(scr), sum(nor));
fprintf('\n--- Model-consistent PATH moments (vs invariant-distribution table) ---\n');
fprintf('CIP  mean = %.1f bps | scr - nor = %.1f bps | rel std (scr/nor) = %.2f\n', ...
    mean(cip_t), mean(cip_t(scr)) - mean(cip_t(nor)), std(cip_t(scr)) / std(cip_t(nor)));
fprintf('FX   scr/nor mean gap = %.0f bps | rel std = %.2f\n', ...
    (mean(FX_t(scr)) / mean(FX_t(nor)) - 1) * 1e4, std(FX_t(scr)) / std(FX_t(nor)));

% Comparison with the DATA cip series over the same dates
D = load('data/LFX_data.mat', 'cip');
cip_d = D.cip * 12e4;                        % annualized bps
ok = isfinite(cip_d) & cip_d ~= 0;
fprintf('\n--- Path vs data CIP (%d valid months) ---\n', sum(ok));
fprintf('corr(model path, data) = %.3f | RMSE = %.1f bps\n', ...
    corr(cip_t(ok), cip_d(ok)), sqrt(mean((cip_t(ok) - cip_d(ok)).^2)));
fprintf('data: mean = %.1f | scr - nor = %.1f | rel std = %.2f\n', ...
    mean(cip_d(ok)), mean(cip_d(ok & scr)) - mean(cip_d(ok & nor)), ...
    std(cip_d(ok & scr)) / std(cip_d(ok & nor)));

out = table((1:T)', sig, p, cip_t, FX_t, 'VariableNames', {'t','sigma_us','p_scr','cip_model_bps','FX_model'});
writetable(out, 'model_path_nopen_lcr.csv');
fprintf('\nwrote model_path_nopen_lcr.csv\n');
