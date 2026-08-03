addpath('functions'); addpath('data');
G = load('data/global_solcalibration_dynare_sigma_us.mat', ...
    'sigma_us_vec','Rm_us_vec','Rm_eu_vec','e_euus_vec','freq','rate_scale','invp');
sv = G.sigma_us_vec(:); cipv = (G.Rm_eu_vec(:).^G.freq - G.Rm_us_vec(:).^G.freq)*G.rate_scale;
invp = G.invp(:); invp = invp/sum(invp);
P = readtable('model_path_nopen_lcr.csv');

fprintf('=== [1] Ergodic vs path reconciliation ===\n');
fprintf('ergodic CIP mean (cip_vec * invp)      = %+.1f bps\n', cipv'*invp);
fprintf('path CIP mean                          = %+.1f bps\n', mean(P.cip_model_bps));
fprintf('sigma: ergodic mean = %.3f  median = %.3f | path mean = %.3f  median = %.3f\n', ...
    sv'*invp, sv(find(cumsum(invp)>=0.5,1)), mean(P.sigma_us), median(P.sigma_us));
fprintf('invariant mass with cip_vec > 0: %.1f%% | path months with model CIP > 0: %.1f%%\n', ...
    100*sum(invp(cipv>0)), 100*mean(P.cip_model_bps>0));
pre = 1:227; post = 240:295;  % pre-2020 vs 2021+
fprintf('path CIP mean pre-2020 = %+.1f | 2021+ = %+.1f  (chain estimated pre-2020 only)\n', ...
    mean(P.cip_model_bps(pre)), mean(P.cip_model_bps(post)));
fprintf('CIP(sigma) floor: min(cip_vec) = %+.1f bps at sigma <= %.2f (normal block)\n', ...
    min(cipv(1:76)), sv(find(cipv(1:76) < min(cipv(1:76))+1, 1, 'last')));

fprintf('\n=== [2] Probability-weighted conditional moments (vs binary split) ===\n');
p = P.p_scr; x = P.cip_model_bps;
wm = @(x,w) sum(x.*w)/sum(w);
fprintf('binary  (p>0.5): scr - nor = %+.1f bps\n', mean(x(p>0.5)) - mean(x(p<=0.5)));
fprintf('weighted (by p): scr - nor = %+.1f bps\n', wm(x,p) - wm(x,1-p));

fprintf('\n=== [3] FX levels & orientation check ===\n');
D = load('data/LFX_data.mat', 'ln_eu_us_t', 'inv_e');
e_data = exp(D.ln_eu_us_t(:));
fprintf('data exp(ln_eu_us): range [%.3f, %.3f] (level ✓)\n', min(e_data), max(e_data));
fprintf('data inv_e:         range [%.3f, %.3f]\n', min(D.inv_e), max(D.inv_e));
fprintf('corr(exp(ln_eu_us), 1./inv_e) = %+.3f   (model e_euus_vec = 1./inv_e_vec)\n', ...
    corr(e_data, 1./D.inv_e(:)));
fprintf('model e_euus_vec: range [%.3f, %.3f] (level; ss approx 1)\n', min(G.e_euus_vec), max(G.e_euus_vec));

T2 = table(sv, cipv, invp, [ones(76,1); 2*ones(76,1)], 'VariableNames', {'sigma','cip_bps','invp','regime'});
writetable(T2, 'cip_profile_nopen_lcr.csv');
