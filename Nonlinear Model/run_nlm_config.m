%% run_nlm_config.m
% Run the nonlinear model (main_LFX.m) under a given calibration, capture the
% figures locally, and rebuild the model-moments table from the workspace.
%
% Invoke e.g.:
%   matlab -batch "nlm_tag='C'; nlm_lambda=1.129; nlm_eta=0.662; nlm_sigma_us=0.314; nlm_iota_ann=0.0855; run_nlm_config"
%
% Any of nlm_lambda / nlm_eta / nlm_sigma_us / nlm_sigma_eu / nlm_iota_ann left
% unset (or NaN) keeps the value from LFX_params_cd_v3.m.
% printit is forced to 0 so the Overleaf folder is never written.

if ~exist('nlm_tag','var'), nlm_tag = 'run'; end
nlm_printit = 0;                    % never touch Overleaf
tag = nlm_tag;                      % local copy (main_LFX clears the workspace)

main_LFX

%% ---- did the script complete? ----
tag  = nlm_tag;                     % restored via the temp-file mechanism
figs = findobj('Type','figure');
fprintf('\n===== NLM run [%s] =====\n', tag);
fprintf('figures created: %d\n', numel(figs));
key = {'E_bp','std_bp','rho_bp_sim','E_cip','std_cip','rho_cip_sim', ...
       'rho_e_sim','std_e','E_e','E_pi_ret_us','std_pi_ret_us'};
present = false(size(key));
for kk = 1:numel(key), present(kk) = evalin('caller', sprintf('exist(''%s'',''var'')==1', key{kk})); end
if all(present)
    fprintf('moment variables: ALL present (script reached the moments section)\n');
else
    fprintf('moment variables MISSING: %s\n', strjoin(key(~present),', '));
end

%% ---- capture figures locally ----
outdir = fullfile('output', ['nlm_' tag]);
if ~exist(outdir,'dir'), mkdir(outdir); end
[~, ord] = sort([figs.Number]);  figs = figs(ord);
nsaved = 0;
for k = 1:numel(figs)
    f = figs(k);
    nm = get(f,'Name');
    if isempty(nm), nm = sprintf('fig%d', f.Number); end
    nm = regexprep(strtrim(nm), '[^\w]+', '_');
    try
        exportgraphics(f, fullfile(outdir, sprintf('%02d_%s.png', f.Number, nm)), 'Resolution',150);
        nsaved = nsaved + 1;
    catch
    end
end
fprintf('figures captured: %d -> %s\n', nsaved, outdir);

%% ---- rebuild the model-moments table from the workspace ----
% (main_LFX writes Mod_Moments.tex to Overleaf but it comes out empty; we
%  reconstruct the same rows here so the numbers are actually visible.)
M = struct('tag', tag);
try
    rows = {
      'FX',        1,              rho_e_sim,       std_e/E_e,        (E_e_r2./E_e_r1-1)*1e4,           rho_e_sim_r2,   rho_e_sim_r1,   std_e_r2./std_e_r1
      'dFX',       E_pi_ret_us,    rho_pi_us_sim,   std_pi_ret_us,    E_pi_ret_us_r2-E_pi_ret_us_r1,    rho_e_sim_r2,   rho_e_sim_r1,   std_pi_ret_us_r2./std_pi_ret_us_r1
      'BP',        E_bp,           rho_bp_sim,      std_bp,           E_bp_r2-E_bp_r1,                  rho_bp_sim_r2,  rho_bp_sim_r1,  std_bp_r2./std_bp_r1
      'CIP',       E_cip,          rho_cip_sim,     std_cip,          E_cip_r2-E_cip_r1,                rho_cip_sim_r2, rho_cip_sim_r1, std_cip_r2./std_cip_r1 };
    fprintf('\n%-6s %10s %10s %10s %12s %10s %10s %10s\n', ...
        'var','mean','autocorr','std','regime diff','rho_r2','rho_r1','std ratio');
    fprintf('%s\n', repmat('-',1,86));
    for ii = 1:size(rows,1)
        fprintf('%-6s %10.3f %10.3f %10.3f %12.2f %10.3f %10.3f %10.3f\n', rows{ii,:});
    end
    M.rows = rows;
catch ME
    fprintf('[moments table unavailable: %s]\n', ME.message);
end

% scalar steady-state moments that main_LFX prints
try
    M.E_bp_ss  = (Rb_us_vec.^(freq)-Rm_us_vec.^(freq))*f_ss(:)*1e4;
    M.E_rd_ss  = (Rm_eu_vec.^(freq)-Rm_us_vec.^(freq))*f_ss(:)*1e4;
    fprintf('\nE[Rb_us-Rm_us] = %.2f bps   E[Rm_eu-Rm_us] = %.2f bps\n', M.E_bp_ss, M.E_rd_ss);
catch
end
M.lambda = lambda_us; M.eta = eta; M.sigma_us = sigma_us;
M.iota_us = iota_us;  M.iota_ann = iota_us*freq*1e4;

if ~exist('output','dir'), mkdir('output'); end
save(fullfile('output',['nlm_moments_' tag '.mat']), '-struct', 'M');
fprintf('[saved] output/nlm_moments_%s.mat\n', tag);
