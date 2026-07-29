%% run_lfx_config.m
% Run the LIVE Cobb-Douglas pipeline (main_LFX.m) for one configuration and
% harvest the simulation moments.  TESTS ONLY:
%   - printit forced to 0
%   - foldername redirected to output/lfx_<cfg>/ (Overleaf never touched;
%     note compute_moments writes .tex tables unguarded by printit)
%
% Invoke:  matlab -batch "cfg='A'; run_lfx_config"
%   'A' : paper calibration + paper MS_params        (committed values)
%   'B' : tv-iota re-estimation, original mu         (iota = mean iota_t)
%   'C' : tv-iota re-estimation, new H_18 mu         (iota = mean iota_t)

if ~exist('cfg','var'), error('set cfg = ''A''|''B''|''C'''); end

% --- install this config's Markov regime parameters ---
src = fullfile('data', ['MS_params_' cfg '.csv']);
if ~isfile(src), error('%s not found', src); end
copyfile(src, fullfile('data','MS_params.csv'));
fprintf('[cfg %s] installed %s -> data/MS_params.csv\n', cfg, src);

% --- install this config's calibration override ---
if exist('_calibration_override.mat','file') && ~exist('_calibration_override_LFXBK.mat','file')
    copyfile('_calibration_override.mat','_calibration_override_LFXBK.mat');
end
switch cfg
    case 'A'
        if exist('_calibration_override_LFXBK.mat','file')
            copyfile('_calibration_override_LFXBK.mat','_calibration_override.mat');
        end
    case {'B','C'}
        if cfg=='B', L = load('../Filtering_iota_tv/_calibration_tv_origliq.mat');
        else,        L = load('../Filtering_iota_tv/_calibration_tv_newliq.mat');
        end
        lambda_us_override_val = L.lambda_opt;
        lambda_eu_override_val = L.lambda_opt;
        eta_override_val       = L.eta_opt;
        iota_ss_override_val   = mean(L.iota_t_opt);   % constant iota for the NLM
        save('_calibration_override.mat','lambda_us_override_val', ...
             'lambda_eu_override_val','eta_override_val','iota_ss_override_val');
        fprintf('[cfg %s] override: lambda=%.4f eta=%.4f iota=%.4f (%.0f bps)\n', ...
            cfg, L.lambda_opt, L.eta_opt, mean(L.iota_t_opt), mean(L.iota_t_opt)*1e4);
end

% --- run the live pipeline, output redirected locally ---
matching_type  = 1;                                    % Cobb-Douglas
lfx_printit    = 0;
lfx_foldername = fullfile(pwd, 'output', ['lfx_' cfg]);
lfx_tag        = cfg;

main_LFX

% --- harvest moments (workspace was cleared; lfx_tag survived the temp file) ---
tag = lfx_tag;
M = struct('cfg', tag);
M.lambda = lambda_us;  M.eta = eta;  M.iota_us = iota_us;
M.iota_ann_bps = iota_us*freq*1e4;
grab = {'E_e','rho_e_sim','std_e','E_e_r1','E_e_r2','rho_e_sim_r1','rho_e_sim_r2','std_e_r1','std_e_r2', ...
        'E_pi_ret_us','rho_pi_us_sim','std_pi_ret_us','E_pi_ret_us_r1','E_pi_ret_us_r2','std_pi_ret_us_r1','std_pi_ret_us_r2', ...
        'E_bp','rho_bp_sim','std_bp','E_bp_r1','E_bp_r2','rho_bp_sim_r1','rho_bp_sim_r2','std_bp_r1','std_bp_r2', ...
        'E_cip','rho_cip_sim','std_cip','E_cip_r1','E_cip_r2','rho_cip_sim_r1','rho_cip_sim_r2','std_cip_r1','std_cip_r2', ...
        'E_dw','rho_dw_sim','std_dw','E_dw_r1','E_dw_r2', ...
        'E_ff','rho_ff_sim','std_ff','E_ff_r1','E_ff_r2'};
for k = 1:numel(grab)
    if exist(grab{k},'var'), M.(grab{k}) = eval(grab{k}); else, M.(grab{k}) = NaN; end
end
if ~exist('output','dir'), mkdir('output'); end
save(fullfile('output',['lfx_metrics_' tag '.mat']), '-struct', 'M');
fprintf('\n[cfg %s DONE] BP: mean=%.2f std=%.2f rdiff=%.2f | CIP: mean=%.2f rdiff=%.2f | FX std/E=%.4f\n', ...
    tag, M.E_bp, M.std_bp, M.E_bp_r2-M.E_bp_r1, M.E_cip, M.E_cip_r2-M.E_cip_r1, M.std_e/M.E_e);
fprintf('[saved] output/lfx_metrics_%s.mat\n', tag);

% --- restore the live calibration override ---
if exist('_calibration_override_LFXBK.mat','file')
    copyfile('_calibration_override_LFXBK.mat','_calibration_override.mat');
end
