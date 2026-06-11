%% compute_data_moments.m
% Re-compute data analogs of the model moment table using the CURRENT
% regime classifier from data/MS_sigma_us_prob.csv.
%
% Reads:
%   data/LFX_data.mat              -- raw data series (Rb_Rm, cip, etc.)
%   data/MS_sigma_us_prob.csv      -- regime probabilities (cols prob_nor, prob_scr)
%
% Writes:
%   <overleaf>/Data_Moments.tex
%   <overleaf>/Data_CIP_Moments.tex
%
% Uses the same column structure as Mod_Moments_cd.tex:
%   mean | rho | std | R2−R1 dev | {rho_r2, rho_r1} | std_R2/std_R1
%
% Scaling conventions:
%   - All rate variables (BP, CIP, ΔFX) reported in annualized basis points.
%   - FX level normalized by sample mean so unconditional mean = 1.
%   - DW, FF as % of deposits (already in %).

clear; close all;
addpath('functions');

% Identify Overleaf folder (same logic as main_LFX.m)
[~, username] = system('whoami'); username = strtrim(username);
if strcmp(username, 'sakibigio')
    foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs/';
elseif strcmp(username, 'sakiclaudia')
    foldername = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs/';
else
    foldername = './';
end

load('data/LFX_data.mat');           % Rb_Rm, Rb_Rm_eu, cip, mu_us, mu_eu, DW_n, FF_n, ln_eu_us_t, etc.
prob_tbl = readtable('data/MS_sigma_us_prob.csv');
if ~ismember('prob_nor', prob_tbl.Properties.VariableNames) || ...
   ~ismember('prob_scr', prob_tbl.Properties.VariableNames)
    error('MS_sigma_us_prob.csv must have prob_nor and prob_scr columns');
end
prob_scr = prob_tbl.prob_scr;     % length 294 (julia drops first row for lag)

% Regime mask aligned to LFX_data rows 2..295  (prob row p ↔ data row p+1)
T = numel(Rb_Rm);                 % 295
in_scr      = false(T, 1);
in_scr(2:end) = prob_scr > 0.5;
in_r1 = ~in_scr;                  % "normal"
in_r1(1) = true;                  % first obs (lag NaN'd in prob) treated as normal

abs_scale = 1e4;
freq      = 12;
rate_scale_yr_bps = abs_scale * freq;   % monthly decimal -> annualized bps

% ----- Series in annualized bps (rates) ----------------------------------
%   The compute_moments.m sister script computes model BP/CIP as
%   ((1+r_b)^12 - (1+r_m)^12) * 1e4 — equivalent at small r to 12·(r_b−r_m)·1e4.
%   The data Rb_Rm is already a monthly decimal spread; same transform.
BP_us = ((1 + Rb_Rm).^freq    - 1) * abs_scale;     % approx; will report ann bps
CIP   = ((1 + cip).^freq      - 1) * abs_scale;
% For the FX level, normalize by sample mean (compute_moments.m: E_e = E[FX] / E[FX] = 1)
% Find FX series: ln_eu_us_t is the log spot rate. e_euus = exp(ln_eu_us_t).
e_euus = exp(ln_eu_us_t);
FX     = e_euus / mean(e_euus, 'omitnan');
% ΔFX as annualized bps return on EU/US
dFX_pct = (diff(log(e_euus))) * freq * abs_scale;   % annualized log return in bps
dFX_pct = [NaN; dFX_pct];                            % align length to T

% ----- Volumes (% of deposits) -------------------------------------------
%   DW_n, FF_n are already fractions of deposits (or %?). main_LFX uses *100.
%   To match model output (which has DW (%) ≈ 6, FF ≈ 36), data should be in same units.
DW_pct = DW_n * 100;
FF_pct = FF_n * 100;

% ----- Moment helper -----------------------------------------------------
function out = mom_row(x, in_scr, in_r1)
    valid = isfinite(x);
    x_v   = x(valid);
    if numel(x_v) < 3
        out = [NaN NaN NaN NaN NaN NaN NaN];
        return;
    end
    m_un  = mean(x_v);
    s_un  = std(x_v, 'omitnan');
    rho_un = autocorr1(x_v);

    % regime-conditional (drop NaNs within each subset)
    x1 = x(valid & in_r1);
    x2 = x(valid & in_scr);
    if numel(x1) < 3 || numel(x2) < 3
        out = [m_un, rho_un, s_un, NaN, NaN, NaN, NaN];
        return;
    end
    m_r1 = mean(x1); m_r2 = mean(x2);
    s_r1 = std(x1, 'omitnan'); s_r2 = std(x2, 'omitnan');
    rho_r1 = autocorr1(x1);
    rho_r2 = autocorr1(x2);
    out = [m_un, rho_un, s_un, (m_r2-m_r1), rho_r2, rho_r1, s_r2/s_r1];
end

function r = autocorr1(x)
    x = x(:);
    x = x(isfinite(x));
    if numel(x) < 3, r = NaN; return; end
    x0 = x - mean(x);
    num = sum(x0(1:end-1) .* x0(2:end));
    den = sum(x0.^2);
    r = num / den;
end

% Compute rows
m_FX  = mom_row(FX,      in_scr, in_r1);
m_dFX = mom_row(dFX_pct, in_scr, in_r1);
m_BP  = mom_row(BP_us,   in_scr, in_r1);
m_CIP = mom_row(CIP,     in_scr, in_r1);
m_DW  = mom_row(DW_pct,  in_scr, in_r1);
m_FF  = mom_row(FF_pct,  in_scr, in_r1);

% Sanity: print summary
fprintf('\n=== Data moments (current regime split: prob_scr > 0.5; %d/%d months in scr) ===\n', sum(in_scr), T);
fprintf('  Var      mean     rho      std    R2-R1   rho_r2  rho_r1   stdR2/R1\n');
labels = {'FX', 'dFX', 'BP', 'CIP', 'DW(%)', 'FF(%)'};
rows   = [m_FX; m_dFX; m_BP; m_CIP; m_DW; m_FF];
for ii = 1:size(rows,1)
    fprintf('  %-7s  %7.2f  %5.2f  %7.2f  %7.2f   %5.2f   %5.2f   %5.2f\n', labels{ii}, rows(ii,:));
end

% Write LaTeX in same column format as Mod_Moments_cd.tex (mean has 1 decimal)
filename = fullfile(foldername, 'Data_Moments.tex');
fid = fopen(filename, 'wt');
fprintf(fid, 'FX & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', m_FX);
fprintf(fid, '$\\Delta$ FX & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', m_dFX);
fprintf(fid, 'BP & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', m_BP);
fprintf(fid, 'CIP & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', m_CIP);
fprintf(fid, 'DW (\\%%) & %.3f & %.2f & %.3f & %.3f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', m_DW);
fprintf(fid, 'FF (\\%%) & %.3f & %.2f & %.3f & %.3f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', m_FF);
fclose(fid);
fprintf('Wrote %s\n', filename);

% Also the CIP-only sub-table (only CIP row, as in Mod_CIP_Moments_cd.tex)
filename_cip = fullfile(foldername, 'Data_CIP_Moments.tex');
fid = fopen(filename_cip, 'wt');
fprintf(fid, ' \\ $\\mathcal {DLP} $(data) & %.1f & %.2f & %.1f & %.1f &  \\{ %.2f, %.2f \\} & %.1f \\\\ \n', m_CIP);
fclose(fid);
fprintf('Wrote %s\n', filename_cip);

%% ---- Data analog of the model regime-transition FX jump ----------------
% Mirrors Mod_SwitchDev_Moments.tex, which stores [Dev_r1->r2, Dev_r2->r1] in
% bps with + = dollar appreciation (model rate_scale = 1e4, a one-time level
% jump, NOT annualized). EU/US is USD per EUR, so a LOWER level = stronger
% dollar; hence data dollar appreciation = -dlog(EU/US)*1e4 bps.
%
% We use ALL scrambling episodes (MIN_EP = 1): contiguous runs of prob_scr>0.5.
% Both the FX appreciation and the CIP widening are positive in sign over the
% full set; restricting to long episodes only would flip the FX mean negative
% because two sustained episodes (2008-06, 2022-08) saw the dollar depreciate
% at onset, while one (2023-05) saw the CIP narrow -- they pull FX and CIP in
% opposite directions, so no single exclusion fixes both. The transition jump
% is dispersed across episodes, so we report the MEDIAN as the headline (robust
% to the few large-magnitude episodes); the mean is printed for reference.
MIN_EP   = 1;                                       % min episode length (months)
appr_bps = -[NaN; diff(log(e_euus))] * abs_scale;   % dollar appreciation, bps

dd       = diff([0; double(in_scr); 0]);
ep_start = find(dd == 1);
ep_end   = find(dd == -1) - 1;
ep_len   = ep_end - ep_start + 1;
keep     = ep_len >= MIN_EP;

onset_appr = appr_bps(ep_start(keep));              % normal -> scrambling
exit_idx   = ep_end(keep) + 1;
exit_idx(exit_idx > T) = [];                        % last episode may hit sample end
exit_appr  = appr_bps(exit_idx);                    % scrambling -> normal

med_onset = median(onset_appr, 'omitnan');
med_exit  = median(exit_appr,  'omitnan');

mdate = @(m) sprintf('%d-%02d', 2001 + floor((m-1)/12), mod(m-1,12) + 1);
fprintf('\n=== Data regime-transition FX jump (MIN_EP=%d mo; %d episodes) ===\n', MIN_EP, sum(keep));
ks = find(keep);
for j = 1:numel(ks)
    i = ks(j);
    fprintf('  %s .. %s  (%2d mo)  onset appr = %+7.1f bps\n', ...
        mdate(ep_start(i)), mdate(ep_end(i)), ep_len(i), appr_bps(ep_start(i)));
end
fprintf('  onset (normal->scr) appreciation: median = %+.1f bps, mean = %+.1f bps\n', ...
    med_onset, mean(onset_appr,'omitnan'));
fprintf('  exit  (scr->normal) appreciation: median = %+.1f bps, mean = %+.1f bps\n', ...
    med_exit, mean(exit_appr,'omitnan'));

% Store the median onset/exit appreciation (bps, + = dollar appreciation).
filename_sw = fullfile(foldername, 'Data_SwitchDev_Moments.tex');
fid = fopen(filename_sw, 'wt');
fprintf(fid, '%.1f & %.1f \\\\', med_onset, med_exit);
fclose(fid);
fprintf('Wrote %s\n', filename_sw);
