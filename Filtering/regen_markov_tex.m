%% Regenerate the Markov-estimates LaTeX table from MS_sigma_us_params.csv
% Stand-alone: reads the CSV and writes the Overleaf .tex file directly,
% without re-running main_filter.m.

cd(fileparts(mfilename('fullpath')));   % this file's own folder

params_file = 'data/MS_sigma_us_params.csv';
foldername = overleaf_dir();
mt_suffix = '_cd';

params_tbl = readtable(params_file);
pnames = params_tbl.param;
pvals  = params_tbl.value;
pstruct = struct();
for kk = 1:numel(pnames)
    pname = strrep(pnames{kk}, '.', '_');
    pstruct.(pname) = pvals(kk);
end

tex_file = fullfile(foldername, ['tab_markov_estimates' mt_suffix '.tex']);
fid = fopen(tex_file, 'wt');
fprintf(fid, '\\begin{tabular}{lcc}\n');
fprintf(fid, '    \\toprule\n');
fprintf(fid, '    \\multicolumn{3}{c}{\\textbf{Within Regime Processes}} \\\\\n');
fprintf(fid, '    \\midrule\n');
fprintf(fid, '    \\textbf{Coefficient} & \\textbf{Scrambling} & \\textbf{Normal} \\\\\n');
fprintf(fid, '    \\midrule\n');
fprintf(fid, '    $\\hat{\\sigma}_{ss}$ & %.3f & %.3f \\\\\n', ...
    pstruct.sigma_ss_scr, pstruct.sigma_ss_nor);
fprintf(fid, '                        & (%.3f) & (%.3f) \\\\\n', ...
    pstruct.sigma_ss_scr_se, pstruct.sigma_ss_nor_se);
fprintf(fid, '    $\\rho^{\\sigma,us}$  & %.3f & %.3f \\\\\n', ...
    pstruct.rho_scr, pstruct.rho_nor);
fprintf(fid, '                        & (%.3f) & (%.3f) \\\\\n', ...
    pstruct.rho_scr_se, pstruct.rho_nor_se);
fprintf(fid, '    $\\Sigma^{\\sigma,us}$ & %.3f & %.3f \\\\\n', ...
    pstruct.Sigma_scr, pstruct.Sigma_nor);
fprintf(fid, '                         & (%.3f) & (%.3f) \\\\\n', ...
    pstruct.Sigma_scr_se, pstruct.Sigma_nor_se);
fprintf(fid, '    \\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\vspace{6pt}\n');
fprintf(fid, '\\begin{tabular}{lcc}\n');
fprintf(fid, '    \\toprule\n');
fprintf(fid, '    \\multicolumn{2}{c}{\\textbf{Transition Matrix}} \\\\\n');
fprintf(fid, '    \\midrule\n');
fprintf(fid, '    & \\textbf{Scrambling} & \\textbf{Normal} \\\\\n');
fprintf(fid, '    \\midrule\n');
fprintf(fid, '    $\\Pr(Z_t = Z_{t-1})$ & %.1f\\%% & %.1f\\%% \\\\\n', ...
    (1 - pstruct.trans_scr) * 100, (1 - pstruct.trans_nor) * 100);
fprintf(fid, '    \\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fclose(fid);
fprintf('LaTeX Markov table written to:\n  %s\n', tex_file);
