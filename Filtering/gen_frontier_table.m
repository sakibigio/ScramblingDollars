%% gen_frontier_table.m — LaTeX table for the feasibility-frontier grid.
% Reads data/robustness_iota.mat (written by estimate_robustness_iota.m,
% frontier redesign 2026-08-06) and writes a tabular fragment to the
% Overleaf quantfigs folder: Tab_robustness_iota_frontier.tex.
% Columns: iota (%, annualized), eta on the frontier, re-estimated lambda,
% criterion value, and the criterion relative to the baseline point.
% The baseline row (middle of the grid) is marked.

load('data/robustness_iota.mat', 'rob_iota', 'iota_grid', 'eta_grid', ...
     'iota_baseline', 'eta_baseline');
N = numel(rob_iota);
[~, i_base] = min(abs(iota_grid - iota_baseline));
fv_base = rob_iota(i_base).fval;

ovl = overleaf_dir();
fid = fopen(fullfile(ovl, 'Tab_robustness_iota_frontier.tex'), 'wt');
for ii = 1:N
    tag = '';
    if ii == i_base, tag = '$^{\dagger}$'; end
    fprintf(fid, '%.2f%s & %.3f & %.3f & %.3f & %.2f \\\\ \n', ...
        rob_iota(ii).iota * 100, tag, rob_iota(ii).eta, ...
        rob_iota(ii).lambda, rob_iota(ii).fval, ...
        rob_iota(ii).fval / fv_base);
end
fclose(fid);
fprintf('Wrote %s\n', fullfile(ovl, 'Tab_robustness_iota_frontier.tex'));

% Console copy for the log / reply letter
fprintf('\n%8s %8s %10s %11s %10s\n', 'iota(%)', 'eta', 'lambda*', 'criterion', 'rel.base');
for ii = 1:N
    if ii == i_base, mark = '  <-- baseline'; else, mark = ''; end
    fprintf('%8.2f %8.3f %10.3f %11.3f %10.2f%s\n', ...
        rob_iota(ii).iota * 100, rob_iota(ii).eta, rob_iota(ii).lambda, ...
        rob_iota(ii).fval, rob_iota(ii).fval / fv_base, mark);
end
