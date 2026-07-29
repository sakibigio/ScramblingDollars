%% show_iota_sweep.m - collect nonlinear-model moments across iota values
d = dir(fullfile('output','nlm_moments_*.mat'));
rows = {};
for k = 1:numel(d)
    M = load(fullfile(d(k).folder, d(k).name));
    if ~isfield(M,'rows'), continue; end
    r = M.rows;
    get1 = @(nm,col) r{strcmp(r(:,1),nm), col};
    rows(end+1,:) = { M.tag, M.iota_ann, M.lambda, M.eta, ...
        get1('BP',2), get1('BP',3), get1('BP',4), ...
        get1('CIP',2), get1('CIP',4), ...
        get1('FX',3), get1('FX',4), ...
        get1('BP',6), get1('BP',7) };  %#ok<AGROW>
end
% sort by iota
[~,ix] = sort(cell2mat(rows(:,2))); rows = rows(ix,:);

fprintf('\n%-8s %8s %7s %6s | %8s %8s %8s | %8s %8s | %8s %8s | %7s %7s\n', ...
 'tag','iota bps','lambda','eta','BP mean','BP acorr','BP std','CIP mean','CIP std','FX acorr','FX std','rho_r2','rho_r1');
fprintf('%s\n', repmat('-',1,120));
for k = 1:size(rows,1)
    fprintf('%-8s %8.0f %7.3f %6.3f | %8.2f %8.3f %8.2f | %8.2f %8.2f | %8.3f %8.3f | %7.3f %7.3f\n', rows{k,:});
end
fprintf('\nrho_r1 = persistence, NORMAL regime;  rho_r2 = persistence, SCRAMBLING regime\n');
