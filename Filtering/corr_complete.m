function r = corr_complete(x, y)
%CORR_COMPLETE Pearson correlation over pairwise-complete observations.
%
%   r = corr_complete(x, y)   ==   corr(x, y, 'rows', 'complete')
%
% Uses only base-MATLAB corrcoef, so it does not depend on the Statistics
% Toolbox signature. That matters because Dynare ships its own two-argument
% corr() at matlab/missing/stats/corr.m; when a Dynare startup file is on the
% path it shadows the toolbox version and every corr(..., 'rows', 'complete')
% call in this repo dies with a "too many input arguments" error.
%
% Rows where either input is NaN are dropped, matching 'rows','complete'.
% Returns NaN when fewer than two complete observations remain.

x = x(:);
y = y(:);

ok = ~isnan(x) & ~isnan(y);
if nnz(ok) < 2
    r = NaN;
    return
end

c = corrcoef(x(ok), y(ok));
r = c(1, 2);
end
