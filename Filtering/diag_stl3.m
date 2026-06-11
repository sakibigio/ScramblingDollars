%% Search every sheet for STLFSI / 0.3587
cd('/Users/sakibigio/Dropbox/Scrabmling for Dollars/Code 2025 ScramblingDollars/Filtering');
fn = 'data/LFX_datainputs.xlsx';
sheets = sheetnames(fn);
target_val = 0.3587;
tol = 1e-3;

for ss = 1:numel(sheets)
    sname = char(sheets(ss));
    try
        M = readmatrix(fn, 'Sheet', sname);
        [r, c] = find(abs(M - target_val) < tol);
        if ~isempty(r)
            for kk = 1:numel(r)
                % column letter
                n = c(kk); letters = '';
                while n > 0
                    rem_ = mod(n-1,26);  letters = [char('A'+rem_), letters]; %#ok<AGROW>
                    n = floor((n-1)/26);
                end
                fprintf('Found 0.3587 in sheet "%s" at row %d, column %s (#%d)\n', sname, r(kk), letters, c(kk));
            end
        end
        % Also check size
        fprintf('Sheet "%-30s": %d x %d  matrix\n', sname, size(M,1), size(M,2));
    catch ME
        fprintf('Sheet "%s": readmatrix failed (%s)\n', sname, ME.message);
    end
end
