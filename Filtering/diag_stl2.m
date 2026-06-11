%% Deeper diagnostic
cd('/Users/sakibigio/Dropbox/Scrabmling for Dollars/Code 2025 ScramblingDollars/Filtering');

fn = 'data/LFX_datainputs.xlsx';
fprintf('Sheet names:\n');
sheets = sheetnames(fn);
for ii = 1:numel(sheets), fprintf('  %2d: "%s"\n', ii, sheets(ii)); end

% Try readmatrix on column EK with different options
fprintf('\n--- readmatrix EK1:EK350 ---\n');
try
    M = readmatrix(fn, 'Sheet', 'DataCounterpart', 'Range', 'EK1:EK350');
    fprintf('readmatrix size: %d x %d, finite: %d, first finite at row %s\n', ...
        size(M,1), size(M,2), sum(isfinite(M)), ...
        num2str(find(isfinite(M), 1, 'first')));
    if sum(isfinite(M)) > 0
        f = M(isfinite(M));
        fprintf('first 5 finite: '); disp(f(1:min(5,end))');
    end
catch ME
    fprintf('readmatrix failed: %s\n', ME.message);
end

% Try detectImportOptions to see how MATLAB perceives the columns
fprintf('\n--- detectImportOptions ---\n');
try
    opts = detectImportOptions(fn, 'Sheet', 'DataCounterpart');
    nvars = numel(opts.VariableNames);
    fprintf('Detected %d variables. Last 5 names:\n', nvars);
    for ii = max(1, nvars-5):nvars
        fprintf('  var %d: %s\n', ii, opts.VariableNames{ii});
    end
catch ME
    fprintf('detectImportOptions failed: %s\n', ME.message);
end

% Search column header containing "STLFSI"
fprintf('\n--- Searching row 1 for "STLFSI" header ---\n');
try
    hdr = readcell(fn, 'Sheet', 'DataCounterpart', 'Range', 'A1:GZ1');
    for jj = 1:size(hdr, 2)
        v = hdr{1, jj};
        if (ischar(v) || isstring(v)) && contains(string(v), 'STL', 'IgnoreCase', true)
            % Build column letter
            n = jj;
            letters = '';
            while n > 0
                r = mod(n-1, 26);
                letters = [char('A' + r), letters]; %#ok<AGROW>
                n = floor((n-1) / 26);
            end
            fprintf('  Found "STL" header at column %s (#%d): "%s"\n', letters, jj, char(string(v)));
        end
    end
catch ME
    fprintf('header search failed: %s\n', ME.message);
end
