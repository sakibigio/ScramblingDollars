%% Diagnostic: where does the STL FSI series actually sit in the spreadsheet?
cd(fileparts(mfilename('fullpath')));   % this file's own folder

fn = 'data/LFX_datainputs.xlsx';
sh = 'DataCounterpart';

% Sanity: read a known-good column (BV = DW_n) so we know xlsread works
fprintf('--- Sanity: known-good column BV (DW_n) ---\n');
num_bv = xlsread(fn, sh, 'BV27:BV35');
fprintf('BV27:BV35 -> '); disp(num_bv');

% Show overall sheet extent via readcell (handles empties better)
% Scan full EK column over full row range using readcell
try
    raw_ek_full = readcell(fn, 'Sheet', sh, 'Range', 'EK1:EK400');
    fprintf('readcell EK1:EK400 -> %d rows.  Non-missing/numeric count:\n', size(raw_ek_full,1));
    n_num = 0;  n_str = 0;  first_num_row = NaN;
    for ii = 1:size(raw_ek_full,1)
        v = raw_ek_full{ii,1};
        if isnumeric(v) && ~ismissing(v)
            n_num = n_num + 1;
            if isnan(first_num_row), first_num_row = ii; end
        elseif (ischar(v) || isstring(v)) && ~ismissing(v)
            n_str = n_str + 1;
            if n_str <= 3, fprintf('  String at row %d: "%s"\n', ii, char(v)); end
        end
    end
    fprintf('  -> numeric cells: %d, string cells: %d, first numeric at row %s\n', ...
        n_num, n_str, num2str(first_num_row));
catch ME
    fprintf('readcell scan failed: %s\n', ME.message);
end

% Check header row -- what columns near EK have headers?
try
    raw_hdr = readcell(fn, 'Sheet', sh, 'Range', 'EI1:EM2');
    fprintf('\nHeader area EI1:EM2 (rows 1 and 2):\n');
    for ii = 1:size(raw_hdr,1)
        for jj = 1:size(raw_hdr,2)
            v = raw_hdr{ii,jj};
            col_letter = char('A' + jj - 1 + 4*26);  % approximate
            if ~ismissing(v)
                if isnumeric(v)
                    fprintf('  (r%d, off%d) = %g\n', ii, jj, v);
                else
                    fprintf('  (r%d, off%d) = "%s"\n', ii, jj, char(string(v)));
                end
            end
        end
    end
catch ME
    fprintf('header readcell failed: %s\n', ME.message);
end

% Try a wide row range on column EK
fprintf('\n--- Column EK, rows 1:350 ---\n');
[num_ek, ~, raw_ek] = xlsread(fn, sh, 'EK1:EK350');
fprintf('num_ek size: %d x %d, non-NaN entries: %d\n', size(num_ek,1), size(num_ek,2), sum(~isnan(num_ek)));
if sum(~isnan(num_ek)) > 0
    first_num = find(~isnan(num_ek), 1, 'first');
    last_num  = find(~isnan(num_ek), 1, 'last');
    fprintf('First numeric value at relative row %d, last at %d (within EK1:EK350)\n', first_num, last_num);
    vals = num_ek(~isnan(num_ek));
    fprintf('First %d numeric values: ', min(5,numel(vals)));
    disp(vals(1:min(5,numel(vals)))');
end

% Header / context
fprintf('\nRaw EK cells (first 5 and around row 27):\n');
for ii = [1:3 25:30]
    if ii <= size(raw_ek,1)
        val = raw_ek{ii,1};
        if ischar(val) || isstring(val)
            fprintf('  row %d: "%s"\n', ii, char(val));
        elseif isnumeric(val)
            fprintf('  row %d: %g\n', ii, val);
        else
            fprintf('  row %d: [empty]\n', ii);
        end
    end
end

% Also check neighbors EJ, EL
fprintf('\n--- Quick scan of EJ and EL header rows ---\n');
[~, ~, raw_ej] = xlsread(fn, sh, 'EJ1:EJ30');
[~, ~, raw_el] = xlsread(fn, sh, 'EL1:EL30');
for ii = [1 2 25 26 27]
    if ii <= size(raw_ej,1)
        fprintf('  EJ%d: ', ii);
        v = raw_ej{ii,1};
        if ischar(v)||isstring(v), fprintf('"%s"\n', char(v));
        elseif isnumeric(v), fprintf('%g\n', v);
        else, fprintf('[empty]\n'); end
    end
    if ii <= size(raw_el,1)
        fprintf('  EL%d: ', ii);
        v = raw_el{ii,1};
        if ischar(v)||isstring(v), fprintf('"%s"\n', char(v));
        elseif isnumeric(v), fprintf('%g\n', v);
        else, fprintf('[empty]\n'); end
    end
end
