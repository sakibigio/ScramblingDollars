function foldername = overleaf_dir()
%OVERLEAF_DIR Folder the figure and table scripts write into.
%
% Resolution order:
%   1. SCRAMBLING_QUANTFIGS environment variable, if set (use this to point a
%      fresh machine at its own Overleaf sync folder). Created if absent.
%   2. The authors' two known machines, if that folder currently exists.
%   3. Filtering/output/ -- created on demand. This is the branch a replicator
%      lands on: everything still runs, output just stays inside the repo
%      instead of writing into somebody else's Dropbox. It is also the safety
%      net when an author's Dropbox sync path has moved.
%
% Always returns a path with a trailing filesep, since callers build paths by
% string concatenation.

foldername = getenv('SCRAMBLING_QUANTFIGS');

if ~isempty(foldername)
    if exist(foldername, 'dir') ~= 7
        mkdir(foldername);
    end
else
    [~, username] = system('whoami');
    username = strtrim(username);
    switch username
        case 'sakibigio'
            foldername = '/Users/sakibigio/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs';
        case 'sakiclaudia'
            foldername = '/Users/sakiclaudia/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs';
        otherwise
            foldername = '';
    end

    if isempty(foldername) || exist(foldername, 'dir') ~= 7
        foldername = fullfile(fileparts(mfilename('fullpath')), 'output');
        if exist(foldername, 'dir') ~= 7
            mkdir(foldername);
        end
    end
end

if foldername(end) ~= filesep
    foldername = [foldername filesep];
end
end
