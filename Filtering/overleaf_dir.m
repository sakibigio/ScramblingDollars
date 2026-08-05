function foldername = overleaf_dir()
%OVERLEAF_DIR Folder the figure and table scripts write into.
%
% Overleaf is the destination ONLY on a machine that already has the authors'
% Overleaf sync folder. Everyone else -- referees, replicators, coauthors
% without access -- gets a folder inside the repo. Nothing ever writes into a
% Dropbox path that doesn't exist locally.
%
% Resolution order:
%   1. SCRAMBLING_QUANTFIGS environment variable, if set. Created if absent.
%      This is the supported way to redirect output on any machine.
%   2. $HOME/Dropbox/Apps/Overleaf/<project>/quantfigs, but only if it already
%      exists -- true on the authors' machines, false everywhere else.
%   3. Filtering/output/, created on demand.
%
% Always returns a path with a trailing filesep, since callers build paths by
% string concatenation.

OVERLEAF_REL = 'Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs';

foldername = getenv('SCRAMBLING_QUANTFIGS');

if ~isempty(foldername)
    if exist(foldername, 'dir') ~= 7
        mkdir(foldername);
    end
else
    foldername = fullfile(getenv('HOME'), OVERLEAF_REL);
    if exist(foldername, 'dir') ~= 7
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
