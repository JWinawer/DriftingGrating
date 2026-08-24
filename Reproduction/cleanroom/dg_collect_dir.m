function root = dg_collect_dir()
% DG_COLLECT_DIR  Locate the extracted-data folder, wherever it is mounted.
%
%   root = dg_collect_dir()
%
% The server extraction (../server_extract/) writes ret_*.mat, glm_*.mat,
% runbetas_*.mat, labels_*.mat and gainSummary.csv into one folder. That folder lives in
% different places on different machines -- a local copy under $HOME on the analysis
% machine, an SMB mount of the same share elsewhere -- so nothing downstream should
% hard-code it. Candidates are tried in order:
%
%   1. $DG_COLLECT, if set        (explicit override, wins over everything)
%   2. ~/dg_collect               (local copy)
%   3. /Volumes/jaw288/dg_collect (SMB mount of the analysis account's home)
%
% A candidate wins only if it holds runbetas_*.mat AND ret_*.mat, because partial
% mirrors exist: ~/dg_collect on the laptop carries gainSummary.csv and the glm_*.mat
% files but none of the run-wise betas, and silently selecting it makes every caller
% that needs runs fail one file at a time. A candidate that exists but is incomplete is
% kept only as a last resort, so the failure is still reported against a real path.
%
% Returns '' if none exists, which callers should treat the way they treat a missing
% file rather than as an error -- OBSERVER_GAIN_WEIGHTS, for instance, falls back to
% unit gain weights with a warning.

    if ispc, home = getenv('USERPROFILE'); else, home = getenv('HOME'); end
    cands = {getenv('DG_COLLECT'), fullfile(home, 'dg_collect'), ...
             '/Volumes/jaw288/dg_collect'};

    root = '';
    for i = 1:numel(cands)
        if isempty(cands{i}) || ~isfolder(cands{i}), continue; end
        if ~isempty(dir(fullfile(cands{i}, 'runbetas_*.mat'))) && ...
           ~isempty(dir(fullfile(cands{i}, 'ret_*.mat')))
            root = cands{i};  return
        end
        if isempty(root), root = cands{i}; end       % incomplete; keep looking
    end
end
