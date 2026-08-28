function assert_same_observers(subsA, subsB, contextA, contextB)
% ASSERT_SAME_OBSERVERS  Refuse to combine two results measured on different people.
%
%   assert_same_observers(subsA, subsB, 'dg', 'da')
%
% Call this before ANY arithmetic that pairs observers across two results -- a
% difference, a paired test, a correlation, a shared bootstrap draw. It errors
% unless the two lists are identical AND in the same order, because every array in
% this pipeline carries observers in the subject dimension positionally: element i
% of a dg array and element i of a da array are only the same person when the two
% lists match exactly.
%
% This exists because the failure it prevents is SILENT. Before the subject lists
% were split, dg and da always shared one cfg.subjects, so pairing by position was
% correct by construction. With dg at 13 observers and da at 7 the same code returns
% a full set of plausible numbers computed from mismatched people, with no error and
% nothing obviously wrong in the output. See CONFIG_REPRO and SUBJECTS_FOR.

    if nargin < 3 || isempty(contextA), contextA = 'first set';  end
    if nargin < 4 || isempty(contextB), contextB = 'second set'; end

    a = cellstr(subsA(:));
    b = cellstr(subsB(:));

    if numel(a) == numel(b) && all(strcmp(a, b))
        return
    end

    onlyA = setdiff(a, b, 'stable');
    onlyB = setdiff(b, a, 'stable');
    msg = sprintf(['cannot pair observers across ''%s'' (n=%d) and ''%s'' (n=%d) -- ' ...
                   'these results are not measured on the same people, in the same order.\n'], ...
                  contextA, numel(a), contextB, numel(b));
    if ~isempty(onlyA)
        msg = [msg sprintf('  only in %s : %s\n', contextA, strjoin(onlyA', ', '))];
    end
    if ~isempty(onlyB)
        msg = [msg sprintf('  only in %s : %s\n', contextB, strjoin(onlyB', ', '))];
    end
    if isempty(onlyA) && isempty(onlyB)
        msg = [msg '  same observers, DIFFERENT ORDER -- positional pairing would be wrong.\n'];
    end
    msg = [msg 'Use subjects_for(cfg, ''matched'') for anything that crosses experiments.'];
    error('assert_same_observers:mismatch', '%s', msg);
end
