function ok = validate_table_rows(subject, csvPath, opts)
% VALIDATE_TABLE_ROWS  Check BUILD_TABLE_ROWS against the shipped table.
%
%   ok = validate_table_rows('sub-0037')
%
% Rebuilds an observer who is ALREADY in allsubjectsTable.csv and compares every
% column against the shipped rows. Run this before trusting BUILD_TABLE_ROWS on a new
% observer: the table is the input to every clean-room result, and a rebuild that is
% subtly wrong would be appended silently.
%
% TWO FAILURE MODES THIS EXISTS TO CATCH, both of which passed a count-based check:
%
%   * Off-by-one in the vertex index. COLLECT_LABELS stores label indices already
%     converted to 1-based, while createTables.m adds the 1 itself. Adding it twice
%     shifts every vertex by one. Every label keeps its exact SIZE, so counts agree
%     and only boundary vertices move into the neighbouring area.
%   * NaN columns compared with max(), which IGNORES NaN. A column that is entirely
%     NaN on one side and real on the other reports max|diff| = 0, i.e. a perfect
%     match, when nothing matched at all. NaN patterns are compared explicitly here.

    if nargin < 2 || isempty(csvPath)
        thisDir = fileparts(mfilename('fullpath'));
        csvPath = fullfile(fileparts(fileparts(thisDir)), 'Support', 'allsubjectsTable.csv');
    end
    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'projects'), opts.projects = {'dg'}; end
    if ~isfield(opts,'tol'),      opts.tol = 1e-9; end

    T = build_table_rows(subject, opts);

    ds = tabularTextDatastore(csvPath, 'Delimiter', ',');
    ds.ReadSize = 100000;
    acc = {}; seen = false;
    while hasdata(ds)
        c = read(ds);
        k = string(c.subject) == string(subject);
        if any(k)
            acc{end+1} = c(k,:); seen = true; %#ok<AGROW>
        elseif seen
            break    % chunk with no matches AFTER we started: we are past this subject
        end
    end
    if isempty(acc)
        error('validate_table_rows:absent', '%s is not in %s.', subject, csvPath);
    end
    R = vertcat(acc{:});

    fprintf('validate_table_rows(%s)\n  shipped %d rows, rebuilt %d rows\n', ...
            subject, height(R), height(T));
    ok = height(R) == height(T);
    if ~ok, fprintf('  ROW COUNT DIFFERS -- stopping.\n'); return; end

    vn = T.Properties.VariableNames;
    if ~isequal(vn, R.Properties.VariableNames)
        fprintf('  COLUMN NAMES/ORDER DIFFER -- stopping.\n'); ok = false; return
    end

    worst = 0; bad = {}; skipped = 0;
    for i = 1:numel(vn)
        a = R.(vn{i}); b = T.(vn{i});
        if ismember(vn{i}, {'subject','visual_area'})
            as = string(a); bs = string(b);
            as(ismissing(as)) = ""; bs(ismissing(bs)) = "";
            if ~all(as == bs)
                bad{end+1} = sprintf('%s: %d of %d differ', vn{i}, sum(as~=bs), numel(as)); %#ok<AGROW>
            end
        elseif strcmp(vn{i}, 'included')
            if ~isequal(logical(double(a)), logical(b))
                bad{end+1} = sprintf('%s: logical mismatch', vn{i}); %#ok<AGROW>
            end
        elseif ~isColumnForProjects(vn{i}, opts.projects)
            skipped = skipped + 1;    % not rebuilt in this call, nothing to compare
        else
            a = double(a); b = double(b);
            nanmis = nnz(xor(isnan(a), isnan(b)));
            d = abs(a - b); d = d(~isnan(d));
            mx = max([d; 0]);
            worst = max(worst, mx);
            if nanmis > 0 || mx >= opts.tol
                bad{end+1} = sprintf('%s: NaN-mismatch %d, max|diff| %.3g', vn{i}, nanmis, mx); %#ok<AGROW>
            end
        end
    end

    fprintf('  columns compared %d, skipped %d (not in opts.projects)\n', numel(vn)-skipped, skipped);
    fprintf('  worst numeric difference: %.3g\n', worst);
    ok = isempty(bad);
    if ok
        fprintf('  RESULT: reproduces the shipped table.\n');
    else
        fprintf('  RESULT: DIFFERS --\n'); fprintf('    %s\n', bad{:});
    end
end

function tf = isColumnForProjects(name, projects)
% Beta columns belong to one experiment; metadata belongs to both.
    if startsWith(name, 'cartexp') || startsWith(name, 'dg_')
        tf = ismember('dg', projects);
    elseif startsWith(name, 'polexp') || startsWith(name, 'da_')
        tf = ismember('da', projects);
    else
        tf = true;
    end
end
