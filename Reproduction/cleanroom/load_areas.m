function T = load_areas(cfg, areas)
% LOAD_AREAS  LOAD_AND_FILTER, for more than one visual area.
%
%   T = load_areas()                       % V1 V2 V3 hV4
%   T = load_areas(cfg, {'V1','V2','V3'})
%
% Same contract as LOAD_AND_FILTER -- read the 1.3 GB CSV once, keep the rows that
% matter, cache to _cache/areas.mat, and load from that cache ever after. The only
% change is that the area filter takes a list instead of cfg.roi, so the extrastriate
% supplement can use it. About 127k of the CSV's 2.22M rows are V1/V2/V3/hV4, so the
% cache is small enough to move between machines, which the CSV is not.
%
% RUN THIS ON A MACHINE THAT HAS THE CSV LOCALLY. Streaming 1.3 GB over an SMB mount
% to build the cache works but is slow; building it where the CSV lives and copying
% the cache is much faster.
%
% ON THE LABELLING. This reads cfg.csvPath, i.e. Support/allsubjectsTable.csv. Do NOT
% substitute Support/summaryTables_wleftV2d/allsubjectsTable.mat, which looks
% interchangeable -- same 2,222,377 x 38 shape -- but is a different parcellation:
% it carries a separate left_V2d area and has 341 fewer V1 vertices, spread over 7 of
% the 8 observers. Every V1 number in this repository, and the exact agreement of
% PRECISION_WEIGHTED_CELLS with PRECISION_WEIGHTED_TABLE, is against the CSV's
% labelling. Switching is a decision about parcellation, not about file format.
%
% Returns one table with a visual_area column; filter downstream.

    if nargin < 1 || isempty(cfg), cfg = config_repro(); end
    if nargin < 2 || isempty(areas), areas = {'V1','V2','V3','hV4'}; end

    cacheFile = fullfile(cfg.cacheDir, 'areas.mat');
    if isfile(cacheFile) && ~cfg.force_reload
        S = load(cacheFile, 'T', 'areas');
        if all(ismember(areas, S.areas))
            T = S.T(ismember(string(S.T.visual_area), string(areas)), :);
            fprintf('load_areas: cached (%d rows, %s) from %s\n', ...
                    height(T), strjoin(areas, ' '), cacheFile);
            return
        end
        fprintf('load_areas: cache lacks %s; rebuilding\n', ...
                strjoin(setdiff(areas, S.areas), ' '));
    end

    assert(isfile(cfg.csvPath), 'load_areas:csv', ...
           'CSV not found: %s\nRun this where the CSV lives, then copy _cache/areas.mat.', ...
           cfg.csvPath);

    vars = [cfg.metaCols, ...
            cfg.dg.oriCols, {cfg.dg.blank, cfg.dg.betaStd, cfg.dg.betaMean}, ...
            cfg.da.oriCols, {cfg.da.blank, cfg.da.betaStd, cfg.da.betaMean}];

    fprintf('load_areas: scanning %s for %s ...\n', cfg.csvPath, strjoin(areas, ' '));
    ds = tabularTextDatastore(cfg.csvPath, 'Delimiter', ',');
    ds.SelectedVariableNames = vars;
    ds.ReadSize = 250000;

    chunks = {};  nSeen = 0;
    while hasdata(ds)
        c = read(ds);
        nSeen = nSeen + height(c);
        m = ismember(string(c.visual_area), string(areas));
        if any(m), chunks{end+1} = c(m, :); end %#ok<AGROW>
        fprintf('  ... %d rows scanned\n', nSeen);
    end
    T = vertcat(chunks{:});

    if ~isfolder(cfg.cacheDir), mkdir(cfg.cacheDir); end
    save(cacheFile, 'T', 'areas', '-v7.3');
    fprintf('load_areas: cached %d rows to %s\n', height(T), cacheFile);
    for a = string(areas)
        fprintf('   %-5s %7d rows\n', a, nnz(string(T.visual_area) == a));
    end
end
