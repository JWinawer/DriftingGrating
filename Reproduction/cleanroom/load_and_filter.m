function T = load_and_filter(cfg)
% LOAD_AND_FILTER  Read the 1.3 GB CSV once, keep only V1 rows, cache to _cache/v1.mat.
%
% The CSV has ~2.22M rows but only ~49k are V1; caching the V1 subset makes every
% later iteration fast. Reads in chunks via a tabularTextDatastore so memory stays
% bounded. Returns a table with the metadata + per-experiment beta columns named in
% CONFIG_REPRO.

    if nargin < 1, cfg = config_repro(); end
    cacheFile = fullfile(cfg.cacheDir, 'v1.mat');

    if isfile(cacheFile) && ~cfg.force_reload
        S = load(cacheFile, 'T');
        T = S.T;
        fprintf('load_and_filter: loaded cached V1 table (%d rows) from %s\n', height(T), cacheFile);
        return;
    end

    vars = [cfg.metaCols, ...
            cfg.dg.oriCols, {cfg.dg.blank, cfg.dg.betaStd, cfg.dg.betaMean}, ...
            cfg.da.oriCols, {cfg.da.blank, cfg.da.betaStd, cfg.da.betaMean}];

    fprintf('load_and_filter: scanning %s ...\n', cfg.csvPath);
    ds = tabularTextDatastore(cfg.csvPath, 'Delimiter', ',');
    ds.SelectedVariableNames = vars;
    ds.ReadSize = 250000;   % rows per chunk

    chunks = {};
    nSeen = 0;
    while hasdata(ds)
        c = read(ds);
        nSeen = nSeen + height(c);
        isV1 = strcmp(string(c.visual_area), cfg.roi);
        if any(isV1)
            chunks{end+1} = c(isV1, :); %#ok<AGROW>
        end
        fprintf('  ... %d rows scanned\n', nSeen);
    end
    T = vertcat(chunks{:});

    if ~isfolder(cfg.cacheDir), mkdir(cfg.cacheDir); end
    save(cacheFile, 'T', '-v7.3');
    fprintf('load_and_filter: cached V1 table (%d rows) to %s\n', height(T), cacheFile);
end
