function T = load_allconditions(cfg)
% LOAD_ALLCONDITIONS  Like LOAD_AND_FILTER but keeps all 13 conditions per experiment.
%
% LOAD_AND_FILTER caches only the 4 stationary orientations + blank, which is all the
% figures need. This variant additionally caches the 8 motion conditions and the
% precomputed beta_mean/beta_std, so the per-vertex normalisation can be interrogated
% (e.g. building a divisor from conditions that are NOT analysed in Figs 5-8).
%
% Cache: _cache/v1_allcond.mat. Column name lists are exposed as cfg.dg.motionCols /
% cfg.da.motionCols via ALLCONDITION_COLS below.

    if nargin < 1, cfg = config_repro(); end
    cacheFile = fullfile(cfg.cacheDir, 'v1_allcond.mat');

    if isfile(cacheFile) && ~cfg.force_reload
        S = load(cacheFile, 'T');  T = S.T;
        fprintf('load_allconditions: loaded cached V1 table (%d rows)\n', height(T));
        return;
    end

    [dgMot, daMot] = allcondition_cols();
    vars = [cfg.metaCols, ...
            cfg.dg.oriCols, dgMot, {cfg.dg.blank, cfg.dg.betaStd, cfg.dg.betaMean}, ...
            cfg.da.oriCols, daMot, {cfg.da.blank, cfg.da.betaStd, cfg.da.betaMean}];

    fprintf('load_allconditions: scanning %s ...\n', cfg.csvPath);
    ds = tabularTextDatastore(cfg.csvPath, 'Delimiter', ',');
    ds.SelectedVariableNames = vars;
    ds.ReadSize = 250000;
    chunks = {}; nSeen = 0;
    while hasdata(ds)
        c = read(ds);  nSeen = nSeen + height(c);
        isV1 = strcmp(string(c.visual_area), cfg.roi);
        if any(isV1), chunks{end+1} = c(isV1, :); end %#ok<AGROW>
        fprintf('  ... %d rows scanned\n', nSeen);
    end
    T = vertcat(chunks{:});
    if ~isfolder(cfg.cacheDir), mkdir(cfg.cacheDir); end
    save(cacheFile, 'T', '-v7.3');
    fprintf('load_allconditions: cached V1 table (%d rows)\n', height(T));
end
