function S = collect_gain_areas(varargin)
% COLLECT_GAIN_AREAS  pRF gain for the eight visual areas, not V1 alone.
%
%   S = collect_gain_areas()
%   S = collect_gain_areas('subjects', {'sub-0037'})
%   S = collect_gain_areas('validate', true)   % V1 check against gainSummary.csv
%
% The existing gainSummary.csv holds ONE gain per observer, computed over V1 only
% (DG_GAININV1: V1_REmanual, 4-8 deg, R2 > 0.1). Applying a V1-derived scalar to V2, V3
% or hV4 special-cases V1 as the source, which is what EXTRASTRIATE.md's decision 1 set
% out to remove. This computes gain over the union of V1 V2 V3 V3a V3b hV4 MT MST so a
% per-observer x map factor is available.
%
% WHAT GAIN IS. Not the fitted beta. vistasoft's rfGaussian2d makes Gaussians of unit
% HEIGHT, so beta shrinks as roughly 1/sigma^2 and is not comparable across vertices of
% different pRF size. Following ../../AnalysisCode/01_calculate_observer_gain/README.md,
% gain is the largest excursion of the model's own predicted time series away from
% baseline -- already percent BOLD, already HRF-convolved:
%
%       gain = max_t | beta(1) * (allstimimages * RF)(t) |
%
% RMMODELGAIN does the work; this only chooses which vertices to ask for. The V1
% restriction in DG_COMPUTEGAIN was for speed, not design -- 'voxels' takes any index
% set -- so widening it needs no new transfer beyond the same results.mat.
%
% COST. results.mat is ~19.6 MB per subject per protocol and both are needed (the
% analysis uses the mov/stat average), so ~313 MB over 16 files, about 23 s each on the
% Abu Dhabi mount. datafiles.mat (659 MB per subject) is NOT needed: the predicted time
% series is already in percent BOLD because the data were converted before fitting.
%
% NO ECCENTRICITY OR R2 FILTER IS APPLIED HERE, and that is deliberate. The selection
% criteria for the extrastriate maps are still open -- 4-8 deg leaves hV4 with 48 empty
% cells of 128, so a relaxed band may be needed -- and gain must be computed over the
% SAME vertices the analysis uses, or the normalisation does not correspond to the data
% being normalised. Storing gain per vertex decouples the two: any band can be
% summarised afterwards, and the saved vector can be indexed by whatever vertex set a
% future analysis settles on, without going back to the server. The summary table
% reports 4-8, 2-10, 1-12 and whole-map for convenience only. Vertices that were unfit
% get 0 from RMMODELGAIN, which is why the summary applies R2 and eccentricity itself.
%
% Writes gain_areas_<subject>.mat into the collect dir, holding vertIndex, areaMask,
% areaNames, gainMov, gainStat and gainAvg (their per-vertex mean), plus a per-area
% summary table.
%
% THE MOV/STAT COMBINATION IS A GEOMETRIC MEAN, not the arithmetic one gainSummary.csv
% used (JW/Rania, superseding it). Gain is a multiplicative scale factor and its
% distribution across vertices is right-skewed, so the geometric mean is the meaningful
% average of two protocol estimates of the same underlying responsiveness; it is also
% what OBSERVER_GAIN_WEIGHTS already uses to form the group summary. The per-protocol
% meanGain_mov and meanGain_stat columns are unchanged and DO reproduce
% gainSummary.csv exactly, which is what the validation checks -- only their
% combination differs
% returned as S and written to gain_areas_summary.csv.
%
% Requires /Volumes/Vision mounted, and the handoff folder on the path.

    AREA_LABELS = { 'V1','V1'; 'V2','V2'; 'V3','V3'; 'V3a','V3a'; ...
                    'V3b','V3b'; 'hV4','hV4'; 'MT','pMT'; 'MST','pMST' };

    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(fileparts(thisDir)), ...
                     'AnalysisCode', '01_calculate_observer_gain'));
    addpath(fullfile(fileparts(thisDir), 'cleanroom'));

    p = inputParser;
    p.addParameter('movRoot',  ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
                                'data_bids/derivatives/prfvista_mov'], @ischar);
    p.addParameter('statRoot', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
                                'data_bids/derivatives/prfvista'], @ischar);
    p.addParameter('root',     dg_collect_dir(), @ischar);
    p.addParameter('subjects', {}, @iscell);
    p.addParameter('force',    false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('validate', true,  @(x) islogical(x) || isnumeric(x));
    p.addParameter('metric',   'maxabs', @ischar);
    p.parse(varargin{:});
    opt = p.Results;

    assert(~isempty(which('rmModelGain')), ...
           'rmModelGain.m not on the path -- expected in AnalysisCode/01_calculate_observer_gain.');
    assert(isfolder(opt.movRoot),  'not mounted: %s', opt.movRoot);
    assert(isfolder(opt.statRoot), 'not mounted: %s', opt.statRoot);

    cfg  = config_repro();
    subs = opt.subjects;  if isempty(subs), subs = cfg.subjects; end
    areaNames = AREA_LABELS(:,1).';
    nA = numel(areaNames);

    rows = {};
    for ii = 1:numel(subs)
        s = subs{ii};
        out = fullfile(opt.root, sprintf('gain_areas_%s.mat', s));
        t0 = tic;

        Rt = load(fullfile(opt.root, sprintf('ret_%s.mat', s)), 'eccen','vexpl','nLH');
        Lb = load(fullfile(opt.root, sprintf('labels_%s.mat', s)));
        nLH = double(Rt.nLH);  nTot = numel(Rt.eccen);

        sets = cell(1, nA);
        for a = 1:nA
            fl = ['lh_' AREA_LABELS{a,2} '_REmanual'];
            fr = ['rh_' AREA_LABELS{a,2} '_REmanual'];
            sets{a} = [double(Lb.(fl)(:)); double(Lb.(fr)(:)) + nLH];
        end
        v = unique(vertcat(sets{:}));
        v = v(v >= 1 & v <= nTot);
        m = false(numel(v), nA);
        for a = 1:nA, m(:,a) = ismember(v, sets{a}); end

        if isfile(out) && ~opt.force
            fprintf('  %-14s skipped (exists)\n', s);
            G = load(out);
        else
            gMov  = gain_for(opt.movRoot,  s, v, nTot, opt.metric);
            gStat = gain_for(opt.statRoot, s, v, nTot, opt.metric);
            G = struct('subject', s, 'vertIndex', v, 'areaMask', m, ...
                       'areaNames', {areaNames}, 'nLH', nLH, ...
                       'gainMov', gMov, 'gainStat', gStat, ...
                       'gainAvg', sqrt(max(gMov,0) .* max(gStat,0)), ...
                       'note', ['gain over the union of the eight areas, percent BOLD, ' ...
                                'from rmModelGain. No ecc/R2 filter applied; apply one ' ...
                                'when summarising. gainAvg is the geometric mean of the ' ...
                                'two protocols.'], ...
                       'collected', datetime('now'));
            save(out, '-struct', 'G', '-v7.3');
            fprintf('  %-14s %5d vertices, %.0f s -> %s\n', s, numel(v), toc(t0), out);
        end

        ecc = double(Rt.eccen(G.vertIndex));  ve = double(Rt.vexpl(G.vertIndex));
        for a = 1:nA
            for band = {[4 8], [2 10], [1 12], [0 inf]}
                b = band{1};
                % No >0 filter, matching DG_COMPUTEGAIN: it averaged over exactly the
                % ROI-and-criteria vertex set, unfit zeros included. And meanGain_avg is
                % the ARITHMETIC mean of the two protocol summaries, as in
                % gainSummary.csv -- (5.4672 + 3.5887)/2 = 4.5279 for sub-0037 V1 --
                % not a per-vertex combination.
                % k masks vertIndex; gainMov/gainStat are FULL-SURFACE vectors, so they
                % must be indexed by the vertex numbers, not by the mask.
                k  = G.areaMask(:,a) & ecc >= b(1) & ecc <= b(2) & ve > cfg.r2min;
                vk = G.vertIndex(k);
                gm = mean(G.gainMov(vk));  gs = mean(G.gainStat(vk));
                rows(end+1,:) = {s, areaNames{a}, sprintf('%g-%g', b(1), b(2)), numel(vk), ...
                                 gm, gs, sqrt(gm*gs), ...
                                 sqrt(median(G.gainMov(vk)) * median(G.gainStat(vk)))};  %#ok<AGROW>
            end
        end
    end

    S = cell2table(rows, 'VariableNames', {'subject','area','eccBand','nVertices', ...
                                           'meanGain_mov','meanGain_stat', ...
                                           'meanGain_avg','medianGain_avg'});
    csv = fullfile(opt.root, 'gain_areas_summary.csv');
    writetable(S, csv);
    fprintf('\nwrote %s (%d rows)\n', csv, height(S));

    if opt.validate, validate_v1(S, opt.root); end
end

% ------------------------------------------------------------------------
function g = gain_for(root, sub, voxels, nTot, metric)
% One protocol, one subject. Returns a full-surface vector, 0 outside `voxels`.
    hits = dir(fullfile(root, sub, 'ses-*', '*results.mat'));
    assert(~isempty(hits), 'no *results.mat under %s', fullfile(root, sub));
    f = fullfile(hits(1).folder, hits(1).name);
    % model/params are nested inside a 'results' struct in these files, so unwrap
    % before handing to rmModelGain (DG_COMPUTEGAIN's findModelParams does the same).
    V = load(f);
    if isfield(V,'model') && isfield(V,'params')
        rm = struct('model', V.model, 'params', V.params);
    else
        rm = [];
        for fn = fieldnames(V).'
            n = V.(fn{1});
            if isstruct(n) && isscalar(n) && isfield(n,'model') && isfield(n,'params')
                rm = struct('model', n.model, 'params', n.params);  break
            end
        end
        assert(~isempty(rm), 'no model/params in %s (top level: %s)', ...
               f, strjoin(fieldnames(V).', ', '));
    end
    clear V
    gv = rmModelGain(rm, 'voxels', voxels, 'metric', metric);
    g  = zeros(nTot, 1);
    if numel(gv) == nTot, g = gv(:);            % returned full-surface
    else,                 g(voxels) = gv(:);    % returned only the requested subset
    end
end

% ------------------------------------------------------------------------
function validate_v1(S, root)
% The V1 4-8 deg rows must reproduce the existing gainSummary.csv, which was computed
% by the independent DG_COMPUTEGAIN path. Any disagreement means the vertex set or the
% metric has drifted, and nothing downstream should be trusted until it is resolved.
    f = fullfile(root, 'gainSummary.csv');
    if ~isfile(f), fprintf('\n(no gainSummary.csv to validate against)\n'); return; end
    G = readtable(f);
    % Only the per-protocol columns are comparable: meanGain_avg now uses a geometric
    % mean where gainSummary.csv used an arithmetic one.
    k = strcmp(S.area, 'V1') & strcmp(S.eccBand, '4-8');
    T = S(k, :);
    fprintf('\nV1 4-8 deg against gainSummary.csv (independent DG_COMPUTEGAIN path)\n');
    fprintf('%-15s %8s %8s %10s %10s %9s\n', 'subject','n new','n old','mov new','mov old','|diff|');
    d = [];
    for i = 1:height(T)
        r = strcmp(G.subject, T.subject{i});
        if ~any(r), continue; end
        d(end+1) = abs(T.meanGain_mov(i) - G.meanGain_mov(r));  %#ok<AGROW>
        fprintf('%-15s %8d %8d %10.4f %10.4f %9.2e\n', T.subject{i}, T.nVertices(i), ...
                G.nVertices(r), T.meanGain_mov(i), G.meanGain_mov(r), d(end));
    end
    if ~isempty(d)
        if max(d) < 1e-6, verdict = 'REPRODUCES'; else, verdict = 'DOES NOT REPRODUCE'; end
        fprintf('max |difference| = %.3e  -->  %s\n', max(d), verdict);
    end
end
