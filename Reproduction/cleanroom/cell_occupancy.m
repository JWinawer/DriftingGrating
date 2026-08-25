function [N, subjects] = cell_occupancy(varargin)
% CELL_OCCUPANCY  Vertices per (observer x polar-angle ROI) for any visual map.
%
%   N = cell_occupancy('area','MT')                  % 4-8 deg, nSubj x nPA
%   N = cell_occupancy('area','hV4','eccRange',[2 10])
%
% Counts only. Occupancy depends on area membership, eccentricity and pRF R2 --
% never on the betas -- so it is computable from labels_<subject>.mat and
% ret_<subject>.mat alone, without the runbetas_areas_* extractions. That matters
% because the coverage question is asked about maps whose beta extraction may not
% be present on the machine asking it.
%
% The inclusion filter, the hemisphere offset and the Benson -> conventional angle
% conversion are the same ones LOAD_RUNBETAS_AREA and COLLECT_RUNWISE_BETAS_AREAS
% apply, so N reproduces the empty-cell and median-vertex columns of
% ../supplement/spec_areas_coverage_spec.csv exactly (checked for V1, V3b, hV4 and
% MT in both bands).
%
% AREA NAMING follows Support/allsubjectsTable.csv, as the extraction does: MT and
% MST are the pMT / pMST labels, V3a and V3b are lower case.

    AREA_LABELS = struct('V1','V1', 'V2','V2', 'V3','V3', 'V3a','V3a', ...
                         'V3b','V3b', 'hV4','hV4', 'MT','pMT', 'MST','pMST');

    p = inputParser;
    p.addParameter('area', 'V1', @ischar);
    p.addParameter('eccRange', [], @(x) isempty(x) || numel(x)==2);
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.parse(varargin{:});
    opt = p.Results;

    cfg = config_repro();
    if ~isempty(opt.eccRange), cfg.eccRange = opt.eccRange; end
    assert(isfield(AREA_LABELS, opt.area), 'cell_occupancy:area', ...
           'unknown area ''%s''; expected one of %s.', opt.area, ...
           strjoin(fieldnames(AREA_LABELS).', ' '));
    lbl = AREA_LABELS.(opt.area);

    subjects = cfg.subjects;
    nS = numel(subjects);  nP = numel(cfg.paBins);
    N  = zeros(nS, nP);

    for si = 1:nS
        s  = subjects{si};
        Lb = load(fullfile(opt.root, sprintf('labels_%s.mat', s)));
        R  = load(fullfile(opt.root, sprintf('ret_%s.mat', s)), ...
                  'eccen','vexpl','angle_adj','nLH');
        fl = ['lh_' lbl '_REmanual'];  fr = ['rh_' lbl '_REmanual'];
        assert(isfield(Lb, fl) && isfield(Lb, fr), 'cell_occupancy:label', ...
               '%s has no %s / %s.', s, fl, fr);
        v = unique([double(Lb.(fl)(:)); double(Lb.(fr)(:)) + double(R.nLH)]);

        good = double(R.eccen(v)) >= cfg.eccRange(1) & double(R.eccen(v)) <= cfg.eccRange(2) ...
             & double(R.vexpl(v)) > cfg.r2min;
        v = v(good);
        if isempty(v), continue; end

        conv = mod(90 - double(R.angle_adj(v)), 360);        % Benson -> conventional
        [~, w] = min(abs(mod(conv - cfg.paBins(:).' + 180, 360) - 180), [], 2);
        N(si,:) = accumarray(w(:), 1, [nP 1]).';
    end
end
