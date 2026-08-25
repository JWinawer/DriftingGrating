function P = spec_pooled(varargin)
% SPEC_POOLED  One harmonic fit to all observers' vertices at once.
%
%   P = spec_pooled()                                  % V1, 4-8 deg
%   P = spec_pooled('weighting','equalcell')
%   P = spec_pooled('data', D, 'dropCells', mask)      % reuse a loaded dataset
%
% THE ALTERNATIVE ESTIMATOR. The settled specification fits the four-term harmonic
% model once per observer and averages the eight fits (../SPECIFICATION.md section 1).
% This stacks every observer's vertices into ONE weighted least-squares fit instead.
% Both are computed on every call and returned side by side, because the whole point
% of the comparison is that they are the same model differing only in where the
% observers are combined.
%
% THE WEIGHTS ARE TWO-LEVEL. Within an observer, equal coverage over the eight
% polar-angle bins, exactly as the specification does it (HARMONIC_WEIGHTS at 45 deg)
% -- so cortical magnification does not tilt the fit toward the horizontal meridian.
% Across observers, one of two choices, and the difference between them is the whole
% missing-data question:
%
%   'equalobserver'  every observer carries total weight 1/nObs, spread equally over
%                    the ROIs they HAVE. An observer with three surviving ROIs still
%                    counts as much as one with eight, and those three ROIs each get
%                    2.67x the weight they would otherwise carry.  <- DEFAULT
%   'equalcell'      every (observer x ROI) cell carries the same weight, 1/(8*nObs).
%                    A missing cell contributes nothing and is not compensated, so an
%                    observer with three surviving ROIs carries 3/8 of a full one.
%
% On complete data the two are identical (asserted below), so they are not two
% analyses of V1 -- they are two ways of handling exactly the loss this exists to
% study.
%
% WITHOUT THE CROSS-OBSERVER NORMALISATION, pooling would weight observers by vertex
% count: HARMONIC_WEIGHTS rescales to mean(w) == 1, hence sum(w) == nVertex, which is
% harmless in a per-observer fit and is a 1.6x observer-weight range here if you stack
% the fits without rescaling. That is the "map size" bias, and it is removed here
% rather than assumed absent.
%
% WHAT POOLING ACTUALLY DOES TO THE OBSERVER WEIGHTS. Weighted least squares is linear
% in y, so with M_i = X_i' W_i X_i the pooled fit is exactly
%
%       b_pool = (sum_i M_i)^{-1} (sum_i M_i b_i)
%
% a MATRIX-weighted average of the same per-observer fits the specification averages
% equally. That is the mechanism behind any robustness it has -- an observer whose
% coverage collapsed has a small M_i and is quietly down-weighted -- and it is also
% the reason the two estimators answer slightly different questions. P.leverage
% reports diag((sum M)^{-1} M_i) per coefficient, which sums to exactly 1 over
% observers and says how much of each coefficient each observer is supplying.
%
% INFERENCE. A single pooled fit has no per-observer estimate, so the specification's
% t test across observers has nothing to run on. The interval here is a DELETE-ONE-
% OBSERVER JACKKNIFE, refitting the pool without each observer in turn, on n-1 df.
% That keeps the observer as the unit of error, which standing fact 6 requires, but it
% is not the same interval as the specification's and should not be quoted as if it
% were.
%
% Returns P with, for each of .dg and .da:
%   .b           1 x 4  pooled coefficients
%   .asym        1 x 4  pooled asymmetries (2*b at the eight ROI centres)
%   .asymAvg     1 x 4  the SPECIFICATION's estimate: mean over per-observer fits
%   .bObs/.asymObs  nSubj x 4  the per-observer fits both estimators are built from
%   .jackSE/.jackLo/.jackHi/.jackP  1 x 4  jackknife interval on the pooled estimate
%   .asymJack/.jackIdx  the delete-one replicates behind it, and which observer each
%                omits, so the caller can jackknife the dg-minus-da difference too
%   .tLo/.tHi/.tP   1 x 4  the specification's t interval on .asymAvg, for comparison
%   .leverage    nSubj x 4, columns summing to 1
%   .nOcc        nSubj x 1  occupied weighting bins per observer
%   .nObs        observers contributing
% and at top level .weighting, .area, .eccRange, .subjects, .data (reusable).

    p = inputParser;
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('area', 'V1', @ischar);
    p.addParameter('eccRange', [], @(x) isempty(x) || numel(x)==2);
    p.addParameter('gain', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('dropCells', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
    p.addParameter('weighting', 'equalobserver', ...
                   @(x) any(strcmpi(x, {'equalobserver','equalcell'})));
    p.addParameter('data', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('jackknife', true, @islogical);
    p.parse(varargin{:});
    opt = p.Results;

    cfg = config_repro();
    if ~isempty(opt.eccRange), cfg.eccRange = opt.eccRange; end
    nS = numel(cfg.subjects);  nP = numel(cfg.paBins);
    drop = opt.dropCells;
    if isempty(drop), drop = false(nS, nP); else, drop = logical(drop); end
    assert(isequal(size(drop), [nS nP]), 'spec_pooled:dropCells', ...
           'dropCells must be %d x %d.', nS, nP);

    % Loading dominates the cost and does not depend on the deletion, so the whole
    % dataset is loaded ONCE and every deletion is applied in memory. DIAGNOSE_POOLED_FIT
    % runs hundreds of deletions off a single load.
    if isempty(opt.data), D = load_pooled_data(cfg, opt); else, D = opt.data; end

    P.weighting = lower(opt.weighting);  P.area = opt.area;
    P.eccRange = cfg.eccRange;  P.subjects = cfg.subjects;  P.data = D;
    P.names = {'horiz-vert','card-obl','rad-tang','polc-polo'};

    for en = {'dg','da'}
        e = en{1};
        F = fit_pool(D, e, cfg, drop, P.weighting, true(nS,1));
        E = struct('b', F.b, 'bObs', F.bObs, 'leverage', F.leverage, ...
                   'nOcc', F.nOcc, 'nObs', F.nObs, 'nVert', F.nVert);
        E.asym    = asym_of(F.b, cfg, cfg.(e));
        E.asymObs = nan(nS, 4);
        for si = 1:nS
            if all(isfinite(F.bObs(si,:)))
                E.asymObs(si,:) = asym_of(F.bObs(si,:), cfg, cfg.(e));
            end
        end
        ok = all(isfinite(E.asymObs), 2);
        E.asymAvg = mean(E.asymObs(ok,:), 1);

        % The specification's own interval on the per-observer average, so the two
        % estimators can be read against each other rather than against nothing.
        [E.tLo, E.tHi, E.tP] = t_interval(E.asymObs(ok,:));

        if opt.jackknife && nnz(ok) >= 3
            idx = find(ok);  n = numel(idx);
            aJ = nan(n, 4);
            for k = 1:n
                sel = ok;  sel(idx(k)) = false;
                aJ(k,:) = asym_of(fit_pool(D, e, cfg, drop, P.weighting, sel).b, ...
                                  cfg, cfg.(e));
            end
            se = sqrt((n-1)/n * sum((aJ - mean(aJ,1)).^2, 1));
            tc = tinv(0.975, n-1);
            E.jackSE = se;  E.jackLo = E.asym - tc*se;  E.jackHi = E.asym + tc*se;
            E.jackP  = 2 * tcdf(-abs(E.asym ./ se), n-1);
            % The REPLICATES are returned, not just the SE, so a caller can jackknife
            % any function of the estimate -- the dg-minus-da context effect above all,
            % which a pooled fit cannot form within observer and so must form here.
            E.asymJack = aJ;  E.jackIdx = idx;
        else
            E.jackSE = nan(1,4); E.jackLo = nan(1,4);
            E.jackHi = nan(1,4); E.jackP = nan(1,4);
            E.asymJack = [];    E.jackIdx = [];
        end
        P.(e) = E;
    end
end

% ------------------------------------------------------------------------
function F = fit_pool(D, en, cfg, drop, weighting, subjSel)
% One weighted least-squares fit over the selected observers' vertices, plus the
% per-observer fits that the same normal equations contain. Returning both from one
% pass is what makes the two estimators exactly comparable: identical vertices,
% identical weights, identical design.
    nS = numel(cfg.subjects);
    M  = zeros(4);  r = zeros(4,1);
    Mi = zeros(4,4,nS);
    F.bObs = nan(nS,4);  F.nOcc = zeros(nS,1);  F.nVert = zeros(nS,1);

    for si = 1:nS
        if ~subjSel(si), continue; end
        A = D.(en){si};
        if isempty(A), continue; end
        m = ~ismember(A.wedge, find(drop(si,:)));
        if ~any(m), continue; end

        tv = A.thetaV(m);
        Y  = A.Y(m,:);
        [wv, nOcc] = cell_weights(tv, weighting, nS);
        if nOcc == 0, continue; end
        F.nOcc(si) = nOcc;  F.nVert(si) = nnz(m);

        X  = harmonic_predictors(tv, cfg.(en), struct('expanded', false));
        sw = repmat(sqrt(wv(:)), size(Y,2), 1);
        Xw = X .* sw;  yw = Y(:) .* sw;
        Mi(:,:,si) = Xw' * Xw;
        ri         = Xw' * yw;
        M = M + Mi(:,:,si);  r = r + ri;
        F.bObs(si,:) = (Mi(:,:,si) \ ri).';
    end

    F.b    = (M \ r).';
    F.nObs = nnz(all(isfinite(F.bObs), 2));
    F.leverage = nan(nS, 4);
    Minv = inv(M);
    for si = 1:nS
        if F.nVert(si) > 0, F.leverage(si,:) = diag(Minv * Mi(:,:,si)).'; end
    end
end

% ------------------------------------------------------------------------
function [w, nOcc] = cell_weights(tv, weighting, nS)
% Two-level weights. HARMONIC_WEIGHTS supplies the within-observer equal-coverage
% weighting the specification already uses, at 45 deg; this only fixes the SCALE, which
% a per-observer fit is free to ignore and a pooled fit is not.
%
% HARMONIC_WEIGHTS returns w with mean 1, i.e. sum(w) == nVertex. Rescaling to c, in
% which every OCCUPIED bin totals exactly 1, removes the vertex-count dependence:
%   c = w * nOcc / nVertex.
    [wh, info] = harmonic_weights(tv, 'equalcoverage', 8);
    nV   = numel(wh);
    nOcc = info.nBin - info.nEmpty;
    if nV == 0 || nOcc == 0, w = wh; return; end
    c = wh * (nOcc / nV);                       % every occupied bin totals 1
    switch weighting
        case 'equalobserver', w = c / (nOcc * nS);      % observer totals 1/nS
        case 'equalcell',     w = c / (info.nBin * nS); % each cell totals 1/(8*nS)
    end
end

% ------------------------------------------------------------------------
function a = asym_of(b, cfg, expCfg)
% Coefficients -> the four asymmetries, through the UNMODIFIED COMPUTE_ASYMMETRIES at
% the eight ROI centres, so the units match every other route in the pipeline.
    if any(~isfinite(b)), a = nan(1,4); return; end
    o   = struct('expanded', false, 'weighting', 'equalcoverage');
    Yh  = predict_harmonic(b, cfg.paBins(:), expCfg, o);
    Asy = compute_asymmetries(Yh.', cfg, expCfg);
    a   = nan(1,4);
    for j = 1:4, a(j) = mean(Asy.(Asy.order{j}).diff, 'omitnan'); end
end

% ------------------------------------------------------------------------
function [lo, hi, pv] = t_interval(A)
    n  = size(A,1);
    m  = mean(A,1);  se = std(A,0,1) / sqrt(n);
    tc = tinv(0.975, n-1);
    lo = m - tc*se;  hi = m + tc*se;
    pv = 2 * tcdf(-abs(m./se), n-1);
end

% ------------------------------------------------------------------------
function D = load_pooled_data(cfg, opt)
% Run-averaged, demeaned, gain-scaled responses per observer, loaded once. The gain
% is per observer x map and is applied HERE, at the observer boundary, exactly where
% SPEC_PROFILES applies it -- it must not move just because the fit did.
    nS = numel(cfg.subjects);
    if opt.gain
        band   = sprintf('%g-%g', cfg.eccRange(1), cfg.eccRange(2));
        gscale = observer_gain_weights(cfg, opt.area, band);
        if ~all(isfinite(gscale)), gscale = observer_gain_weights(cfg); end
    else
        gscale = ones(nS,1);
    end
    D.subjects = cfg.subjects;  D.gscale = gscale;
    for en = {'dg','da'}
        e = en{1};
        D.(e) = cell(nS,1);
        for si = 1:nS
            fA = fullfile(opt.root, sprintf('runbetas_areas_%s_%s.mat', cfg.subjects{si}, e));
            f1 = fullfile(opt.root, sprintf('runbetas_%s_%s.mat',       cfg.subjects{si}, e));
            if isfile(fA), f = fA; elseif isfile(f1), f = f1; else, continue; end
            [A, ok] = load_runbetas_area(load(f), cfg, e, opt.root, opt.area);
            if ~ok, continue; end
            B   = mean(A.runBeta, 3, 'omitnan');
            col = cfg.(e).oriIdx - 25 + 8;               % contrasts 26..29 -> cols 9..12
            C   = double(B(:, col));                     % blank cancels in the demeaning
            D.(e){si} = struct('Y', (C - mean(C,2)) * gscale(si), ...
                               'thetaV', A.thetaV, 'wedge', A.wedge);
        end
    end
end
