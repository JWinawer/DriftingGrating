function res = fit_harmonic_vertex(D, cfg, expCfg, thetaVsource, opts)
% FIT_HARMONIC_VERTEX  Fit the four-term harmonic model to per-vertex responses.
%
%   res = fit_harmonic_vertex(D, cfg, expCfg, thetaVsource, opts)
%
% Fits, per subject, the model
%   y_vk = b1*cos(2*theta) + b2*cos(4*theta)
%        + b3*cos(2*(theta-thetaV)) + b4*cos(4*(theta-thetaV))
% by weighted least squares over that subject's (nVertex * 4) demeaned responses,
% then summarises across the 8 observers with a bootstrap (cfg.nBoot, cfg.ciLevel).
%
% The fit is PER SUBJECT and then averaged, not pooled over all vertices at once:
% V1 vertex counts range about 3.4k-8.1k across observers, so a pooled fit would
% silently weight observers by V1 surface area. This also matches the repo's
% resample-the-8-observers convention (BOOTSTRAP_CI, FIT_LME_FIG7:40-54).
%
% VERTEX WEIGHTING (opts.weighting, see HARMONIC_WEIGHTS). The default
% 'equalcoverage' gives every 15-deg polar-angle bin the same total weight, which
% reproduces the equal-weighting-of-wedges of the manuscript ROI analysis and makes the
% b1/b3 design correlation ~0 instead of the ~+0.35 that natural V1 vertex density
% imposes. 'natural' recovers the old unweighted behaviour. Weights enter the fit, the
% R2 and the collinearity diagnostics alike, via b = (X.*sqrt(w))\(y.*sqrt(w)); each
% vertex contributes nOri stacked rows and all of them carry that vertex's weight.
%
% Inputs
%   D            : from HARMONIC_VERTEX_DATA.
%   thetaVsource : 'continuous' -> D.tvCont, each vertex's own pRF polar angle (Fit B)
%                  'binned'     -> D.tvBin,  its 22.5-deg-quantised wedge centre (Fit A)
%                  Fit A reproduces the ROI-based analysis; Fit A minus Fit B is the
%                  within-wedge local-orientation artifact, isolated (same vertices,
%                  same weighting, only thetaV differs).
%   opts         : .expanded  (default false) fit the full 8-column harmonic basis
%                  .weighting (default 'equalcoverage'; 'natural' for unweighted)
%                  .weightSource (default '' = same as thetaVsource) which thetaV the
%                     equal-coverage weights are BINNED from, independently of which
%                     one the design uses. This matters, and the two callers want
%                     different things:
%                     - Leaving it '' makes Fit A weight the eight wedge centres
%                       equally, which is EXACTLY the manuscript ROI weighting, so
%                       Fit A then reproduces FIT_LME_FIG7 to the printed precision.
%                       That is what RUN_HARMONIC_MODEL section 1 wants.
%                     - Pinning it to 'continuous' for BOTH fits makes Fit A minus
%                       Fit B a pure thetaV manipulation, with the weights held
%                       fixed. That is what section 2 wants, because otherwise part
%                       of A - B is the weighting moving from 8 wedges to 24 bins
%                       (worth about a third of the dg horiz-vert shift).
%                  .nBoot, .ciLevel, .seed  (default cfg.nBoot, cfg.ciLevel, 0)
%
% Output res
%   .b        nSubj x nPred per-subject coefficients (NaN for dropped columns)
%   .bMean    1 x nPred  mean across subjects
%   .ci       nPred x 2  bootstrap CI of the mean
%   .r2       nSubj x 1  weighted variance of the demeaned data explained, per subject
%   .n        nSubj x 1  vertices per subject
%   .effN     nSubj x 1  Kish effective vertex count under the weights
%   .condNum  nSubj x 1  cond(X'WX)
%   .vif      nSubj x nPred  no-intercept variance inflation factors
%   .corr     nPred x nPred uncentered predictor correlation (last subject's design)
%   .wMeanCos2 nSubj x 1 weighted mean of cos(2*thetaV). This EQUALS the b1/b3 entry
%             of .corr for that subject -- the analytic identity asserted in
%             RUN_HARMONIC_MODEL section 2.
%   .names, .thetaVsource, .expName, .dropped, .weighting

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'expanded'),     opts.expanded     = false;           end
    if ~isfield(opts,'weighting'),    opts.weighting    = 'equalcoverage'; end
    if ~isfield(opts,'weightSource'), opts.weightSource = '';              end
    if ~isfield(opts,'nBoot'),        opts.nBoot        = cfg.nBoot;       end
    if ~isfield(opts,'ciLevel'),      opts.ciLevel      = cfg.ciLevel;     end
    if ~isfield(opts,'seed'),         opts.seed         = 0;               end

    tv = pick_thetaV(D, thetaVsource);
    if isempty(opts.weightSource)
        tvW = tv;
        wSrc = thetaVsource;
    else
        tvW = pick_thetaV(D, opts.weightSource);
        wSrc = opts.weightSource;
    end
    nS = numel(cfg.subjects);

    % Predictor names/count from a one-vertex probe.
    [~, probe] = harmonic_predictors(0, expCfg, opts);
    names = probe.names;
    nP    = numel(names);

    res = struct('names', {names}, 'thetaVsource', thetaVsource, ...
                 'expName', expCfg.name, 'weighting', lower(opts.weighting), ...
                 'weightSource', wSrc);
    res.b         = nan(nS, nP);
    res.r2        = nan(nS, 1);
    res.n         = zeros(nS, 1);
    res.effN      = nan(nS, 1);
    res.condNum   = nan(nS, 1);
    res.vif       = nan(nS, nP);
    res.wMeanCos2 = nan(nS, 1);
    res.dropped   = false(1, nP);

    for si = 1:nS
        m = D.subj == si;
        res.n(si) = nnz(m);
        if res.n(si) < nP, continue; end

        Y  = D.Y(m, :);
        tvS = tv(m);
        X  = harmonic_predictors(tvS, expCfg, opts);
        y  = Y(:);

        % Vertex weights, one per vertex, replicated across that vertex's nOri stacked
        % rows. sqrt(w) is what multiplies the design: minimising sum((X.*sqrt(w))*b -
        % y.*sqrt(w)).^2 is minimising sum(w.*(y - X*b).^2), which is what we want.
        % Multiplying by w itself would minimise sum(w.^2 .* r.^2) instead.
        % Weights are binned from tvW (see opts.weightSource); wMeanCos2 uses the
        % DESIGN angle tvS, because it is the design's b1/b3 correlation it predicts.
        [wV, wInfo] = harmonic_weights(tvW(m), opts.weighting);
        res.effN(si)      = wInfo.effN;
        res.wMeanCos2(si) = sum(wV .* cosd(2*tvS)) / sum(wV);
        sw = repmat(sqrt(wV), size(Y, 2), 1);          % (nVertex*nOri) x 1

        % Drop columns that are identically zero. sin(4*theta) is zero at all four
        % 45-deg-spaced orientations of the Cartesian experiment, and sin(4*d) is
        % likewise zero for the polar one -- a property of the sampling, not the data.
        keepCol = ~all(abs(X) < 1e-12, 1);
        res.dropped = res.dropped | ~keepCol;

        Xw = X(:,keepCol) .* sw;
        yw = y .* sw;

        b = Xw \ yw;
        res.b(si, keepCol) = b.';

        sst = sum(yw.^2);                      % y is already per-vertex mean-zero
        ssr = sum((yw - Xw*b).^2);
        res.r2(si) = 1 - ssr/sst;

        [res.condNum(si), v, Rc] = design_diagnostics(Xw);
        res.vif(si, keepCol) = v;
        res.corr = nan(nP);
        res.corr(keepCol, keepCol) = Rc;
    end

    res.bMean = mean(res.b, 1, 'omitnan');
    res.ci    = nan(nP, 2);
    for j = 1:nP
        xj = res.b(:, j).';
        if all(isfinite(xj))
            res.ci(j, :) = bootstrap_ci(xj, opts.nBoot, opts.ciLevel, opts.seed);
        end
    end
end

% ------------------------------------------------------------------------
function tv = pick_thetaV(D, src)
    switch lower(src)
        case 'continuous', tv = D.tvCont;
        case 'binned',     tv = D.tvBin;
        otherwise
            error('fit_harmonic_vertex:badSource', ...
                  'thetaVsource must be ''continuous'' or ''binned'', got ''%s''.', src);
    end
end

% ------------------------------------------------------------------------
function [cn, vif, Rc] = design_diagnostics(X)
% Condition number and no-intercept VIFs. The model has no intercept and the
% predictors are mean-zero across conditions by construction, so the uncentered
% correlation is the right normalisation here. X is the ALREADY-WEIGHTED design
% (X.*sqrt(w)), so X'*X is the weighted Gram X'WX and Rc describes the design that
% was actually fitted.
    S  = X.' * X;
    cn = cond(S);
    dn = sqrt(diag(S));
    Rc = S ./ (dn * dn.');
    vif = diag(inv(Rc)).';
end
