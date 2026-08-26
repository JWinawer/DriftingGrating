function S = spec_profiles(varargin)
% SPEC_PROFILES  Everything Figures 5/6 need, under the settled specification.
%
%   S = spec_profiles()                         % V1, 4-8 deg, the specification
%   S = spec_profiles('area','V2')
%   S = spec_profiles('area','V3a','eccRange',[2 10])
%   S = spec_profiles('route','roi')            % the ALTERNATIVE route, for comparison
%
% THE SPECIFICATION (../SPECIFICATION.md section 1), applied unchanged to every map:
%   gain        on, per observer x MAP, applied at the observer boundary
%   model       four-term harmonic, fitted per vertex, CONTINUOUS thetaV
%   weighting   equal coverage at 45 deg -- the polar-angle ROIs themselves
%   fitting     per observer, then averaged across observers
%
% TWO ROUTES, and both are always computed. 'harmonic' (the specification) fits the
% four-term model per vertex with continuous thetaV; 'roi' bins vertices into the eight
% polar-angle wedges and contrasts within each, which is the manuscript route. They are
% the same estimator re-parameterised and agree exactly on complete data; they diverge
% in proportion to empty cells, because the ROI route averages over whichever wedges an
% observer happens to have and the harmonic route never bins. 'route' selects which one
% .asym carries; .asymHarmonic and .asymROI always carry both, so the comparison costs
% nothing extra.
%
% Returns per-observer quantities only. Nothing here averages across observers or
% forms an interval; that is SPEC_TABLES' job, and keeping the boundary sharp is what
% lets the gain rescaling stay where the specification puts it (it does NOT commute
% with precision weighting, so it can never be applied to a finished group estimate).
%
% WHY THIS EXISTS ALONGSIDE DIAGNOSE_WITHIN_OBSERVER_ERROR. That function returns the
% four asymmetries per observer and their measurement error, which is what the tables
% need. Figures 5/6 also need the polar-angle-RESOLVED profile behind each asymmetry,
% both as the model predicts it and as the data show it. This computes both, from the
% same files and the same vertex set (LOAD_RUNBETAS_AREA), and asserts that its
% asymmetries match that function's to 1e-12 -- so the figure and the table cannot
% drift apart.
%
% Output S, with S.dg and S.da each carrying
%   .b        nSubj x 4   fitted coefficients, gain-scaled
%   .asym     nSubj x 4   the four asymmetries (2*b read off at the 8 ROI centres)
%   .mPro/.mCon/.mDiff   nSubj x nPA x 4  MODEL profile at the 8 ROI centres
%   .oPro/.oCon/.oDiff   nSubj x nPA x 4  OBSERVED wedge profile, same classification
%   .oAsym    nSubj x 4   wedge-mean of .oDiff -- the ROI route, for comparison
%   .asymHarmonic/.asymROI  nSubj x 4  both routes, always; .asym is whichever
%             'route' selected
%   .sigma    nSubj x 4   within-observer SE of each asymmetry, from bootstrapping
%             RUNS on the selected route -- the input precision weighting needs
%   .fine     struct: .centres, .obs, .mdl  (nSubj x nFine x 3) the three contrasts
%             that are measurable at a single polar angle, data and model, plus
%             .denseCentres / .mdlDense, the model on a 0.5-deg grid for plotting
%   .dense    struct: .centres (1 x nDense) and .mPro/.mCon (nSubj x nDense x 4),
%             the SAME model profile as .mPro/.mCon but on a 0.5-deg polar-angle
%             grid, so Figures 5/6 can draw the fit as a curve rather than as
%             straight segments between the eight wedge centres. It passes exactly
%             through those centres, which is asserted per observer.
%   .fineLbl  1x3 labels for those contrasts
%   .nVert    nSubj x nPA vertices per wedge;  .nTot nSubj x 1
% and at top level .subjects, .paBins, .area, .eccRange, .gainScale, .names.

    p = inputParser;
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('area', 'V1', @ischar);
    p.addParameter('eccRange', [], @(x) isempty(x) || numel(x)==2);
    p.addParameter('nFine', 24, @isscalar);      % 15 deg, nests inside the 45 deg bins
    p.addParameter('gain', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('route', 'harmonic', @(x) any(strcmpi(x,{'harmonic','roi'})));
    p.addParameter('verify', true, @islogical);
    % SIMULATED CELL LOSS, for DIAGNOSE_CELL_LOSS only: nSubj x nPA logical marking
    % (observer x ROI) cells to empty out. Empty by default, and empty means the
    % settled specification, unchanged. nBoot is exposed for the same reason -- the
    % simulation runs this pipeline dozens of times and the run bootstrap is ~75% of
    % the cost -- and it stays at 500 for anything reported.
    p.addParameter('dropCells', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
    p.addParameter('nBoot', 500, @isscalar);
    p.parse(varargin{:});
    opt = p.Results;

    cfg = config_repro();
    if ~isempty(opt.eccRange), cfg.eccRange = opt.eccRange; end
    wBins = 8;                                   % 45 deg -- the specification
    expn  = {'dg','da'};
    nS    = numel(cfg.subjects);
    nP    = numel(cfg.paBins);
    nF    = opt.nFine;
    drop  = opt.dropCells;
    if isempty(drop), drop = false(nS, nP); else, drop = logical(drop); end
    assert(isequal(size(drop), [nS nP]), 'spec_profiles:dropCells', ...
           'dropCells must be %d x %d.', nS, nP);

    % --- gain, per observer x map, exactly as DIAGNOSE_WITHIN_OBSERVER_ERROR does it
    if opt.gain
        band   = sprintf('%g-%g', cfg.eccRange(1), cfg.eccRange(2));
        gscale = observer_gain_weights(cfg, opt.area, band);
        if ~all(isfinite(gscale))
            warning('spec_profiles:gainFallback', ...
                    'no per-map gain for %s %s; falling back to the V1 scalar.', ...
                    opt.area, band);
            gscale = observer_gain_weights(cfg);
        end
    else
        gscale = ones(nS,1);
    end

    S.subjects = cfg.subjects;  S.paBins = cfg.paBins;  S.area = opt.area;
    S.eccRange = cfg.eccRange;  S.gainScale = gscale;
    S.names = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    S.route = lower(opt.route);
    S.hasModel = strcmp(S.route, 'harmonic');   % the ROI route fits no per-vertex model

    fineCtr = (0:nF-1) * (360/nF) + (360/nF)/2;      % bin centres, conventional deg

    for ei = 1:2
        en = expn{ei};  expCfg = cfg.(en);
        E = struct('b', nan(nS,4), 'asym', nan(nS,4), 'oAsym', nan(nS,4), ...
                   'asymHarmonic', nan(nS,4), 'asymROI', nan(nS,4), 'sigma', nan(nS,4), ...
                   'mPro', nan(nS,nP,4), 'mCon', nan(nS,nP,4), 'mDiff', nan(nS,nP,4), ...
                   'oPro', nan(nS,nP,4), 'oCon', nan(nS,nP,4), 'oDiff', nan(nS,nP,4), ...
                   'nVert', zeros(nS,nP), 'nTot', zeros(nS,1));
        obsF = nan(nS,nF,3);  mdlF = nan(nS,nF,3);  wedF = nan(nS,nP,3);
        % The model is evaluated on a DENSE grid for display. Reading it off the same
        % 15-deg bin centres as the data aliases the fourth harmonic: cos(4*thetaV)
        % sampled every 15 deg advances 60 deg per step and so revisits the same three
        % values, which plots as a step function rather than a cosine. The fit is
        % unaffected -- this is display sampling only.
        denseCtr = 0:0.5:360;
        mdlD = nan(nS, numel(denseCtr), 3);
        dPro = nan(nS, numel(denseCtr), 4);   % the pro/con classes on that same grid
        dCon = nan(nS, numel(denseCtr), 4);

        for si = 1:nS
            [A, ok] = load_one(cfg, opt.root, cfg.subjects{si}, en, opt.area, drop(si,:));
            if ~ok, continue; end

            % --- run-averaged demeaned per-vertex responses -----------------
            B   = mean(A.runBeta, 3, 'omitnan');
            col = cfg.(en).oriIdx - 25 + 8;              % CONTRASTS 26..29 -> cols 9..12
            C   = double(B(:, col));                     % blank cancels in the demeaning
            Y   = (C - mean(C,2)) * gscale(si);          % gain at the observer boundary
            tv  = A.thetaV;

            % --- the fit: continuous thetaV, equal coverage at 45 deg -------
            o    = struct('expanded', false, 'weighting', 'equalcoverage');
            X    = harmonic_predictors(tv, expCfg, o);
            wV   = harmonic_weights(tv, 'equalcoverage', wBins);
            sw   = repmat(sqrt(wV(:)), size(Y,2), 1);
            keep = ~all(abs(X) < 1e-12, 1);              % sin columns vanish by design
            b    = nan(1, size(X,2));
            b(keep) = ((X(:,keep) .* sw) \ (Y(:) .* sw)).';
            E.b(si,:) = b;

            % --- MODEL profile at the eight ROI centres --------------------
            Yh  = predict_harmonic(b, cfg.paBins(:), expCfg, o);      % nPA x nOri
            Am  = compute_asymmetries(Yh.', cfg, expCfg);
            % --- OBSERVED wedge profile, identical classification ----------
            % Equal-coverage weighting at 45 deg makes every vertex in a wedge carry
            % the same weight, so the weighted wedge aggregate IS the plain mean.
            Mo = nan(size(Y,2), nP);
            for pIdx = 1:nP
                m = A.wedge == pIdx;
                E.nVert(si,pIdx) = nnz(m);
                if any(m), Mo(:,pIdx) = mean(Y(m,:), 1).'; end
            end
            Ao = compute_asymmetries(Mo, cfg, expCfg);

            for j = 1:4
                f = Am.order{j};
                E.mPro(si,:,j)  = Am.(f).pro(:).';
                E.mCon(si,:,j)  = Am.(f).con(:).';
                E.mDiff(si,:,j) = Am.(f).diff(:).';
                E.oPro(si,:,j)  = Ao.(f).pro(:).';
                E.oCon(si,:,j)  = Ao.(f).con(:).';
                E.oDiff(si,:,j) = Ao.(f).diff(:).';
                E.asymHarmonic(si,j) = mean(Am.(f).diff, 'omitnan');
                E.asymROI(si,j)      = mean(Ao.(f).diff, 'omitnan');
                E.oAsym(si,j)        = E.asymROI(si,j);
            end
            E.nTot(si) = size(Y,1);

            % --- the SAME model profile, on a dense polar-angle grid --------
            % Read off at the eight wedge centres, mPro/mCon plot as straight
            % segments between eight points; the model itself is continuous in
            % thetaV, so it can be drawn as the curve it is. CLASS_PROFILES is
            % that curve, and the assertion below is what ties it to the points.
            [pD, cD] = class_profiles(b, denseCtr, cfg, o);
            dPro(si,:,:) = reshape(pD, [1 size(pD)]);
            dCon(si,:,:) = reshape(cD, [1 size(cD)]);
            if opt.verify
                [p8, c8] = class_profiles(b, cfg.paBins, cfg, o);
                dChk = max(abs([p8(:) - reshape(E.mPro(si,:,:), [], 1); ...
                                c8(:) - reshape(E.mCon(si,:,:), [], 1)]));
                assert(dChk < 1e-12, 'spec_profiles:denseMismatch', ...
                    ['%s %s: the dense model curve misses the plotted wedge ' ...
                     'centres by %.3g. The curve and the markers would disagree.'], ...
                    cfg.subjects{si}, en, dChk);
            end

            % --- the three contrasts measurable at a single polar angle ----
            % Only 3 of the 4 harmonic degrees of freedom are observable per vertex
            % (../supplement/SUPPLEMENT_harmonic_model.md section S2.3), so a profile
            % against continuous thetaV can show these three and no more. The fourth
            % coefficient is identified ACROSS vertices, by the thetaV modulation of
            % the first and third -- which is exactly what the curves below display.
            edges = linspace(0, 360, nF+1);
            [~, fb] = histc(mod(tv,360), edges); %#ok<HISTC>
            fb(fb > nF) = nF;
            Yf = predict_harmonic(b, fineCtr(:), expCfg, o);          % nF x nOri
            for fi = 1:nF
                m = fb == fi;
                if nnz(m) >= 3, obsF(si,fi,:) = abc(mean(Y(m,:),1)); end
            end
            mdlF(si,:,:) = abc(Yf);
            mdlD(si,:,:) = abc(predict_harmonic(b, denseCtr(:), expCfg, o));
            % The same three contrasts at WEDGE resolution: what the ROI route sees.
            % Piecewise constant within each 45 deg wedge, which is the assumption the
            % ROI route makes about polar-angle structure and the harmonic route does not.
            wedF(si,:,:) = abc(Mo.');
        end

        if S.hasModel, E.asym = E.asymHarmonic; else, E.asym = E.asymROI; end

        if strcmp(en,'dg')
            E.fineLbl = {'horizontal - vertical','45\circ - 135\circ','cardinal - oblique'};
        else
            E.fineLbl = {'radial - tangential','ccw spiral - cw spiral', ...
                         'polar cardinal - polar oblique'};
        end
        E.fine = struct('centres', fineCtr, 'obs', obsF, 'mdl', mdlF, ...
                        'denseCentres', denseCtr, 'mdlDense', mdlD, ...
                        'wedgeCentres', cfg.paBins, 'wedge', wedF);
        E.dense = struct('centres', denseCtr, 'mPro', dPro, 'mCon', dCon);
        S.(en) = E;
    end

    % --- within-observer sigma, on the selected route --------------------------
    % Measured by resampling RUNS (DIAGNOSE_WITHIN_OBSERVER_ERROR), not inferred: the
    % condition-wise design has no replication, so a model fitted to it cannot separate
    % a reliable observer from a noisy one. Precision weighting needs these; equal
    % weighting does not, but they are cheap here and reporting both is the point.
    W = diagnose_within_observer_error('root', opt.root, 'area', opt.area, ...
            'eccRange', opt.eccRange, 'route', S.route, 'thetaV', tvSrc(S.route), ...
            'gain', opt.gain, 'weighting', 'equalcoverage', 'quiet', true, ...
            'dropCells', drop, 'nBoot', opt.nBoot);
    for ei = 1:2
        en = expn{ei};
        S.(en).sigma = W.seBoot(:,:,ei);
        if opt.verify, checkRoute(S.(en).asym, W.full(:,:,ei), opt.area, en, S.route); end
    end
end

% ------------------------------------------------------------------------
function [pro, con] = class_profiles(b, tv, cfg, o)
% CLASS_PROFILES  The four asymmetries' pro and con classes, continuous in polar angle.
%
%   [pro, con] = class_profiles(b, tv, cfg, o)   % pro, con are numel(tv) x 4
%
% This is the same quantity COMPUTE_ASYMMETRIES returns as .pro and .con, but readable
% at any polar angle instead of only at the eight wedge centres. It exists so Figures
% 5/6 can draw the fitted model as a curve.
%
% WHY COMPUTE_ASYMMETRIES CANNOT DO THIS. It labels each presented stimulus by exact
% equality -- horizontal is local orientation == 0, radial is offset from the radius
% == 0 -- so a class is populated only where a stimulus lands exactly on it. That is
% every 45 deg and nowhere in between, and at an intermediate polar angle the label
% would find no stimulus and return NaN. The classes are not undefined there; the
% stimulus set simply does not sample them.
%
% WHAT IS EVALUATED INSTEAD. The fitted model is a continuous function of the bar
% orientation theta and the polar angle thetaV,
%     y = b1*cos(2*theta) + b2*cos(4*theta) + b3*cos(2*(theta-thetaV)) + b4*cos(4*(theta-thetaV))
% so each class can be read off at the orientation that DEFINES it, for any thetaV:
%     horizontal theta=0     vertical theta=90    oblique   mean of 45 and 135
%     radial     theta=thetaV  tangential thetaV+90  polar-oblique mean of +/-45
% Those two frames are exactly the geometry cfg.dg and cfg.da already carry, so they
% are used here as FRAMES rather than as experiments: cfg.dg supplies the four fixed
% Cartesian orientations and cfg.da the four that rotate with the radius, whichever
% experiment b was fitted to. No new convention is introduced, and none is restated.
%
% At the eight wedge centres both frames coincide with the presented stimuli, so this
% reproduces COMPUTE_ASYMMETRIES exactly -- the curve passes through the plotted
% markers rather than near them. SPEC_PROFILES asserts that, per observer, to 1e-12.
%
% This is display sampling only. Nothing fitted, tabled or tested reads it.
    Yc = predict_harmonic(b, tv(:), cfg.dg, o);   % Cartesian frame: [H  V  ob ob]
    Yp = predict_harmonic(b, tv(:), cfg.da, o);   % radial frame:    [R  T  po po]
    % Columns in COMPUTE_ASYMMETRIES' order: HV, cardObl, radTan, polcardPolobl.
    pro = [Yc(:,1), mean(Yc(:,1:2),2), Yp(:,1), mean(Yp(:,1:2),2)];
    con = [Yc(:,2), mean(Yc(:,3:4),2), Yp(:,2), mean(Yp(:,3:4),2)];
end

% ------------------------------------------------------------------------
function t = tvSrc(route)
% The ROI route bins thetaV to the wedge centres by construction; the harmonic route
% is the specification's continuous one. Keeping these paired is what makes sigma the
% measurement error OF the estimate being reported, rather than of a different one.
    if strcmp(route,'roi'), t = 'binned'; else, t = 'continuous'; end
end

% ------------------------------------------------------------------------
function checkRoute(a, w, area, en, route)
% The figure and the table must not be able to drift apart: assert that the
% asymmetries computed here equal the ones PRECISION_WEIGHTED_TABLE consumes.
    m = isfinite(a) & isfinite(w);
    if ~any(m(:)), return; end
    d = max(abs(a(m) - w(m)));
    assert(d < 1e-9, 'spec_profiles:mismatch', ...
        ['%s %s (%s): asymmetries differ from DIAGNOSE_WITHIN_OBSERVER_ERROR by %.3g. ' ...
         'The figure and the tables would show different numbers.'], area, en, route, d);
    fprintf('spec_profiles(%s, %s, %s): agrees with diagnose_within_observer_error to %.1e\n', ...
            area, en, route, d);
end

% ------------------------------------------------------------------------
function v = abc(Y)
% The three single-polar-angle contrasts, from responses in oriCols order.
% Both experiments order oriCols as [cardinal-ish pair, oblique pair]:
%   dg {horizontal, vertical, rightleaning, leftleaning}
%   da {pinwheel(radial), annulus(tangential), ccspiral, cspiral}
% so one expression serves both, and reads as the matched-frame contrast in each.
    v = [Y(:,1) - Y(:,2), Y(:,3) - Y(:,4), ...
         mean(Y(:,1:2),2) - mean(Y(:,3:4),2)];
    v = reshape(v, [size(Y,1), 1, 3]);
    if size(Y,1) == 1, v = reshape(v, [1 1 3]); end
end

% ------------------------------------------------------------------------
function [A, ok] = load_one(cfg, root, subj, en, area, dropWedges)
    fA = fullfile(root, sprintf('runbetas_areas_%s_%s.mat', subj, en));
    f1 = fullfile(root, sprintf('runbetas_%s_%s.mat',       subj, en));
    if isfile(fA), f = fA; elseif isfile(f1), f = f1;
    else, A = struct(); ok = false; return
    end
    [A, ok] = load_runbetas_area(load(f), cfg, en, root, area, dropWedges);
end

