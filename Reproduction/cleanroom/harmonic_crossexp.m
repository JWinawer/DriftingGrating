function res = harmonic_crossexp(Ddg, Dda, cfg, thetaVsource, opts)
% HARMONIC_CROSSEXP  Cross-experiment prediction and the nested context-effect test.
%
%   res = harmonic_crossexp(Ddg, Dda, cfg, thetaVsource, opts)
%
% The context-free null is that ONE set of four weights drives both experiments, the
% two differing only in the geometry of their local orientations. This function tests
% that null two ways, per subject, then bootstraps across the 8 observers.
%
% (1) CROSS-PREDICTION, both directions with equal billing:
%       fit dg -> predict da, and fit da -> predict dg.
%     Reported as variance explained in the held-out experiment, against that
%     experiment's own-fit R2 (the ceiling). Also with a single free multiplicative
%     gain, which absorbs "the polar responses are simply weaker overall" without
%     letting the four reference-frame weights differ.
%
% (2) NESTED MODELS on the concatenated data, three levels:
%       (i)   shared b1..b4 across experiments
%       (ii)  shared b1..b4 + one free gain on da
%       (iii) fully separate b1..b4 per experiment
%     (i)->(ii) isolates overall response magnitude. (ii)->(iii) is the actual
%     reference-frame context effect, and is what the manuscript's claim rests on:
%     if (iii) buys nothing over (ii), the Cartesian-vs-polar difference is geometry
%     plus a scale factor.
%
% R2 is computed against the concatenated demeaned data with a common denominator
% (SST = sum of squares of both experiments' responses), so the levels are comparable.
% No F-tests: the "n" here is millions of vertices that are not independent, so the
% honest statistics are the cross-prediction and the across-observer bootstrap.
%
% VERTEX WEIGHTING (opts.weighting, default 'equalcoverage'; see HARMONIC_WEIGHTS).
% Each subject's design and data are scaled by sqrt(w) ONCE, in SUBJECT_DESIGN, before
% anything else touches them. Every quantity below is a ratio of sums of squares or a
% least-squares solve, so that single scaling makes the fits, the R2s, the free gains
% and the design diagnostics all weighted consistently -- there is no second place
% where the weights have to be remembered.
%
% Output res
%   .r2      struct of nSubj x 1 vectors:
%            .dgSelf .daSelf                own-experiment fits (ceiling)
%            .dg2da  .da2dg                 cross-predicted, no gain
%            .dg2daG .da2dgG                cross-predicted, free gain
%            .shared .sharedGain .separate  nested levels (i)/(ii)/(iii)
%   .gain    .dg2da .da2dg .joint           nSubj x 1 fitted gains
%   .bShared      nSubj x 4     level (i) coefficients
%   .bSharedGain  nSubj x 4     level (ii) coefficients
%   .dBeta        nSubj x 4     b_dg - b_da (raw)
%   .dBetaGain    nSubj x 4     b_dg - g*b_da, gain-equalised
%   .ci           struct of bootstrap CIs on the means of the above
%   .corrDG .corrDA   4 x 4 uncentered correlation of the WEIGHTED designs (last
%                 subject), and .condDG .condDA their condition numbers
%   .names, .thetaVsource, .weighting

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'nBoot'),     opts.nBoot     = cfg.nBoot;       end
    if ~isfield(opts,'ciLevel'),   opts.ciLevel   = cfg.ciLevel;     end
    if ~isfield(opts,'seed'),      opts.seed      = 0;               end
    if ~isfield(opts,'weighting'),    opts.weighting    = 'equalcoverage'; end
    if ~isfield(opts,'weightSource'), opts.weightSource = thetaVsource;    end

    % THIS FUNCTION PAIRS THE TWO EXPERIMENTS BY POSITION: the same si indexes Ddg
    % and Dda below, so element si of every output is only meaningful if it is the
    % same person in both. That used to be guaranteed, because dg and da shared one
    % cfg.subjects. It is not guaranteed now -- dg can carry 13 observers and da 7 --
    % and the failure is silent, producing a full set of plausible numbers from
    % mismatched people. So refuse rather than assume.
    if isfield(Ddg, 'subjects') && isfield(Dda, 'subjects')
        assert_same_observers(Ddg.subjects, Dda.subjects, 'dg', 'da');
        subs = Ddg.subjects;
    else
        subs = cfg.subjects;   % older callers; no list to check against
    end
    res.subjects = subs;
    nS = numel(subs);
    res.names = {'hVv','cVo','rVt','pcVpo'};
    res.thetaVsource = thetaVsource;
    res.weighting    = lower(opts.weighting);
    res.corrDG = nan(4); res.corrDA = nan(4);
    res.condDG = NaN;    res.condDA = NaN;

    z = @() nan(nS, 1);
    res.r2 = struct('dgSelf',z(),'daSelf',z(),'dg2da',z(),'da2dg',z(), ...
                    'dg2daG',z(),'da2dgG',z(),'shared',z(),'sharedGain',z(),'separate',z());
    res.gain = struct('dg2da',z(),'da2dg',z(),'joint',z());
    res.bShared     = nan(nS, 4);
    res.bSharedGain = nan(nS, 4);
    res.bDG         = nan(nS, 4);
    res.bDA         = nan(nS, 4);

    for si = 1:nS
        [Xg, yg] = subject_design(Ddg, cfg.dg, si, thetaVsource, opts);
        [Xa, ya] = subject_design(Dda, cfg.da, si, thetaVsource, opts);
        if isempty(yg) || isempty(ya), continue; end

        [res.condDG, res.corrDG] = design_corr(Xg);
        [res.condDA, res.corrDA] = design_corr(Xa);

        sstg = sum(yg.^2);  ssta = sum(ya.^2);  sstBoth = sstg + ssta;

        % --- own-experiment fits (the ceiling for cross-prediction) ---
        bg = Xg \ yg;   res.bDG(si,:) = bg.';
        ba = Xa \ ya;   res.bDA(si,:) = ba.';
        res.r2.dgSelf(si) = 1 - sum((yg - Xg*bg).^2)/sstg;
        res.r2.daSelf(si) = 1 - sum((ya - Xa*ba).^2)/ssta;

        % --- cross-prediction, no gain and with a free gain ---
        [res.r2.dg2da(si), res.r2.dg2daG(si), res.gain.dg2da(si)] = cross_r2(Xa, ya, bg, ssta);
        [res.r2.da2dg(si), res.r2.da2dgG(si), res.gain.da2dg(si)] = cross_r2(Xg, yg, ba, sstg);

        % --- nested level (i): shared coefficients ---
        bs = [Xg; Xa] \ [yg; ya];
        res.bShared(si,:) = bs.';
        res.r2.shared(si) = 1 - (sum((yg - Xg*bs).^2) + sum((ya - Xa*bs).^2))/sstBoth;

        % --- nested level (ii): shared coefficients + one free gain on da ---
        [bsg, g] = fit_shared_gain(Xg, yg, Xa, ya);
        res.bSharedGain(si,:) = bsg.';
        res.gain.joint(si)    = g;
        res.r2.sharedGain(si) = 1 - (sum((yg - Xg*bsg).^2) + sum((ya - g*(Xa*bsg)).^2))/sstBoth;

        % --- nested level (iii): fully separate ---
        res.r2.separate(si) = 1 - (sum((yg - Xg*bg).^2) + sum((ya - Xa*ba).^2))/sstBoth;
    end

    % Coefficient differences: raw, and after equalising overall magnitude.
    res.dBeta     = res.bDG - res.bDA;
    res.dBetaGain = res.bDG - res.gain.joint .* res.bDA;

    % --- bootstrap CIs over the 8 observers ---
    res.ci = struct();
    fn = fieldnames(res.r2);
    for k = 1:numel(fn)
        res.ci.(fn{k}) = boot_row(res.r2.(fn{k}).', opts);
    end
    for nm = {'bDG','bDA','bShared','bSharedGain','dBeta','dBetaGain'}
        M = res.(nm{1});
        C = nan(size(M,2), 2);
        for j = 1:size(M,2), C(j,:) = boot_row(M(:,j).', opts); end
        res.ci.(nm{1}) = C;
    end
    res.ci.gainJoint = boot_row(res.gain.joint.', opts);
end

% ------------------------------------------------------------------------
function [X, y] = subject_design(D, expCfg, si, src, opts)
% Returns the sqrt(w)-SCALED design and data, so that every plain least-squares solve
% and every sum-of-squares ratio downstream is automatically the weighted one.
    m = D.subj == si;
    if ~any(m), X = []; y = []; return; end
    tv  = pick(D, src, m);
    tvW = pick(D, opts.weightSource, m);      % see FIT_HARMONIC_VERTEX opts.weightSource

    X = harmonic_predictors(tv, expCfg);
    Y = D.Y(m, :);
    y = Y(:);

    w  = harmonic_weights(tvW, opts.weighting);
    sw = repmat(sqrt(w), size(Y, 2), 1);
    X  = X .* sw;
    y  = y .* sw;
end

% ------------------------------------------------------------------------
function tv = pick(D, src, m)
    switch lower(src)
        case 'continuous', tv = D.tvCont(m);
        case 'binned',     tv = D.tvBin(m);
        otherwise, error('harmonic_crossexp:src', 'bad thetaV source ''%s''.', src);
    end
end

% ------------------------------------------------------------------------
function [cn, Rc] = design_corr(X)
% Condition number and uncentered correlation of the already-weighted design.
    S  = X.' * X;
    cn = cond(S);
    dn = sqrt(diag(S));
    Rc = S ./ (dn * dn.');
end

% ------------------------------------------------------------------------
function [r2, r2g, g] = cross_r2(X, y, b, sst)
% Variance of y explained by coefficients fitted on the OTHER experiment, without
% and with a single free multiplicative gain.
    yhat = X * b;
    r2   = 1 - sum((y - yhat).^2)/sst;
    den  = sum(yhat.^2);
    if den > 0, g = (yhat.' * y)/den; else, g = NaN; end
    r2g  = 1 - sum((y - g*yhat).^2)/sst;
end

% ------------------------------------------------------------------------
function [b, g] = fit_shared_gain(Xg, yg, Xa, ya)
% y_dg = X_dg*b ; y_da = g * X_da*b. Bilinear in (b,g), so alternate. dg's gain is
% fixed at 1, which is what identifies the scale of b.
    g = 1;
    for it = 1:100
        b = [Xg; g*Xa] \ [yg; ya];
        p = Xa * b;
        den = sum(p.^2);
        if den <= 0, break; end
        gNew = (p.' * ya)/den;
        if abs(gNew - g) < 1e-10, g = gNew; break; end
        g = gNew;
    end
    b = [Xg; g*Xa] \ [yg; ya];
end

% ------------------------------------------------------------------------
function ci = boot_row(x, opts)
    if all(isfinite(x))
        ci = bootstrap_ci(x, opts.nBoot, opts.ciLevel, opts.seed);
    else
        ci = [NaN NaN];
    end
end
