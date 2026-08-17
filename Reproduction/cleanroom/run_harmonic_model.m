function R = run_harmonic_model(doZscore)
% RUN_HARMONIC_MODEL  Per-vertex harmonic model: fit, cross-predict, report.
%
%   R = run_harmonic_model()          % raw betas (the adopted variant)
%   R = run_harmonic_model(true)      % z-scored sensitivity variant
%
% Replaces the wedge binning of Figs 5-8 with a per-vertex regression on each
% vertex's own pRF polar angle, to separate two explanations of the
% Cartesian-vs-polar difference in the orientation asymmetries:
%   (a) within-ROI geometry -- a "horizontal" grating is only exactly radial for a
%       vertex on the horizontal meridian, and the wedges are 45 deg wide;
%   (b) genuine longer-range context effects.
%
% See ../HARMONIC_MODEL.md for the model, the conventions and the conclusions.
% TEST_HARMONIC_MODEL must pass before any of this is interpreted.

    if nargin < 1 || isempty(doZscore), doZscore = false; end
    variant = 'raw'; if doZscore, variant = 'zscored'; end

    cfg = config_repro();
    T   = load_allconditions(cfg);
    R   = struct('variant', variant, 'cfg', cfg);

    Ddg = harmonic_vertex_data(T, cfg, cfg.dg, doZscore);
    Dda = harmonic_vertex_data(T, cfg, cfg.da, doZscore);
    R.nVertex = size(Ddg.Y, 1);

    banner(sprintf('PER-VERTEX HARMONIC MODEL  (%s betas, %d V1 vertices)', variant, R.nVertex));
    fprintf(['Convention: theta = orientation of the BARS, conventional visual-field deg\n' ...
             '(0 = horizontal, 90 = vertical). +b1 horizontal>vertical, +b2 cardinal>oblique,\n' ...
             '+b3 radial>tangential, +b4 polar-cardinal>polar-oblique. Reported as 2*b, the\n' ...
             'pro-minus-con difference, matching the manuscript''s deltas.\n' ...
             'Vertex weighting: EQUAL POLAR-ANGLE COVERAGE is primary (24 bins of 15 deg,\n' ...
             'w = 1/count, mean 1), matching the equal weighting of the eight wedges in the\n' ...
             'published ROI analysis and making the four predictors orthogonal. Natural\n' ...
             'vertex density is reported alongside in section 2. See HARMONIC_WEIGHTS.\n']);

    R.validate = section_validate(T, cfg, doZscore);
    R.fits     = section_fitAB(Ddg, Dda, cfg);
    R.cross    = section_cross(Ddg, Dda, cfg);
    R.roi      = section_roi(Ddg, Dda, cfg, R.fits);
    R.spec     = section_spec(Ddg, Dda, cfg);
    R.sens     = section_sensitivity(T, cfg, doZscore);

    write_csv(R, cfg, variant);
    R.figs = plot_harmonic(R, Ddg, Dda, cfg, variant);

    outMat = fullfile(cfg.cacheDir, sprintf('harmonic_%s.mat', variant));
    if ~isfolder(cfg.cacheDir), mkdir(cfg.cacheDir); end
    save(outMat, 'R', '-v7.3');
    fprintf('\nrun_harmonic_model: wrote %s\n', outMat);
end

% ========================================================================
function V = section_validate(T, cfg, doZscore)
% The per-vertex machinery must reproduce the published ROI pipeline before any of
% its departures mean anything.
%
% REFERENCE COLUMN. Z-scoring has been dropped (local_qc/REPORT.md section 4; the
% revised manuscript, confirmed 2026-08-17, reports percent signal change). The four
% constants per experiment transcribed in VALIDATE_AGAINST_MANUSCRIPT are reproduced
% to 3 decimals by the Z-SCORED pipeline and not by the raw one, so as transcribed
% they are in sigma units and are SUPERSEDED; they are kept only as the reference for
% the z-scored sensitivity run.
%
% The live reference is the raw one, in PERCENT SIGNAL CHANGE, computed here from the
% published route (bin_and_aggregate + compute_asymmetries on blank-subtracted
% contrasts) rather than transcribed -- these are the values the revised manuscript
% should carry.
    banner('1. VALIDATION -- does the per-vertex machinery reproduce the ROI pipeline?');
    if doZscore
        ref.dg = [-1.155 -0.40 0.23 0.06];      % SUPERSEDED sigma-unit manuscript values
        ref.da = [-0.45  -0.06 0.60 0.17];
        refLbl = 'old(sigma)';
    else
        ref.dg = [-0.4798 -0.2044 0.1159 0.0292];   % published route, % signal change
        ref.da = [-0.2114 -0.0296 0.1501 0.0336];
        refLbl = 'ref(%BOLD)';
    end
    nmA = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    V = struct();
    for e = {'dg','da'}
        en = e{1}; expCfg = cfg.(en);
        D = harmonic_vertex_data(T, cfg, expCfg, doZscore);

        Mpub   = bin_and_aggregate(T, cfg, expCfg, doZscore);
        Apub   = compute_asymmetries(Mpub, cfg, expCfg);
        indPub = cellfun(@(a) mean(mean(Apub.(a).diff,1)), Apub.order);

        [Amed, ~]     = harmonic_roi_roundtrip(D.Y, D, cfg, expCfg, 'median');
        indMed        = cellfun(@(a) mean(mean(Amed.(a).diff,1)), Amed.order);
        [Amean, Mmn]  = harmonic_roi_roundtrip(D.Y, D, cfg, expCfg, 'mean');
        indMean       = cellfun(@(a) mean(mean(Amean.(a).diff,1)), Amean.order);
        rl            = fit_lme_fig7(Mmn, cfg, expCfg, [], false);
        rA            = fit_harmonic_vertex(D, cfg, expCfg, 'binned');

        fprintf('\n--- %s ---\n%-11s %11s %8s %8s %8s %8s %8s\n', en, ...
                'asymmetry',refLbl,'pubMed','rtMed','rtMean','LMEmean','FitA*2');
        for j = 1:4
            fprintf('%-11s % 11.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f\n', nmA{j}, ...
                    ref.(en)(j), indPub(j), indMed(j), indMean(j), rl.delta(j), 2*rA.bMean(j));
        end
        V.(en) = struct('ref',ref.(en),'refLabel',refLbl,'pubMed',indPub,'rtMed',indMed, ...
                        'rtMean',indMean,'lmeMean',rl.delta,'fitA',2*rA.bMean);
    end
    fprintf(['\npubMed = published route (wedge medians of blank-subtracted contrasts).\n' ...
             'rtMed/rtMean = the same wedges applied to the per-vertex DEMEANED data. The four\n' ...
             'asymmetries are zero-sum contrasts, so demeaning cancels exactly under the mean\n' ...
             'but not under the median -- hence rtMean == LMEmean exactly, while rtMed differs\n' ...
             'slightly from pubMed. FitA*2 is the per-vertex fit with thetaV quantised to the\n' ...
             'wedge centre. Under the primary equal-polar-angle-coverage weighting its weights\n' ...
             'are binned from those same wedge centres, so every wedge carries equal total\n' ...
             'weight -- which is exactly what the ROI route does. FitA*2 therefore now equals\n' ...
             'LMEmean to the printed precision, where the old natural-density fit differed by\n' ...
             'up to 0.03 (it weighted observers'' wedges by vertex count).\n']);
end

% ========================================================================
function F = section_fitAB(Ddg, Dda, cfg)
% Fit A (wedge-centre thetaV) vs Fit B (true per-vertex thetaV). Same vertices, same
% weighting, ONLY thetaV differs -- so A minus B is the within-wedge local-orientation
% artifact, isolated.
%
% Reported under BOTH vertex weightings, so the two are like-for-like:
%   PRIMARY   'equalcoverage' -- every 15-deg polar-angle bin carries equal total
%             weight. This is what the published ROI analysis does (it weights the
%             eight wedges equally), and it makes the four predictors orthogonal.
%   SECONDARY 'natural'       -- one vertex, one vote. V1 over-represents the
%             horizontal meridian, which correlates the b1 and b3 columns at ~+0.35.
% See HARMONIC_WEIGHTS and ../HARMONIC_MODEL.md.
%
% weightSource is pinned to 'continuous' for BOTH fits here, so the weights are held
% fixed and only thetaV moves. Without the pin Fit A would bin its weights from the
% wedge centres and Fit B from the true angles, and roughly a third of the dg
% horiz-vert A->B shift would be that change of weighting rather than geometry. (The
% unpinned Fit A is still the right one for section 1, where the point is precisely
% that equal weight per wedge reproduces the published ROI analysis.)
    banner('2. FIT A vs FIT B -- how much is within-ROI local-orientation geometry?');
    nm = {'horiz-vert','card-obl','rad-tang','polc-polo'};

    modes  = {'equalcoverage', 'natural'};
    suffix = {'', 'nat'};
    lbl    = {'PRIMARY -- equal polar-angle coverage', ...
              'SECONDARY -- natural vertex density (the old behaviour)'};

    F = struct('gap', struct());
    for wi = 1:2
        o  = struct('weighting', modes{wi}, 'weightSource', 'continuous');
        sx = suffix{wi};
        fprintf('\n\n***** %s *****\n', lbl{wi});
        check_orthogonality(Ddg, Dda, cfg, modes{wi});

        for e = {'dg','da'}
            en = e{1};
            D = Ddg; if strcmp(en,'da'), D = Dda; end
            rA = fit_harmonic_vertex(D, cfg, cfg.(en), 'binned',     o);
            rB = fit_harmonic_vertex(D, cfg, cfg.(en), 'continuous', o);
            F.(en).(['A' sx]) = rA;  F.(en).(['B' sx]) = rB;

            fprintf('\n--- %s ---   (2*b, with %d%% bootstrap CI over %d observers;\n', ...
                    en, cfg.ciLevel, numel(cfg.subjects));
            fprintf('               weights held fixed across A and B, binned from the true pRF angle)\n');
            fprintf('%-11s %22s %22s %9s\n', 'asymmetry', 'Fit A (wedge centre)', ...
                    'Fit B (true pRF angle)', 'B - A');
            for j = 1:4
                fprintf('%-11s  % 7.3f [% 6.3f % 6.3f]  % 7.3f [% 6.3f % 6.3f]  % 8.3f\n', nm{j}, ...
                        2*rA.bMean(j), 2*rA.ci(j,1), 2*rA.ci(j,2), ...
                        2*rB.bMean(j), 2*rB.ci(j,1), 2*rB.ci(j,2), ...
                        2*(rB.bMean(j)-rA.bMean(j)));
            end
            fprintf('   R2  A %.4f  B %.4f     cond(X''X) A %.2f  B %.2f   maxVIF B %.2f\n', ...
                    mean(rA.r2), mean(rB.r2), mean(rA.condNum), mean(rB.condNum), ...
                    max(mean(rB.vif,1,'omitnan')));
        end

        F.gap.(modes{wi}) = report_gap(F, sx, nm);
    end

    fprintf(['\nThe PRIMARY (equal-coverage) percentages are the headline numbers. They are\n' ...
             'not interchangeable with the natural-density ones: at natural density the b1\n' ...
             'and b3 columns are correlated ~+0.35, so part of what Fit B moves is the\n' ...
             'collinearity re-splitting rather than the geometry correction.\n']);
end

% ------------------------------------------------------------------------
function G = report_gap(F, sx, nm)
% The quantity the manuscript's claim rests on: the Cartesian-vs-polar GAP in each
% asymmetry, and how much of it Fit B's exact per-vertex local orientation removes.
    G = struct('gapA', nan(1,4), 'gapB', nan(1,4), 'explained', nan(1,4));
    for j = 1:4
        G.gapA(j) = 2*(F.dg.(['A' sx]).bMean(j) - F.da.(['A' sx]).bMean(j));
        G.gapB(j) = 2*(F.dg.(['B' sx]).bMean(j) - F.da.(['B' sx]).bMean(j));
        if abs(G.gapA(j)) > 0
            G.explained(j) = 100*(1 - abs(G.gapB(j))/abs(G.gapA(j)));
        end
    end
    fprintf('\nCartesian-vs-polar gap (2*b_dg - 2*b_da) and what geometry explains:\n');
    fprintf('%-11s %10s %10s %22s\n', 'asymmetry', 'Fit A gap', 'Fit B gap', ...
            'explained by geometry');
    for j = 1:4
        fprintf('%-11s % 10.3f % 10.3f %20.1f%%\n', nm{j}, G.gapA(j), G.gapB(j), G.explained(j));
    end
end

% ------------------------------------------------------------------------
function check_orthogonality(Ddg, Dda, cfg, weighting)
% ASSERTED, not eyeballed. Two things must hold for the weights to be trustworthy.
%
% (1) The uncentered correlation between the b1 column, cos(2*theta), and the b3
%     column, cos(2*(theta-thetaV)), is ANALYTICALLY the weighted mean of
%     cos(2*thetaV):  r(b1,b3) = sum(w.*cosd(2*tv))/sum(w).  Checking the design
%     diagnostic against that closed form is a strong test that the same weights
%     reached both the design and the summary -- a sqrt(w)-vs-w slip, or weights
%     applied to the fit but not the diagnostics, breaks it immediately.
% (2) The four cross-harmonic blocks (b1-b2, b1-b4, b2-b3, b3-b4) pair a 2nd with a
%     4th harmonic and are zero at EVERY vertex individually, hence zero under any
%     weighting. If one of these drifts, the sampling angles are wrong, not the weights.
    cross = [1 2; 1 4; 2 3; 3 4];
    for e = {'dg','da'}
        en = e{1};
        D = Ddg; if strcmp(en,'da'), D = Dda; end
        r = fit_harmonic_vertex(D, cfg, cfg.(en), 'continuous', ...
                                struct('weighting', weighting, 'nBoot', 1));
        got = r.corr(1,3);
        want = r.wMeanCos2(find(isfinite(r.wMeanCos2), 1, 'last'));
        assert(abs(got - want) < 1e-12, 'run_harmonic_model:orthIdentity', ...
               '%s/%s: corr(b1,b3)=%.15g but weighted mean cos(2 thetaV)=%.15g.', ...
               en, weighting, got, want);
        for k = 1:size(cross,1)
            v = r.corr(cross(k,1), cross(k,2));
            assert(abs(v) < 1e-12, 'run_harmonic_model:orthCross', ...
                   '%s/%s: corr(b%d,b%d) = %.3e, expected 0.', ...
                   en, weighting, cross(k,1), cross(k,2), v);
        end
        fprintf(['   %s design (%s): corr(b1,b3) = % .4f == weighted mean cos(2*thetaV)\n' ...
                 '        for the last observer; across observers % .4f. b1-b2, b1-b4, b2-b3,\n' ...
                 '        b3-b4 all 0 to 1e-12. cond(X''WX) %.2f, effective n %.0f%% of vertices. PASS\n'], ...
                en, weighting, got, mean(r.wMeanCos2, 'omitnan'), mean(r.condNum), ...
                100*mean(r.effN./r.n, 'omitnan'));
    end
end

% ========================================================================
function X = section_cross(Ddg, Dda, cfg)
    banner('3. CROSS-EXPERIMENT PREDICTION AND THE NESTED CONTEXT TEST');
    X = harmonic_crossexp(Ddg, Dda, cfg, 'continuous', struct('weighting','equalcoverage'));
    fprintf('(equal polar-angle coverage weighting; Fit B geometry)\n');
    f = @(s) sprintf('% .4f [% .4f % .4f]', mean(X.r2.(s)), X.ci.(s)(1), X.ci.(s)(2));

    fprintf('\nVariance explained in the held-out experiment (%d%% CI):\n', cfg.ciLevel);
    fprintf('   dg fitted on itself  (ceiling)   %s\n', f('dgSelf'));
    fprintf('   da fitted on itself  (ceiling)   %s\n', f('daSelf'));
    fprintf('   dg -> da                         %s\n', f('dg2da'));
    fprintf('   dg -> da, free gain              %s   gain %.3f\n', f('dg2daG'), mean(X.gain.dg2da));
    fprintf('   da -> dg                         %s\n', f('da2dg'));
    fprintf('   da -> dg, free gain              %s   gain %.3f\n', f('da2dgG'), mean(X.gain.da2dg));

    fprintf('\nNested models on the concatenated data (common denominator):\n');
    fprintf('   (i)   shared b1..b4              %s\n', f('shared'));
    fprintf('   (ii)  shared b1..b4 + da gain    %s   gain %.3f\n', f('sharedGain'), mean(X.gain.joint));
    fprintf('   (iii) separate b1..b4            %s\n', f('separate'));
    d1 = mean(X.r2.sharedGain) - mean(X.r2.shared);
    d2 = mean(X.r2.separate)   - mean(X.r2.sharedGain);
    fprintf('   (i)->(ii)  overall gain          +%.4f  (%.0f%% of the shared->separate gap)\n', ...
            d1, 100*d1/(d1+d2));
    fprintf('   (ii)->(iii) reference-frame      +%.4f  (%.0f%%)  <-- the context effect\n', ...
            d2, 100*d2/(d1+d2));

    nm = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    fprintf('\nCoefficient differences b_dg - b_da (* = %d%% CI excludes 0):\n', cfg.ciLevel);
    fprintf('%-11s %26s %26s\n', 'asymmetry', 'raw', 'gain-equalised');
    for j = 1:4
        s1 = star(X.ci.dBeta(j,:));  s2 = star(X.ci.dBetaGain(j,:));
        fprintf('%-11s  % 7.4f [% 7.4f % 7.4f]%s  % 7.4f [% 7.4f % 7.4f]%s\n', nm{j}, ...
                mean(X.dBeta(:,j)),     X.ci.dBeta(j,1),     X.ci.dBeta(j,2),     s1, ...
                mean(X.dBetaGain(:,j)), X.ci.dBetaGain(j,1), X.ci.dBetaGain(j,2), s2);
    end
end

% ========================================================================
function O = section_roi(Ddg, Dda, cfg, F)
% Push observed and cross-predicted per-vertex responses through the published ROI
% pipeline, so the model's claim lands in the manuscript's own units.
    banner('4. ROI ROUND-TRIP -- observed vs cross-predicted asymmetries');
    nm = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    O = struct();
    pairs = {'da','dg'; 'dg','da'};      % {target, source}
    for k = 1:2
        tgt = pairs{k,1};  src = pairs{k,2};
        Dt = Ddg; if strcmp(tgt,'da'), Dt = Dda; end
        bSrc = F.(src).B.bMean;

        Yobs  = Dt.Y;
        Ypred = predict_harmonic(bSrc, Dt.tvCont, cfg.(tgt));
        g     = (Ypred(:).' * Yobs(:)) / sum(Ypred(:).^2);

        Aobs  = harmonic_roi_roundtrip(Yobs,    Dt, cfg, cfg.(tgt), 'mean');
        Apred = harmonic_roi_roundtrip(Ypred,   Dt, cfg, cfg.(tgt), 'mean');
        Apg   = harmonic_roi_roundtrip(g*Ypred, Dt, cfg, cfg.(tgt), 'mean');
        gv = @(A,a) mean(mean(A.(a).diff,1));

        fprintf('\n--- %s observed vs predicted from %s coefficients (gain %.3f) ---\n', tgt, src, g);
        fprintf('%-11s %10s %12s %12s\n', 'asymmetry', 'observed', 'predicted', 'pred*gain');
        for j = 1:4
            a = Aobs.order{j};
            fprintf('%-11s % 10.3f % 12.3f % 12.3f\n', nm{j}, gv(Aobs,a), gv(Apred,a), gv(Apg,a));
        end
        O.(tgt) = struct('obs',Aobs,'pred',Apred,'predGain',Apg,'gain',g,'src',src);
    end
    fprintf(['\nThese are wedge MEANS of demeaned responses, so they are directly comparable\n' ...
             'between observed and predicted but sit slightly below the published medians.\n']);
end

% ========================================================================
function S = section_spec(Ddg, Dda, cfg)
% Is the 4-term model adequate? The complete harmonic basis at harmonics 2 and 4 adds
% four sin columns; they should vanish under left-right visual-field symmetry.
    banner('5. SPECIFICATION TEST -- expanded harmonic basis');
    o.expanded = true;
    S = struct();
    for e = {'dg','da'}
        en = e{1};
        D = Ddg; if strcmp(en,'da'), D = Dda; end
        r  = fit_harmonic_vertex(D, cfg, cfg.(en), 'continuous', o);
        r4 = fit_harmonic_vertex(D, cfg, cfg.(en), 'continuous');
        S.(en) = r;
        fprintf('\n--- %s (identically zero, dropped: %s) ---\n', en, strjoin(r.names(r.dropped),','));
        for j = 1:numel(r.names)
            if r.dropped(j), continue; end
            core = ''; if j <= 4, core = sprintf('   4-term: % .4f', r4.bMean(j)); end
            fprintf('   %-8s b=% .4f  CI[% .4f % .4f]%s%s\n', r.names{j}, r.bMean(j), ...
                    r.ci(j,1), r.ci(j,2), star(r.ci(j,:)), core);
        end
        fprintf('   R2 expanded %.4f vs 4-term %.4f\n', mean(r.r2), mean(r4.r2));
    end
end

% ========================================================================
function S = section_sensitivity(T, cfg, doZscore)
% pRF polar-angle error attenuates the thetaV-dependent terms (regression dilution)
% and, for da, leaks b1 into b3. Raising the pRF R2 floor bounds it.
    banner('6. SENSITIVITY -- pRF R2 floor (bounds the regression-dilution caveat)');
    fprintf('(Fit B geometry, equal polar-angle coverage weighting)\n');
    thr = [0.1 0.2 0.3 0.5];
    nm  = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    S = struct('thr', thr);
    for e = {'dg','da'}
        en = e{1};
        fprintf('\n--- %s (2*b) ---\n%7s %8s %11s %11s %11s %11s\n', en, ...
                'r2min','nVertex', nm{:});
        B = nan(numel(thr), 4);
        for i = 1:numel(thr)
            c2 = cfg; c2.r2min = thr(i);
            D2 = harmonic_vertex_data(T, c2, cfg.(en), doZscore);
            r  = fit_harmonic_vertex(D2, c2, cfg.(en), 'continuous');
            B(i,:) = 2*r.bMean;
            fprintf('%7.2f %8d % 11.4f % 11.4f % 11.4f % 11.4f\n', thr(i), size(D2.Y,1), B(i,:));
        end
        S.(en) = B;
    end

    % does the cross-experiment difference survive a stricter floor?
    c2 = cfg; c2.r2min = 0.3;
    Dg = harmonic_vertex_data(T, c2, cfg.dg, doZscore);
    Da = harmonic_vertex_data(T, c2, cfg.da, doZscore);
    X  = harmonic_crossexp(Dg, Da, c2, 'continuous');
    S.cross03 = X;
    fprintf('\nCross-experiment differences b_dg - b_da at pRF R2 > 0.3:\n');
    for j = 1:4
        fprintf('   %-11s raw % 7.4f [% 7.4f % 7.4f]%s   gain-eq % 7.4f [% 7.4f % 7.4f]%s\n', ...
                nm{j}, mean(X.dBeta(:,j)), X.ci.dBeta(j,1), X.ci.dBeta(j,2), star(X.ci.dBeta(j,:)), ...
                mean(X.dBetaGain(:,j)), X.ci.dBetaGain(j,1), X.ci.dBetaGain(j,2), star(X.ci.dBetaGain(j,:)));
    end
end

% ========================================================================
function write_csv(R, cfg, variant)
    if ~isfolder(cfg.figDir), mkdir(cfg.figDir); end
    nm = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    flds = {'A','B','Anat','Bnat'};
    fitL = {'A','B','A','B'};
    wL   = {'equalcoverage','equalcoverage','natural','natural'};
    rows = {};
    for e = {'dg','da'}
        en = e{1};
        for k = 1:numel(flds)
            r = R.fits.(en).(flds{k});
            for j = 1:4
                rows(end+1,:) = {variant, wL{k}, en, fitL{k}, nm{j}, r.bMean(j), ...
                                 2*r.bMean(j), 2*r.ci(j,1), 2*r.ci(j,2), mean(r.r2)}; %#ok<AGROW>
            end
        end
    end
    Tc = cell2table(rows, 'VariableNames', {'variant','weighting','experiment','fit', ...
                    'asymmetry','beta','delta','ciLo','ciHi','r2'});
    outCsv = fullfile(cfg.figDir, sprintf('harmonic_coefficients_%s.csv', variant));
    writetable(Tc, outCsv);
    fprintf('\nrun_harmonic_model: wrote %s\n', outCsv);
end

% ========================================================================
function s = star(ci)
    s = '';
    if all(isfinite(ci)) && ~(ci(1) <= 0 && 0 <= ci(2)), s = ' *'; end
end

function banner(txt)
    fprintf('\n\n%s\n%s\n%s\n', repmat('=',1,78), txt, repmat('=',1,78));
end
