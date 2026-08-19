function test_harmonic_model()
% TEST_HARMONIC_MODEL  Assertions guarding the per-vertex harmonic model.
%
%   test_harmonic_model()
%
% Runs before any interpretation. Sections 1-3 are pure geometry (no data); section 4
% is a synthetic-recovery check against the real polar-angle distribution.
%
% Section 1 is the important one: it asserts that the four harmonic predictors reduce
% EXACTLY to the repo's +1/0/-1 asymmetry codes at the eight wedge centres, for both
% experiments. Sign/frame slips of precisely this kind produced both of the retracted
% bugs recorded in ../FINDINGS.md and ../AUDIT.md, so this is an assert, not an
% eyeball check.

    cfg = config_repro();
    fprintf('\n================ test_harmonic_model ================\n');

    section1_codes(cfg);
    section2_identities(cfg);
    section3_decomposition(cfg);
    section4_recovery(cfg);

    fprintf('\nAll harmonic-model checks PASSED.\n');
end

% ------------------------------------------------------------------------
function section1_codes(cfg)
% The four predictors, at the wedge centres, must equal LME_CODES exactly.
    fprintf('\n-- 1. predictor codes vs lme_codes at the 8 wedge centres --\n');
    for e = {'dg','da'}
        en = e{1}; expCfg = cfg.(en);
        X  = lme_codes(cfg, expCfg);                 % 4 fields, each nOri x nPA
        nO = numel(expCfg.oriAngle);
        nP = numel(cfg.paBins);
        got = zeros(nO, nP, 4);
        for p = 1:nP
            Xp = harmonic_predictors(cfg.paBins(p), expCfg);   % (1*nOri) x 4
            got(:, p, :) = reshape(Xp, [nO, 1, 4]);
        end
        ref = cat(3, X.hVv, X.cVo, X.rVt, X.pcVpo);
        err = max(abs(got(:) - ref(:)));
        fprintf('   %s: max |harmonic - lme_codes| = %.3e\n', en, err);
        assert(err < 1e-12, 'test_harmonic_model:codes', ...
               '%s predictors do not match lme_codes (err %.3e).', en, err);
    end
    fprintf('   -> cos(2t)=hVv, cos(4t)=cVo, cos(2(t-tV))=rVt, cos(4(t-tV))=pcVpo. PASS\n');
end

% ------------------------------------------------------------------------
function section2_identities(cfg)
% Analytic properties the construction must have.
    fprintf('\n-- 2. analytic identities --\n');
    tv = (0:1:359).';                       % dense sweep of polar angles
    opts.expanded = true;

    for e = {'dg','da'}
        en = e{1}; expCfg = cfg.(en);
        [X, info] = harmonic_predictors(tv, expCfg, opts);
        nV = info.nVertex; nO = info.nOri;
        R  = @(j) reshape(X(:,j), [nV nO]);

        % (a) every predictor sums to zero across the four conditions, so the model
        %     needs no intercept and per-vertex demeaning does not distort it.
        for j = 1:4
            s = max(abs(sum(R(j), 2)));
            assert(s < 1e-12, 'test_harmonic_model:sum0', ...
                   '%s predictor %d does not sum to zero across conditions (%.3e).', en, j, s);
        end
        fprintf('   %s: all 4 predictors sum to zero across conditions. PASS\n', en);

        % (b) the 4th-harmonic collinearity. The two 4th-harmonic columns are
        %     proportional AT EACH VERTEX, with ratio cos(4 tV); which one is the
        %     vertex-independent factor swaps between experiments.
        if expCfg.isPolar
            lhs = R(2); rhs = cosd(4*tv) .* R(4);   % cos(4 t) = cos(4 tV)*cos(4(t-tV))
            lbl = 'cos(4 t) == cos(4 tV)*cos(4(t-tV))';
        else
            lhs = R(4); rhs = cosd(4*tv) .* R(2);   % cos(4(t-tV)) = cos(4 tV)*cos(4 t)
            lbl = 'cos(4(t-tV)) == cos(4 tV)*cos(4 t)';
        end
        assert(max(abs(lhs(:)-rhs(:))) < 1e-12, 'test_harmonic_model:h4', ...
               '%s: 4th-harmonic proportionality failed.', en);
        fprintf('   %s: %s  (b2/b4 near-collinear). PASS\n', en, lbl);

        % (c) the vertex-independent term, and the identically-zero sin column
        spread = @(M) max(max(M,[],1) - min(M,[],1));
        if expCfg.isPolar
            v = R(3);                              % cos(2(t-tV)) for da
            assert(spread(v) < 1e-12, 'test_harmonic_model:daConst', ...
                   'da: cos(2(t-tV)) varies across vertices.');
            fprintf('   da: cos(2(t-tV)) is vertex-independent -> under the model, da''s\n');
            fprintf('       radial-tangential asymmetry is predicted to EQUAL dg''s. PASS\n');
            zc = 8;  zname = 'sin(4(t-tV))';
        else
            v = R(1);                              % cos(2 t) for dg
            assert(spread(v) < 1e-12, 'test_harmonic_model:dgConst', ...
                   'dg: cos(2 t) varies across vertices.');
            fprintf('   dg: cos(2 t) is vertex-independent. PASS\n');
            zc = 6;  zname = 'sin(4 t)';
        end
        assert(info.isZeroCol(zc), 'test_harmonic_model:zerocol', ...
               '%s: expected %s to be identically zero.', en, zname);
        fprintf('   %s: %s is identically zero -> only 3 of 4 harmonic dof are\n', en, zname);
        fprintf('       observable per vertex. PASS\n');
    end
end

% ------------------------------------------------------------------------
function section3_decomposition(cfg)
% The A/B/C decomposition must reconstruct the data exactly, and must agree with the
% model curves when the data ARE the model.
    fprintf('\n-- 3. A/B/C decomposition --\n');
    rng(0);
    tv = 360*rand(500,1);
    b  = [0.7 -0.3 0.45 0.12];
    for e = {'dg','da'}
        en = e{1}; expCfg = cfg.(en);

        Y = predict_harmonic(b, tv, expCfg);          % noiseless model data
        H = harmonic_decompose(Y, tv, expCfg, b);

        assert(H.reconErr < 1e-12, 'test_harmonic_model:recon', ...
               '%s: 3-term reconstruction error %.3e.', en, H.reconErr);
        eA = max(abs(H.A - H.Ahat));
        eB = max(abs(H.B - H.Bhat));
        eC = max(abs(H.C - H.Chat));
        fprintf('   %s (%s frame): recon %.2e | A %.2e  B %.2e  C %.2e\n', ...
                en, H.frame, H.reconErr, eA, eB, eC);
        assert(max([eA eB eC]) < 1e-12, 'test_harmonic_model:curves', ...
               '%s: A/B/C model curves disagree with the decomposition.', en);
    end
    fprintf('   -> the model is exactly three scalar regressions on thetaV. PASS\n');
end

% ------------------------------------------------------------------------
function section4_recovery(cfg)
% Recover known coefficients from synthetic data at the REAL polar-angle
% distribution, and measure how much pRF angle error attenuates b3/b4.
    fprintf('\n-- 4. synthetic recovery at the real thetaV distribution --\n');
    T  = load_allconditions(cfg);
    bTrue = [-0.55 -0.20 0.12 0.03];

    for e = {'dg','da'}
        en = e{1}; expCfg = cfg.(en);
        D  = harmonic_vertex_data(T, cfg, expCfg, false);

        % noise scaled to the real per-vertex spread, so the SNR is realistic
        sd = std(D.Y(:));
        rng(1);
        Dsim = D;
        Dsim.Y = predict_harmonic(bTrue, D.tvCont, expCfg) + sd*randn(size(D.Y));

        r = fit_harmonic_vertex(Dsim, cfg, expCfg, 'continuous');
        fprintf('   %s clean : true [% .3f % .3f % .3f % .3f]\n', en, bTrue);
        fprintf('   %s        : fit  [% .3f % .3f % .3f % .3f]\n', en, r.bMean);
        assert(max(abs(r.bMean - bTrue)) < 0.02, 'test_harmonic_model:recovery', ...
               '%s: coefficients not recovered (max err %.3f).', en, max(abs(r.bMean-bTrue)));

        % pRF polar-angle error attenuates the thetaV-dependent terms toward zero.
        for jit = [5 10 20]
            Djit = Dsim;
            Djit.tvCont = mod(D.tvCont + jit*randn(size(D.tvCont)), 360);
            rj = fit_harmonic_vertex(Djit, cfg, expCfg, 'continuous');
            fprintf('   %s jitter %2d deg: [% .3f % .3f % .3f % .3f]\n', en, jit, rj.bMean);
        end
    end
    fprintf('   -> recovery exact without jitter; the jitter rows quantify the\n');
    fprintf('      regression-dilution caveat on the thetaV-dependent terms.\n');
end
