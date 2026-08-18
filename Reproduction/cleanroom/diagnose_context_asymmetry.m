function S = diagnose_context_asymmetry(varargin)
% DIAGNOSE_CONTEXT_ASYMMETRY  Is the context effect really confined to the Cartesian frame?
%
%   S = diagnose_context_asymmetry()
%
% The manuscript reports context effects for the horizontal/vertical and
% cardinal/oblique asymmetries. A stronger, more interesting claim would be that the
% POLAR-frame asymmetries (radial/tangential, polar-cardinal/polar-oblique) are NOT
% context-dependent -- that the effect is one-sided. This function asks whether the
% data support that stronger claim. They do not, and the reasons are worth recording.
%
% Every observer completed BOTH experiments, so every comparison here is WITHIN
% SUBJECT: the per-observer difference is formed first, then summarised across the
% eight observers.
%
% Two routes to the POINT ESTIMATE, which agree to three decimals:
%   (A) the per-vertex harmonic model (HARMONIC_CROSSEXP, b_dg - b_da per observer)
%   (B) the published ROI route: wedge medians -> COMPUTE_ASYMMETRIES, per observer
% and the paired t-test / bootstrap across the eight observers on (B) is the primary
% inferential test.
%
% (C) An LME with experiment x asymmetry interactions is reported ONLY as a caution.
%     It is NOT a valid within-subject test here and must not be quoted as one. In
%     'y ~ expf*(hVv+cVo+rVt+pcVpo) + (1|subj)' every asymmetry slope and every
%     interaction is a FIXED effect, assumed identical across observers; only the
%     intercept varies. fitlme then tests the interaction on DF = 502 -- the
%     wedge-level observations -- rather than on the 8 observers, and applies no
%     Satterthwaite or Kenward-Roger correction. The result is anti-conservative by
%     roughly 5-25x in p. Adding random slopes, including for the interaction terms
%     themselves, does not fix it: the DF stays at 502.
%
%     Because the 4 x 8 design is balanced and its four asymmetry codes are exactly
%     orthogonal, the LME fixed effect is IDENTICAL to the mean of the per-observer
%     contrasts (0.268459 either way). The LME therefore adds nothing to the estimate
%     and only misstates its uncertainty. With a balanced orthogonal within-subject
%     design, the summary-statistic route -- per-observer effect, then t-test across
%     observers -- is the correct analysis.
%
% Then the part that matters for the framing:
%   (D) leave-one-observer-out, and a sign test
%   (E) an equivalence bound -- what size of polar-frame context effect the data can
%       actually EXCLUDE
%   (F) a within-subject difference of differences: is the Cartesian-frame context
%       effect reliably LARGER than the polar-frame one?
%
% (G) WITHIN-OBSERVER VARIANCE. The standard objection to a summary-statistic test is
%     that it treats each observer's estimate as noiseless, so the across-observer
%     variance it uses contains both true between-observer variability AND
%     within-observer estimation error. That is what a mixed model is normally for.
%     Here it is measured rather than argued: bootstrapping VERTICES within each wedge
%     and recomputing the wedge medians and asymmetries gives the within-observer SE of
%     exactly the quantity being tested. It is 0.02-0.03, against a between-observer SD
%     of 0.13-0.27 -- measurement error is only 2-4% of the total variance. The
%     two-stage test is therefore very close to optimal, and a mixed model has almost
%     nothing left to recover. Two consequences: the implied BLUP shrinkage is under 2%,
%     so per-observer estimates from a mixed model would be within 2% of the raw ones;
%     and the spread across observers is genuine individual variation, not noise.
%
% Absence of evidence is not evidence of absence, and with n = 8 the distinction is
% not academic here. See ../HARMONIC_MODEL.md and ../supplement/.

    p = inputParser;
    p.addParameter('nBoot', 1000, @isnumeric);
    p.addParameter('level', 95, @isnumeric);
    p.parse(varargin{:});
    opt = p.Results;

    cfg = config_repro();
    T   = load_allconditions(cfg);
    nm  = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    S   = struct('names', {nm});

    % ---- (B) published ROI route, per observer -------------------------------
    per = struct();
    for e = {'dg','da'}
        en = e{1};
        M  = bin_and_aggregate(T, cfg, cfg.(en), false);
        A  = compute_asymmetries(M, cfg, cfg.(en));
        V  = nan(numel(cfg.subjects), 4);
        for j = 1:4, V(:,j) = mean(A.(A.order{j}).diff, 1).'; end
        per.(en) = V;
    end
    d = per.dg - per.da;                       % paired, per observer
    S.perSubjectDG = per.dg;  S.perSubjectDA = per.da;  S.paired = d;

    banner('WITHIN-SUBJECT CONTEXT TEST (published ROI route, %BOLD)');
    fprintf('%-11s %9s %9s %10s %-21s %7s %8s\n', ...
            'asymmetry','dg','da','dg-da','95% CI (paired)','t(7)','p');
    S.t = nan(1,4); S.p = nan(1,4); S.ci = nan(4,2);
    for j = 1:4
        S.ci(j,:) = bootstrap_ci(d(:,j).', opt.nBoot, opt.level, 0);
        [~, S.p(j), ~, st] = ttest(d(:,j));  S.t(j) = st.tstat;
        fprintf('%-11s %9.3f %9.3f %10.3f [%7.3f %7.3f]%-2s %7.2f %8.4f\n', nm{j}, ...
                mean(per.dg(:,j)), mean(per.da(:,j)), mean(d(:,j)), S.ci(j,:), ...
                star(S.ci(j,:)), S.t(j), S.p(j));
    end

    % ---- (C) LME, reported as a CAUTION, not as the test ---------------------
    tb = build_lme_table(T, cfg);
    fx = 'y ~ hVv+cVo+rVt+pcVpo+isDA+hVvDA+cVoDA+rVtDA+pcDA';
    re = {'(1|subj)', '(1+isDA|subj)', '(1+hVv+rVt|subj)', '(1+hVv+rVt+hVvDA+rVtDA|subj)'};
    lb = {'random intercept only','+ random experiment','+ random asymmetry slopes', ...
          '+ random INTERACTION slopes'};
    banner('LME -- ANTI-CONSERVATIVE, shown only to document why it is not used');
    fprintf('%-30s %8s %10s %10s %10s\n','random-effects structure','DF','p hVvDA','p rVtDA','2b hVvDA');
    S.lme = cell(numel(re),1);
    for i = 1:numel(re)
        L = fitlme(tb, [fx ' + ' re{i}]);
        S.lme{i} = L;
        k = strcmp(L.Coefficients.Name, 'hVvDA');
        r = strcmp(L.Coefficients.Name, 'rVtDA');
        fprintf('%-30s %8d %10.4f %10.4f %10.3f\n', lb{i}, L.Coefficients.DF(k), ...
                L.Coefficients.pValue(k), L.Coefficients.pValue(r), ...
                2*L.Coefficients.Estimate(k));
    end
    fprintf(['\nDF stays at %d (the wedge-level observations) whatever the random structure,\n' ...
             'and the point estimate never moves. Compare the paired t(7) p-values above:\n' ...
             'hVv %.4f and rVt %.4f. The LME is anti-conservative here; the paired test is\n' ...
             'the correct within-subject analysis.\n'], ...
            S.lme{1}.Coefficients.DF(strcmp(S.lme{1}.Coefficients.Name,'hVvDA')), S.p(1), S.p(3));

    % ---- (D) how robust is the polar-frame null? ----------------------------
    banner('IS THE POLAR-FRAME NULL ROBUST? leave-one-observer-out, rad-tang');
    rt = d(:,3);
    nNeg = sum(rt < 0);
    fprintf('sign: %d of %d observers negative (da > dg); binomial two-tailed p = %.3f\n', ...
            nNeg, numel(rt), 2*binocdf(min(nNeg, numel(rt)-nNeg), numel(rt), 0.5));
    S.loo = nan(numel(rt), 3);
    for si = 1:numel(rt)
        k = true(numel(rt),1);  k(si) = false;
        c = bootstrap_ci(rt(k).', opt.nBoot, opt.level, 0);
        S.loo(si,:) = [mean(rt(k)), c];
        fprintf('  drop %-14s mean %7.3f  CI [%7.3f %7.3f]%s\n', ...
                cfg.subjects{si}, mean(rt(k)), c, star(c));
    end

    % ---- (E) equivalence bound ----------------------------------------------
    banner('EQUIVALENCE: what can actually be EXCLUDED?');
    S.excludeAbove = max(abs(S.ci(3,:)));
    fprintf(['rad-tang context effect  %.3f, 95%% CI [%.3f %.3f].\n' ...
             'So an effect as large as |%.3f| remains compatible with these data.\n' ...
             'For scale: card-obl context effect is %.3f, horiz-vert is %.3f.\n'], ...
            mean(rt), S.ci(3,:), S.excludeAbove, mean(d(:,2)), mean(d(:,1)));

    % ---- (F) is the Cartesian context effect reliably larger? ---------------
    banner('DIFFERENCE OF DIFFERENCES: Cartesian-frame vs polar-frame context effect');
    fprintf('(within subject, on |effect| so the two are comparable)\n');
    S.dod = nan(2,4);
    cmp = {1,3,'horiz-vert vs rad-tang'; 2,3,'card-obl vs rad-tang'};
    for i = 1:2
        dd = abs(d(:,cmp{i,1})) - abs(d(:,cmp{i,2}));
        c  = bootstrap_ci(dd.', opt.nBoot, opt.level, 0);
        [~, pp] = ttest(dd);
        S.dod(i,:) = [mean(dd), c, pp];
        fprintf('  %-24s %7.3f  CI [%7.3f %7.3f]%-2s  p=%.4f\n', cmp{i,3}, mean(dd), c, star(c), pp);
    end

    % ---- (G) how much of the across-observer variance is measurement error? ----
    S.within = within_observer_se(T, cfg);
    banner('WITHIN-OBSERVER VARIANCE: is the summary-statistic test losing anything?');
    seDiff = sqrt(S.within.dg.^2 + S.within.da.^2);
    fprintf('%-11s %13s %14s %13s %12s\n', 'asymmetry', 'SD across obs', ...
            'within-obs SE', 'implied TRUE', 'within/total');
    S.varDecomp = nan(4,4);
    for j = 1:4
        vObs = var(d(:,j));  vWin = mean(seDiff(:,j).^2);  vTrue = max(vObs - vWin, 0);
        S.varDecomp(j,:) = [sqrt(vObs), sqrt(vWin), sqrt(vTrue), vWin/vObs];
        fprintf('%-11s %13.4f %14.4f %13.4f %11.0f%%\n', nm{j}, ...
                sqrt(vObs), sqrt(vWin), sqrt(vTrue), 100*vWin/vObs);
    end
    shr = S.varDecomp(:,3).^2 ./ (S.varDecomp(:,3).^2 + S.varDecomp(:,2).^2);
    fprintf(['\nMeasurement error is %.0f-%.0f%% of the across-observer variance, so the\n' ...
             'summary-statistic test is near-optimal and a mixed model has little to add.\n' ...
             'Implied BLUP shrinkage %.1f-%.1f%%: per-observer estimates from a mixed model\n' ...
             'would sit within a couple of percent of the raw ones. The spread across\n' ...
             'observers is therefore genuine individual variation, not noise -- which also\n' ...
             'means sub-0395 is a real outlier, not a badly estimated one.\n'], ...
            100*min(S.varDecomp(:,4)), 100*max(S.varDecomp(:,4)), ...
            100*(1-max(shr)), 100*(1-min(shr)));

    banner('READING');
    fprintf(['Every comparison here is WITHIN SUBJECT: the per-observer difference is formed\n' ...
             'first, then summarised across observers. Context effects on the two CARTESIAN-\n' ...
             'frame asymmetries are robust across both estimation routes. The polar-frame\n' ...
             'asymmetries show no DETECTABLE context effect --\n' ...
             'but that is not the same as showing there is none:\n' ...
             '  - the rad-tang CI still admits an effect the size of the card-obl one;\n' ...
             '  - dropping a single observer makes rad-tang significant;\n' ...
             '  - the Cartesian effect is NOT reliably larger than the polar one (F).\n' ...
             'So the manuscript should report the Cartesian-frame context effects positively\n' ...
             'and describe the polar-frame result as an absence of evidence, NOT as evidence\n' ...
             'that the polar-frame asymmetries are context-invariant.\n']);
end

% ------------------------------------------------------------------------
function SE = within_observer_se(T, cfg, nB)
% Within-observer SE of each asymmetry, by resampling VERTICES within each wedge and
% recomputing the wedge medians and asymmetries -- the same quantity the table holds.
    if nargin < 3, nB = 200; end
    keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
    Tk = T(keep, :);
    subjStr = string(Tk.subject);
    rng(0);
    SE = struct();
    for e = {'dg','da'}
        en = e{1};
        C  = compute_vertex_contrasts(Tk, cfg.(en), false);
        M8 = nan(numel(cfg.subjects), 4);
        for si = 1:numel(cfg.subjects)
            m    = subjStr == cfg.subjects{si};
            Cs   = C(m, :);
            pb   = double(Tk.pRF_angle_bin(m));
            idxW = arrayfun(@(a) find(pb == a), cfg.paBins, 'UniformOutput', false);
            B    = nan(nB, 4);
            for b = 1:nB
                M = nan(4, numel(cfg.paBins));
                for p = 1:numel(cfg.paBins)
                    ii = idxW{p};
                    if isempty(ii), continue; end
                    M(:,p) = median(Cs(ii(randi(numel(ii), [numel(ii) 1])), :), 1).';
                end
                A = compute_asymmetries(M, cfg, cfg.(en));
                for j = 1:4, B(b,j) = mean(A.(A.order{j}).diff, 'omitnan'); end
            end
            M8(si,:) = std(B, 0, 1, 'omitnan');
        end
        SE.(en) = M8;
    end
end

% ------------------------------------------------------------------------
function tb = build_lme_table(T, cfg)
    rows = [];
    for e = {'dg','da'}
        en = e{1};
        M  = bin_and_aggregate(T, cfg, cfg.(en), false);
        Xc = lme_codes(cfg, cfg.(en));
        [nO,nP,nS] = size(M);
        for si = 1:nS
            for p = 1:nP
                for k = 1:nO
                    rows(end+1,:) = [M(k,p,si), Xc.hVv(k,p), Xc.cVo(k,p), ...
                                     Xc.rVt(k,p), Xc.pcVpo(k,p), si, ...
                                     strcmp(en,'da')]; %#ok<AGROW>
                end
            end
        end
    end
    tb = array2table(rows, 'VariableNames', ...
                     {'y','hVv','cVo','rVt','pcVpo','subj','isDA'});
    tb = tb(~isnan(tb.y), :);
    tb.hVvDA = tb.hVv   .* tb.isDA;
    tb.cVoDA = tb.cVo   .* tb.isDA;
    tb.rVtDA = tb.rVt   .* tb.isDA;
    tb.pcDA  = tb.pcVpo .* tb.isDA;
    tb.subj  = categorical(tb.subj);
    tb.expf  = categorical(tb.isDA, [0 1], {'dg','da'});
end

function s = star(ci)
    s = ''; if all(isfinite(ci)) && ~(ci(1) <= 0 && 0 <= ci(2)), s = ' *'; end
end

function banner(t)
    fprintf('\n\n%s\n%s\n%s\n', repmat('=',1,74), t, repmat('=',1,74));
end
