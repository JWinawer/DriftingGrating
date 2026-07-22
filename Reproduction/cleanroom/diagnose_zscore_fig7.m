% DIAGNOSE_ZSCORE_FIG7  Why z-scoring reverses the rank of radTan and HV in Fig 7B (da).
%
% In the z-scored Fig 7 the largest polar-grating asymmetry is radial-vs-tangential
% (0.603) and horizontal-vs-vertical is second (-0.446); without z-scoring the order
% reverses (0.150 vs -0.211). This script produces the four pieces of evidence that
% account for that, documented in ../ZSCORE_FIG7.md:
%
%   1. Divisor ladder     - at which level of beta_std variation the differential
%                           amplification appears (global / subject / wedge / vertex).
%   2. Amplification identity - amp = mean(1/s) + cov(a,1/s)/mean(a), evaluated per
%                           asymmetry, showing one covariance with two opposite effects.
%   3. Harmonic identity  - for da, radTan is the polar-angle MEAN of pinwheel-annulus
%                           and HV is its cos(2*theta) MODULATION, so they trade off.
%   4. Robustness         - leave-one-subject-out and a subject bootstrap on the ordering.
%
% Run from this folder (needs config_repro, load_and_filter, compute_asymmetries).

cfg = config_repro();
T   = load_and_filter(cfg);
keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
T = T(keep, :);
subjV = string(T.subject);  pabV = T.pRF_angle_bin;
nS = numel(cfg.subjects);   nP = numel(cfg.paBins);
asymNames = {'HV','cardObl','radTan','polcardPolobl'};

%% ---------------------------------------------------------------- 1. ladder
rungNames = {'D0 raw','D1 global scalar','D2 per subject','D3 per subj x wedge','D4 per vertex'};
for E = {cfg.dg, cfg.da}
    E = E{1}; %#ok<FXSET>
    X = zeros(height(T),4);
    for k = 1:4, X(:,k) = T.(E.oriCols{k}) - T.(E.blank); end
    s = T.(E.betaStd);

    D = {ones(height(T),1), median(s)*ones(height(T),1), nan(height(T),1), nan(height(T),1), s};
    for si = 1:nS
        inS = subjV == cfg.subjects{si};
        D{3}(inS) = median(s(inS));
        for pi = 1:nP
            idx = inS & (pabV == cfg.paBins(pi));
            if any(idx), D{4}(idx) = median(s(idx)); end
        end
    end

    fprintf('\n================ %s : divisor ladder ================\n', E.name);
    fprintf('median beta_std = %.4f\n\n%-22s %9s %9s %9s %9s\n', median(s), 'rung', asymNames{:});
    A0 = [];
    for r = 1:5
        A = compute_asymmetries(wedgeMedian(X./D{r}, subjV, pabV, cfg), cfg, E);
        v = cellfun(@(n) mean(A.(n).diff(:),'omitnan'), asymNames);
        if r == 1, A0 = v; end
        fprintf('%-22s %9.4f %9.4f %9.4f %9.4f\n', rungNames{r}, v);
        if r > 1, fprintf('%-22s %9.2f %9.2f %9.2f %9.2f   <- amp vs raw\n', '  (ratio)', v./A0); end
    end
end

%% --------------------------------------------- 2-4. da: the reversal in detail
E = cfg.da;
X = zeros(height(T),4);
for k = 1:4, X(:,k) = T.(E.oriCols{k}) - T.(E.blank); end
s = T.(E.betaStd);

Mraw = wedgeMedian(X, subjV, pabV, cfg);
Mz   = wedgeMedian(X./s, subjV, pabV, cfg);
Araw = compute_asymmetries(Mraw, cfg, E);
Az   = compute_asymmetries(Mz,   cfg, E);
rawRT = mean(Araw.radTan.diff,1,'omitnan');  rawHV = mean(Araw.HV.diff,1,'omitnan');
zRT   = mean(Az.radTan.diff,1,'omitnan');    zHV   = mean(Az.HV.diff,1,'omitnan');

sSubj = zeros(1,nS); respSubj = zeros(1,nS); r2Subj = zeros(1,nS); sdgSubj = zeros(1,nS);
for si = 1:nS
    inS = subjV == cfg.subjects{si};
    sSubj(si)    = median(s(inS));
    sdgSubj(si)  = median(T.(cfg.dg.betaStd)(inS));
    r2Subj(si)   = median(T.pRF_r2(inS));
    respSubj(si) = median(mean(X(inS,:),2));
end

fprintf('\n================ da : per-subject structure ================\n');
fprintf(['%-14s %7s %8s %8s | %8s %8s | %8s %8s | %7s\n'], ...
    'subject','betaSD','meanResp','medR2','rawRT','zRT','rawHV','zHV','wt%%');
w = (1./sSubj)/sum(1./sSubj)*100;
for si = 1:nS
    fprintf('%-14s %7.3f %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f | %6.1f%%\n', ...
        cfg.subjects{si}, sSubj(si), respSubj(si), r2Subj(si), ...
        rawRT(si), zRT(si), rawHV(si), zHV(si), w(si));
end
fprintf('%-14s %7.3f %8s %8s | %8.3f %8.3f | %8.3f %8.3f | %6.1f%%\n', 'GROUP', median(s), '', '', ...
    mean(rawRT), mean(zRT), mean(rawHV), mean(zHV), 12.5);
fprintf('("wt%%" = each subject''s effective weight under z-scoring; uniform would be 12.5%%)\n');

fprintf('\n---- amplification identity: amp = mean(1/s) + cov(a,1/s)/mean(a) ----\n');
inv_s = 1./sSubj;
fprintf('mean_i(1/betaSD_i) = %.3f  <- generic amplification if a and 1/s were unrelated\n', mean(inv_s));
for c = 1:2
    if c==1, nm='radTan'; a=rawRT; else, nm='HV'; a=rawHV; end
    cv = mean((a-mean(a)).*(inv_s-mean(inv_s)));
    fprintf('%-7s mean(a)=%+.4f  cov(a,1/s)=%+.4f  cov/mean(a)=%+.3f  => amp %.2f (actual %.2f)\n', ...
        nm, mean(a), cv, cv/mean(a), mean(inv_s)+cv/mean(a), mean(a.*inv_s)/mean(a));
end

fprintf('\n---- harmonic identity: radTan = DC, HV = cos(2*theta) of the same pattern ----\n');
iP = @(d) find(cfg.paBins == d);
for c = 1:2
    if c==1, nm='raw'; M=Mraw; else, nm='zscored'; M=Mz; end
    A  = compute_asymmetries(M, cfg, E);
    rp = mean(A.radTan.diff, 2, 'omitnan');
    q  = squeeze(mean(M(4,:,:) - M(3,:,:), 3));            % cspiral - ccspiral per wedge
    c2r = (rp(iP(0))+rp(iP(180))-rp(iP(90))-rp(iP(270)))/8;
    c2s = (q(iP(45))+q(iP(225))-q(iP(135))-q(iP(315)))/8;
    fprintf('%-8s radTan = mean_theta(rad-tang)          = %+.4f\n', nm, mean(rp));
    fprintf('%-8s HV     = cos2th(rad-tang) %+.4f + cos2th(spirals) %+.4f = %+.4f (direct %+.4f)\n', ...
        nm, c2r, c2s, c2r+c2s, mean(mean(A.HV.diff,1,'omitnan')));
end

fprintf('\n---- per-wedge rad-tang ----\n%5s %10s %10s\n','PA','raw','zscored');
rpRaw = mean(Araw.radTan.diff,2,'omitnan'); rpZ = mean(Az.radTan.diff,2,'omitnan');
for pi = 1:nP, fprintf('%5d %10.3f %10.3f\n', cfg.paBins(pi), rpRaw(pi), rpZ(pi)); end

fprintf('\n---- robustness of the ordering ----\n');
fprintf('all 8 subjects: raw  radTan %.3f vs |HV| %.3f -> %s\n', mean(rawRT), abs(mean(rawHV)), ...
    pick(mean(rawRT)>abs(mean(rawHV)),'radTan larger','HV larger'));
fprintf('all 8 subjects: z    radTan %.3f vs |HV| %.3f -> %s\n', mean(zRT), abs(mean(zHV)), ...
    pick(mean(zRT)>abs(mean(zHV)),'radTan larger','HV larger'));
for si = 1:nS
    k = setdiff(1:nS, si);
    fprintf('  drop %-14s z: radTan %.3f vs |HV| %.3f -> %s\n', cfg.subjects{si}, ...
        mean(zRT(k)), abs(mean(zHV(k))), pick(mean(zRT(k))>abs(mean(zHV(k))),'radTan larger','HV LARGER'));
end
rng(0); B = 10000; winZ = 0; winR = 0;
for b = 1:B
    k = randi(nS,[1 nS]);
    winZ = winZ + (mean(zRT(k))   > abs(mean(zHV(k))));
    winR = winR + (mean(rawRT(k)) > abs(mean(rawHV(k))));
end
fprintf('subject bootstrap (B=%d): P(radTan > |HV|)  z-scored %.2f , raw %.2f\n', B, winZ/B, winR/B);

fprintf('\n---- what subject-level beta_std tracks ----\n');
fprintf('corr(da betaSD, median mean-response) = %+.3f\n', corr(sSubj', respSubj'));
fprintf('corr(da betaSD, median pRF r2)        = %+.3f\n', corr(sSubj', r2Subj'));
fprintf('corr(da betaSD, dg betaSD)            = %+.3f (Spearman %+.3f)\n', ...
    corr(sSubj', sdgSubj'), corr(sSubj', sdgSubj', 'type','Spearman'));

%% ------------------------------------------------------------------ helpers
function M = wedgeMedian(C, subjV, pabV, cfg)
    nS = numel(cfg.subjects); nP = numel(cfg.paBins);
    M = nan(4, nP, nS);
    for si = 1:nS
        inS = subjV == cfg.subjects{si};
        for pi = 1:nP
            idx = inS & (pabV == cfg.paBins(pi));
            if any(idx), M(:,pi,si) = median(C(idx,:), 1); end
        end
    end
end
function o = pick(c,a,b), if c, o = a; else, o = b; end, end
