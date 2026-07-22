% COMPARE_SUBJECT_WEIGHTING  Three ways to weight observers in the group mean.
%
% Every choice in this analysis implies a subject weighting, and the per-vertex
% z-scoring makes one implicitly. This script makes the three candidates explicit and
% evaluates them side by side (see ../ZSCORE_FIG7.md section 6):
%
%   W1 equal weight        - one observer, one vote (what the raw analysis does)
%   W2 precision weight    - w ∝ 1/SE^2, SE = within-subject bootstrap standard error of
%                            that subject's own asymmetry estimate (resampling vertices).
%                            This is the standard inverse-variance / random-effects weight.
%   W3 amplitude weight    - w ∝ 1/beta_std, which is what per-vertex z-scoring produces.
%
% W2 and W3 answer different questions: W2 asks "how well measured is this subject's
% effect?", W3 asks "how big is this subject's overall response?". They are not
% interchangeable, and here they point in opposite directions.

cfg = config_repro();
T   = load_and_filter(cfg);
keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
T = T(keep,:);
subjV = string(T.subject); pabV = T.pRF_angle_bin;
nS = numel(cfg.subjects); nP = numel(cfg.paBins);
asymNames = {'HV','cardObl','radTan','polcardPolobl'};
B = 500;

for E = {cfg.dg, cfg.da}
    E = E{1}; %#ok<FXSET>
    X = zeros(height(T),4);
    for k = 1:4, X(:,k) = T.(E.oriCols{k}) - T.(E.blank); end
    s = T.(E.betaStd);

    % cell index lists (subject x wedge)
    idxCell = cell(nP, nS);
    for si = 1:nS
        i0 = subjV == cfg.subjects{si};
        for pi = 1:nP, idxCell{pi,si} = find(i0 & (pabV == cfg.paBins(pi))); end
    end

    % point estimate, per subject
    M = nan(4,nP,nS);
    for si=1:nS, for pi=1:nP
        ix = idxCell{pi,si}; if ~isempty(ix), M(:,pi,si) = median(X(ix,:),1); end
    end, end
    A = compute_asymmetries(M, cfg, E);
    aSubj = zeros(4,nS);
    for k=1:4, aSubj(k,:) = mean(A.(asymNames{k}).diff,1,'omitnan'); end

    % within-subject bootstrap over vertices -> SE of each subject's asymmetry
    rng(0);
    boot = nan(4,nS,B);
    for b = 1:B
        Mb = nan(4,nP,nS);
        for si=1:nS, for pi=1:nP
            ix = idxCell{pi,si};
            if ~isempty(ix), Mb(:,pi,si) = median(X(ix(randi(numel(ix),numel(ix),1)),:),1); end
        end, end
        Ab = compute_asymmetries(Mb, cfg, E);
        for k=1:4, boot(k,:,b) = mean(Ab.(asymNames{k}).diff,1,'omitnan'); end
    end
    SE = std(boot,0,3);                       % 4 x nS

    sSubj = zeros(1,nS);
    for si=1:nS, sSubj(si) = median(s(subjV==cfg.subjects{si})); end

    fprintf('\n================ %s ================\n', E.name);
    fprintf('per-subject radTan and H-V with within-subject bootstrap SE (raw units):\n');
    fprintf('%-14s %18s %18s | %8s %8s %8s\n','subject','radTan +/- SE','H-V +/- SE', ...
            'w equal','w prec','w ampl');
    kRT = 3; kHV = 1;
    wEq = ones(1,nS)/nS;
    wPr = (1./SE(kRT,:).^2); wPr = wPr/sum(wPr);         % shown for radTan
    wAm = (1./sSubj);        wAm = wAm/sum(wAm);
    for si=1:nS
        fprintf('%-14s %9.3f +/- %-6.3f %9.3f +/- %-6.3f | %7.1f%% %7.1f%% %7.1f%%\n', ...
            cfg.subjects{si}, aSubj(kRT,si), SE(kRT,si), aSubj(kHV,si), SE(kHV,si), ...
            100*wEq(si), 100*wPr(si), 100*wAm(si));
    end
    fprintf('corr(precision weight, amplitude weight) across subjects = %+.3f\n', ...
        corr(wPr(:), wAm(:)));

    % The decisive quantity: is the spread across subjects measurement noise, or real?
    fprintf('\nvariance decomposition (is between-subject spread just noise?):\n');
    for k = [3 1]
        sdBetween = std(aSubj(k,:));
        seWithin  = mean(SE(k,:));
        fprintf('  %-14s between-subject SD %.3f vs mean within-subject SE %.3f  -> var ratio %.0fx\n', ...
            asymNames{k}, sdBetween, seWithin, (sdBetween/seWithin)^2);
    end
    fprintf(['  (ratio >> 1 means observers genuinely differ; inverse-variance weights then\n' ...
             '   converge on equal weights, and down-weighting an observer discards signal,\n' ...
             '   not noise.)\n']);

    fprintf('\ngroup means under each weighting (raw units):\n');
    fprintf('%-22s %9s %9s %9s %9s   %s\n','weighting',asymNames{:},'ordering');
    schemes = {'W1 equal', wEq; 'W2 precision', [] ; 'W3 amplitude (1/betaSD)', wAm};
    for w = 1:3
        v = zeros(1,4);
        for k = 1:4
            if w==2, ww = 1./SE(k,:).^2; ww = ww/sum(ww);   % per-asymmetry precision
            else,    ww = schemes{w,2}; end
            v(k) = sum(ww .* aSubj(k,:));
        end
        ord = 'H-V larger'; if v(3) > abs(v(1)), ord = 'radTan larger'; end
        fprintf('%-22s %9.3f %9.3f %9.3f %9.3f   %s\n', schemes{w,1}, v, ord);
    end
end

fprintf(['\nNote: W3 is applied here as a per-subject scalar for comparability. The published\n' ...
         'per-vertex z-scoring is W3 plus a units change; the ladder in diagnose_zscore_fig7.m\n' ...
         'shows the per-subject scalar reproduces its whole effect on Fig 7.\n']);
