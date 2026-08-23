function K = calibrate_run_mismatch(varargin)
% CALIBRATE_RUN_MISMATCH  Measure what RUN_MISMATCH_LOCAL has to assume.
%
%   K = calibrate_run_mismatch()
%
% Uses the observed time series fetched by ../server_extract/collect_timeseries.m for one
% session to measure three things the local estimator cannot get:
%
%   1. MEASURED R2_m for every ordered run pair, against the local prediction
%      R2_m = R2_c * (2*rho*k - k^2). The gap IS 2*<e,yhat_m>/||y||^2, the term the local
%      estimate drops -- so this measures it over 56 pairs instead of the 2 available
%      from Fig 4A.
%   2. The PER-VERTEX -> ROI R2 relationship, i.e. how much averaging over vertices buys.
%      This is set by the cross-vertex noise correlation and is what blocks converting the
%      scale-free ratio into an absolute ROI-level R2.
%   3. Measured per-vertex R2_c against GLMsingle's own R2run, as a sanity check that the
%      prediction is being built correctly.
%
% NOTE ON BETAS. Predictions use the OVERALL condition betas (mean of the per-run betas
% across all runs), not the run's own betas. Using a run's own betas to predict that run
% would overfit and inflate R2_c. This matches how Fig 4A's prediction is built, from
% GLMsingle's betas which are fit across all runs.
%
% NOTE ON LEVEL. Fig 4A averages over vertices AND observers; only one observer is
% fetched here, so the ROI-level R2 below is single-observer and will be lower than the
% figure's 0.42/0.28. The ratio statistics are unaffected.

    p = inputParser;
    p.addParameter('root', '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('tsdir', '/Users/jaw288/dg_collect/timeseries', @ischar);
    p.addParameter('glmsingle', '/Users/jaw288/repos/Code/Toolboxes/GLMsingle', @ischar);
    p.addParameter('subject', 'sub-0037', @ischar);
    p.addParameter('experiment', 'dg', @ischar);
    p.parse(varargin{:});
    opt = p.Results;

    addpath(genpath(opt.glmsingle));
    cfg = config_repro();

    B  = load(fullfile(opt.root, sprintf('runbetas_%s_%s.mat', opt.subject, opt.experiment)));
    G  = load(fullfile(opt.root, sprintf('glm_%s_%s.mat', opt.subject, opt.experiment)), 'HRFindex','R2run');
    Rt = load(fullfile(opt.root, sprintf('ret_%s.mat', opt.subject)), 'eccen','vexpl');

    v    = B.v1Index;
    good = double(Rt.eccen(v)) >= cfg.eccRange(1) & double(Rt.eccen(v)) <= cfg.eccRange(2) ...
         & double(Rt.vexpl(v)) > cfg.r2min;
    n    = B.nRun;

    betaAll = permute(mean(double(B.runBeta(good,:,:)), 3), [2 1 3]);   % 13 x nVert, OVERALL
    hrfIdx  = double(G.HRFindex(v(good)));
    hrfLib  = getcanonicalhrflibrary(3, 1).';
    design  = load_designs(opt.root, opt.experiment, n);

    % observed + predicted, per run
    Y = cell(1,n); P = cell(1,n);
    for r = 1:n
        S = load(fullfile(opt.tsdir, sprintf('ts_%s_%s_run%d.mat', opt.subject, opt.experiment, r)), 'ts');
        y = double(S.ts(good,:)).';                    % nTR x nVert
        Y{r} = y - mean(y,1);
        P{r} = predict_run(design{r}, betaAll, hrfLib, hrfIdx);
    end
    nV = size(Y{1},2);

    % ---- 1. measured vs predicted, per run pair ----
    rowsV = []; rowsR = [];
    r2cV = nan(1,n); r2cR = nan(1,n);
    for r = 1:n
        y = Y{r}; yc = P{r};
        r2cV(r) = median(1 - sum((y-yc).^2,1)./sum(y.^2,1));
        ym_ = mean(y,2); yc_ = mean(yc,2); ym_ = ym_-mean(ym_); yc_ = yc_-mean(yc_);
        r2cR(r) = 1 - sum((ym_-yc_).^2)/sum(ym_.^2);
        for m = 1:n
            if m==r, continue; end
            yy = P{m};
            % per vertex
            meas = median(1 - sum((y-yy).^2,1)./sum(y.^2,1));
            rho  = median(sum(yc.*yy,1)./sqrt(sum(yc.^2,1).*sum(yy.^2,1)));
            k    = median(sqrt(sum(yy.^2,1)./sum(yc.^2,1)));
            rowsV(end+1,:) = [r m meas r2cV(r)*(2*rho*k-k^2) rho k]; %#ok<AGROW>
            % ROI mean
            yr_ = mean(yy,2); yr_ = yr_-mean(yr_);
            measR = 1 - sum((ym_-yr_).^2)/sum(ym_.^2);
            rhoR  = corr(yc_, yr_);  kR = norm(yr_)/norm(yc_);
            rowsR(end+1,:) = [r m measR r2cR(r)*(2*rhoR*kR-kR^2) rhoR kR]; %#ok<AGROW>
        end
    end

    K.perVertex = rowsV; K.roi = rowsR; K.r2cVert = r2cV; K.r2cROI = r2cR;
    bar = repmat('=',1,88);

    fprintf('\n%s\nCALIBRATION: %s %s, %d runs, %d analysed vertices\n%s\n', ...
            bar, opt.subject, opt.experiment, n, nV, bar);
    fprintf('measured R2_correct  per-vertex median %.3f   ROI-mean %.3f\n', mean(r2cV), mean(r2cR));
    fprintf('GLMsingle R2run      per-vertex median %.3f   (sanity check)\n', ...
            median(median(double(squeeze(G.R2run(v(good),1,1,1:n))),1))/100);

    fprintf('\n-- 2. what averaging over vertices buys --\n');
    fprintf('ROI-mean R2 / per-vertex R2 = %.2f x\n', mean(r2cR)/mean(r2cV));
    E = cell2mat(cellfun(@(a,b) a-b, Y, P, 'UniformOutput', false));
    idx = randperm(nV, min(400,nV));
    C = corr(E(:,idx)); C = C(triu(true(size(C)),1));
    fprintf('mean cross-vertex residual correlation = %.3f (400-vertex sample)\n', mean(C));

    % ---- 3. the two inner products the formula assumes are zero ----
    ec = []; em = [];
    for r = 1:n
        y = Y{r}; yc = P{r}; e = y - yc; ny2 = sum(y.^2,1);
        ec(end+1) = median(2*sum(e.*yc,1)./ny2); %#ok<AGROW>
        for m = 1:n
            if m==r, continue; end
            em(end+1) = median(2*sum(e.*P{m},1)./ny2); %#ok<AGROW>
        end
    end
    K.eDotYc = ec; K.eDotYm = em;
    fprintf('\n-- 3. the terms the formula assumes are zero (per vertex, /||y||^2) --\n');
    fprintf('   2*<e,yhat_c>/||y||^2 : mean %+.4f   <- shrinkage of an out-of-sample prediction\n', mean(ec));
    fprintf('   2*<e,yhat_m>/||y||^2 : mean %+.4f  SD %.4f\n', mean(em), std(em));

    for lev = 1:2
        if lev==1, T = rowsV; nm='PER VERTEX (median over vertices)';
        else,      T = rowsR; nm='ROI MEAN (single observer)'; end
        gap = T(:,3) - T(:,4);
        fprintf('\n-- 1. measured vs predicted R2_m, %s --\n', nm);
        fprintf('   %d run pairs | measured mean %+.3f  predicted mean %+.3f\n', ...
                size(T,1), mean(T(:,3)), mean(T(:,4)));
        fprintf('   gap (measured - predicted): mean %+.3f  SD %.3f  range [%+.3f %+.3f]\n', ...
                mean(gap), std(gap), min(gap), max(gap));
        fprintf('   measured R2_m negative in %d of %d pairs\n', sum(T(:,3)<0), size(T,1));
        fprintf('   corr(measured, predicted) across pairs = %.3f\n', corr(T(:,3), T(:,4)));
    end
    fprintf(['\nThe gap is 2*<e,yhat_m>/||y||^2. A mean near zero means the local estimator is\n' ...
             'unbiased; the SD is its per-pair accuracy.\n']);
end

function D = load_designs(root, en, nRun)
    D = cell(1, nRun);
    for r = 1:nRun
        f = dir(fullfile(root, 'design', en, sprintf('*_Run%d_*design_Run%d.mat', r, r)));
        e = load(fullfile(f(1).folder, f(1).name)); e = e.expDes;
        nTR = round(e.total_s); X = zeros(nTR,13); per = e.stimDur_s + e.itiDur_s;
        for k = 1:size(e.trialMat,1)
            onsetTR = round(e.runPadding_s + (k-1)*per) + 1;
            typ = e.trialMat(k,2);
            if     typ==0, ci=13;
            elseif typ==1, ci = 9 + find([0 45 90 135]==e.trialMat(k,3),1) - 1;
            else,          ci = find((0:45:315)==e.trialMat(k,4),1);
            end
            if ~isempty(ci) && onsetTR<=nTR, X(onsetTR,ci)=1; end
        end
        D{r} = X;
    end
end

function Y = predict_run(X, beta, hrfLib, hrfIdx)
    nTR = size(X,1); Y = zeros(nTR, size(beta,2));
    for h = unique(hrfIdx(:)).'
        sel = hrfIdx==h; if ~any(sel), continue; end
        Xc = zeros(nTR,13);
        for c = 1:13, vv = conv(X(:,c), hrfLib(:,h)); Xc(:,c) = vv(1:nTR); end
        Y(:,sel) = Xc * beta(:,sel);
    end
    Y = Y - mean(Y,1);
end
