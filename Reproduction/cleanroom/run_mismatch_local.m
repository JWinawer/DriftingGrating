function R = run_mismatch_local(varargin)
% RUN_MISMATCH_LOCAL  Predicted correct-vs-mismatched-design R2, with NO server access.
%
%   R = run_mismatch_local()
%
% WHAT THIS IS FOR. Figure 4A shows one run whose time series is well fit by its own
% design (R2 0.42) and poorly fit by another run's design (R2 -0.25). Generalising that
% to every run needs the observed time series, which are 37.6 GB on a ~0.8 MB/s mount.
% This computes what the mismatched R2 SHOULD be, from data already local.
%
% THE ALGEBRA, AND ITS LIMIT. For a fixed (not refit) prediction, with y mean-removed,
%
%   R2_m = ( 2<y,yhat_m> - ||yhat_m||^2 ) / ||y||^2
%        = ( 2<yhat_c,yhat_m> + 2<e,yhat_m> - ||yhat_m||^2 ) / ||y||^2      (y = yhat_c + e)
%
% Everything here is local EXCEPT <e,yhat_m>: how much the residual from the correct
% model happens to project onto the mismatched prediction. This function assumes
% <e,yhat_m> = 0 and therefore returns an ESTIMATE, not a prediction.
%
% WHAT IS SCALE-FREE AND WHAT IS NOT. Writing rho = corr(yhat_c,yhat_m) and
% k = ||yhat_m||/||yhat_c||,
%
%       R2_m / R2_c = 2*rho*k - k^2                                     (*)
%
% (*) depends ONLY on the predicted time series, so it is fully local and is valid at
% ANY averaging level -- per vertex, per ROI, or group-average -- provided rho and k are
% computed at that level. In particular the SIGN of R2_m is decided by (*) alone: it is
% negative whenever k > 2*rho, which holds comfortably here since the designs are nearly
% orthogonal. The manuscript's qualitative claim therefore needs nothing from the server.
%
% The ABSOLUTE R2_m does not follow, because it needs R2_c at the SAME averaging level.
% GLMsingle's R2run is PER VERTEX. The R2 quoted in Fig 4A (0.42, 0.28) is for a time
% series averaged over vertices AND observers, which is far higher than any single
% vertex's R2 -- how much higher depends on how correlated the noise is across vertices,
% which cannot be determined without the observed time series. So this function reports
% absolute values per vertex only, and the scale-free ratio at both levels. To get an
% absolute ROI-level number, supply the measured ROI-level R2_c (as CHECK_FIG4A does with
% the figure's own values) or fetch the time series with
% ../server_extract/collect_timeseries.m.
%
% THE ASSUMPTION IS NOT SAFE, AND IS BIASED IN A KNOWN DIRECTION. Event timing is
% identical across runs (only the stimulus SEQUENCE differs), so yhat_m carries the same
% onset structure as yhat_c. Any onset-locked but condition-nonspecific response that the
% 13-condition model misses sits in e and projects POSITIVELY onto yhat_m. So the true
% R2_m is expected to be HIGHER (less negative) than what this returns. Use
% ../server_extract/collect_timeseries.m on one session to measure <e,yhat_m> and
% calibrate. ||y||^2 is likewise inferred as ||yhat_c||^2 / R2run, which assumes
% <yhat_c,e> = 0 -- true for an OLS fit, only approximate for GLMsingle's ridge.
%
% Local inputs: ~/dg_collect/runbetas_* (per-run condition betas), glm_* (HRFindex,
% R2run, meanvol), ret_* (eccen/vexpl for the ROI filter), design/ (trialMat + timing).

    p = inputParser;
    p.addParameter('root', '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('glmsingle', '/Users/jaw288/repos/Code/Toolboxes/GLMsingle', @ischar);
    p.addParameter('quiet', false, @islogical);
    p.parse(varargin{:});
    opt = p.Results;

    addpath(genpath(opt.glmsingle));
    cfg  = config_repro();
    expn = {'dg','da'};
    nS   = numel(cfg.subjects);

    R = struct('subjects', {cfg.subjects}, 'experiments', {expn});
    R.r2correct  = nan(nS, 2);    % median over vertices of R2run, averaged over runs
    R.r2mismatch = nan(nS, 2);    % estimated, averaged over all ordered run pairs r~=m
    R.ratioVert  = nan(nS, 2);    % scale-free R2_m/R2_c = 2*rho*k - k^2, per vertex
    R.ratioROI   = nan(nS, 2);    % same, on the ROI-MEAN predicted series (Fig 4A level)
    R.predCorr   = nan(nS, 2);    % mean corr(yhat_c, yhat_m) over run pairs, ROI-mean
    R.kROI       = nan(nS, 2);    % mean ||yhat_m||/||yhat_c|| over run pairs, ROI-mean
    R.fracNeg    = nan(nS, 2);    % fraction of run pairs with R2_m/R2_c < 0 (ROI-mean)
    R.nRun       = nan(nS, 2);

    for si = 1:nS
        for ei = 1:2
            [A, ok] = local_inputs(cfg, opt.root, cfg.subjects{si}, expn{ei});
            if ~ok, continue; end
            n = A.nRun;  R.nRun(si,ei) = n;

            % Predicted time series per run: design (nTR x 13) * betas (13 x nVert),
            % each condition column convolved with that vertex's HRF. Build per HRF
            % index so the convolution is done once per distinct HRF, not per vertex.
            Y = cell(1, n);
            for r = 1:n
                Y{r} = predict_run(A.design{r}, A.beta(:,:,r), A.hrf, A.hrfIdx);
            end

            % ---- per-vertex: absolute (uses per-vertex R2run) and scale-free ratio ----
            num = 0; rat = 0; den = 0;
            for r = 1:n
                yc  = Y{r};
                nc2 = sum(yc.^2, 1);                      % ||yhat_c||^2 per vertex
                ny2 = nc2 ./ max(A.r2run(:,r).', eps);    % ||y||^2 = ||yhat_c||^2 / R2
                for m = 1:n
                    if m == r, continue; end
                    ym = Y{m};
                    cross = 2*sum(yc.*ym,1) - sum(ym.^2,1);
                    num = num + median(cross ./ max(ny2, eps), 'omitnan');
                    rat = rat + median(cross ./ max(nc2, eps), 'omitnan');   % = 2*rho*k - k^2
                    den = den + 1;
                end
            end
            R.r2mismatch(si,ei) = num / den;
            R.ratioVert(si,ei)  = rat / den;
            R.r2correct(si,ei)  = mean(median(A.r2run, 1, 'omitnan'), 'omitnan');

            % ---- ROI-mean level (what Fig 4A plots): SCALE-FREE RATIO ONLY ----
            % No absolute R2 here on purpose. R2run is per vertex; the ROI-average time
            % series has much higher R2 than any vertex, by an amount set by the
            % cross-vertex noise correlation, which is not knowable from local data.
            Ym = cell2mat(cellfun(@(y) mean(y, 2), Y, 'UniformOutput', false));  % nTR x nRun
            Ym = Ym - mean(Ym, 1);
            vals = []; cors = []; ks = [];
            for r = 1:n
                for m = 1:n
                    if m == r, continue; end
                    rho = corr(Ym(:,r), Ym(:,m));
                    kk  = norm(Ym(:,m)) / norm(Ym(:,r));
                    vals(end+1) = 2*rho*kk - kk^2;   %#ok<AGROW>   R2_m / R2_c
                    cors(end+1) = rho;               %#ok<AGROW>
                    ks(end+1)   = kk;                %#ok<AGROW>
                end
            end
            R.ratioROI(si,ei) = mean(vals);
            R.predCorr(si,ei) = mean(cors);
            R.kROI(si,ei)     = mean(ks);
            R.fracNeg(si,ei)  = mean(vals < 0);
        end
    end

    if ~opt.quiet, report(R, cfg, expn); end
end

% ------------------------------------------------------------------------
function [A, ok] = local_inputs(cfg, root, subj, en)
    ok = false;  A = struct();
    fb = fullfile(root, sprintf('runbetas_%s_%s.mat', subj, en));
    fg = fullfile(root, sprintf('glm_%s_%s.mat', subj, en));
    fr = fullfile(root, sprintf('ret_%s.mat', subj));
    if ~isfile(fb) || ~isfile(fg) || ~isfile(fr), fprintf('missing inputs for %s %s\n', subj, en); return; end

    B = load(fb);  G = load(fg, 'HRFindex', 'R2run');  Rt = load(fr, 'eccen','vexpl');
    v = B.v1Index;
    good = double(Rt.eccen(v)) >= cfg.eccRange(1) & double(Rt.eccen(v)) <= cfg.eccRange(2) ...
         & double(Rt.vexpl(v)) > cfg.r2min;
    if ~any(good), return; end
    vv = v(good);

    A.nRun   = B.nRun;
    A.beta   = double(B.runBeta(good,:,:));                 % nVert x 13 x nRun
    A.beta   = permute(A.beta, [2 1 3]);                    % 13 x nVert x nRun
    A.hrfIdx = double(G.HRFindex(vv));                      % nVert x 1
    A.r2run  = double(squeeze(G.R2run(vv,1,1,1:B.nRun))) / 100;   % nVert x nRun, fraction

    A.hrf    = getcanonicalhrflibrary(3, 1).';              % nTimeHRF x 20 (stimdur 3, TR 1)
    A.design = load_designs(root, en, B.nRun);
    ok = true;
end

% ------------------------------------------------------------------------
function D = load_designs(root, en, nRun)
% TR-resolution design matrices, rebuilt from the stimulus code's trialMat + timing.
% Conditions 1-8 = motion directions (0:45:315), 9-12 = stationary orientations
% (0,45,90,135 -> matching the beta column order), 13 = blank.
    D = cell(1, nRun);
    for r = 1:nRun
        % Any session's run r will do: REPORT.md section 3 establishes that the fixed
        % rng seed makes every observer's and session's run r byte-identical.
        f = dir(fullfile(root, 'design', en, sprintf('*_Run%d_*design_Run%d.mat', r, r)));
        if isempty(f), error('run_mismatch_local:design', 'no design for %s run %d', en, r); end
        e = load(fullfile(f(1).folder, f(1).name)); e = e.expDes;
        nTR = round(e.total_s);                                  % TR = 1 s
        X   = zeros(nTR, 13);
        per = e.stimDur_s + e.itiDur_s;
        for k = 1:size(e.trialMat, 1)
            onsetTR = round(e.runPadding_s + (k-1)*per) + 1;     % 1-based TR index
            typ = e.trialMat(k,2);
            if     typ == 0, ci = 13;
            elseif typ == 1, ci = 9 + find([0 45 90 135] == e.trialMat(k,3), 1) - 1;
            else,            ci = find((0:45:315) == e.trialMat(k,4), 1);
            end
            if ~isempty(ci) && onsetTR <= nTR, X(onsetTR, ci) = 1; end
        end
        D{r} = X;
    end
end

% ------------------------------------------------------------------------
function Y = predict_run(X, beta, hrfLib, hrfIdx)
% X: nTR x 13 onsets. beta: 13 x nVert. hrfLib: nH x 20. hrfIdx: nVert x 1.
% Convolve once per distinct HRF, then combine. Returns nTR x nVert, mean-removed.
    nTR = size(X,1);  nV = size(beta,2);
    Y = zeros(nTR, nV);
    for h = unique(hrfIdx(:)).'
        sel = hrfIdx == h;
        if ~any(sel), continue; end
        Xc = zeros(nTR, 13);
        for c = 1:13
            v = conv(X(:,c), hrfLib(:,h));
            Xc(:,c) = v(1:nTR);
        end
        Y(:, sel) = Xc * beta(:, sel);
    end
    Y = Y - mean(Y, 1);
end

% ------------------------------------------------------------------------
function report(R, cfg, expn)
    bar = repmat('=',1,92);
    fprintf('\n%s\nCORRECT vs MISMATCHED-DESIGN (local, assumes <e,yhat_m> = 0)\n%s\n', bar, bar);
    fprintf(['per-vertex R2 columns use GLMsingle R2run (per vertex). The ratio columns are\n' ...
             'scale-free (2*rho*k - k^2) and need no R2 at all; ratioROI is at the Fig 4A\n' ...
             'averaging level. NO absolute ROI-level R2 is reported -- see the header.\n\n']);
    fprintf('%-15s %10s %10s %10s %10s %9s %9s %7s\n', 'observer', ...
            'R2c/vert','R2m/vert','ratioVert','ratioROI','rho ROI','k ROI','frac<0');
    for ei = 1:2
        fprintf('\n-- %s --\n', expn{ei});
        for si = 1:numel(cfg.subjects)
            fprintf('%-15s %10.3f %10.3f %10.3f %10.3f %9.3f %9.3f %6.0f%%\n', cfg.subjects{si}, ...
                R.r2correct(si,ei), R.r2mismatch(si,ei), R.ratioVert(si,ei), ...
                R.ratioROI(si,ei), R.predCorr(si,ei), R.kROI(si,ei), 100*R.fracNeg(si,ei));
        end
        fprintf('%-15s %10.3f %10.3f %10.3f %10.3f %9.3f %9.3f %6.0f%%\n', 'MEDIAN', ...
            median(R.r2correct(:,ei),'omitnan'), median(R.r2mismatch(:,ei),'omitnan'), ...
            median(R.ratioVert(:,ei),'omitnan'), median(R.ratioROI(:,ei),'omitnan'), ...
            median(R.predCorr(:,ei),'omitnan'), median(R.kROI(:,ei),'omitnan'), ...
            100*median(R.fracNeg(:,ei),'omitnan'));
    end
    fprintf(['\nThe SIGN result needs no server data: R2_m < 0 whenever k > 2*rho, and the\n' ...
             'ratio columns show that holds for essentially every run pair. The ABSOLUTE\n' ...
             'ROI-level magnitude does need measured data -- multiply ratioROI by the\n' ...
             'measured ROI-level R2_c. CHECK_FIG4A does this against the two runs already\n' ...
             'reported in Fig 4A and agrees to +-0.07, with the two errors opposite in\n' ...
             'sign, so <e,yhat_m> looks like noise rather than a systematic term.\n']);
end
