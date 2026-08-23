function C = check_fig4a(varargin)
% CHECK_FIG4A  Validate RUN_MISMATCH_LOCAL against the two runs already reported in Fig 4A.
%
%   C = check_fig4a()
%
% Fig 4A reports measured values for the group-average time series of one ROI (right
% horizontal V1, Cartesian experiment): correct-design R2 = 0.42 (run 2) and 0.28
% (run 5), and opposite-run R2 = -0.25 and -0.53. That is enough to test the assumption
% RUN_MISMATCH_LOCAL rests on, with no server access.
%
% For a fixed (not refit) prediction, writing rho = corr(yhat_c, yhat_m) and
% k = ||yhat_m|| / ||yhat_c||,
%
%       R2_m = R2_c * (2*rho*k - k^2)                                   (*)
%
% provided <e, yhat_m> = 0 and <e, yhat_c> = 0. rho and k are computable locally from
% the predicted time series; R2_c is read off the figure. So (*) turns each panel into a
% falsifiable prediction of the reported mismatched R2. A discrepancy is a direct
% measurement of the neglected <e, yhat_m> term.

    p = inputParser;
    p.addParameter('root', '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('glmsingle', '/Users/jaw288/repos/Code/Toolboxes/GLMsingle', @ischar);
    p.addParameter('runs', [2 5], @isnumeric);       % the two runs shown in Fig 4A
    p.addParameter('wedge', 1, @isnumeric);          % cfg.paBins(1) = 0 deg = right horizontal
    p.addParameter('r2correct', [0.42 0.28], @isnumeric);
    p.addParameter('r2mismatch', [-0.25 -0.53], @isnumeric);
    p.parse(varargin{:});
    opt = p.Results;

    addpath(genpath(opt.glmsingle));
    cfg = config_repro();

    % Group-average predicted ROI time series for each of the two runs.
    Yg = zeros(280, numel(opt.runs));  nUsed = 0;
    for si = 1:numel(cfg.subjects)
        [A, ok] = inputs_for(cfg, opt.root, cfg.subjects{si}, 'dg', opt.wedge);
        if ~ok, continue; end
        for j = 1:numel(opt.runs)
            r = opt.runs(j);
            if r > A.nRun, continue; end
            y = predict_run(A.design{r}, A.beta(:,:,r), A.hrf, A.hrfIdx);
            Yg(:,j) = Yg(:,j) + mean(y, 2);
        end
        nUsed = nUsed + 1;
    end
    Yg = Yg / nUsed;
    Yg = Yg - mean(Yg, 1);

    fprintf('\n%s\nFIG 4A CHECK  (group average over %d observers, ROI = %d deg wedge, dg)\n%s\n', ...
            repmat('=',1,84), nUsed, cfg.paBins(opt.wedge), repmat('=',1,84));
    fprintf('%-28s %10s %10s %12s %12s %10s\n', 'panel','rho','k','R2c (fig)','R2m pred','R2m (fig)');

    C = struct('runs', opt.runs, 'rho', nan(1,2), 'k', nan(1,2), ...
               'predicted', nan(1,2), 'reported', opt.r2mismatch);
    for j = 1:numel(opt.runs)
        c = j;  m = 3 - j;                       % "the opposite run"
        yc = Yg(:,c);  ym = Yg(:,m);
        rho = corr(yc, ym);
        k   = norm(ym) / norm(yc);
        pred = opt.r2correct(j) * (2*rho*k - k^2);
        C.rho(j) = rho;  C.k(j) = k;  C.predicted(j) = pred;
        fprintf('run %d observed, run %d design %10.3f %10.3f %12.2f %12.3f %10.2f\n', ...
                opt.runs(c), opt.runs(m), rho, k, opt.r2correct(j), pred, opt.r2mismatch(j));
    end

    % The gap IS the neglected term: R2m_reported - R2m_predicted = 2<e,yhat_m>/||y||^2.
    C.gap = C.reported - C.predicted;
    fprintf('\ngap (reported - predicted) = %s\n', sprintf('%+.3f  ', C.gap));
    fprintf(['This gap is exactly 2*<e,yhat_m>/||y||^2, the term the local estimate drops.\n' ...
             'Positive gap => the correct model''s residual projects POSITIVELY onto the\n' ...
             'mismatched prediction (shared event timing), so the local estimate is too\n' ...
             'pessimistic. Negative gap => the opposite.\n']);
end

% ------------------------------------------------------------------------
function [A, ok] = inputs_for(cfg, root, subj, en, wedge)
    ok = false;  A = struct();
    fb = fullfile(root, sprintf('runbetas_%s_%s.mat', subj, en));
    fg = fullfile(root, sprintf('glm_%s_%s.mat', subj, en));
    fr = fullfile(root, sprintf('ret_%s.mat', subj));
    if ~isfile(fb) || ~isfile(fg) || ~isfile(fr), return; end
    B = load(fb);  G = load(fg, 'HRFindex');  Rt = load(fr, 'eccen','vexpl','angle_adj');
    v = B.v1Index;
    good = double(Rt.eccen(v)) >= cfg.eccRange(1) & double(Rt.eccen(v)) <= cfg.eccRange(2) ...
         & double(Rt.vexpl(v)) > cfg.r2min;
    ang  = mod(90 - double(Rt.angle_adj(v)), 360);
    [~, wi] = min(abs(mod(ang - cfg.paBins(:).' + 180, 360) - 180), [], 2);
    good = good & (wi == wedge);
    if ~any(good), return; end
    A.nRun   = B.nRun;
    A.beta   = permute(double(B.runBeta(good,:,:)), [2 1 3]);
    A.hrfIdx = double(G.HRFindex(v(good)));
    A.hrf    = getcanonicalhrflibrary(3, 1).';
    A.design = load_designs(root, en, B.nRun);
    ok = true;
end

function D = load_designs(root, en, nRun)
    D = cell(1, nRun);
    for r = 1:nRun
        f = dir(fullfile(root, 'design', en, sprintf('*_Run%d_*design_Run%d.mat', r, r)));
        e = load(fullfile(f(1).folder, f(1).name)); e = e.expDes;
        nTR = round(e.total_s);  X = zeros(nTR, 13);  per = e.stimDur_s + e.itiDur_s;
        for k = 1:size(e.trialMat,1)
            onsetTR = round(e.runPadding_s + (k-1)*per) + 1;
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

function Y = predict_run(X, beta, hrfLib, hrfIdx)
    nTR = size(X,1);  Y = zeros(nTR, size(beta,2));
    for h = unique(hrfIdx(:)).'
        sel = hrfIdx == h;  if ~any(sel), continue; end
        Xc = zeros(nTR, 13);
        for c = 1:13, v = conv(X(:,c), hrfLib(:,h)); Xc(:,c) = v(1:nTR); end
        Y(:, sel) = Xc * beta(:, sel);
    end
    Y = Y - mean(Y, 1);
end
