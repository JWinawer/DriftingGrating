function S = splithalf_reliability(varargin)
% SPLITHALF_RELIABILITY  Split-half (over RUNS) reliability of the analysed quantities.
%
%   S = splithalf_reliability()
%
% WHY. The GLM-summary section of the manuscript currently reports fit quality (R2)
% from single examples. R2 is the wrong currency here: the reference condition is
% full-field pink noise, so R2 measures cross-condition DIFFERENTIATION, not
% responsiveness (see ../local_qc/REPORT.md section 1). What Figures 5-8 actually
% rest on is that the analysed profile -- 4 orientations x 8 polar-angle wedges per
% observer -- is REPRODUCIBLE. That is measured by resampling the measurement, i.e.
% splitting the runs.
%
% Uses ~/dg_collect/runbetas_<subject>_<exp>.mat (mean GLMsingle single-trial beta per
% V1 vertex x condition x run), written by ../server_extract/collect_runwise_betas.m.
% The design is balanced -- 4 trials per condition per run -- so a k-vs-k run split is
% balanced by construction.
%
% FOUR RELIABILITIES, at two spatial levels x two contrast levels. Each is the
% correlation between the two halves, Fisher-z averaged over all balanced splits and
% Spearman-Brown corrected to the full run count:
%
%   vertexFull  : vertex x orientation contrast map (nVert x 4). Includes each
%                 vertex's overall response level, which is reliable but carries no
%                 orientation information.
%   vertexDiff  : same, after removing each vertex's mean over the 4 orientations --
%                 the purely orientation-differential fine-grained pattern.
%   roiFull     : the analysed profile, 4 orientations x 8 wedges (32 values).
%   roiDiff     : the analysed profile after removing each wedge's mean over the 4
%                 orientations. This strips wedge-level gain differences (reliable,
%                 but not orientation effects) and leaves exactly the subspace that
%                 carries all four asymmetries. THIS IS THE STRICT TEST.
%
% Plus, per experiment, the within-observer SE of each of the four asymmetries from
% the same splits, against the between-observer SD -- the statistic that says whether
% the n=8 inference is limited by measurement noise or by individual differences.
%
% Odd run counts (sub-0255 dg has 9, sub-0395 da has 6): for an odd count each run in
% turn is dropped and the even procedure run on the rest, then averaged. The
% Spearman-Brown step then corrects to the even count, which slightly understates the
% reliability of the full 9-run estimate -- conservative.

    p = inputParser;
    p.addParameter('root', '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('quiet', false, @islogical);
    p.addParameter('aggregator', '', @ischar);   % override cfg.aggregator ('mean'|'median')
    p.parse(varargin{:});
    opt = p.Results;

    cfg  = config_repro();
    if ~isempty(opt.aggregator), cfg.aggregator = opt.aggregator; end
    opt.aggregator = cfg.aggregator;
    nm   = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    expn = {'dg','da'};
    nS   = numel(cfg.subjects);

    S = struct('subjects', {cfg.subjects}, 'names', {nm}, 'experiments', {expn}, ...
               'aggregator', opt.aggregator);
    fields = {'vertexFull','vertexDiff','roiFull','roiDiff'};
    for f = 1:numel(fields), S.(fields{f}) = nan(nS, 2); end
    S.nRun    = nan(nS, 2);
    S.nVert   = nan(nS, 2);
    S.nSplit  = nan(nS, 2);
    S.asym    = nan(nS, 4, 2);   % full-data asymmetry estimate
    S.asymSE  = nan(nS, 4, 2);   % within-observer SE from the splits

    for si = 1:nS
        for ei = 1:2
            f = fullfile(opt.root, sprintf('runbetas_%s_%s.mat', cfg.subjects{si}, expn{ei}));
            if ~isfile(f)
                fprintf('missing %s -- run ../server_extract/collect_runwise_betas.m\n', f);
                return
            end
            D = load(f);
            [A, ok] = prep_vertices(D, cfg, opt.root);
            if ~ok, continue; end
            A.expn = expn{ei};
            A.agg  = cfg.aggregator;

            n = D.nRun;
            S.nRun(si,ei)  = n;
            S.nVert(si,ei) = size(A.runBeta, 1);

            % full-data asymmetry estimate
            S.asym(si,:,ei) = asym_from_profile(profile_from_runs(A, 1:n, cfg), cfg, expn{ei});

            % ---- balanced splits ------------------------------------------------
            splits = balanced_splits(n);
            S.nSplit(si,ei) = size(splits, 1);
            nsp = size(splits, 1);
            z   = nan(nsp, 4);           % Fisher z per split, one column per field
            dA  = nan(nsp, 4);           % half1 - half2 asymmetry difference

            for k = 1:nsp
                h1 = splits{k,1};  h2 = splits{k,2};
                [C1, M1] = halves(A, h1, cfg, expn{ei});
                [C2, M2] = halves(A, h2, cfg, expn{ei});

                z(k,1) = fz(corr_vec(C1, C2));
                z(k,2) = fz(corr_vec(C1 - mean(C1,2), C2 - mean(C2,2)));
                z(k,3) = fz(corr_vec(M1, M2));
                z(k,4) = fz(corr_vec(M1 - mean(M1,1), M2 - mean(M2,1)));

                dA(k,:) = asym_from_profile(M1, cfg, expn{ei}) ...
                        - asym_from_profile(M2, cfg, expn{ei});
            end

            r = tanh(mean(z, 1, 'omitnan'));            % average correlation of a HALF
            r = 2*r ./ (1 + r);                          % Spearman-Brown -> full length
            for jj = 1:4, S.(fields{jj})(si,ei) = r(jj); end

            % SD of a half-estimate = SD(diff)/sqrt(2); full-data SE = that/sqrt(2)
            S.asymSE(si,:,ei) = std(dA, 0, 1, 'omitnan') / 2;
        end
    end

    % ---- observer gain weighting (the published route) ----------------------
    [S.scale, S.gain] = observer_gain_weights(cfg);
    S.asymW   = S.asym   .* S.scale;
    S.asymSEW = S.asymSE .* S.scale;

    if ~opt.quiet, report(S, cfg, nm, expn, fields); end
end

% ------------------------------------------------------------------------
function [A, ok] = prep_vertices(D, cfg, root)
% Restrict to the analysed vertices (4-8 deg, pRF R2 > 0.1) and attach wedge labels.
% Same filter and wedge convention as DIAGNOSE_WITHIN_OBSERVER_ERROR.
    ok = false;  A = struct();
    R = load(fullfile(root, sprintf('ret_%s.mat', D.subject)), 'eccen','vexpl','angle_adj');
    v = D.v1Index;
    good = double(R.eccen(v)) >= cfg.eccRange(1) & double(R.eccen(v)) <= cfg.eccRange(2) ...
         & double(R.vexpl(v)) > cfg.r2min;
    if ~any(good), return; end
    A.agg = 'mean';
    A.runBeta = D.runBeta(good, :, :);
    ang  = double(R.angle_adj(v(good)));
    conv = mod(90 - ang, 360);
    [~, A.wedge] = min(abs(mod(conv - cfg.paBins(:).' + 180, 360) - 180), [], 2);
    ok = true;
end

% ------------------------------------------------------------------------
function [C, M] = halves(A, runs, cfg, en)
% C : nVert x 4 orientation-minus-blank contrasts, averaged over the given runs.
% M : 4 x nPA wedge medians of C -- the analysed profile.
    B   = mean(A.runBeta(:, :, runs), 3, 'omitnan');
    col = cfg.(en).oriIdx - 25 + 8;                 % CONTRASTS.json 26..29 -> cols 9..12
    C   = double(B(:, col)) - double(B(:, 13));     % orientation minus blank
    nP  = numel(cfg.paBins);
    M   = nan(numel(col), nP);
    for p = 1:nP
        m = A.wedge == p;
        if any(m)
            if strcmp(A.agg, 'median'), M(:,p) = median(C(m,:), 1).';
            else,                       M(:,p) = mean(C(m,:), 1).';
            end
        end
    end
end

function M = profile_from_runs(A, runs, cfg)
    [~, M] = halves(A, runs, cfg, A.expn);
end

% ------------------------------------------------------------------------
function a = asym_from_profile(M, cfg, en)
    Asy = compute_asymmetries(M, cfg, cfg.(en));
    a = nan(1,4);
    for j = 1:4, a(j) = mean(Asy.(Asy.order{j}).diff, 'omitnan'); end
end

% ------------------------------------------------------------------------
function splits = balanced_splits(n)
% All disjoint k-vs-k run splits, each counted once. For odd n, drop each run in turn.
    if mod(n,2) == 0
        C = nchoosek(1:n, n/2);
        C = C(1:size(C,1)/2, :);                    % complement pairs counted once
        splits = cell(size(C,1), 2);
        for k = 1:size(C,1)
            splits{k,1} = C(k,:);
            splits{k,2} = setdiff(1:n, C(k,:));
        end
    else
        splits = {};
        for d = 1:n
            keep = setdiff(1:n, d);
            sub  = balanced_splits(n-1);
            for k = 1:size(sub,1)
                splits(end+1, :) = {keep(sub{k,1}), keep(sub{k,2})}; %#ok<AGROW>
            end
        end
    end
end

% ------------------------------------------------------------------------
function r = corr_vec(X, Y)
    x = X(:); y = Y(:);
    g = isfinite(x) & isfinite(y);
    if nnz(g) < 3, r = NaN; return; end
    r = corr(x(g), y(g));
end

function z = fz(r)
    r = max(min(r, 0.9999), -0.9999);
    z = atanh(r);
end

% ------------------------------------------------------------------------
function report(S, cfg, nm, expn, fields)
    bar = repmat('=',1,86);
    labels = {'vertex, full', 'vertex, orientation-differential', ...
              'ROI profile, full (32 values)', 'ROI profile, orientation-differential'};

    fprintf('\n%s\nSPLIT-HALF RELIABILITY OVER RUNS (Spearman-Brown corrected to full run count)\n%s\n', bar, bar);
    fprintf('%-15s %5s %5s %6s   %12s %12s %12s %12s\n', 'observer','nRun','nRun','nVert', ...
            'vertFull','vertDiff','roiFull','roiDiff');
    fprintf('%-15s %5s %5s %6s   %6s %5s %6s %5s %6s %5s %6s %5s\n', '', 'dg','da','dg', ...
            'dg','da','dg','da','dg','da','dg','da');
    fprintf('%s\n', repmat('-',1,86));
    for si = 1:numel(cfg.subjects)
        fprintf('%-15s %5d %5d %6d   ', cfg.subjects{si}, S.nRun(si,1), S.nRun(si,2), S.nVert(si,1));
        for jj = 1:4
            fprintf('%6.3f %5.3f ', S.(fields{jj})(si,1), S.(fields{jj})(si,2));
        end
        fprintf('\n');
    end
    fprintf('%s\n', repmat('-',1,86));
    fprintf('%-15s %5s %5s %6s   ', 'MEDIAN','','','');
    for jj = 1:4
        fprintf('%6.3f %5.3f ', median(S.(fields{jj})(:,1),'omitnan'), ...
                                 median(S.(fields{jj})(:,2),'omitnan'));
    end
    fprintf('\n\n');
    for jj = 1:4
        fprintf('  %-40s dg %.3f [%.3f %.3f]   da %.3f [%.3f %.3f]\n', labels{jj}, ...
            median(S.(fields{jj})(:,1),'omitnan'), min(S.(fields{jj})(:,1)), max(S.(fields{jj})(:,1)), ...
            median(S.(fields{jj})(:,2),'omitnan'), min(S.(fields{jj})(:,2)), max(S.(fields{jj})(:,2)));
    end
    fprintf('\n  splits per observer-experiment: %s\n', mat2str(unique(S.nSplit(isfinite(S.nSplit))).'));

    wtd = ~all(S.scale == 1);
    for ei = 1:2
        fprintf('\n%s\nASYMMETRY RELIABILITY -- %s experiment (gain-weighted: %d)\n%s\n', ...
                bar, upper(expn{ei}), wtd, bar);
        fprintf('%-12s %10s %12s %14s %12s %8s %9s\n', 'asymmetry','group mean', ...
                'SD across obs','within-obs SE','within/total','SNR','p disatt');
        n = size(S.asymW, 1);
        for j = 1:4
            a    = S.asymW(:,j,ei);
            m    = mean(a, 'omitnan');
            sdO  = std(a, 'omitnan');
            seW  = mean(S.asymSEW(:,j,ei), 'omitnan');
            frac = min((seW^2)/(sdO^2), 1);
            sdT  = sqrt(max(sdO^2 - seW^2, 0));
            pT   = 2*tcdf(-abs(m/(sdT/sqrt(n))), n-1);
            fprintf('%-12s %10.3f %12.3f %14.3f %11.0f%% %8.1f %9.4f\n', nm{j}, m, sdO, seW, ...
                    100*frac, abs(m)/seW, pT);
        end
    end

    outdir = fullfile(cfg.reproDir, 'local_qc');
    if ~strcmp(S.aggregator, 'mean'), fprintf('\n(aggregator = %s; not saved)\n', S.aggregator); return; end
    save(fullfile(outdir, 'splithalf_reliability.mat'), '-struct', 'S');
    % Both the raw half-vs-half correlation and the Spearman-Brown corrected value, so
    % the table says which is which (they tell rather different stories: raw ROI-diff
    % medians are 0.77/0.46, corrected 0.87/0.62). raw = r/(2-r) inverts 2r/(1+r).
    raw = @(r) r ./ (2 - r);
    T = table(cfg.subjects(:), S.nRun(:,1), S.nRun(:,2), S.nVert(:,1), ...
              raw(S.vertexFull(:,1)), raw(S.vertexFull(:,2)), ...
              raw(S.vertexDiff(:,1)), raw(S.vertexDiff(:,2)), ...
              raw(S.roiFull(:,1)),    raw(S.roiFull(:,2)), ...
              raw(S.roiDiff(:,1)),    raw(S.roiDiff(:,2)), ...
              S.vertexFull(:,1), S.vertexFull(:,2), S.vertexDiff(:,1), S.vertexDiff(:,2), ...
              S.roiFull(:,1), S.roiFull(:,2), S.roiDiff(:,1), S.roiDiff(:,2), ...
        'VariableNames', {'subject','nRun_dg','nRun_da','nVert_dg', ...
              'rawHalf_vertexFull_dg','rawHalf_vertexFull_da', ...
              'rawHalf_vertexDiff_dg','rawHalf_vertexDiff_da', ...
              'rawHalf_roiFull_dg','rawHalf_roiFull_da', ...
              'rawHalf_roiDiff_dg','rawHalf_roiDiff_da', ...
              'sbCorr_vertexFull_dg','sbCorr_vertexFull_da', ...
              'sbCorr_vertexDiff_dg','sbCorr_vertexDiff_da', ...
              'sbCorr_roiFull_dg','sbCorr_roiFull_da', ...
              'sbCorr_roiDiff_dg','sbCorr_roiDiff_da'});
    writetable(T, fullfile(outdir, 'splithalf_reliability.csv'));
    fprintf('\nwrote %s/splithalf_reliability.{mat,csv}\n', outdir);
end
