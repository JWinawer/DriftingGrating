function P = precision_weighted_cells(varargin)
% PRECISION_WEIGHTED_CELLS  PRECISION_WEIGHTED_TABLE, made safe against empty cells.
%
%   P = precision_weighted_cells()
%   P = precision_weighted_cells('csv', 'out.csv')
%
% Same estimand, same weighting, same intervals as PRECISION_WEIGHTED_TABLE. The one
% change is what each observer's number is a mean OF.
%
% THE PROBLEM. The published route takes each observer's asymmetry as the mean over
% their eight polar-angle wedges. In V1 every (observer, wedge) cell is populated, so
% that mean is over the same eight wedges for everyone and the route is sound. In the
% extrastriate maps cells go empty -- and not at random: map boundaries lie on the
% meridians, which is exactly the axis the asymmetry varies along. Averaging "whichever
% wedges survived" then means something different for each observer, and the group mean
% mixes estimands. Note this is NOT a broken contrast: the inclusion mask is
% orientation-independent (BIN_AND_AGGREGATE:16), so all four orientations of a cell
% vanish together and the within-wedge subtraction is never partially missing. What
% breaks is only the average over wedges.
%
% THE FIX. Estimate the wedge profile alpha_p from the incomplete design with
% FIT_CELL_META, subtract it, and then average whichever cells are present:
%
%       y_i = mean over present p of ( d_ip - alpha_p )
%
% After removing the polar-angle profile, which wedges an observer happens to have is
% no longer confounded with where the asymmetry is large, so the means are comparable
% again. The alpha_p sum to zero over the wedge set, so WITH COMPLETE DATA y_i is the
% plain wedge mean and this function returns PRECISION_WEIGHTED_TABLE's numbers
% exactly -- V1 does not move. VALIDATE_CELL_META checks that as an identity.
%
% WHY NOT JUST READ mu OFF FIT_CELL_META. The cell-level GLS estimates the same mu and
% recovers it about equally well under cell loss (RMSE 0.016 vs 0.016 in simulation,
% against 0.032 for the published route). But on COMPLETE V1 data its answer moves by
% up to 0.05 depending on how the within-observer sampling covariance is modelled --
% full bootstrap covariance, diagonal only, or compound symmetric all give different
% mu, because GLS re-weights wedges WITHIN an observer using a covariance estimated
% from 8 runs. That is LME.md section 5's own warning ("at n = 8, estimated-weight GLS
% can add variance rather than remove it") showing up in practice. De-trending uses the
% model only for the wedge profile, which is what the missing cells actually require,
% and leaves the weighting at the observer level where it is already validated.
%
% CAVEAT. alpha_p is treated as a known constant when forming y_i, so its own
% uncertainty is not propagated into sigma_i. With complete data that is exact (nothing
% is being adjusted); as cells go missing it makes the intervals slightly too narrow.
%
% Returns P.asym, P.context (as PRECISION_WEIGHTED_TABLE) plus P.alpha, the fitted
% wedge profile, and P.nMissing. Optional 'csv' writes the two tables out.

    p = inputParser;
    p.addParameter('csv', '', @ischar);
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('area', 'V1', @ischar);
    p.addParameter('route', 'roi', @ischar);
    p.addParameter('thetaV', 'binned', @ischar);
    p.addParameter('gain', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('eccRange', [], @(x) isempty(x) || numel(x)==2);
    p.parse(varargin{:});
    opt = p.Results;

    W    = diagnose_within_observer_error('root', opt.root, 'area', opt.area, ...
                                          'eccRange', opt.eccRange, 'route', opt.route, ...
                                          'thetaV', opt.thetaV, 'gain', opt.gain, ...
                                          'quiet', true);
    expn = {'dg','da'};
    nS   = numel(W.subjects);
    nP   = size(W.cell, 2);

    y  = nan(nS, 4, 2);   s2 = nan(nS, 4, 2);
    P.alpha    = nan(nP, 4, 2);
    P.nMissing = squeeze(sum(~isfinite(W.cell(:,:,1,:)), 2));   % nSubj x nExp

    for ei = 1:2
        for j = 1:4
            [d, oi, wj, S, keep] = pack(W, j, ei);
            R = fit_cell_meta(d, oi, wj, S);

            aW = zeros(1, nP);  aW(R.wedgeLevels) = R.alpha;
            P.alpha(:,j,ei) = aW.';

            for si = 1:nS
                kp = find(keep(si,:));
                if isempty(kp), continue; end
                y(si,j,ei)  = mean(reshape(W.cell(si,kp,j,ei), 1, []) - aW(kp));
                Cv = W.cellCov(kp,kp,j,ei,si);
                s2(si,j,ei) = sum(Cv(:)) / numel(kp)^2;
            end
        end
    end

    rows = {};
    for ei = 1:2
        for j = 1:4
            rows(end+1,:) = est(y(:,j,ei), s2(:,j,ei), expn{ei}, W.names{j});  %#ok<AGROW>
        end
    end
    P.asym = mk(rows);

    rows = {};
    for j = 1:4
        rows(end+1,:) = est(y(:,j,1)-y(:,j,2), s2(:,j,1)+s2(:,j,2), 'dg-da', W.names{j}); %#ok<AGROW>
    end
    P.context = mk(rows);

    show('ASYMMETRIES, per experiment', P.asym);
    show('CONTEXT EFFECTS (dg - da)',   P.context);
    if all(P.nMissing(:) == 0)
        fprintf(['\nEvery cell populated, so no wedge profile was subtracted and these are\n' ...
                 'PRECISION_WEIGHTED_TABLE''s numbers exactly.\n']);
    else
        fprintf('\n%d of %d cells empty; wedge profile subtracted before averaging.\n', ...
                sum(P.nMissing(:)), nS*nP*2);
    end

    if ~isempty(opt.csv)
        writetable([P.asym; P.context], opt.csv);
        fprintf('wrote %s\n', opt.csv);
    end
end

% ------------------------------------------------------------------------
function [d, oi, wj, S, keep] = pack(W, j, ei)
% Flatten the populated cells, with the block-diagonal within-observer covariance.
    nS   = size(W.cell, 1);
    keep = isfinite(squeeze(W.cell(:,:,j,ei)));
    d = []; oi = []; wj = []; blocks = {};
    for si = 1:nS
        kp = find(keep(si,:));
        if isempty(kp), continue; end
        d  = [d;  reshape(W.cell(si,kp,j,ei), [], 1)];     %#ok<AGROW>
        oi = [oi; repmat(si, numel(kp), 1)];               %#ok<AGROW>
        wj = [wj; kp(:)];                                  %#ok<AGROW>
        blocks{end+1} = W.cellCov(kp,kp,j,ei,si);          %#ok<AGROW>
    end
    S = blkdiag(blocks{:});  S = (S + S.')/2;
end

% ------------------------------------------------------------------------
function r = est(y, s2, expLabel, name)
% Identical to PRECISION_WEIGHTED_TABLE's estimator, on observers that have data.
    ok = isfinite(y) & isfinite(s2);
    y  = y(ok);  s2 = s2(ok);  n = numel(y);

    tau2 = max(0, var(y,0) - mean(s2));
    w    = 1./(tau2 + s2);  w = w/sum(w);
    tcrit = tinv(0.975, n-1);

    mEq = mean(y);     seEq = std(y)/sqrt(n);
    mW  = sum(w.*y);   seW  = sqrt(sum(w.^2 .* (tau2 + s2)));
    [~, pEq] = ttest(y);
    pW  = 2*(1 - tcdf(abs(mW/seW), n-1));

    r = {expLabel, name, n, mEq, mEq-tcrit*seEq, mEq+tcrit*seEq, pEq, ...
         mW, mW-tcrit*seW, mW+tcrit*seW, pW, sqrt(tau2), sqrt(mean(s2)), max(w)/min(w)};
end

function T = mk(rows)
    T = cell2table(rows, 'VariableNames', {'experiment','asymmetry','n', ...
        'equal','equal_lo','equal_hi','equal_p', ...
        'weighted','weighted_lo','weighted_hi','weighted_p', ...
        'tau','mean_sigma','weight_ratio'});
end

function show(ttl, T)
    fprintf('\n%s\n%s\n', ttl, repmat('=',1,92));
    fprintf('%-6s %-11s %3s %22s %22s %7s %7s %6s\n', ...
            'exp','asymmetry','n','equal-weighted [95% CI]','precision-w [95% CI]', ...
            'tau','sigma','w max/min');
    for i = 1:height(T)
        fprintf('%-6s %-11s %3d %7.3f [%6.3f %6.3f] %7.3f [%6.3f %6.3f] %7.3f %7.3f %6.2fx\n', ...
            T.experiment{i}, T.asymmetry{i}, T.n(i), ...
            T.equal(i), T.equal_lo(i), T.equal_hi(i), ...
            T.weighted(i), T.weighted_lo(i), T.weighted_hi(i), ...
            T.tau(i), T.mean_sigma(i), T.weight_ratio(i));
    end
end
