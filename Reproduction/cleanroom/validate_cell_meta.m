function V = validate_cell_meta(varargin)
% VALIDATE_CELL_META  Does the cell-level machinery do what it claims, on V1?
%
%   V = validate_cell_meta()
%   V = validate_cell_meta('nRep', 200)
%
% V1 is the control: every (observer, wedge) cell is populated (minimum 22 vertices),
% so the published observer-level route in PRECISION_WEIGHTED_TABLE is valid there and
% anything new must agree with it before being trusted on the extrastriate maps, where
% cells go empty and no such check is available.
%
% TEST 1 -- exact reduction of the fitter. Collapse each observer to one cell (d = the
%   wedge average, S = that observer's bootstrap variance, tau2cell = 0, tau2obs fixed
%   at the method-of-moments value) and FIT_CELL_META must return the published
%   weighted mean and its SE as an algebraic identity. Tests the GLS, the weighting and
%   the interval, and nothing else. Passes at ~1e-17.
%
% TEST 2 -- complete-data agreement of the two candidate estimators.
%   PRECISION_WEIGHTED_CELLS must reproduce the published table EXACTLY, because with
%   no missing cells the fitted wedge profile sums to zero over all eight wedges and
%   cancels out of every observer mean. It does, at ~1e-16.
%
%   Reading mu straight off the cell-level GLS does NOT reproduce it, and the size of
%   the discrepancy depends on how the within-observer sampling covariance is modelled:
%
%       asymmetry        published   full S    diag S   compound-symmetric S
%       dg horiz-vert      -0.549    -0.505    -0.535        -0.549
%       dg card-obl        -0.213    -0.171    -0.204        -0.213
%       dg rad-tang         0.112     0.124     0.132         0.111
%
%   Only the compound-symmetric case recovers the published route, because that is the
%   one case where the intercept separates from the wedge effects and GLS reduces to
%   weighting observer means. With any heterogeneous S the GLS re-weights wedges WITHIN
%   an observer using a covariance estimated from 8 runs, and moves mu by up to 0.05 --
%   larger than several of the effects being estimated. This is LME.md section 5's own
%   caveat, that estimated-weight GLS at n = 8 can add variance rather than remove it,
%   appearing in practice. It is why PRECISION_WEIGHTED_CELLS uses the fit only for the
%   wedge profile and keeps the weighting at the observer level.
%
% TEST 3 -- recovery under cell loss, which is the whole point. Delete cells from V1 in
%   the pattern the extrastriate maps show (loss on the meridian wedges, where the map
%   boundaries are), and measure how far each route moves from its own full-data answer.
%   The published route inherits the deletion pattern; both model-based routes largely
%   do not.
%
% Returns V with .test1, .test2, .test3 tables. Prints all three.

    p = inputParser;
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('nRep', 200, @isnumeric);
    p.addParameter('seed', 7, @isnumeric);
    p.parse(varargin{:});
    opt = p.Results;

    W    = diagnose_within_observer_error('root', opt.root, 'quiet', true);
    cfg  = config_repro();
    expn = {'dg','da'};
    nS   = numel(W.subjects);
    nP   = size(W.cell, 2);

    % ================================================================= TEST 1
    rows = {};
    for ei = 1:2
        for j = 1:4
            yv = W.full(:,j,ei);  s2 = W.seBoot(:,j,ei).^2;
            tau2 = max(0, var(yv,0) - mean(s2));
            R = fit_cell_meta(yv, (1:nS).', ones(nS,1), s2, ...
                              'tau2', [tau2 0], 'wedgeEffects', false);
            w  = 1./(tau2 + s2);  w = w/sum(w);
            mW = sum(w.*yv);      seW = sqrt(sum(w.^2 .* (tau2 + s2)));
            rows(end+1,:) = {expn{ei}, W.names{j}, mW, R.mu, abs(R.mu-mW), ...
                             abs(R.se-seW)};                       %#ok<AGROW>
        end
    end
    V.test1 = cell2table(rows, 'VariableNames', ...
        {'exp','asymmetry','closedForm','fitCellMeta','muAbsDiff','seAbsDiff'});
    banner('TEST 1  exact reduction of FIT_CELL_META to the published estimator');
    disp(V.test1);
    fprintf('max |mu diff| = %.3e   max |se diff| = %.3e\n', ...
            max(V.test1.muAbsDiff), max(V.test1.seAbsDiff));

    % ================================================================= TEST 2
    Pub  = precision_weighted_table('root', opt.root);
    Cell = precision_weighted_cells('root', opt.root);
    dAsym = max(abs(Cell.asym.weighted    - Pub.asym.weighted));
    dCtx  = max(abs(Cell.context.weighted - Pub.context.weighted));
    dCiA  = max(abs(Cell.asym.weighted_hi - Pub.asym.weighted_hi));

    V.test2 = table(dAsym, dCtx, dCiA, 'VariableNames', ...
        {'maxDiffAsym','maxDiffContext','maxDiffCI'});
    banner('TEST 2  PRECISION_WEIGHTED_CELLS reproduces the published table on V1');
    fprintf('max |asymmetry diff| = %.3e\nmax |context diff|   = %.3e\nmax |CI diff|        = %.3e\n', ...
            dAsym, dCtx, dCiA);
    assert(max([dAsym dCtx dCiA]) < 1e-12, 'validate_cell_meta:test2', ...
           'complete-data reduction is not exact -- the wedge profile is not cancelling');

    % ================================================================= TEST 3
    merid = find(ismember(cfg.paBins, [0 90 180 270]));
    rng(opt.seed);
    rows = {};
    for ei = 1:2
        for j = 1:4
            [dF, oF, wF, SF] = pack(W, j, ei, true(nS,nP));
            refGls = fit_cell_meta(dF, oF, wF, SF).mu;
            yF = W.full(:,j,ei);  s2F = W.seBoot(:,j,ei).^2;
            refPub = wmean(yF, s2F);

            eOld = nan(opt.nRep,1); eGls = nan(opt.nRep,1); eDet = nan(opt.nRep,1);
            for r = 1:opt.nRep
                keep = true(nS, nP);
                for si = 1:nS
                    keep(si, merid(randperm(numel(merid), randi([1 3])))) = false;
                end

                yk = nan(nS,1); s2k = nan(nS,1);
                for si = 1:nS
                    kp = find(keep(si,:));
                    yk(si)  = mean(reshape(W.cell(si,kp,j,ei), 1, []));
                    Cv = W.cellCov(kp,kp,j,ei,si);
                    s2k(si) = sum(Cv(:))/numel(kp)^2;
                end
                eOld(r) = wmean(yk, s2k) - refPub;             % the current route

                [d2, o2, w2, S2] = pack(W, j, ei, keep);
                Rg = fit_cell_meta(d2, o2, w2, S2);
                eGls(r) = Rg.mu - refGls;                      % cell-level GLS

                aW = zeros(1,nP);  aW(Rg.wedgeLevels) = Rg.alpha;
                yd = nan(nS,1);
                for si = 1:nS
                    kp = find(keep(si,:));
                    yd(si) = mean(reshape(W.cell(si,kp,j,ei), 1, []) - aW(kp));
                end
                eDet(r) = wmean(yd, s2k) - refPub;             % de-trend, as shipped
            end

            rows(end+1,:) = {expn{ei}, W.names{j}, refPub, ...
                             mean(eOld), rms_(eOld), mean(eGls), rms_(eGls), ...
                             mean(eDet), rms_(eDet)};           %#ok<AGROW>
        end
    end
    V.test3 = cell2table(rows, 'VariableNames', ...
        {'exp','asymmetry','fullData','biasOld','rmseOld', ...
         'biasGls','rmseGls','biasDetrend','rmseDetrend'});
    banner(sprintf('TEST 3  recovery under meridian cell loss (%d draws)', opt.nRep));
    disp(V.test3);
    fprintf(['old = wedge mean over surviving cells (the current route)\n' ...
             'gls = mu straight off FIT_CELL_META\n' ...
             'det = de-trend then average (PRECISION_WEIGHTED_CELLS, as shipped)\n\n' ...
             'mean RMSE   old %.4f   gls %.4f   detrend %.4f\n'], ...
            mean(V.test3.rmseOld), mean(V.test3.rmseGls), mean(V.test3.rmseDetrend));
end

% ------------------------------------------------------------------------
function [d, oi, wj, S] = pack(W, j, ei, keep)
    nS = size(W.cell,1);
    d = []; oi = []; wj = []; blocks = {};
    for si = 1:nS
        kp = find(keep(si,:) & isfinite(reshape(W.cell(si,:,j,ei),1,[])));
        if isempty(kp), continue; end
        d  = [d;  reshape(W.cell(si,kp,j,ei), [], 1)];   %#ok<AGROW>
        oi = [oi; repmat(si, numel(kp), 1)];             %#ok<AGROW>
        wj = [wj; kp(:)];                                %#ok<AGROW>
        blocks{end+1} = W.cellCov(kp,kp,j,ei,si);        %#ok<AGROW>
    end
    S = blkdiag(blocks{:});  S = (S + S.')/2;
end

function m = wmean(y, s2)
    tau2 = max(0, var(y,0) - mean(s2));
    w = 1./(tau2 + s2);  m = sum(w.*y)/sum(w);
end

function r = rms_(x), r = sqrt(mean(x.^2)); end

function banner(t)
    fprintf('\n%s\n%s\n%s\n', repmat('=',1,78), t, repmat('=',1,78));
end
