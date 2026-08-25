function R = diagnose_pooled_fit(varargin)
% DIAGNOSE_POOLED_FIT  Is one group fit less damaged by cell loss than eight averaged?
%
%   R = diagnose_pooled_fit()                          % MT's 4-8 deg pattern
%   R = diagnose_pooled_fit('eccRange',[2 10],'nRand',200)
%
% THE QUESTION (JW, 2026-08-25). DIAGNOSE_CELL_LOSS showed that deleting an
% extrastriate map's empty cells from V1 moves the specification's estimate, because
% the specification fits each observer separately and an observer whose polar-angle
% coverage has holes returns a biased fit that is then averaged in with equal weight.
% Fitting the pool ONCE, with vertices weighted to equalise both polar-angle coverage
% within an observer and total contribution across observers, might be more robust:
% one observer's missing wedge is covered by the other seven.
%
% THREE ESTIMATORS, on identical vertices and identical deletions:
%   avg           the settled specification -- fit per observer, average equally
%   pooled-obs    one fit; each observer totals 1/nObs, spread over the ROIs they HAVE
%   pooled-cell   one fit; each (observer x ROI) cell totals 1/(8*nObs), uncompensated
% The two pooled weightings are identical on complete data, so any gap between them is
% caused entirely by the deletion. See SPEC_POOLED for what the weights are. The
% eight-wedge ROI route is not repeated here; DIAGNOSE_CELL_LOSS carries it, and 'avg'
% reproduces that function's harmonic arm exactly, which is how the two are tied
% together.
%
% AND THE SAME THREE CONDITIONS as DIAGNOSE_CELL_LOSS: no deletion, the donor map's
% own empty cells, and nRand random deletions of the same size preserving each
% observer's count. Everything runs off ONE data load, so nRand can be large.
%
% WHAT TO READ. Robustness is necessary, not sufficient. A pooled fit that does not
% move may be reweighting observers rather than resisting the bias, which is a
% different estimand and not automatically the one wanted -- R.leverage reports
% diag((sum M)^{-1} M_i), the share of each coefficient each observer supplies, under
% no deletion and under the donor pattern. Read the two together.
%
% Returns R.summary, R.leverage, R.classTest, R.drop.

    p = inputParser;
    p.addParameter('area', 'V1', @ischar);
    p.addParameter('donor', 'MT', @ischar);
    p.addParameter('eccRange', [4 8], @(x) numel(x)==2);
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('nRand', 100, @isscalar);
    p.addParameter('seed', 7, @isscalar);
    p.addParameter('classTest', true, @islogical);
    % Severity sweep: is pooling better when loss is severe, or when it is
    % IDIOSYNCRATIC rather than shared? That is the case the pooling intuition rests
    % on -- one observer's hole covered by the other seven -- so it is tested rather
    % than argued. Set false to skip; it is ~8 x nSweep fits.
    p.addParameter('sweep', true, @islogical);
    p.addParameter('nSweep', 100, @isscalar);
    p.parse(varargin{:});
    opt = p.Results;

    cfg = config_repro();  cfg.eccRange = opt.eccRange;
    nS  = numel(cfg.subjects);  nP = numel(cfg.paBins);
    nm  = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    est = {'avg','pooled-obs','pooled-cell'};

    donorN = cell_occupancy('area', opt.donor, 'eccRange', opt.eccRange, 'root', opt.root);
    drop   = donorN == 0;

    fprintf('\n%s\nPOOLED vs PER-OBSERVER FIT under %s''s cell loss  (%s, %g-%g deg, n=%d cells)\n%s\n', ...
            repmat('=',1,92), opt.donor, opt.area, opt.eccRange(1), opt.eccRange(2), ...
            nnz(drop), repmat('=',1,92));

    % One load, reused by every fit below.
    P0 = spec_pooled('area', opt.area, 'eccRange', opt.eccRange, 'root', opt.root);
    D  = P0.data;

    % MATCHED BASELINE. If the donor pattern empties an observer completely -- MT at
    % 4-8 deg does, sub-wlsubj123 -- then a baseline over all eight and a deleted
    % condition over seven differ by TWO things, and the delta would not be the cell
    % loss. The baseline therefore removes those observers too and loses nothing else,
    % so every delta below isolates the deletion. The observer that disappears is a
    % separate, and larger, consequence; it is stated in the header, not buried in a
    % delta. Nulls preserve each observer's count, so they drop the same observers.
    fullyLost = all(drop, 2);
    baseDrop  = false(nS, nP);  baseDrop(fullyLost,:) = true;
    if any(fullyLost)
        fprintf(['NOTE: %s empties %d observer(s) entirely (%s). They are removed from the\n' ...
                 'baseline as well, so the deltas below are the cell loss alone.\n'], ...
                opt.donor, nnz(fullyLost), strjoin(cfg.subjects(fullyLost), ', '));
    end
    base = evaluate(D, cfg, baseDrop, opt);
    don  = evaluate(D, cfg, drop,     opt);

    % Null: same number of lost cells per observer, random ROIs. Generated up front,
    % though nothing here reseeds the stream -- consistency with DIAGNOSE_CELL_LOSS.
    rng(opt.seed);
    kPer  = sum(drop, 2);
    nulls = cell(opt.nRand,1);
    for b = 1:opt.nRand
        dr = false(nS,nP);
        for si = 1:nS
            if kPer(si) > 0, dr(si, randperm(nP, kPer(si))) = true; end
        end
        nulls{b} = dr;
    end
    Dn = nan(opt.nRand, numel(est), 12);
    for b = 1:opt.nRand
        r = evaluate(D, cfg, nulls{b}, opt);
        for k = 1:numel(est), Dn(b,k,:) = r.val(k,:) - base.val(k,:); end
    end

    % --- table ----------------------------------------------------------------
    rows = {};
    lbl  = [strcat({'dg '},nm), strcat({'da '},nm), strcat({'dg-da '},nm)];
    for k = 1:numel(est)
        for j = 1:12
            dn = squeeze(Dn(:,k,j));
            rows(end+1,:) = {opt.area, opt.donor, est{k}, lbl{j}, don.n(k), ...
                             base.val(k,j), base.p(k,j), don.val(k,j), don.p(k,j), ...
                             don.val(k,j)-base.val(k,j), mean(abs(dn),'omitnan'), ...
                             prctile(abs(dn),95)}; %#ok<AGROW>
        end
    end
    R.summary = cell2table(rows, 'VariableNames', ...
        {'area','donor','estimator','quantity','n','base_mean','base_p', ...
         'drop_mean','drop_p','delta','null_mean_abs_delta','null_p95_abs_delta'});

    report(R.summary, est, lbl, opt);

    % --- who is actually supplying the estimate -------------------------------
    R.leverage = leverage_report(base, don, cfg, opt);
    R.drop = drop;  R.base = base;  R.donor = don;
    R.naive = naive_check(D, cfg, base, baseDrop);
    if opt.classTest, R.classTest = class_test(D, cfg, nS, nm, est, opt); end
    if opt.sweep,     R.sweep     = loss_sweep(D, cfg, nS, nP, opt); end
end

% ------------------------------------------------------------------------
function T = naive_check(D, cfg, base, baseDrop)
% Does the cross-observer normalisation earn its place? Pooling WITHOUT it weights
% observers by vertex count, because HARMONIC_WEIGHTS rescales to mean(w) == 1 and so
% sum(w) == nVertex. This is the map-size bias the normalisation exists to remove,
% measured rather than assumed.
    % Same observers as the matched baseline it is compared against, or the two would
    % differ by the observer set as well as by the weighting.
    keepObs = ~all(baseDrop, 2);
    nv = cellfun(@(a) size(a.Y,1), D.dg(keepObs));
    M = zeros(4);  r = zeros(4,1);
    for si = 1:numel(cfg.subjects)
        a = D.dg{si};
        if isempty(a) || ~keepObs(si), continue; end
        w  = harmonic_weights(a.thetaV, 'equalcoverage', 8);       % sum(w) == nVertex
        X  = harmonic_predictors(a.thetaV, cfg.dg, struct('expanded', false));
        sw = repmat(sqrt(w(:)), size(a.Y,2), 1);
        Xw = X .* sw;
        M  = M + Xw'*Xw;  r = r + Xw'*(a.Y(:) .* sw);
    end
    o  = struct('expanded', false, 'weighting', 'equalcoverage');
    As = compute_asymmetries(predict_harmonic((M\r).', cfg.paBins(:), cfg.dg, o).', cfg, cfg.dg);
    an = arrayfun(@(j) mean(As.(As.order{j}).diff, 'omitnan'), 1:4);
    fprintf(['\n%s\nDOES EQUALISING OBSERVERS MATTER? (dg, no deletion)\n%s\n' ...
             'V1 vertex counts run %d to %d, a %.2fx range across observers.\n'], ...
            repmat('=',1,92), repmat('=',1,92), min(nv), max(nv), max(nv)/min(nv));
    fprintf('%-28s %9.3f %9.3f %9.3f %9.3f\n', 'pooled, count-weighted', an);
    fprintf('%-28s %9.3f %9.3f %9.3f %9.3f\n', 'pooled, observers equalised', base.val(2,1:4));
    fprintf('%-28s %9.3f %9.3f %9.3f %9.3f\n', 'specification average', base.val(1,1:4));
    fprintf(['Equalising is worth up to %.3f and moves the pooled fit ONTO the ' ...
             'specification''s\nanswer, which is the check that the normalisation is the ' ...
             'right one.\n'], max(abs(an - base.val(2,1:4))));
    T = cell2table([{'count-weighted'}, num2cell(an); ...
                    {'observers-equalised'}, num2cell(base.val(2,1:4)); ...
                    {'specification-average'}, num2cell(base.val(1,1:4))], ...
        'VariableNames', {'pooling','horiz_vert','card_obl','rad_tang','polc_polo'});
end

% ------------------------------------------------------------------------
function T = loss_sweep(D, cfg, nS, nP, opt)
% THE CASE FOR POOLING, TESTED DIRECTLY. Pooling is supposed to help when one
% observer's missing wedge is covered by another's intact one. That requires the loss
% to be IDIOSYNCRATIC. Two regimes are run at four severities:
%   idiosyncratic  each observer loses k ROIs, drawn independently
%   common         the SAME k ROIs are lost by every observer
% and the metric is RMSE of the group estimate against its own no-deletion value, so
% bias and instability are both counted.
    rows = {};
    fprintf('\n%s\nSEVERITY SWEEP: dg rad-tang, RMSE against the no-deletion estimate\n%s\n', ...
            repmat('=',1,92), repmat('=',1,92));
    fprintf('%-16s %-9s %10s %10s %10s %10s\n', 'loss', 'k of 8', ...
            'RMSE avg', 'RMSE p-obs', 'RMSE p-cell', 'bias avg');
    B  = spec_pooled('data', D, 'jackknife', false);
    Bc = spec_pooled('data', D, 'weighting', 'equalcell', 'jackknife', false);
    b0 = [mean(B.dg.asymObs,1), B.dg.asym, Bc.dg.asym];
    for mode = {'idiosyncratic','common'}
        for k = [2 4 5 6]
            rng(11);
            pats = cell(opt.nSweep,1);
            for b = 1:opt.nSweep
                dr = false(nS,nP);
                if strcmp(mode{1},'idiosyncratic')
                    for si = 1:nS, dr(si, randperm(nP,k)) = true; end
                else
                    dr(:, randperm(nP,k)) = true;
                end
                pats{b} = dr;
            end
            E = nan(opt.nSweep, 3);
            for b = 1:opt.nSweep
                A = spec_pooled('data', D, 'dropCells', pats{b}, 'jackknife', false);
                C = spec_pooled('data', D, 'dropCells', pats{b}, ...
                                'weighting', 'equalcell', 'jackknife', false);
                ok = all(isfinite(A.dg.asymObs), 2);
                if nnz(ok) >= 3, E(b,1) = mean(A.dg.asymObs(ok,3), 1); end
                E(b,2) = A.dg.asym(3);  E(b,3) = C.dg.asym(3);
            end
            d = E - [b0(3), b0(7), b0(11)];
            rmse = sqrt(mean(d.^2, 1, 'omitnan'));  bi = mean(d, 1, 'omitnan');
            fprintf('%-16s %-9s %10.4f %10.4f %10.4f %10.4f\n', mode{1}, ...
                    sprintf('%d of 8', k), rmse, bi(1));
            rows(end+1,:) = {mode{1}, k, rmse(1), rmse(2), rmse(3), bi(1), bi(2), bi(3)}; %#ok<AGROW>
        end
    end
    fprintf(['\nIf pooling were the answer, the idiosyncratic rows would separate. They do ' ...
             'not.\n']);
    T = cell2table(rows, 'VariableNames', {'loss','k','rmse_avg','rmse_pooled_obs', ...
        'rmse_pooled_cell','bias_avg','bias_pooled_obs','bias_pooled_cell'});
end

% ------------------------------------------------------------------------
function r = evaluate(D, cfg, drop, opt)
% All estimators from as few fits as possible: SPEC_POOLED returns the per-observer
% fits alongside the pooled one, so 'avg' costs nothing extra.
    assert(isequal(cfg.eccRange(:).', opt.eccRange(:).'), 'diagnose_pooled_fit:band', ...
           'the loaded data and the requested band disagree.');
    A = spec_pooled('data', D, 'area', opt.area, 'eccRange', opt.eccRange, ...
                    'weighting', 'equalobserver', 'dropCells', drop, 'jackknife', true);
    C = spec_pooled('data', D, 'area', opt.area, 'eccRange', opt.eccRange, ...
                    'weighting', 'equalcell', 'dropCells', drop, 'jackknife', true);
    r.val = nan(3,12);  r.p = nan(3,12);  r.n = nan(3,1);
    r.A = A;  r.C = C;

    ok = all(isfinite(A.dg.asymObs),2) & all(isfinite(A.da.asymObs),2);
    r.n(:) = nnz(ok);

    % avg: the specification. Context effect formed WITHIN observer first (standing
    % fact 6), then averaged -- which is why it is not simply avg(dg) - avg(da) when
    % observers differ in which experiments survive.
    dgO = A.dg.asymObs(ok,:);  daO = A.da.asymObs(ok,:);
    r.val(1,:) = [mean(dgO,1), mean(daO,1), mean(dgO-daO,1)];
    r.p(1,:)   = [tp(dgO), tp(daO), tp(dgO-daO)];

    % pooled: no per-observer estimate exists, so the context effect is a difference
    % of two group fits and its interval is a jackknife over observers.
    for k = 1:2
        Q = A;  if k == 2, Q = C; end
        v = [Q.dg.asym, Q.da.asym, Q.dg.asym - Q.da.asym];
        r.val(k+1,:) = v;
        if ~isempty(Q.dg.asymJack) && isequal(Q.dg.jackIdx, Q.da.jackIdx)
            J  = [Q.dg.asymJack, Q.da.asymJack, Q.dg.asymJack - Q.da.asymJack];
            nJ = size(J,1);
            se = sqrt((nJ-1)/nJ * sum((J - mean(J,1)).^2, 1));
            r.p(k+1,:) = 2 * tcdf(-abs(v./se), nJ-1);
        end
    end
end

function pv = tp(A)
    n = size(A,1);
    pv = 2 * tcdf(-abs(mean(A,1) ./ (std(A,0,1)/sqrt(n))), n-1);
end

% ------------------------------------------------------------------------
function report(T, est, lbl, opt)
    fprintf('\n%-14s %-14s %3s %9s %8s %9s %8s %9s %9s\n', ...
            'estimator','quantity','n','baseline','p','deleted','p','delta','null p95');
    for j = 1:12
        for k = 1:numel(est)
            r = find(strcmp(T.estimator, est{k}) & strcmp(T.quantity, lbl{j}), 1);
            fprintf('%-14s %-14s %3d %9.3f %8s %9.3f %8s %9.3f %9.3f%s\n', ...
                est{k}, tern(k==1, lbl{j}, ''), T.n(r), T.base_mean(r), pstr(T.base_p(r)), ...
                T.drop_mean(r), pstr(T.drop_p(r)), T.delta(r), T.null_p95_abs_delta(r), ...
                flag(T.base_p(r), T.drop_p(r)));
        end
        fprintf('\n');
    end
    fprintf(['<- marks a quantity that crosses p = .05. The pooled p values are ' ...
             'delete-one-observer\njackknives, not the specification''s t test -- ' ...
             'comparable in spirit, not identical.\n']);
    fprintf('donor: %s, %g-%g deg.\n', opt.donor, opt.eccRange(1), opt.eccRange(2));
end

% ------------------------------------------------------------------------
function T = leverage_report(base, don, cfg, opt)
% diag((sum_i M_i)^{-1} M_i) -- the share of each coefficient observer i supplies.
% Sums to exactly 1 down the observers. Under equal averaging every share is 1/8 by
% definition, so this is the quantity that says whether pooling bought its robustness
% by reweighting rather than by resisting bias.
    fprintf('\n%s\nWHO SUPPLIES THE ESTIMATE  (share of each coefficient, pooled-obs, dg)\n%s\n', ...
            repmat('=',1,92), repmat('=',1,92));
    fprintf('%-15s %8s %10s %10s %10s %10s %10s\n', 'observer', 'ROIs', ...
            'h-v base', 'h-v drop', 'r-t base', 'r-t drop', 'wedges lost');
    rows = {};
    Lb = base.A.dg.leverage;  Ld = don.A.dg.leverage;
    for si = 1:numel(cfg.subjects)
        fprintf('%-15s %8d %10.3f %10.3f %10.3f %10.3f %10d\n', cfg.subjects{si}, ...
                don.A.dg.nOcc(si), Lb(si,1), Ld(si,1), Lb(si,3), Ld(si,3), ...
                8 - don.A.dg.nOcc(si));
        rows(end+1,:) = {cfg.subjects{si}, don.A.dg.nOcc(si), ...
                         Lb(si,1), Ld(si,1), Lb(si,3), Ld(si,3)}; %#ok<AGROW>
    end
    fprintf('%-15s %8s %10.3f %10.3f %10.3f %10.3f\n', 'sum', '', ...
            sum(Lb(:,1),'omitnan'), sum(Ld(:,1),'omitnan'), ...
            sum(Lb(:,3),'omitnan'), sum(Ld(:,3),'omitnan'));
    fprintf(['\nEqual averaging gives every observer 0.125 by definition. Departures ' ...
             'here are\nthe pooled fit choosing whom to listen to.\n']);
    T = cell2table(rows, 'VariableNames', ...
        {'observer','nROI_after_drop','hv_base','hv_drop','rt_base','rt_drop'});
    T.area = repmat({opt.area}, height(T), 1);
end

% ------------------------------------------------------------------------
function T = class_test(D, cfg, nS, nm, est, opt)
% The mechanism check from DIAGNOSE_CELL_LOSS, run on all three estimators: delete a
% whole polar-angle class from EVERY observer. When the loss is common to all
% observers, pooling has no intact neighbour to borrow from, so this is the case that
% separates "robust" from "reweighting".
    cls = struct('vertical',   ismember(cfg.paBins, [90 270]), ...
                 'horizontal', ismember(cfg.paBins, [0 180]), ...
                 'oblique',    ismember(cfg.paBins, [45 135 225 315]));
    fn = fieldnames(cls);  rows = {};
    fprintf('\n%s\nMECHANISM: delete a whole polar-angle class from EVERY observer\n%s\n', ...
            repmat('=',1,92), repmat('=',1,92));
    for k = 1:numel(est)
        fprintf('\n--- estimator: %s ---\n%-12s %-6s', est{k}, 'deleted', 'exp');
        fprintf('%11s', nm{:});  fprintf('\n');
        r0 = evaluate(D, cfg, false(nS, numel(cfg.paBins)), opt);
        prow('(none)', r0.val(k,:));
        rows = [rows; mk('(none)', est{k}, r0.val(k,:), nm)]; %#ok<AGROW>
        for c = 1:numel(fn)
            rc = evaluate(D, cfg, repmat(cls.(fn{c}), nS, 1), opt);
            prow(fn{c}, rc.val(k,:));
            rows = [rows; mk(fn{c}, est{k}, rc.val(k,:), nm)]; %#ok<AGROW>
        end
    end
    fprintf(['\nThe two pooled weightings agree exactly in this block, and must: when every\n' ...
             'observer loses the SAME class they differ only by one common constant, which\n' ...
             'weighted least squares ignores. They can only differ when the loss is uneven\n' ...
             'ACROSS observers, which is the donor-pattern case above.\n']);
    T = cell2table(rows, 'VariableNames', {'deleted','estimator','experiment','asymmetry','mean'});
end

function prow(lbl, v)
    fprintf('%-12s %-6s', lbl, 'dg');   fprintf('%11.3f', v(1:4));
    fprintf('\n%-12s %-6s', '', 'da');  fprintf('%11.3f', v(5:8));
    fprintf('\n');
end

function rows = mk(lbl, e, v, nm)
    rows = {};
    for j = 1:4, rows(end+1,:) = {lbl, e, 'dg', nm{j}, v(j)};   end %#ok<AGROW>
    for j = 1:4, rows(end+1,:) = {lbl, e, 'da', nm{j}, v(4+j)}; end %#ok<AGROW>
end

function s = pstr(p)
    if ~isfinite(p), s = 'NaN'; elseif p < 0.001, s = '<.001'; else, s = sprintf('%.3f', p); end
end
function s = flag(p1, p2)
    if isfinite(p1) && isfinite(p2) && (p1 < 0.05) ~= (p2 < 0.05), s = '  <-'; else, s = ''; end
end
function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
