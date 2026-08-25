function R = diagnose_loss_structure(varargin)
% DIAGNOSE_LOSS_STRUCTURE  Is the shift from cell loss systematic, and is loss the
% reason MT is not reported by polar angle?
%
%   R = diagnose_loss_structure()
%   R = diagnose_loss_structure('nRand',200,'nSub',200)
%
% TWO QUESTIONS, both raised by JW on 2026-08-25 after DIAGNOSE_CELL_LOSS.
%
% PART A -- SYSTEMATIC OR IDIOSYNCRATIC? Deleting MT's empty cells from V1 moved some
% asymmetries up and some down, which looks like noise. It is not: an extrastriate
% map's loss has a COMMON part and an IDIOSYNCRATIC part, and only the common part can
% bias the group. MT at 4-8 deg empties the vertical meridian (90 and 270 deg) for all
% eight observers; the other 14 empty cells are scattered over particular observers.
% This holds the vertical meridian fixed and RANDOMISES the scattered part, two ways:
%
%   permute-observers   the same 14-cell loss profile, reassigned to observers at random
%   permute-and-rois    each observer keeps a count, but the ROIs are redrawn from the
%                       six non-vertical ones
%
% If the delta under randomisation centres on the vertical-meridian-only delta, the
% shift is systematic and the scatter is just which observer happened to be missing.
% If it centres on zero, the shift really was idiosyncratic.
%
% PART B -- WHY MT IS NOT REPORTED. DIAGNOSE_CELL_LOSS deleted MT's empty cells from V1
% and left V1's density intact in every cell that survived, so it tested the HOLES
% only. That is one of the two things wrong with MT and not the bigger one: V1 with
% MT's holes still holds ~7000 vertices, MT holds ~545, a median of ONE per cell. This
% subsamples V1 down to each map's actual per-cell counts -- holes and sparsity
% together -- and reports the spread of the estimate over repeated draws. The spread is
% a LOWER bound on the real map's, since V1's vertices are the better-fit ones.
%
% Run over several maps at once, the two parts say which side of the section 7
% criterion each map falls on and why.
%
% Returns R.partA, R.partB, R.pattern.

    p = inputParser;
    p.addParameter('area', 'V1', @ischar);
    p.addParameter('donor', 'MT', @ischar);
    p.addParameter('eccRange', [4 8], @(x) numel(x)==2);
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('nRand', 200, @isscalar);
    p.addParameter('nSub', 200, @isscalar);
    p.addParameter('subMaps', {'V3','V3a','hV4','MT'}, @iscell);
    p.addParameter('seed', 23, @isscalar);
    p.parse(varargin{:});
    opt = p.Results;

    cfg = config_repro();  cfg.eccRange = opt.eccRange;
    nS  = numel(cfg.subjects);  nP = numel(cfg.paBins);
    nm  = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    lbl = [strcat({'dg '},nm), strcat({'da '},nm), strcat({'dg-da '},nm)];

    P0 = spec_pooled('area', opt.area, 'eccRange', opt.eccRange, 'root', opt.root);
    D  = P0.data;

    donorN = cell_occupancy('area', opt.donor, 'eccRange', opt.eccRange, 'root', opt.root);
    drop   = donorN == 0;

    % Split the donor's loss into the part every observer shares and the rest.
    common = all(drop, 1);                       % 1 x nPA, true where ALL are empty
    extraK = sum(drop & ~common, 2);             % per observer, the idiosyncratic count
    freeP  = find(~common);                      % ROIs still available to lose

    R.pattern = struct('drop', drop, 'common', common, 'extraK', extraK);

    fprintf('\n%s\nPART A -- IS THE SHIFT SYSTEMATIC?  %s, %g-%g deg, donor %s\n%s\n', ...
            repmat('=',1,100), opt.area, opt.eccRange(1), opt.eccRange(2), opt.donor, ...
            repmat('=',1,100));
    fprintf(['%s empties %s deg for ALL %d observers (%d cells). A further %d cells are ' ...
             'scattered\nover particular observers (%s). The scattered part is what gets ' ...
             'randomised below.\n'], opt.donor, ...
            strjoin(arrayfun(@(x) sprintf('%g',x), cfg.paBins(common), 'uni', 0), '/'), ...
            nS, nS*nnz(common), sum(extraK), ...
            strjoin(arrayfun(@(k) sprintf('%d',k), extraK.', 'uni', 0), ' '));

    % --- the three deterministic conditions -----------------------------------
    vertOnly = repmat(common, nS, 1);
    dActual  = delta_of(D, cfg, drop,     opt);
    dVert    = delta_of(D, cfg, vertOnly, opt);

    % --- the randomisations ----------------------------------------------------
    rs = RandStream('twister', 'Seed', opt.seed);
    modes = {'permute-observers','permute-and-rois'};
    DR = nan(opt.nRand, numel(modes), 12);
    for mi = 1:numel(modes)
        for b = 1:opt.nRand
            dr = vertOnly;
            switch modes{mi}
                case 'permute-observers'
                    % Reassign the exact per-observer loss profiles to other observers.
                    perm = randperm(rs, nS);
                    src  = drop & ~common;
                    dr   = dr | src(perm, :);
                case 'permute-and-rois'
                    % Keep the counts (permuted), redraw which non-vertical ROIs.
                    k = extraK(randperm(rs, nS));
                    for si = 1:nS
                        if k(si) > 0
                            dr(si, freeP(randperm(rs, numel(freeP), k(si)))) = true;
                        end
                    end
            end
            DR(b,mi,:) = delta_of(D, cfg, dr, opt);
        end
    end

    rows = {};
    fprintf(['\n%-16s %9s %9s   %9s %9s %9s   %9s %9s\n'], 'quantity', 'MT actual', ...
            'vert only', 'perm-obs', 'sd', '90%% range', 'perm-roi', 'sd');
    for j = 1:12
        a = squeeze(DR(:,1,j));  b = squeeze(DR(:,2,j));
        q = prctile(a, [5 95]);
        fprintf('%-16s %9.3f %9.3f   %9.3f %9.3f  [%5.3f %5.3f]  %9.3f %9.3f\n', ...
                lbl{j}, dActual(j), dVert(j), mean(a,'omitnan'), std(a,'omitnan'), ...
                q(1), q(2), mean(b,'omitnan'), std(b,'omitnan'));
        rows(end+1,:) = {lbl{j}, dActual(j), dVert(j), mean(a,'omitnan'), ...
                         std(a,'omitnan'), q(1), q(2), mean(b,'omitnan'), ...
                         std(b,'omitnan')}; %#ok<AGROW>
    end
    R.partA = cell2table(rows, 'VariableNames', {'quantity','delta_actual', ...
        'delta_vertical_only','permobs_mean','permobs_sd','permobs_p5','permobs_p95', ...
        'permroi_mean','permroi_sd'});
    fprintf(['\nThe "vert only" column is the loss EVERY observer shares. Where the ' ...
             'randomised means\nsit on it, the shift is systematic and the scattered ' ...
             'cells only add scatter.\n']);

    % --- PART B ---------------------------------------------------------------
    R.partB = density_floor(D, cfg, opt, nm);
end

% ------------------------------------------------------------------------
function d = delta_of(D, cfg, drop, opt)
% The specification's estimate under a deletion, minus its estimate with nothing
% deleted, over the SAME observers. Randomised patterns lose a different observer on
% each draw, so the baseline has to be recomputed per draw or the delta would carry the
% observer swap as well as the cell loss.
    nS = numel(cfg.subjects);
    baseDrop = false(nS, numel(cfg.paBins));
    baseDrop(all(drop,2), :) = true;
    d = est_of(D, cfg, drop, opt) - est_of(D, cfg, baseDrop, opt);
end

function v = est_of(D, cfg, drop, opt)
% The settled specification: fit per observer, average equally. Context effects formed
% WITHIN observer first (standing fact 6). SPEC_POOLED returns the per-observer fits,
% so this is the specification's estimator and not a re-implementation of it.
    A = spec_pooled('data', D, 'area', opt.area, 'eccRange', opt.eccRange, ...
                    'dropCells', drop, 'jackknife', false);
    ok = all(isfinite(A.dg.asymObs),2) & all(isfinite(A.da.asymObs),2);
    if nnz(ok) < 3, v = nan(1,12); return; end
    dgO = A.dg.asymObs(ok,:);  daO = A.da.asymObs(ok,:);
    v = [mean(dgO,1), mean(daO,1), mean(dgO-daO,1)];
end

% ------------------------------------------------------------------------
function T = density_floor(D, cfg, opt, nm)
% PART B. Reproduce each map's WHOLE coverage profile in V1 -- the empty cells and the
% vertex counts in the cells that survive -- and see how much the estimate moves from
% draw to draw. This is the quantity the section 7 criterion is really about.
    nS = numel(cfg.subjects);
    truth = est_of(D, cfg, false(nS, numel(cfg.paBins)), opt);
    fprintf('\n%s\nPART B -- HOW MUCH DATA IS LEFT?  V1 subsampled to each map''s coverage\n%s\n', ...
            repmat('=',1,100), repmat('=',1,100));
    fprintf('%-8s %7s %8s %8s %8s   %s\n', 'map', 'empty', 'median', 'vertices', 'obs', ...
            'dg rad-tang / dg horiz-vert:  mean [5th, 95th] over draws');
    rows = {};
    for mi = 1:numel(opt.subMaps)
        m = opt.subMaps{mi};
        N = cell_occupancy('area', m, 'eccRange', opt.eccRange, 'root', opt.root);
        rs = RandStream('twister', 'Seed', 101 + mi);
        E  = nan(opt.nSub, 12);
        for b = 1:opt.nSub
            Ds = subsample_cells(D, cfg, N, rs);
            E(b,:) = est_of(Ds, cfg, false(nS, numel(cfg.paBins)), opt);
        end
        nObsAlive = nnz(sum(N,2) > 0);
        fprintf('%-8s %7d %8g %8d %8d   rt %6.3f [%6.3f %6.3f]   hv %6.3f [%6.3f %6.3f]\n', ...
                m, nnz(N==0), median(N(:)), sum(N(:)), nObsAlive, ...
                mean(E(:,3),'omitnan'), prctile(E(:,3),5), prctile(E(:,3),95), ...
                mean(E(:,1),'omitnan'), prctile(E(:,1),5), prctile(E(:,1),95));
        for j = 1:4
            rows(end+1,:) = {m, nnz(N==0), median(N(:)), sum(N(:)), nObsAlive, ...
                             nm{j}, truth(j), mean(E(:,j),'omitnan'), ...
                             std(E(:,j),'omitnan'), prctile(E(:,j),5), ...
                             prctile(E(:,j),95)}; %#ok<AGROW>
        end
    end
    fprintf(['\nV1 itself, undiluted: dg rad-tang %.3f, dg horiz-vert %.3f, from %d ' ...
             'vertices.\n'], truth(3), truth(1), sum(cellfun(@(a) size(a.Y,1), D.dg)));
    fprintf(['A map whose 90%% band spans more than the effect being measured cannot ' ...
             'resolve it,\nhowever the cells are combined. That is what the section 7 ' ...
             'criterion screens for --\nnot the holes, which PART A shows the ' ...
             'specification largely absorbs.\n']);
    T = cell2table(rows, 'VariableNames', {'map','empty_cells','median_per_cell', ...
        'total_vertices','observers','asymmetry','v1_truth','sub_mean','sub_sd', ...
        'sub_p5','sub_p95'});
end
