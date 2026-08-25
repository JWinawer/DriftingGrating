function [Dsub, capped] = subsample_cells(D, cfg, keepCounts, rs)
% SUBSAMPLE_CELLS  Thin a loaded dataset down to a target vertex count per cell.
%
%   [Dsub, capped] = subsample_cells(D, cfg, keepCounts, rs)
%
% D is a SPEC_POOLED data struct. keepCounts is nSubj x nPA, the number of vertices to
% keep in each (observer x polar-angle ROI) cell; a cell with 0 is emptied, and a count
% at or above what the cell holds keeps all of it. rs is a RandStream, so the draw is
% controllable without touching the global stream the fitting path reseeds.
%
% WHY THIS EXISTS, alongside the dropCells option. Deleting cells reproduces an
% extrastriate map's HOLES while leaving V1's density intact everywhere else, which
% answers "does a structured gap in polar-angle coverage bias the fit". It does not
% answer "is there enough data left to fit at all" -- V1 with MT's holes still holds
% thousands of vertices, and MT holds a few hundred. This reproduces the map's whole
% coverage profile, holes and sparsity together, so the two failure modes can be told
% apart instead of being confounded.
%
% The subsample is WITHOUT replacement and uniform within a cell. It therefore
% reproduces the count but not the reliability: V1's vertices have better pRF fits and
% larger responses than an extrastriate map's, so whatever spread this produces is a
% LOWER BOUND on the real map's.
%
% `capped` counts cells where the target exceeded what was available, so a target that
% could not be met is visible rather than silent.

    nS = numel(cfg.subjects);  nP = numel(cfg.paBins);
    assert(isequal(size(keepCounts), [nS nP]), 'subsample_cells:size', ...
           'keepCounts must be %d x %d.', nS, nP);
    Dsub = D;  capped = 0;

    % ONE mask per observer, applied to BOTH experiments. The inclusion filter is
    % built from eccentricity and pRF vexpl, neither of which depends on the
    % experiment, so dg and da index the same vertices in the same order -- asserted
    % here, because drawing the two independently would break that and
    % SPEC_AREAS_SUMMARY's nVert assertion downstream.
    for si = 1:nS
        A = D.dg{si};  B = D.da{si};
        if isempty(A), continue; end
        assert(~isempty(B) && isequal(A.wedge, B.wedge), 'subsample_cells:pair', ...
               '%s: dg and da index different vertices.', cfg.subjects{si});
        keep = false(numel(A.wedge), 1);
        for pIdx = 1:nP
            idx  = find(A.wedge == pIdx);
            want = keepCounts(si, pIdx);
            if want <= 0 || isempty(idx), continue; end
            if want >= numel(idx)
                keep(idx) = true;
                if want > numel(idx), capped = capped + 1; end
            else
                keep(idx(randperm(rs, numel(idx), want))) = true;
            end
        end
        for en = {'dg','da'}
            e = en{1};
            Dsub.(e){si} = struct('Y', D.(e){si}.Y(keep,:), ...
                                  'thetaV', D.(e){si}.thetaV(keep), ...
                                  'wedge',  D.(e){si}.wedge(keep));
        end
    end
end
