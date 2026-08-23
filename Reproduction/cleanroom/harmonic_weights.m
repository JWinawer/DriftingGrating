function [w, info] = harmonic_weights(thetaV, mode)
% HARMONIC_WEIGHTS  Per-vertex weights giving equal polar-angle coverage.
%
%   [w, info] = harmonic_weights(thetaV)                  % 'equalcoverage'
%   [w, info] = harmonic_weights(thetaV, 'natural')       % all ones
%
% WHY THIS EXISTS. The published ROI analysis (BIN_AND_AGGREGATE + FIT_LME_FIG7)
% aggregates each of the eight polar-angle wedges first and then weights the eight
% wedges EQUALLY. Under that weighting the four asymmetry predictors of LME_CODES are
% EXACTLY orthogonal -- their Gram matrix is diag(16,32,16,32) with zero off-diagonals.
%
% A per-vertex fit weights by actual vertex density instead, and V1 over-represents the
% horizontal meridian. The consequence is not cosmetic: for the four-term harmonic model
% the uncentered correlation between the cos(2*theta) column (b1, horizontal-vs-vertical)
% and the cos(2*(theta-thetaV)) column (b3, radial-vs-tangential) is exactly the weighted
% mean of cos(2*thetaV),
%
%       r(b1,b3) = sum(w.*cosd(2*thetaV)) / sum(w),
%
% which at natural density comes out near +0.35. b1 and b3 are precisely the pair whose
% separation the model exists to adjudicate, so a 0.35 design correlation is an unwanted
% confound imported from cortical magnification.
%
% Re-weighting the vertices so every polar-angle bin contributes the same total weight
% drives that correlation to ~0 and restores the published design's orthogonality, while
% keeping thetaV CONTINUOUS. The only remaining difference from the published analysis is
% then the thing actually under study: each vertex's true pRF polar angle versus its
% 45-deg-wide wedge centre.
%
% THE WEIGHTING. Bin thetaV into 24 bins of 15 deg over [0,360). Each vertex gets
% w = 1/count(its bin), so every occupied bin carries the same total weight regardless of
% how many vertices fall in it. Weights are then rescaled to mean(w) == 1, which leaves
% weighted least squares unchanged but keeps the Gram matrix on the same scale as the
% unweighted one (sum(w) == nVertex), so cond() and VIFs stay comparable across modes.
%
% 15 deg is three bins per 45-deg wedge -- fine enough that the residual correlation is
% ~0.02 rather than 0.35, coarse enough that every bin holds many vertices at the
% per-subject counts here (3.4k-8.1k V1 vertices).
%
% HOW TO APPLY IT. Weighted least squares must minimise sum(w.*(y - X*b).^2), i.e.
%
%       b = (X .* sqrt(w)) \ (y .* sqrt(w));
%
% NOT (X.*w)\(y.*w), which minimises sum(w.^2 .* r.^2) and is a different estimator.
% In the stacked design of HARMONIC_PREDICTORS each vertex contributes nOri rows, all
% carrying the SAME weight, so the vector multiplying the design is
% repmat(sqrt(w), nOri, 1). The same sqrt(w) scaling must be applied when forming R2 and
% the collinearity diagnostics, or the reported design correlation will not describe the
% design that was actually fitted.
%
% Inputs
%   thetaV : nVertex x 1 pRF polar angles in degrees (conventional; any range, mod 360
%            is applied here). Call it PER SUBJECT -- the weights equalise coverage
%            within a subject, matching the per-subject fits of FIT_HARMONIC_VERTEX.
%   mode   : 'equalcoverage' (default) or 'natural' (returns all ones, i.e. the
%            unweighted per-vertex fit, kept available as the secondary specification).
%
% Outputs
%   w    : nVertex x 1 weights, mean 1.
%   info : .mode, .edges (1 x 25), .nBin (24), .count (1 x 24 vertices per bin),
%          .nEmpty (bins with no vertices), .nUnbinned (vertices with non-finite
%          thetaV, which keep weight 1), .effN (Kish effective sample size
%          sum(w)^2/sum(w.^2), <= nVertex; how much the re-weighting costs in
%          precision).

    if nargin < 2 || isempty(mode), mode = 'equalcoverage'; end

    tv = double(thetaV(:));
    nV = numel(tv);

    edges = 0:15:360;
    nBin  = numel(edges) - 1;
    [cnt, ~, ibin] = histcounts(mod(tv, 360), edges);
    cnt = cnt(:);

    switch lower(mode)
        case 'natural'
            w = ones(nV, 1);
        case 'equalcoverage'
            w  = ones(nV, 1);
            ok = ibin > 0;
            w(ok) = 1 ./ cnt(ibin(ok));
            if nV > 0
                w = w / mean(w);            % mean(w) == 1, hence sum(w) == nVertex
            end
        otherwise
            error('harmonic_weights:badMode', ...
                  'mode must be ''equalcoverage'' or ''natural'', got ''%s''.', mode);
    end

    info.mode      = lower(mode);
    info.edges     = edges;
    info.nBin      = nBin;
    info.count     = cnt.';
    info.nEmpty    = nnz(cnt == 0);
    info.nUnbinned = nnz(ibin == 0);
    if nV > 0
        info.effN = sum(w)^2 / sum(w.^2);
    else
        info.effN = 0;
    end
end
