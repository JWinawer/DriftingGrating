function [w, info] = harmonic_weights(thetaV, mode, nBins)
% HARMONIC_WEIGHTS  Per-vertex weights giving equal polar-angle coverage.
%
%   [w, info] = harmonic_weights(thetaV)                  % 'equalcoverage', 24 bins
%   [w, info] = harmonic_weights(thetaV, 'natural')       % all ones
%   [w, info] = harmonic_weights(thetaV, 'equalcoverage', 8)   % 45-deg bins
%
% nBins sets the bin count (default 24, i.e. 15 deg). 8 gives 45-deg bins, which are
% the eight polar-angle ROIs themselves -- see the note on bin width below.
%
% WHY THIS EXISTS. The manuscript ROI analysis (BIN_AND_AGGREGATE + FIT_LME_FIG7)
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
% drives that correlation to ~0 and restores the manuscript design's orthogonality, while
% keeping thetaV CONTINUOUS. The only remaining difference from the manuscript analysis is
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

    if nargin < 2 || isempty(mode),   mode   = 'equalcoverage'; end
    if nargin < 3 || isempty(nBins),  nBins  = 24; end

    tv = double(thetaV(:));
    nV = numel(tv);

    % EMPTY BINS ARE A NON-EVENT. Only vertices that exist are indexed below
    % (w(ok) = 1./cnt(ibin(ok)) touches occupied bins only), so an empty bin never
    % divides, never produces a weight and never creates a missing value. The bin
    % width therefore has nothing to do with the missing-data problem: bins are not
    % analysis units, they only partition the vertices that are present.
    %
    % What the width DOES control is leverage. w = 1/count means a bin holding one
    % vertex gives that vertex the same total weight as a bin holding two hundred. At
    % 15 deg the sparsest bin in V1 holds 2 vertices and the weight ratio between
    % individual vertices reaches ~100x; at 45 deg it holds 33 and the ratio is 11x,
    % for a residual r(b1,b3) of 0.075 instead of 0.016 -- a variance inflation of
    % 1.006, i.e. nothing. 45 deg is also the width of the polar-angle ROIs
    % themselves, so it keeps ONE binning in the whole pipeline and makes the
    % estimand exactly "equal weight per ROI", matching the ROI route.
    edges = linspace(0, 360, nBins + 1);
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
