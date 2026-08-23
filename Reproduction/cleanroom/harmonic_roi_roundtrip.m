function [A, M] = harmonic_roi_roundtrip(Y, D, cfg, expCfg, aggName)
% HARMONIC_ROI_ROUNDTRIP  Push per-vertex responses through the published ROI pipeline.
%
%   [A, M] = harmonic_roi_roundtrip(Y, D, cfg, expCfg, aggName)
%
% Aggregates per-vertex responses into the 4 x 8 x nSubj wedge array and computes the
% four asymmetries with the UNMODIFIED COMPUTE_ASYMMETRIES, so model predictions and
% observations come out in exactly the units and layout of the manuscript's numbers.
% Wedges are taken from D.tvBin (pRF_angle_bin), the same binning as
% BIN_AND_AGGREGATE:29.
%
% aggName : 'mean' (default) or 'median'.
%   Use 'mean' when comparing observed with model-predicted. The four asymmetries are
%   zero-sum contrasts across the four conditions, so with a linear aggregator the
%   per-vertex demeaning cancels exactly and the comparison is clean. The median is
%   NOT linear, so demeaned-median asymmetries do not equal blank-subtracted-median
%   asymmetries; 'median' is provided to reproduce the Figs 5/6 convention, and should
%   be read against the manuscript rather than against the model.
%
% Returns A (see COMPUTE_ASYMMETRIES) and M, nOri x nPA x nSubj.

    if nargin < 5 || isempty(aggName), aggName = 'mean'; end
    switch lower(aggName)
        case 'mean',   agg = @(x) mean(x, 1);
        case 'median', agg = @(x) median(x, 1);
        otherwise, error('harmonic_roi_roundtrip:agg', 'aggName must be mean or median.');
    end

    nO = size(Y, 2);
    nP = numel(cfg.paBins);
    nS = numel(cfg.subjects);

    M = nan(nO, nP, nS);
    for si = 1:nS
        inSubj = D.subj == si;
        for pi = 1:nP
            idx = inSubj & (D.tvBin == cfg.paBins(pi));
            if any(idx)
                M(:, pi, si) = agg(Y(idx, :)).';
            end
        end
    end

    A = compute_asymmetries(M, cfg, expCfg);
end
