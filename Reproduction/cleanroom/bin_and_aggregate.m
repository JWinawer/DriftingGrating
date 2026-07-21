function [M, counts] = bin_and_aggregate(T, cfg, expCfg, doZscore)
% BIN_AND_AGGREGATE  Median orientation contrast within each V1 polar-angle wedge.
%
%   [M, counts] = bin_and_aggregate(T, cfg, expCfg, doZscore)
%
% Applies the analysis inclusion filter (ecc in cfg.eccRange, pRF_r2 > cfg.r2min;
% V1 already enforced by LOAD_AND_FILTER), computes per-vertex contrasts, then takes
% the median across vertices within each polar-angle wedge, separately per subject.
%
% M      : nOri x nPA x nSubj, ordering follows expCfg.oriCols, cfg.paBins, cfg.subjects.
% counts : nPA x nSubj, number of vertices contributing to each wedge/subject.

    keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
    T = T(keep, :);

    C    = compute_vertex_contrasts(T, expCfg, doZscore);
    subj = string(T.subject);
    pab  = T.pRF_angle_bin;

    nO = numel(expCfg.oriCols);
    nP = numel(cfg.paBins);
    nS = numel(cfg.subjects);

    M      = nan(nO, nP, nS);
    counts = zeros(nP, nS);
    for si = 1:nS
        inSubj = subj == cfg.subjects{si};
        for pi = 1:nP
            idx = inSubj & (pab == cfg.paBins(pi));
            counts(pi, si) = nnz(idx);
            if any(idx)
                M(:, pi, si) = median(C(idx, :), 1);
            end
        end
    end
end
