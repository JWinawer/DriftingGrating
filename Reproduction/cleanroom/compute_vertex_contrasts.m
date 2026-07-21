function C = compute_vertex_contrasts(T, expCfg, doZscore)
% COMPUTE_VERTEX_CONTRASTS  Per-vertex orientation-minus-blank contrasts.
%
%   C = compute_vertex_contrasts(T, expCfg, doZscore)
%
% T        : table of vertices (V1 subset), see LOAD_AND_FILTER.
% expCfg   : cfg.dg or cfg.da (defines oriCols, blank, betaStd).
% doZscore : if true, divide each contrast by the per-vertex beta std
%            (reproduces the original per-vertex z-scoring, because the
%            per-vertex mean cancels in the orientation-minus-blank subtraction).
%
% Returns C, nVertex x 4, one column per stationary orientation in expCfg.oriCols
% order, each = (orientation beta - blank beta), optionally / beta_std.

    blank = T.(expCfg.blank);
    nOri  = numel(expCfg.oriCols);
    C = zeros(height(T), nOri);
    for k = 1:nOri
        C(:,k) = T.(expCfg.oriCols{k}) - blank;
    end
    if doZscore
        C = C ./ T.(expCfg.betaStd);
    end
end
