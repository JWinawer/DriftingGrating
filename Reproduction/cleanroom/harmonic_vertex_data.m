function D = harmonic_vertex_data(T, cfg, expCfg, doZscore)
% HARMONIC_VERTEX_DATA  Per-vertex demeaned orientation responses + geometry.
%
%   D = harmonic_vertex_data(T, cfg, expCfg, doZscore)
%
% Applies the manuscript inclusion filter (V1 already enforced by LOAD_AND_FILTER /
% LOAD_ALLCONDITIONS; here ecc in cfg.eccRange and pRF_r2 > cfg.r2min -- the exact
% line from BIN_AND_AGGREGATE:13, so the vertex set matches Figs 5-8), then returns
% each vertex's four stationary responses with the vertex mean removed.
%
% Removing the mean across the four orientations also removes the blank, which is
% why the pink-noise-baseline problem (local_qc/REPORT.md section 1) does not touch
% this model: only orientation DIFFERENCES survive. All four harmonic predictors
% likewise sum to zero across the four conditions, for both experiments, so the
% model needs no intercept and demeaning does not distort the design.
%
% Inputs
%   T        : per-vertex V1 table from LOAD_AND_FILTER or LOAD_ALLCONDITIONS.
%   expCfg   : cfg.dg or cfg.da.
%   doZscore : divide by the per-vertex beta std (sensitivity variant only; the raw
%              analysis is the one adopted -- local_qc/REPORT.md section 4).
%
% Output D
%   .Y       nVertex x nOri, demeaned per-vertex responses (the regression target).
%   .tvCont  nVertex x 1, continuous pRF polar angle  (Fit B).
%   .tvBin   nVertex x 1, wedge-centre polar angle    (Fit A).
%   .subj    nVertex x 1, index into cfg.subjects.
%   .rowIdx  nVertex x 1, row index into the FILTERED table (for later joins).
%   .T       the filtered table itself.
%   .nOri, .expName

    if nargin < 4, doZscore = false; end

    keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & ...
           T.pRF_r2 > cfg.r2min;
    T = T(keep, :);

    C = compute_vertex_contrasts(T, expCfg, doZscore);   % nVertex x nOri
    Y = C - mean(C, 2);                                  % remove each vertex's mean

    tvCont = double(T.pRF_angle);
    tvBin  = double(T.pRF_angle_bin);

    % Subject index in cfg.subjects order.
    subjStr = string(T.subject);
    subj    = zeros(height(T), 1);
    for si = 1:numel(cfg.subjects)
        subj(subjStr == cfg.subjects{si}) = si;
    end

    ok = all(isfinite(Y), 2) & isfinite(tvCont) & isfinite(tvBin) & subj > 0;
    if ~all(ok)
        fprintf('harmonic_vertex_data(%s): dropping %d of %d vertices (non-finite or unknown subject)\n', ...
                expCfg.name, nnz(~ok), numel(ok));
    end

    D.Y       = Y(ok, :);
    D.tvCont  = tvCont(ok);
    D.tvBin   = tvBin(ok);
    D.subj    = subj(ok);
    D.rowIdx  = find(ok);
    D.T       = T(ok, :);
    D.nOri    = size(Y, 2);
    D.expName = expCfg.name;
end
