function G = build_group_matrices_fromCSV(cfg, expCfg, doZscore, T)
% BUILD_GROUP_MATRICES_FROMCSV  Bridge: CSV -> the group arrays stage-04 expects.
%
%   G = build_group_matrices_fromCSV(cfg, expCfg, doZscore, T)
%
% Reproduces what stages 01+03 (main_singlesub + meanWithinLabel) would have
% produced, but starting from allsubjectsTable.csv. Places the four
% orientation-minus-blank contrasts at their CONTRASTS.json indices
% (expCfg.oriIdx, i.e. 26-29) so the existing analysis functions index them
% correctly. Rows 1-25 are left NaN (not needed for Figures 5-8).
%
% G fields (matching the .mat variables the repo loads):
%   medianBOLDpa, meanBOLDpa : nC x nPA x 1 x nSubj   (nC=29, ROI=V1)
%   medianBOLD,   meanBOLD   : nC x 1 x nSubj
% Ordering: contrasts by CONTRASTS.json index; PA by cfg.paBins
% [0 45 90 135 180 225 270 315] (= meanWithinLabel order); subjects by cfg.subjects.

    if nargin < 4 || isempty(T), T = load_and_filter(cfg); end
    keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
    T = T(keep, :);

    C    = compute_vertex_contrasts(T, expCfg, doZscore);   % nVertex x 4, oriCols order
    subj = string(T.subject);
    pab  = T.pRF_angle_bin;

    nC = 29; nP = numel(cfg.paBins); nS = numel(cfg.subjects);
    medianBOLDpa = nan(nC, nP, 1, nS); meanBOLDpa = nan(nC, nP, 1, nS);
    medianBOLD   = nan(nC, 1, nS);     meanBOLD   = nan(nC, 1, nS);

    for s = 1:nS
        inS = subj == cfg.subjects{s};
        for k = 1:numel(expCfg.oriIdx)
            ci = expCfg.oriIdx(k);
            medianBOLD(ci,1,s) = median(C(inS, k));
            meanBOLD(ci,1,s)   = mean(C(inS, k));
        end
        for p = 1:nP
            idx = inS & (pab == cfg.paBins(p));
            if any(idx)
                for k = 1:numel(expCfg.oriIdx)
                    ci = expCfg.oriIdx(k);
                    medianBOLDpa(ci,p,1,s) = median(C(idx, k));
                    meanBOLDpa(ci,p,1,s)   = mean(C(idx, k));
                end
            end
        end
    end

    G.medianBOLDpa = medianBOLDpa; G.meanBOLDpa = meanBOLDpa;
    G.medianBOLD   = medianBOLD;   G.meanBOLD   = meanBOLD;
    G.contrastIdx  = expCfg.oriIdx;
end
