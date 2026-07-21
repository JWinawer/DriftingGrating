function res = fit_lme_fig7(M, cfg, expCfg, refBins, doBoot)
% FIT_LME_FIG7  Joint linear mixed-effects model of the four asymmetries (Fig 7).
%
%   res = fit_lme_fig7(M, cfg, expCfg, refBins)
%
% Fits, over 256 data points (nOri=4 x nPA=8 x nSubj=8),
%   y ~ hVv + cVo + rVt + pcVpo + (1|subject)
% where the predictors code each orientation-at-a-location by asymmetry:
%   hVv, rVt   (1st harmonic) : +1 pro / -1 con / 0 not-in-condition
%   cVo, pcVpo (2nd harmonic) : +1 pro / -1 con
% "pro" = horizontal / cardinal / radial / polar-cardinal (sign matches Figs 5/6).
%
% y is the per-(orientation, PA, subject) wedge-median contrast in M.
% refBins (optional) sets the reference polar angle per column (default cfg.paBins).
%
% Returns res with:
%   .beta   (1x4) fixed-effect coefficients [hVv cVo rVt pcVpo]
%   .delta  (1x4) = 2*beta  (pro-minus-con reconstruction, as reported by the paper)
%   .ci68   (4x2) 68% bootstrap CI of delta (resampling subjects)
%   .names  (1x4) asymmetry labels
%   .lme    the fitted model on the full data

    if nargin < 4 || isempty(refBins), refBins = cfg.paBins; end
    if nargin < 5 || isempty(doBoot),  doBoot  = true; end
    res.names = {'hVv','cVo','rVt','pcVpo'};

    X = lme_codes(cfg, expCfg, refBins);          % 4 code matrices, each nOri x nPA
    tbl = build_table(M, X);
    res.lme = fitlme(tbl, 'y ~ hVv + cVo + rVt + pcVpo + (1|subject)');
    res.intercept = res.lme.Coefficients.Estimate(1);
    beta = res.lme.Coefficients.Estimate(2:5).';  % skip intercept
    res.beta  = beta;
    res.delta = 2 * beta;

    % fixed-effect prediction per (orientation, PA) for Fig 8 (no random effect)
    res.pred = res.intercept + beta(1)*X.hVv + beta(2)*X.cVo + beta(3)*X.rVt + beta(4)*X.pcVpo;

    % --- 68% bootstrap CI on delta, resampling the 8 subjects ---
    if ~doBoot, res.ci68 = nan(4,2); return; end
    nS = numel(cfg.subjects);
    rng(0);
    B = cfg.nBoot;
    bootDelta = nan(B, 4);
    for b = 1:B
        sidx = randi(nS, [1 nS]);
        tb = build_table(M(:,:,sidx), X, sidx);   % relabel subjects to keep them distinct
        try
            lb = fitlme(tb, 'y ~ hVv + cVo + rVt + pcVpo + (1|subject)');
            bootDelta(b,:) = 2 * lb.Coefficients.Estimate(2:5).';
        catch
            % singular resample (e.g. all-same subject) -> skip
        end
    end
    res.ci68 = [prctile(bootDelta, 16).' , prctile(bootDelta, 84).'];
end

% ------------------------------------------------------------------------
function tbl = build_table(M, X, subjLabels)
    [nO, nP, nS] = size(M);
    if nargin < 3, subjLabels = 1:nS; end
    y = M(:); % column-major: ori fastest, then PA, then subj
    [kk, pp, ss] = ndgrid(1:nO, 1:nP, 1:nS);
    kk = kk(:); pp = pp(:); ss = ss(:);
    lin = sub2ind([nO nP], kk, pp);
    hVv   = X.hVv(lin);   cVo = X.cVo(lin);
    rVt   = X.rVt(lin);   pcVpo = X.pcVpo(lin);
    subj  = categorical(subjLabels(ss).');
    tbl = table(y, hVv, cVo, rVt, pcVpo, subj, 'VariableNames', ...
                {'y','hVv','cVo','rVt','pcVpo','subject'});
    tbl = tbl(~isnan(tbl.y), :);
end
