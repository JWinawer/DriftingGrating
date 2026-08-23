function X = lme_codes(cfg, expCfg, refBins)
% LME_CODES  Per (orientation, polar-angle) asymmetry predictor codes.
%
%   X = lme_codes(cfg, expCfg, refBins)
%
% For each stimulus k and wedge p, computes its local orientation
%   loc = mod(oriAngle(k) + isPolar*(refBins(p)-90), 180)
% and the four asymmetry codes:
%   X.hVv   : +1 horizontal-local, -1 vertical-local,   0 oblique-local  (1st harmonic)
%   X.cVo   : +1 cardinal-local (H or V), -1 oblique-local              (2nd harmonic)
%   X.rVt   : +1 radial-local,   -1 tangential-local,   0 polar-oblique (1st harmonic)
%   X.pcVpo : +1 polar-cardinal-local (R or T), -1 polar-oblique-local  (2nd harmonic)
% Each field is nOri x nPA. refBins defaults to cfg.paBins (true wedge angles).

    if nargin < 3 || isempty(refBins), refBins = cfg.paBins; end
    nO = numel(expCfg.oriAngle);
    nP = numel(cfg.paBins);
    [hVv, cVo, rVt, pcVpo] = deal(zeros(nO, nP));
    for p = 1:nP
        th  = refBins(p);
        loc = mod(expCfg.oriAngle(:) + expCfg.isPolar*(th - 90), 180);
        for k = 1:nO
            L = loc(k);
            isH = (L == 0); isV = (L == 90); isOb = ~(isH || isV);
            d = mod(L - th, 180);
            isR = (d == 0); isT = (d == 90); isPO = ~(isR || isT);
            hVv(k,p)   = isH - isV;
            cVo(k,p)   = (isH || isV) - isOb;
            rVt(k,p)   = isR - isT;
            pcVpo(k,p) = (isR || isT) - isPO;
        end
    end
    X = struct('hVv',hVv,'cVo',cVo,'rVt',rVt,'pcVpo',pcVpo);
end
