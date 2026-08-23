function A = compute_asymmetries(M, cfg, expCfg, refBins)
% COMPUTE_ASYMMETRIES  The four orientation asymmetries from wedge-median contrasts.
%
%   A = compute_asymmetries(M, cfg, expCfg, refBins)
%
% M      : nOri x nPA x nSubj wedge-median contrasts (see BIN_AND_AGGREGATE).
% expCfg : cfg.dg or cfg.da (uses .oriAngle and .isPolar).
% refBins: (optional) per-column reference polar angle used to compute each
%          stimulus's local orientation/frame. Defaults to cfg.paBins (the true
%          polar angle of each data column). Passing the repo's ordering
%          [90 45 0 315 270 225 180 135] reproduces the original pipeline's
%          derived-direction convention, for A-vs-original comparison.
%
% Each stimulus's local orientation at polar angle theta is
%     loc = mod(oriAngle + isPolar*(theta-90), 180).
% Cartesian-frame label:  loc==0 -> horizontal, loc==90 -> vertical, else oblique.
% Polar-frame label:      d=mod(loc-theta,180): d==0 -> radial, d==90 -> tangential, else polar-oblique.
% This reduces to the direct computation when the asymmetry's frame matches the
% stimulus (e.g. horizontal/vertical for dg, radial/tangential for da) and derives
% the PA-dependent cases otherwise.
%
% Returns A with:
%   A.cat.{H,V,Ob,R,T,PO} : nPA x nSubj category means (per polar angle, per subject)
%   A.<asym>.pro/.con/.diff : nPA x nSubj (diff = pro - con), for asym in
%                             {HV, cardObl, radTan, polcardPolobl}
%   A.<asym>.proName/.conName : labels
%   A.order : cell list of the 4 asymmetry field names (fixed order)

    [nOri, nPA, nSubj] = size(M); %#ok<ASGLU>
    if nargin < 4 || isempty(refBins), refBins = cfg.paBins; end

    isH = false(nOri, nPA); isV = false(nOri, nPA); isOb = false(nOri, nPA);
    isR = false(nOri, nPA); isT = false(nOri, nPA); isPO = false(nOri, nPA);
    for p = 1:nPA
        th  = refBins(p);
        loc = mod(expCfg.oriAngle(:) + expCfg.isPolar*(th - 90), 180);   % nOri x 1
        isH(:,p) = (loc == 0);
        isV(:,p) = (loc == 90);
        isOb(:,p) = ~(isH(:,p) | isV(:,p));
        d = mod(loc - th, 180);
        isR(:,p) = (d == 0);
        isT(:,p) = (d == 90);
        isPO(:,p) = ~(isR(:,p) | isT(:,p));
    end

    H  = catmean(M, isH);   V  = catmean(M, isV);   Ob = catmean(M, isOb);
    R  = catmean(M, isR);   T  = catmean(M, isT);   PO = catmean(M, isPO);
    A.cat = struct('H',H,'V',V,'Ob',Ob,'R',R,'T',T,'PO',PO);

    A.HV            = mkAsym(H,          V,  'horizontal', 'vertical');
    A.cardObl       = mkAsym((H+V)/2,    Ob, 'cardinal',   'oblique');
    A.radTan        = mkAsym(R,          T,  'radial',     'tangential');
    A.polcardPolobl = mkAsym((R+T)/2,    PO, 'polarCardinal','polarOblique');
    A.order = {'HV','cardObl','radTan','polcardPolobl'};
end

function X = catmean(M, mask)
% Mean over stimuli in each per-theta category. mask is nOri x nPA.
    [~, nPA, nSubj] = size(M);
    X = nan(nPA, nSubj);
    for p = 1:nPA
        idx = mask(:, p);
        if any(idx)
            X(p, :) = squeeze(mean(M(idx, p, :), 1)).';
        end
    end
end

function s = mkAsym(pro, con, proName, conName)
    s.pro = pro; s.con = con; s.diff = pro - con;
    s.proName = proName; s.conName = conName;
end
