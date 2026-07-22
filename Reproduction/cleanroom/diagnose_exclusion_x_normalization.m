% DIAGNOSE_EXCLUSION_X_NORMALIZATION  Normalisation and exclusion are one decision.
%
% Accepting that observers should be brought into commensurate units has a consequence:
% an observer whose gain cannot be measured cannot be normalised. In the polar
% experiment, blank-referenced gain is <= 0 for sub-0037 (-0.03) and sub-0201 (-0.40) --
% their V1 did not respond to the stimuli at all. Standard-deviation divisors paper over
% this by returning a positive number regardless (the noise level, in place of a gain).
%
% So the question is not "normalise or not" but "normalise, having first removed the
% observers for whom normalisation is undefined". This script crosses the two decisions
% and reports whether the Fig 7B ordering is stable once they are made together.
%
% See ../ZSCORE_FIG7.md section 6 and section 8.

cfg = config_repro();
T   = load_allconditions(cfg);
keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
T = T(keep,:);
subjV = string(T.subject); pabV = T.pRF_angle_bin;
nS = numel(cfg.subjects); nP = numel(cfg.paBins);
[~, daMot] = allcondition_cols();
E = cfg.da;

S4 = zeros(height(T),4); M8 = zeros(height(T),8);
for k=1:4, S4(:,k) = T.(E.oriCols{k}); end
for k=1:8, M8(:,k) = T.(daMot{k}); end
blank = T.(E.blank);
X = S4 - blank;

g = struct('std13',zeros(1,nS),'std8mot',zeros(1,nS),'mean8mot',zeros(1,nS),'maxResp',zeros(1,nS));
for si = 1:nS
    i = subjV == cfg.subjects{si};
    g.std13(si)    = median(T.(E.betaStd)(i));
    g.std8mot(si)  = median(std(M8(i,:),0,2));
    g.mean8mot(si) = median(mean(M8(i,:),2) - blank(i));
    g.maxResp(si)  = median(max([S4(i,:) M8(i,:)] - blank(i), [], 2));
end

excluded  = {'sub-0037','sub-0201'};
keepSubj  = ~ismember(cfg.subjects, excluded);
normNames = {'raw (none)','std13 (published)','std8mot','mean8mot','maxResp (ephys-style)'};
normVals  = {ones(1,nS), g.std13, g.std8mot, g.mean8mot, g.maxResp};

fprintf('\nExcluded from the "n=6" rows: %s\n', strjoin(excluded, ', '));
fprintf('(both have non-positive blank-referenced gain in the polar experiment)\n');

for setIdx = 1:2
    if setIdx == 1, use = true(1,nS);  lbl = 'ALL 8 OBSERVERS';
    else,           use = keepSubj;    lbl = 'n = 6 (gain estimable for every observer)'; end

    fprintf('\n================ %s ================\n', lbl);
    fprintf('%-24s %9s %9s %9s %10s   %s\n', 'normaliser','H-V','radTan','ratio','P(rT>|HV|)','ordering');

    for d = 1:numel(normNames)
        gv = normVals{d};
        if any(gv(use) <= 0)
            fprintf('%-24s %9s %9s %9s %10s   INVALID: divisor <= 0\n', normNames{d},'-','-','-','-');
            continue
        end
        Xs = X;
        for si = 1:nS
            i = subjV == cfg.subjects{si};
            Xs(i,:) = X(i,:) / gv(si);
        end
        A = compute_asymmetries(wm(Xs, subjV, pabV, cfg), cfg, E);
        rt = mean(A.radTan.diff,1,'omitnan'); hv = mean(A.HV.diff,1,'omitnan');
        rt = rt(use); hv = hv(use); n = numel(rt);

        rng(0); B = 10000; win = 0;
        for b = 1:B
            k = randi(n,[1 n]);
            win = win + (mean(rt(k)) > abs(mean(hv(k))));
        end
        ord = 'H-V larger'; if mean(rt) > abs(mean(hv)), ord = 'radTan larger'; end
        fprintf('%-24s %9.3f %9.3f %9.2f %10.2f   %s\n', normNames{d}, ...
            mean(hv), mean(rt), mean(rt)/abs(mean(hv)), win/B, ord);
    end
end

fprintf(['\nRead the P column as the evidence for the ordering, not the point estimate:\n' ...
         '0.5 is a coin flip, and anything inside roughly 0.2-0.8 means the dataset does\n' ...
         'not determine which polar-grating asymmetry is larger.\n']);

function M = wm(C,subjV,pabV,cfg)
    nS=numel(cfg.subjects); nP=numel(cfg.paBins); M=nan(4,nP,nS);
    for si=1:nS
        i0=subjV==cfg.subjects{si};
        for pi=1:nP, idx=i0&(pabV==cfg.paBins(pi)); if any(idx), M(:,pi,si)=median(C(idx,:),1); end, end
    end
end
