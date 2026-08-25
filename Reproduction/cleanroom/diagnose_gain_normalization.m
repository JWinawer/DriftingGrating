% DIAGNOSE_GAIN_NORMALIZATION  Is there a stable per-observer BOLD gain to normalise by?
%
% Percent BOLD is not a well-defined biological quantity -- it depends on pulse sequence,
% field strength and physiology -- so raw averaging silently weights observers by an
% arbitrary gain, and normalising per observer (as single-unit work does when it scales
% each neuron's peak response to 1) is a reasonable response. That argument, however,
% presupposes something checkable: that each observer HAS a gain, stable enough to be
% estimated from data other than the effect of interest.
%
% This script tests that presupposition three ways:
%   1. Do different gain definitions agree with each other, within an experiment?
%   2. Does an observer's gain TRANSFER across the two experiments? A gain that is a
%      property of the person should; one that is a property of the session should not.
%   3. Is the gain estimable at all for every observer?
%
% It also demonstrates the identity that makes the "weighting" framing slippery:
% averaging normalised per-subject values with EQUAL weight is arithmetically identical
% to a 1/gain-weighted average of the raw values. Normalisation and reweighting are the
% same operation seen from two sides, so calling one of them unfair proves nothing.

cfg = config_repro();
T   = load_allconditions(cfg);
keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
T = T(keep,:);
subjV = string(T.subject); pabV = T.pRF_angle_bin;
nS = numel(cfg.subjects); nP = numel(cfg.paBins);
[dgMot, daMot] = allcondition_cols();
asymNames = {'HV','cardObl','radTan','polcardPolobl'};

G = struct();   % gain estimates, per experiment
for E = {cfg.dg, cfg.da}
    E = E{1}; %#ok<FXSET>
    if strcmp(E.name,'dg'), mot = dgMot; else, mot = daMot; end
    S4 = zeros(height(T),4); M8 = zeros(height(T),8);
    for k=1:4, S4(:,k) = T.(E.oriCols{k}); end
    for k=1:8, M8(:,k) = T.(mot{k}); end
    blank = T.(E.blank);
    X = S4 - blank;                      % the analysed contrasts

    g = struct();
    g.std13     = zeros(1,nS);   % manuscript divisor: std over all 13 conditions
    g.std8mot   = zeros(1,nS);   % std over the 8 motion conditions (effect-independent)
    g.mean8mot  = zeros(1,nS);   % mean motion drive over blank (effect-independent)
    g.maxResp   = zeros(1,nS);   % ephys-style: peak stimulus response over blank
    for si = 1:nS
        i = subjV == cfg.subjects{si};
        g.std13(si)    = median(T.(E.betaStd)(i));
        g.std8mot(si)  = median(std(M8(i,:),0,2));
        g.mean8mot(si) = median(mean(M8(i,:),2) - blank(i));
        g.maxResp(si)  = median(max([S4(i,:) M8(i,:)] - blank(i), [], 2));
    end
    G.(E.name) = g;
    G.(E.name).X = X;
end

%% 1. do the gain definitions agree within an experiment?
fprintf('\n==== 1. per-observer gain estimates ====\n');
for E = {'dg','da'}
    g = G.(E{1});
    fprintf('\n%s:\n%-14s %9s %9s %10s %9s\n', E{1}, 'subject','std13','std8mot','mean8mot','maxResp');
    for si=1:nS
        fprintf('%-14s %9.3f %9.3f %10.3f %9.3f\n', cfg.subjects{si}, ...
            g.std13(si), g.std8mot(si), g.mean8mot(si), g.maxResp(si));
    end
    fprintf('  corr(std8mot, mean8mot) = %+.3f ; corr(std8mot, maxResp) = %+.3f\n', ...
        corr(g.std8mot', g.mean8mot'), corr(g.std8mot', g.maxResp'));
    neg = cfg.subjects(g.mean8mot <= 0.05);
    if ~isempty(neg)
        fprintf('  *** gain NOT ESTIMABLE (mean motion drive <= 0.05) for: %s\n', strjoin(neg,', '));
    end
end

%% 2. does gain transfer across experiments?
fprintf('\n==== 2. is gain a property of the OBSERVER or of the SESSION? ====\n');
fprintf('If gain were a stable individual trait, these should correlate strongly.\n');
fprintf('%-12s %10s %10s\n','estimate','Pearson','Spearman');
for f = {'std13','std8mot','mean8mot','maxResp'}
    a = G.dg.(f{1})'; b = G.da.(f{1})';
    fprintf('%-12s %10.3f %10.3f\n', f{1}, corr(a,b), corr(a,b,'type','Spearman'));
end

%% 3. what each normalisation does to the polar-grating result
fprintf('\n==== 3. group asymmetries under per-observer gain normalisation (da) ====\n');
E = cfg.da; X = G.da.X;
schemes = {'raw (no normalisation)', ones(1,nS); ...
           'std13 (manuscript)',       G.da.std13; ...
           'std8mot',                 G.da.std8mot; ...
           'mean8mot',                G.da.mean8mot; ...
           'maxResp (ephys-style)',   G.da.maxResp};
fprintf('%-24s %9s %9s %9s %9s   %s\n','normaliser',asymNames{:},'ordering');
perSubj = struct();
for d = 1:size(schemes,1)
    Xs = X;
    for si=1:nS
        i = subjV==cfg.subjects{si};
        Xs(i,:) = X(i,:) / schemes{d,2}(si);
    end
    A = compute_asymmetries(wm(Xs,subjV,pabV,cfg), cfg, E);
    v = cellfun(@(n) mean(A.(n).diff(:),'omitnan'), asymNames);
    ord = 'H-V larger'; if v(3)>abs(v(1)), ord='radTan larger'; end
    bad = '';
    if any(schemes{d,2} <= 0), bad = '  <== INVALID: divisor <= 0 for some observer'; end
    fprintf('%-24s %9.3f %9.3f %9.3f %9.3f   %s%s\n', schemes{d,1}, v, ord, bad);
    perSubj.(sprintf('s%d',d)) = mean(A.radTan.diff,1,'omitnan');
end
fprintf(['\nA divisor <= 0 sign-flips that observer, so those rows are uninterpretable, not\n' ...
         'merely different. Note WHICH divisors stay positive: std-based ones always do,\n' ...
         'because a standard deviation is positive by construction -- including for an\n' ...
         'observer whose V1 did not respond, where it returns the NOISE level in place of a\n' ...
         'gain and then divides by it. The failure is hidden, not absent.\n']);

%% 4. normalisation IS equal weighting, in the normalised units
fprintf('\n==== 4. the identity that makes "reweighting" a matter of viewpoint ====\n');
rt_raw = perSubj.s1; w = 1./G.da.std13;
eqNorm  = mean(rt_raw ./ G.da.std13);        % equal weight, normalised units
wtdRaw  = sum(w.*rt_raw)/sum(w);             % 1/gain weighted, raw units
fprintf('equal-weighted mean of the normalised per-subject radTan = %.4f\n', eqNorm);
fprintf('1/gain-weighted mean of the RAW per-subject radTan       = %.4f\n', wtdRaw);
fprintf('ratio = %.4f, which is exactly mean(1/gain) = %.4f\n', eqNorm/wtdRaw, mean(w));
fprintf(['\nThe two differ only by that global constant, so they rank the four asymmetries\n' ...
         'identically. "Normalising each observer" and "weighting observers by 1/gain" are\n' ...
         'the same operation described from two sides -- so the fact that z-scoring gives\n' ...
         'sub-0037 20%% of the raw-units weight is NOT by itself an objection to it. The\n' ...
         'objection has to be to the divisor, not to the reweighting.\n']);

%% 5. can normalisation reconcile the dissenting observer?
fprintf('\n==== 5. sub-0395 under each normalisation (gain is positive, so signs survive) ====\n');
i395 = find(strcmp(cfg.subjects,'sub-0395'));
fprintf('%-24s %12s %12s %10s\n','normaliser','sub-0395','group mean','ratio');
for d = 1:size(schemes,1)
    v = perSubj.(sprintf('s%d',d));
    fprintf('%-24s %12.3f %12.3f %10.2f\n', schemes{d,1}, v(i395), mean(v), v(i395)/mean(v));
end
fprintf(['\nNo positive scalar can flip a sign: sub-0395 dissents in EVERY set of units, and\n' ...
         'proportionally no less after normalisation. Normalisation does not reconcile that\n' ...
         'observer with the group -- it only reduces the observer''s leverage on the mean.\n']);

function M = wm(C,subjV,pabV,cfg)
    nS=numel(cfg.subjects); nP=numel(cfg.paBins); M=nan(4,nP,nS);
    for si=1:nS
        i0=subjV==cfg.subjects{si};
        for pi=1:nP, idx=i0&(pabV==cfg.paBins(pi)); if any(idx), M(:,pi,si)=median(C(idx,:),1); end, end
    end
end
