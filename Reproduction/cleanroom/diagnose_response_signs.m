% DIAGNOSE_RESPONSE_SIGNS  Are the negative per-subject responses real, or a blank artefact?
%
% Two subjects show a negative median (mean stationary orientation - blank) in V1 4-8 deg
% for the polar experiment, which would mean V1 responded LESS to a grating in its
% receptive field than to blank. This script separates the possible causes by pulling all
% 13 GLMsingle conditions per experiment (LOAD_ALLCONDITIONS):
%
%   - is the blank beta an outlier relative to the 12 stimulus betas?
%   - do the 8 (unanalysed) motion conditions show the same sign as the 4 stationary ones?
%   - what is the beta zero-point, i.e. what does a beta of 0 mean here?
%   - does an effect-independent divisor (motion conditions only) change Fig 7?

cfg = config_repro();
T   = load_allconditions(cfg);
keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
T = T(keep,:);
subjV = string(T.subject); pabV = T.pRF_angle_bin;
nS = numel(cfg.subjects); nP = numel(cfg.paBins);
[dgMot, daMot] = allcondition_cols();

%% 0. what are beta_mean / beta_std computed over?
fprintf('\n==== 0. definition check: beta_mean / beta_std vs the 13 conditions ====\n');
for E = {cfg.dg, cfg.da}
    E = E{1}; %#ok<FXSET>
    if strcmp(E.name,'dg'), mot = dgMot; else, mot = daMot; end
    cols = [E.oriCols, mot, {E.blank}];
    B = zeros(height(T), 13);
    for k = 1:13, B(:,k) = T.(cols{k}); end
    fprintf('%s: max|beta_mean - mean(13)| = %.3e ; max|beta_std - std(13)| = %.3e  (n=13 cols)\n', ...
        E.name, max(abs(T.(E.betaMean) - mean(B,2))), max(abs(T.(E.betaStd) - std(B,0,2))));
end

%% 1. per-subject condition breakdown
for E = {cfg.dg, cfg.da}
    E = E{1}; %#ok<FXSET>
    if strcmp(E.name,'dg'), mot = dgMot; else, mot = daMot; end
    S4 = zeros(height(T),4); M8 = zeros(height(T),8);
    for k = 1:4, S4(:,k) = T.(E.oriCols{k}); end
    for k = 1:8, M8(:,k) = T.(mot{k}); end
    blank = T.(E.blank);
    all12 = [S4 M8];

    fprintf('\n==== 1. %s : per-subject medians of the raw betas ====\n', E.name);
    fprintf('%-14s %8s %8s %8s | %9s %9s %9s\n', 'subject', ...
        'stat4', 'mot8', 'BLANK', 'stat4-blk', 'mot8-blk', 'all12-blk');
    for si = 1:nS
        i = subjV == cfg.subjects{si};
        fprintf('%-14s %8.3f %8.3f %8.3f | %9.3f %9.3f %9.3f\n', cfg.subjects{si}, ...
            median(mean(S4(i,:),2)), median(mean(M8(i,:),2)), median(blank(i)), ...
            median(mean(S4(i,:),2)-blank(i)), median(mean(M8(i,:),2)-blank(i)), ...
            median(mean(all12(i,:),2)-blank(i)));
    end

    % where does the blank sit among the 13 conditions?
    fprintf('\n  blank rank among the 13 conditions (1 = smallest beta), median over vertices:\n  ');
    for si = 1:nS
        i = find(subjV == cfg.subjects{si});
        r = sum([S4(i,:) M8(i,:)] < blank(i), 2) + 1;
        fprintf('%s=%.1f  ', erase(cfg.subjects{si},'sub-'), median(r));
    end
    fprintf('\n  (13 = blank is the LARGEST beta; 1 = the smallest. Positive stimulus drive => low rank)\n');

    % beta zero point: what fraction of all 13 betas are negative?
    fprintf('  fraction of all 13 betas < 0 (pooled vertices) = %.2f ; median of all 13 = %.3f\n', ...
        mean([S4(:); M8(:); blank(:)] < 0), median([S4(:); M8(:); blank(:)]));
end

%% 2. does an effect-independent divisor change Fig 7B?
fprintf('\n==== 2. da Fig 7: alternative per-vertex divisors ====\n');
E = cfg.da;
X = zeros(height(T),4);
for k = 1:4, X(:,k) = T.(E.oriCols{k}) - T.(E.blank); end
M8 = zeros(height(T),8);
for k = 1:8, M8(:,k) = T.(daMot{k}); end

divs = { ...
  'none (raw)',                    ones(height(T),1); ...
  'beta_std, all 13 (manuscript)',  T.(E.betaStd); ...
  'std of the 8 motion conds',     std(M8,0,2); ...
  'mean(8 motion) - blank',        mean(M8,2) - T.(E.blank)};

asymNames = {'HV','cardObl','radTan','polcardPolobl'};
fprintf('%-30s %9s %9s %9s %9s   %s\n', 'divisor', asymNames{:}, 'ordering');
for d = 1:size(divs,1)
    A = compute_asymmetries(wedgeMed(X./divs{d,2}, subjV, pabV, cfg), cfg, E);
    v = cellfun(@(n) mean(A.(n).diff(:),'omitnan'), asymNames);
    ord = 'H-V larger'; if v(3) > abs(v(1)), ord = 'radTan larger'; end
    fprintf('%-30s %9.3f %9.3f %9.3f %9.3f   %s\n', divs{d,1}, v, ord);
end
fprintf(['\nNote: "mean(8 motion) - blank" still contains the blank and goes near zero for\n' ...
         'some vertices, so it is a fragile divisor; the motion-only std does not.\n']);

function M = wedgeMed(C, subjV, pabV, cfg)
    nS=numel(cfg.subjects); nP=numel(cfg.paBins); M=nan(4,nP,nS);
    for si=1:nS
        i0 = subjV==cfg.subjects{si};
        for pi=1:nP
            idx = i0 & (pabV==cfg.paBins(pi));
            if any(idx), M(:,pi,si)=median(C(idx,:),1); end
        end
    end
end
