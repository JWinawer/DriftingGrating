% VALIDATE_AGAINST_MANUSCRIPT  Path A numeric check (z-scored) vs the manuscript.
% Prints, per experiment, the four asymmetries computed two ways -- the independent
% group-mean difference (Figs 5/6) and the joint-LME delta=2*beta (Fig 7) -- next to
% the manuscript's reported values. LME point estimates only (bootstrap off) for speed.

cfg = config_repro();
T   = load_and_filter(cfg);

% manuscript reported values (sigma units); LME entries are |reported| Fig 7
man.dg.ind = [-1.155 -0.40  0.23  0.06];
man.dg.lme = [ 1.10   0.36  0.19  0.00];
man.da.ind = [-0.45  -0.06  0.60  0.17];
man.da.lme = [ 0.02   0.06  0.60  0.18];
names = {'horiz-vert','card-obl','rad-tang','polc-polo'};

fprintf('\n================ Path A validation (z-scored) ================\n');
for e = {'dg','da'}
    en = e{1}; expCfg = cfg.(en);
    M = bin_and_aggregate(T, cfg, expCfg, true);
    A = compute_asymmetries(M, cfg, expCfg);
    ind = cellfun(@(a) mean(mean(A.(a).diff,1)), A.order);      % independent means
    res = fit_lme_fig7(M, cfg, expCfg, [], false);             % LME (no bootstrap)

    fprintf('\n--- %s experiment ---\n', en);
    fprintf('%-11s | independent  (manu)  | LME delta   (|manu|)\n', 'asymmetry');
    for j = 1:4
        flag = '';
        if abs(ind(j) - man.(en).ind(j)) > 0.1, flag = '   <-- DIVERGES'; end
        fprintf('%-11s |   % .3f   (% .3f) |  % .3f    ( %.2f)%s\n', ...
                names{j}, ind(j), man.(en).ind(j), res.delta(j), man.(en).lme(j), flag);
    end
end
fprintf('\nNote: da horiz-vert independent diverges (clean ~0 vs manuscript -0.45); the\n');
fprintf('manuscript''s own Fig 7 LME reports ~0.02, agreeing with the clean-room. See FINDINGS.md.\n');
