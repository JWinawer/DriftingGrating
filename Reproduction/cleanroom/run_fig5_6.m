% RUN_FIG5_6  Generate Figure 5 (dg) and Figure 6 (da) for both variants.
% Clean-room Path A uses the physically-correct polar-angle convention
% (reference = true wedge angle cfg.paBins).

cfg = config_repro();
T   = load_and_filter(cfg);

variants = {'zscored', true; 'raw', false};
figNums  = struct('dg', 5, 'da', 6);

for expName = {'dg','da'}
    e = expName{1}; expCfg = cfg.(e);
    for vi = 1:size(variants,1)
        vName = variants{vi,1}; doZ = variants{vi,2};
        M = bin_and_aggregate(T, cfg, expCfg, doZ);
        A = compute_asymmetries(M, cfg, expCfg);          % correct convention (default ref)
        plot_fig5_6(A, cfg, expCfg, vName, figNums.(e));
    end
end
fprintf('\nrun_fig5_6: done.\n');
