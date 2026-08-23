% RUN_FIG7  Fit the joint LME and draw Figure 7 for both variants.
cfg = config_repro();
T   = load_and_filter(cfg);

variants = {'zscored', true; 'raw', false};
for vi = 1:size(variants,1)
    vName = variants{vi,1}; doZ = variants{vi,2};
    Mdg = bin_and_aggregate(T, cfg, cfg.dg, doZ);
    Mda = bin_and_aggregate(T, cfg, cfg.da, doZ);
    resDG = fit_lme_fig7(Mdg, cfg, cfg.dg);
    resDA = fit_lme_fig7(Mda, cfg, cfg.da);
    plot_fig7(resDG, resDA, cfg, vName);
end
fprintf('\nrun_fig7: done.\n');
