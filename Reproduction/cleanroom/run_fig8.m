% RUN_FIG8  Per-location data + model polar plots (Figure 8), both variants.
% Uses the fixed-effect LME prediction (no bootstrap needed here).
cfg = config_repro();
T   = load_and_filter(cfg);

variants = {'zscored', true; 'raw', false};
labels   = struct('dg','A', 'da','B');
for vi = 1:size(variants,1)
    vName = variants{vi,1}; doZ = variants{vi,2};
    for e = {'dg','da'}
        en = e{1}; expCfg = cfg.(en);
        M   = bin_and_aggregate(T, cfg, expCfg, doZ);
        res = fit_lme_fig7(M, cfg, expCfg, [], false);   % skip bootstrap
        plot_fig8(M, res, cfg, expCfg, vName, labels.(en));
    end
end
fprintf('\nrun_fig8: done.\n');
