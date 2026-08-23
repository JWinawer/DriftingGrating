% RUN_ALL_REPRO  Master driver for the clean-room reproduction (Path A).
% Builds the V1 cache if needed, then regenerates Figures 5, 6, 7, 8 (both the
% z-scored and raw variants) and prints the numeric validation table.
%
% Note: Fig 7 refits the LME 1000x per experiment for bootstrap CIs (~20 s each),
% so a full run takes a few minutes. Run individual run_figN scripts for speed.

cfg = config_repro();
load_and_filter(cfg);         % ensure cache exists

fprintf('\n[1/5] Figures 5 & 6 ...\n'); run_fig5_6;
fprintf('\n[2/5] Figure 7 ...\n');     run_fig7;
fprintf('\n[3/5] Figure 8 ...\n');     run_fig8;
fprintf('\n[4/5] Validation ...\n');   validate_against_manuscript;
fprintf('\n[5/5] Done. Figures in %s\n', cfg.figDir);
