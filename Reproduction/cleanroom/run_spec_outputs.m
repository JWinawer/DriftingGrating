% RUN_SPEC_OUTPUTS  Figures 5/6 and every table, in all three analysis variants.
%
% Route and weighting are orthogonal, so the three variants change one thing at a time
% (SPEC_VARIANTS):
%
%   spec    harmonic model, continuous thetaV, equal weighting     <- PRIMARY
%   roi     eight polar-angle wedges,          equal weighting     spec vs roi   = ROUTE
%   roipw   eight polar-angle wedges,          precision weighting roi vs roipw  = WEIGHTING
%
% Every figure in every variant is drawn on ONE set of axis limits, computed across all
% three, so the variants can be laid side by side and compared by eye. That is the whole
% reason for producing them.
%
%   Reproduction/supplement/figures/   figures, PNG + PDF, TRACKED. Written DIRECTLY,
%     Figure_5_<tag>_dg, Figure_6_<tag>_da, Figure_5_6_<tag>_profile,
%     Figure_S5_<tag>_hierarchy          for tag in spec | roi | roipw
%   Reproduction/supplement/           tables, CSV, TRACKED, one set per variant
%
% Nothing is copied into Reproduction/figures/, which stays the scratch destination for
% the older reproduction and diagnostic plots. One image, one location.
%
% COST, measured rather than guessed. About 2-3 minutes. The dominant term is the
% 500-resample bootstrap over runs inside DIAGNOSE_WITHIN_OBSERVER_ERROR, which measures
% each observer's within-observer sigma: it is ~75% of every call (4.7 s per map x band
% in V1, against 1.2 s at nBoot = 50). File I/O is negligible -- 0.69 s for a pass over
% all 16 runbetas files. If this ever needs to be faster, nBoot is the only real lever,
% and lowering it degrades sigma, which is what the precision weighting rests on.
%
% SPEC_TABLES is handed the profiles SPEC_AREAS_SUMMARY already computed, rather than
% recomputing them. That halved the sweep and changes no number (verified: delta = 0).

cfg = config_repro();
if ~isfolder(cfg.suppFigDir), mkdir(cfg.suppFigDir); end
V = spec_variants();

% --- V1 profiles for every variant first, so the axis limits span all of them -------
fprintf('\n================ V1 4-8 deg, THREE VARIANTS ================\n');
S = cell(1, numel(V));
for k = 1:numel(V)
    fprintf('\n--- %s: %s ---\n', V(k).tag, V(k).label);
    S{k} = spec_profiles('area','V1', 'route', V(k).route);
end
L = spec_axis_limits(S);
fprintf(['\nshared scales across all three variants: polar +/-%.2f, ' ...
         'difference panels [%.2f %.2f], profile [%.2f %.2f]\n'], L.rmax, L.dot, L.prof);

for k = 1:numel(V)
    plot_fig5_6_spec(S{k}, 'dg', 5, cfg.suppFigDir, L, V(k));
    plot_fig5_6_spec(S{k}, 'da', 6, cfg.suppFigDir, L, V(k));
    plot_spec_profile(S{k},        cfg.suppFigDir, L, V(k));
end

fprintf('\n================ V1 TABLES ================\n');
for k = 1:numel(V), spec_tables('area','V1', 'variant', V(k).tag); end

fprintf('\n================ VISUAL-HIERARCHY SWEEP ================\n');
for k = 1:numel(V)
    fprintf('\n########## variant %s ##########\n', V(k).tag);
    A = spec_areas_summary('variant', V(k).tag);
    plot_spec_hierarchy(A, cfg.suppFigDir);
end

fprintf('\nrun_spec_outputs: done.\n');
