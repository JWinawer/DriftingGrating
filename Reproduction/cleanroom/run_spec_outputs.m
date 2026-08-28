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
%
% TWO OBSERVER SETS, and which figure gets which is deliberate:
%
%   Figure 5 (dg)        ALL 13 observers. Within-experiment, so it needs no da.
%   Figure 6 (da)        the 7 valid observers (sub-0395 saw a pilot stimulus).
%   profile figure       the matched 7 -- it puts dg and da side by side.
%   every table          the matched 7 -- the context effect is a paired difference.
%
% dg is therefore computed TWICE, and the second computation is not redundant. The
% gain rescaling multiplies by a group gain formed over whichever observers are in the
% set, so the same observer's dg values are a uniform 6.7% larger inside the 13-set
% than inside the 7-set. Differencing a 13-based dg against a 7-based da would push
% that factor into the contrast on one side only. SPEC_TABLES refuses a mismatched
% pair outright; this is how the right pair is supplied.
fprintf('\n================ V1 4-8 deg, THREE VARIANTS ================\n');
S     = cell(1, numel(V));      % matched 7/7: Figure 6, the profile figure, all tables
Sdg13 = cell(1, numel(V));      % dg on 13:    Figure 5 only
for k = 1:numel(V)
    fprintf('\n--- %s: %s ---\n', V(k).tag, V(k).label);
    S{k}     = spec_profiles('area','V1', 'route', V(k).route);
    Sdg13{k} = spec_profiles('area','V1', 'route', V(k).route, 'dgSubjectMode','all');
end

% Limits span BOTH subject sets, or Figure 5 (13 observers) and Figure 6 (7) would sit
% on different scales -- and dg-versus-da is the comparison the reader most needs to
% make by eye.
L = spec_axis_limits([S, Sdg13]);
fprintf(['\nshared scales across all three variants and both observer sets: polar +/-%.2f, ' ...
         'difference panels [%.2f %.2f], profile [%.2f %.2f]\n'], L.rmax, L.dot, L.prof);

for k = 1:numel(V)
    plot_fig5_6_spec(Sdg13{k}, 'dg', 5, cfg.suppFigDir, L, V(k));   % 13 observers
    plot_fig5_6_spec(S{k},     'da', 6, cfg.suppFigDir, L, V(k));   %  7 observers
    plot_spec_profile(S{k},           cfg.suppFigDir, L, V(k));     %  7, both experiments
end

fprintf('\n================ V1 TABLES ================\n');
% Hand the MATCHED profiles in explicitly rather than letting SPEC_TABLES recompute
% them from the config default. Two reasons: it saves recomputing what was just built,
% and it makes the pairing a property of this call instead of a property of whatever
% cfg.dgSubjectMode happens to be. SPEC_TABLES refuses a mismatched pair either way.
for k = 1:numel(V)
    spec_tables('area','V1', 'variant', V(k).tag, 'profiles', S{k});
end

% The 13-observer dg asymmetries -- the numbers BEHIND FIGURE 5, which the paired
% tables above do not contain because they are the matched 7. 'context' false because
% dg-13 and da-7 are different people: the per-experiment asymmetries are well defined,
% a paired contrast is not, and SPEC_TABLES still refuses one if asked.
fprintf('\n================ V1 dg, ALL 13 OBSERVERS (Figure 5) ================\n');
for k = 1:numel(V)
    spec_tables('area','V1', 'variant', V(k).tag, 'profiles', Sdg13{k}, ...
                'context', false, 'fileSuffix', 'dg13');
end

fprintf('\n================ VISUAL-HIERARCHY SWEEP ================\n');
for k = 1:numel(V)
    fprintf('\n########## variant %s ##########\n', V(k).tag);
    A = spec_areas_summary('variant', V(k).tag);
    plot_spec_hierarchy(A, cfg.suppFigDir);
end

fprintf('\nrun_spec_outputs: done.\n');
