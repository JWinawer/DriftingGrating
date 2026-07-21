% RUN_PATHB_FIGURES  Regenerate Figures 5 & 6 through the ACTUAL repo functions.
%
% Replicates the calls in AnalysisCode/04_plot_betaAsymmetries/plot_NeuralAsymmetries.m
% (plot1_experimentalCond / plot2_experimentalCond / compute_derivativeDirections),
% but fed by the CSV bridge instead of the lab pipeline, and writing to
% Reproduction/figures/bridge. Canonical AnalysisCode files are used unedited.
%
% This is the z-scored variant (the manuscript default). Because it uses the real
% derived-direction routine, the da H-V panel here carries the PA-ordering artifact
% (see FINDINGS.md); the clean-room figures show the corrected version.

clear; clc;
here     = fileparts(mfilename('fullpath'));
reproDir = fileparts(here);
repoDir  = fileparts(reproDir);
addpath(fullfile(reproDir,'cleanroom'));
addpath(genpath(fullfile(repoDir,'AnalysisCode')));

cfg = config_repro();
T   = load_and_filter(cfg);
ps  = loadConfig(fileparts(repoDir));
ps.comparisonName = 'orientation_minus_baseline';
ps.subjects = cfg.subjects;
ps.figureDir = fullfile(reproDir,'figures','bridge');
if ~isfolder(ps.figureDir), mkdir(ps.figureDir); end

nDerived   = {{1:2,3},{1,2}};   % {derived cardinal/oblique}, {derived subset}
mainSubset = [0 1];

for e = {'dg','da'}
    en = e{1}; ps.projectName = en;
    G = build_group_matrices_fromCSV(cfg, cfg.(en), true, T);
    mbpa = G.medianBOLDpa;

    % --- MAIN (direct) conditions ---
    for s = mainSubset
        plot1_experimentalCond(mbpa, 'mainCardinalVsMainOblique', ps, s);
        plot2_experimentalCond(mbpa, 'mainCardinalVsMainOblique', ps, s);
    end

    % --- DERIVED conditions (via the real compute_derivativeDirections) ---
    for ci = 1:numel(nDerived)
        s = mainSubset(ci);
        nm = compute_derivativeDirections(mbpa, ps, s);
        avg = mean(nm(nDerived{ci}{1}, :, :, :), 1);
        proconMatrix = cat(1, avg, nm(nDerived{ci}{2}, :, :, :));
        plot1_experimentalCond(proconMatrix, 'derivedCardinalVsDerivedOblique', ps, s);
        plot2_experimentalCond(proconMatrix, 'derivedCardinalVsDerivedOblique', ps, s);
    end
    fprintf('run_pathB_figures: %s done\n', en);
end
fprintf('\nFigures written to %s\n', ps.figureDir);
