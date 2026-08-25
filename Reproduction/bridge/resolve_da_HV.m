% RESOLVE_DA_HV  Cross-check: repo derived-direction code vs the clean-room.
%
% Builds the bridged da medianBOLDpa from the CSV in the repo's PA order, runs the
% *actual* repo function compute_derivativeDirections.m, and compares the resulting
% da horizontal-vertical asymmetry -- per polar angle and as a group mean -- against
% the clean-room computation and the manuscript.
%
% History: this script previously reported a cardinal-meridian "bug" in
% compute_derivativeDirections.m. That was wrong. The repo function is correct; the
% discrepancy came from two errors in this reproduction (swapped spirals in
% config_repro.m, and this bridge emitting the PA dimension in ascending conventional
% order instead of the repo's Benson order). Both are fixed. See ../STIMULUS_CONVENTIONS.md.

clear; clc;
here      = fileparts(mfilename('fullpath'));
reproDir  = fileparts(here);
repoDir   = fileparts(reproDir);
cleanroom = fullfile(reproDir, 'cleanroom');
analysis  = fullfile(repoDir, 'AnalysisCode');
addpath(cleanroom);
addpath(genpath(analysis));          % repo analysis functions (unedited)

cfg = config_repro();
T   = load_and_filter(cfg);

% projectSettings the repo functions expect
githubDir = fileparts(repoDir);                      % .../Projects  (contains DriftingGrating)
projectSettings = loadConfig(githubDir);
projectSettings.projectName    = 'da';
projectSettings.comparisonName = 'orientation_minus_baseline';

% --- bridged da medianBOLDpa (z-scored), PA dimension in cfg.paBinsRepoOrder ---
G = build_group_matrices_fromCSV(cfg, cfg.da, true, T);
medianBOLDpa = G.medianBOLDpa;                       % 29 x 8 x 1 x 8

% --- run the ACTUAL repo derived-direction code ---
% newMatrix rows for da: 1 = horizontal(local), 2 = vertical(local), 3 = oblique
newMatrix = compute_derivativeDirections(medianBOLDpa, projectSettings);
HV_repo_subj = squeeze(newMatrix(1,:,1,:) - newMatrix(2,:,1,:));   % 8 PA(repo order) x 8 subj
HV_repo_repoOrder = mean(HV_repo_subj, 2).';
HV_repo_group     = mean(mean(HV_repo_subj, 1));

% re-express the repo result on the clean-room's ascending conventional axis
[~, toConv] = ismember(cfg.paBins, cfg.paBinsRepoOrder);
HV_repo_perPA = HV_repo_repoOrder(toConv);

% --- clean-room (independent implementation, conventional order) ---
Mda = bin_and_aggregate(T, cfg, cfg.da, true);
Acr = compute_asymmetries(Mda, cfg, cfg.da);           % default ref = true wedge angle
HV_clean_perPA = mean(Acr.HV.diff, 2).';
HV_clean_group = mean(mean(Acr.HV.diff, 1));

fprintf('\n================= da horizontal - vertical =================\n');
fprintf('polar angle (conventional) : %s\n', mat2str(cfg.paBins));
fprintf('repo  per-PA               : %s\n', mat2str(round(HV_repo_perPA, 3)));
fprintf('clean per-PA               : %s\n', mat2str(round(HV_clean_perPA, 3)));
fprintf('max abs per-PA difference  : %.3e\n', max(abs(HV_repo_perPA - HV_clean_perPA)));
fprintf('-----------------------------------------------------------\n');
fprintf('repo  code (compute_derivativeDirections) group mean : % .3f\n', HV_repo_group);
fprintf('clean-room (independent implementation)   group mean : % .3f\n', HV_clean_group);
fprintf('manuscript Fig 6A reported                           :  -0.450\n');

if max(abs(HV_repo_perPA - HV_clean_perPA)) < 1e-9
    fprintf('\nPASS: repo code and clean-room agree exactly, per polar angle.\n');
else
    fprintf('\nFAIL: paths disagree -- check cfg.paBinsRepoOrder and cfg.da.oriAngle/oriIdx.\n');
end
