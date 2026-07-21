% RESOLVE_DA_HV  Definitive test of the da horizontal-vertical discrepancy.
%
% Builds the bridged da medianBOLDpa from the CSV, then runs the *actual* repo
% function compute_derivativeDirections.m (which the manuscript's Fig 6A used) and
% computes the da H-V asymmetry it produces. Compares, per polar angle, against the
% clean-room (physically-correct convention). This settles whether the manuscript's
% -0.45 is an artifact of the derived-direction polar-angle ordering.

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

% --- bridged da medianBOLDpa (z-scored) ---
G = build_group_matrices_fromCSV(cfg, cfg.da, true, T);
medianBOLDpa = G.medianBOLDpa;                       % 29 x 8 x 1 x 8

% --- run the ACTUAL repo derived-direction code ---
% newMatrix rows for da: 1 = horizontal(local), 2 = vertical(local), 3 = oblique
newMatrix = compute_derivativeDirections(medianBOLDpa, projectSettings);
HV_repo_perPA_subj = squeeze(newMatrix(1,:,1,:) - newMatrix(2,:,1,:));   % 8 PA x 8 subj
HV_repo_perPA = mean(HV_repo_perPA_subj, 2).';                          % avg over subjects
HV_repo_group = mean(mean(HV_repo_perPA_subj, 1));

% --- clean-room (physically-correct convention) ---
Mda = bin_and_aggregate(T, cfg, cfg.da, true);
Acr = compute_asymmetries(Mda, cfg, cfg.da);           % default ref = true wedge angle
HV_clean_perPA = mean(Acr.HV.diff, 2).';
HV_clean_group = mean(mean(Acr.HV.diff, 1));

fprintf('\n================= da horizontal - vertical =================\n');
fprintf('polar-angle order : %s\n', mat2str(cfg.paBins));
fprintf('repo  per-PA      : %s\n', mat2str(round(HV_repo_perPA, 3)));
fprintf('clean per-PA      : %s\n', mat2str(round(HV_clean_perPA, 3)));
fprintf('-----------------------------------------------------------\n');
fprintf('repo  code (compute_derivativeDirections) group mean : % .3f\n', HV_repo_group);
fprintf('clean-room (correct convention)           group mean : % .3f\n', HV_clean_group);
fprintf('manuscript Fig 6A reported                           :  -0.450\n');
fprintf('manuscript Fig 7B (LME) reported                     :   0.020\n');
