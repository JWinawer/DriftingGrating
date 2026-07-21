% RUN_PATHB_VALUES  Reproduce the Fig 5-8 asymmetry VALUES through the existing code.
%
% For both experiments and both variants, builds the bridged group array and computes
% all four asymmetries using the *actual* repo functions:
%   - direct asymmetries   : retrieveProConIdx.m + medianBOLDpa indexing
%   - derived asymmetries  : compute_derivativeDirections.m
% aggregated equally over the 8 wedges (as plot2_experimentalCond does). Prints a
% table comparing existing-code vs clean-room (Path A) vs manuscript, exposing exactly
% which asymmetries the existing pipeline reproduces and which carry the PA-ordering
% artifact.

clear; clc;
here      = fileparts(mfilename('fullpath'));
reproDir  = fileparts(here);
repoDir   = fileparts(reproDir);
addpath(fullfile(reproDir, 'cleanroom'));
addpath(genpath(fullfile(repoDir, 'AnalysisCode')));   % repo functions, unedited

cfg = config_repro();
T   = load_and_filter(cfg);
githubDir = fileparts(repoDir);
ps = loadConfig(githubDir);
ps.comparisonName = 'orientation_minus_baseline';

% manuscript independent (Fig 5/6) values, in asymmetry order [HV cardObl radTan polcardPolobl]
man.dg = [-1.155 -0.40 0.23 0.06];
man.da = [-0.45  -0.06 0.60 0.17];
asymNames = {'HV','cardObl','radTan','polcardPolobl'};

for variant = {'zscored', 'raw'}
  vName = variant{1}; doZ = strcmp(vName,'zscored');
  for e = {'dg','da'}
    en = e{1}; expCfg = cfg.(en);
    ps.projectName = en;

    % --- bridged group array + clean-room array (from same CSV/filter) ---
    G = build_group_matrices_fromCSV(cfg, expCfg, doZ, T);
    mbpa = G.medianBOLDpa;
    Mcr  = bin_and_aggregate(T, cfg, expCfg, doZ);
    Acr  = compute_asymmetries(Mcr, cfg, expCfg);
    clean = cellfun(@(a) mean(mean(Acr.(a).diff,1)), asymNames);

    % --- existing-code values ---
    % direct (main) asymmetries via retrieveProConIdx
    [p1,c1] = retrieveProConIdx(en,'orientation_minus_baseline',1);  % subset 1
    [p0,c0] = retrieveProConIdx(en,'orientation_minus_baseline',0);  % subset 0
    mainSub = directAsym(mbpa, p1, c1);   % dg: HV ; da: radTan
    mainCar = directAsym(mbpa, p0, c0);   % dg: cardObl ; da: polcardPolobl

    % derived asymmetries via the real compute_derivativeDirections
    nm = compute_derivativeDirections(mbpa, ps);          % 3 x 8 x 1 x 8
    derSub = groupmean(squeeze(nm(1,:,1,:) - nm(2,:,1,:)));               % row1 - row2
    derCar = groupmean(squeeze(mean(nm(1:2,:,1,:),1) - nm(3,:,1,:)));      % mean(1,2) - row3

    % assemble into asymmetry order [HV cardObl radTan polcardPolobl]
    if strcmp(en,'dg')
      exist = [mainSub, mainCar, derSub, derCar];   % HV,cardObl direct; radTan,polcard derived
    else
      exist = [derSub, derCar, mainSub, mainCar];   % HV,cardObl derived; radTan,polcard direct
    end

    fprintf('\n=== %s / %s :  existing-code | clean-room | manuscript ===\n', en, vName);
    for j = 1:4
      tag = '';
      if abs(exist(j) - clean(j)) > 0.05, tag = '   <-- existing code diverges from correct'; end
      fprintf('  %-14s % .3f    % .3f    % .3f%s\n', asymNames{j}, exist(j), clean(j), man.(en)(j), tag);
    end
  end
end

% ---- helpers ----
function g = directAsym(mbpa, proIdx, conIdx)
    pro = squeeze(mean(mbpa(proIdx,:,1,:),1));   % nP x nS
    con = squeeze(mean(mbpa(conIdx,:,1,:),1));
    g = groupmean(pro - con);
end
function g = groupmean(d)   % d: nP x nS -> equal-weight over PA, then over subjects
    g = mean(mean(d,1,'omitnan'),2,'omitnan');
end
