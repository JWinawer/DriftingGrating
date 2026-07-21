function cfg = config_repro()
% CONFIG_REPRO  Central constants/paths for the clean-room reproduction (Path A).
%
% Everything downstream reads from this struct so there are no magic numbers or
% hard-coded machine paths scattered across scripts. Paths are resolved relative
% to this file, so the code is portable across machines.

    thisDir = fileparts(mfilename('fullpath'));          % .../Reproduction/cleanroom
    reproDir = fileparts(thisDir);                        % .../Reproduction
    repoDir  = fileparts(reproDir);                       % .../DriftingGrating

    cfg.repoDir      = repoDir;
    cfg.reproDir     = reproDir;
    cfg.cleanroomDir = thisDir;
    cfg.csvPath      = fullfile(repoDir, 'Support', 'allsubjectsTable.csv');
    cfg.cacheDir     = fullfile(thisDir, '_cache');
    cfg.figDir       = fullfile(reproDir, 'figures', 'cleanroom');
    cfg.force_reload = false;   % set true to rebuild the V1 cache from the CSV

    % --- analysis inclusion filter (V1 patch) ---
    cfg.roi       = 'V1';
    cfg.eccRange  = [4 8];      % degrees, inclusive
    cfg.r2min     = 0.1;        % pRF variance explained (fraction), strictly greater than

    % --- subjects (fixed order; those who completed both experiments) ---
    cfg.subjects = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123', ...
                    'sub-wlsubj124','sub-0395','sub-0426','sub-0250'};

    % --- polar-angle wedges: order of the PA dimension in all aggregated arrays ---
    % Matches the data order used by the original meanWithinLabel.m.
    cfg.paBins = [0 45 90 135 180 225 270 315];   % degrees (conventional, 0..360)

    % --- bootstrap ---
    cfg.nBoot   = 1000;
    cfg.ciLevel = 95;           % for Figs 5/6 pairwise difference CIs

    % --- metadata columns to pull from the CSV ---
    cfg.metaCols = {'subject','visual_area','pRF_angle_bin','pRF_angle', ...
                    'pRF_ecc','pRF_r2','pRF_sigma','included'};

    % --- per-experiment column maps -----------------------------------------
    % 'oriCols' are the 4 stationary-orientation beta columns, in the fixed
    % order [horizontal-like, vertical-like, oblique1, oblique2], with 'oriAngle'
    % giving each stimulus's local orientation in degrees (used for the
    % PA-dependent radial/tangential derivation). 'blank' and 'betaStd' name the
    % blank-condition and the per-vertex std columns.
    %
    % 'oriAngle' is each stimulus's local orientation at the upper vertical
    % meridian (UVM, theta=90), in degrees mod 180. 'isPolar' says whether that
    % local orientation rotates with polar angle: local orientation at theta is
    %   mod(oriAngle + isPolar*(theta - 90), 180).
    % For Cartesian gratings it is constant; for polar gratings it rotates.
    % These match the repo's index->UVM mapping (dg names s0/s90/s45/s135 -> 0/90/45/135;
    % da annulus/pinwheel/ccspiral/cspiral share those indices 26/27/28/29).

    % Cartesian experiment (dg): orientations are literal grating orientations.
    cfg.dg.name     = 'dg';
    cfg.dg.isPolar  = false;
    cfg.dg.oriNames = {'horizontal','vertical','rightleaning','leftleaning'};
    cfg.dg.oriCols  = {'cartexp_horizontal_stationary', ...
                       'cartexp_vertical_stationary', ...
                       'cartexp_rightleaning_grating_stationary', ...
                       'cartexp_leftleaning_grating_stationary'};
    cfg.dg.oriAngle = [0 90 45 135];        % horizontal=0, vertical=90, right-leaning=45, left-leaning=135
    cfg.dg.blank    = 'cartexp_blank';
    cfg.dg.betaStd  = 'dg_beta_std';
    cfg.dg.betaMean = 'dg_beta_mean';

    % Polar experiment (da): at UVM, pinwheel(radial)=vertical=90, annulus(tangential)
    % =horizontal=0, ccspiral=45, cspiral=135; these rotate with polar angle.
    cfg.da.name     = 'da';
    cfg.da.isPolar  = true;
    cfg.da.oriNames = {'pinwheel','annulus','ccspiral','cspiral'};
    cfg.da.oriCols  = {'polexp_pinwheel_grating_stationary', ...
                       'polexp_annulus_grating_stationary', ...
                       'polexp_ccspiral_grating_stationary', ...
                       'polexp_cspiral_grating_stationary'};
    cfg.da.oriAngle = [90 0 45 135];        % UVM local orientation of each stimulus
    cfg.da.blank    = 'polexp_blank';
    cfg.da.betaStd  = 'da_beta_std';
    cfg.da.betaMean = 'da_beta_mean';
end
