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
    % TWO figure destinations, with different jobs.
    %   figDir     scratch. Reproduction/figures/ is git-ignored, and the older
    %              reproduction and diagnostic plots (RUN_ALL_REPRO, PLOT_HARMONIC)
    %              write here because they are working output, not deliverables.
    %   suppFigDir deliverables. Tracked despite the blanket *.png / *.pdf ignore
    %              rules, because these are the figures RESULTS.md and the
    %              manuscript refer to. The specification figures write here
    %              DIRECTLY -- they are not copied from figDir, because keeping the
    %              same image in two places is how the two drift apart.
    cfg.figDir       = fullfile(reproDir, 'figures', 'cleanroom');
    cfg.suppFigDir   = fullfile(reproDir, 'supplement', 'figures');
    cfg.force_reload = false;   % set true to rebuild the V1 cache from the CSV

    % --- analysis inclusion filter (V1 patch) ---
    cfg.roi       = 'V1';
    cfg.eccRange  = [4 8];      % degrees, inclusive
    cfg.r2min     = 0.1;        % pRF variance explained (fraction), strictly greater than

    % --- subjects (fixed order; those who completed both experiments) ---
    cfg.subjects = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123', ...
                    'sub-wlsubj124','sub-0395','sub-0426','sub-0250'};

    % --- polar-angle wedges --------------------------------------------------
    % Order of the PA dimension in the CLEAN-ROOM arrays. The CSV's pRF_angle has
    % already been through map_theta (createTables.m:75), so these are CONVENTIONAL
    % angles: 0 = right horizontal meridian, counter-clockwise.
    cfg.paBins = [0 45 90 135 180 225 270 315];   % degrees (conventional, 0..360)

    % The ORIGINAL pipeline uses a different order. meanWithinLabel.m bins the raw
    % Benson angle_adj (it does NOT call map_theta) in ascending Benson order, so the
    % PA dimension of its medianBOLDpa carries the conventional angles below.
    % compute_derivativeDirections.m:48, plot1_experimentalCond.m:121 and lme1_fit.m:88
    % all assume exactly this. ANY array handed to a repo stage-04 function must be in
    % this order -- getting it wrong reflects the wedges about 45 deg, which swaps the
    % four cardinal meridians. See ../STIMULUS_CONVENTIONS.md sections 3 and 5.
    cfg.paBinsRepoOrder = [90 45 0 315 270 225 180 135];

    % --- across-vertex aggregator, and observer gain rescaling ------------------
    % MEAN, not median (JW, 2026-08-19). The repo's meanWithinLabel.m saves both
    % meanBOLDpa and medianBOLDpa; lme1_fit.m (Fig 7) already reads meanBOLDpa, while
    % plot_NeuralAsymmetries.m (Figs 5/6) passed medianBOLDpa. Mean + gain rescaling
    % reproduces all eight manuscript asymmetries to +-0.003; the median route does not
    % (it misses dg horiz-vert by 0.07). See ../local_qc/RELIABILITY.md.
    cfg.aggregator = 'mean';        % 'mean' | 'median'

    % Per-observer pRF gain rescaling, as applied in lme1_fit.m/plot1_/plot2_.
    % Set cfg.gainFile = '' to disable. See OBSERVER_GAIN_WEIGHTS.
    cfg.collectDir = dg_collect_dir();   % see DG_COLLECT_DIR; '' if not mounted
    cfg.gainFile = fullfile(cfg.collectDir, 'gainSummary.csv');
    cfg.gainMean = 'geometric';     % 'geometric' (manuscript) | 'arithmetic' (repo code)

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
    cfg.dg.oriIdx   = [26 27 28 29];        % CONTRASTS.json index per oriCols entry (s0/s90/s45/s135 _v_b)
    cfg.dg.blank    = 'cartexp_blank';
    cfg.dg.betaStd  = 'dg_beta_std';
    cfg.dg.betaMean = 'dg_beta_mean';

    % Polar experiment (da): at UVM, pinwheel(radial)=vertical=90, annulus(tangential)
    % =horizontal=0, cspiral=45, ccspiral=135; these rotate with polar angle.
    % The spiral identities are fixed by CONTRASTS.json (idx 28 = scspiral_v_b = s45,
    % idx 29 = sccspiral_v_b = s135) and by createTables.m:146-149, and were confirmed
    % against sub-0255's raw betas to 5e-15 (STIMULUS_CONVENTIONS.md section 6).
    % They were swapped here
    % until 2026-07-22, which sign-flipped the four OBLIQUE wedges of every derived
    % asymmetry and produced a spurious da horiz-vert of -0.041 instead of -0.446.
    cfg.da.name     = 'da';
    cfg.da.isPolar  = true;
    cfg.da.oriNames = {'pinwheel','annulus','ccspiral','cspiral'};
    cfg.da.oriCols  = {'polexp_pinwheel_grating_stationary', ...
                       'polexp_annulus_grating_stationary', ...
                       'polexp_ccspiral_grating_stationary', ...
                       'polexp_cspiral_grating_stationary'};
    cfg.da.oriAngle = [90 0 135 45];        % UVM local orientation: pinwheel/annulus/ccspiral/cspiral
    cfg.da.oriIdx   = [27 26 29 28];        % CONTRASTS.json index: pinwheel=27, annulus=26, ccspiral=29, cspiral=28
    cfg.da.blank    = 'polexp_blank';
    cfg.da.betaStd  = 'da_beta_std';
    cfg.da.betaMean = 'da_beta_mean';
end
