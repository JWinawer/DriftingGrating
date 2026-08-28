function T = build_table_rows(subject, opts)
% BUILD_TABLE_ROWS  One subject's allsubjectsTable rows, from ~/dg_collect.
%
%   T = build_table_rows('sub-0442')
%   T = build_table_rows('sub-0442', struct('projects', {{'dg'}}))
%
% Reproduces AnalysisCode/01_process_singlesubjectGLM/CREATETABLES.M for a single
% observer, reading the ALREADY-COLLECTED local files rather than the mounted volume:
%
%   ret_<subj>.mat        retinotopy   (angle_adj, eccen, vexpl, sigma)
%   labels_<subj>.mat     FreeSurfer labels
%   betamaps_<subj>_<p>.mat   the betas_nonzscored.mat betamaps for project <p>
%
% WHY THIS EXISTS. createTables.m has its own hardcoded list of the 8 observers who did
% both experiments, and it assumes every observer has BOTH a dg and a da session -- it
% would fail looking for the missing da. dg was run on 13. This builds the five dg-only
% observers' rows so they can be appended, filling the polexp_* columns with NaN, and
% it needs no mount because everything it reads is in ~/dg_collect.
%
% FAITHFULNESS IS CHECKED, NOT ASSUMED. VALIDATE_TABLE_ROWS rebuilds an observer who is
% already in the table and compares column by column. Do not trust this function on a
% new observer without running that first.
%
% See ../../AGENTS.md standing fact 7.

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'root'),     opts.root = dg_collect_dir(); end
    if ~isfield(opts,'projects'), opts.projects = {'dg'}; end

    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(fileparts(thisDir)), 'AnalysisCode', 'general'));  % map_theta

    % Column layout, in the CSV's own order.
    dgStim = {'cartexp_vertical_grating_rightwards_motion','cartexp_horizontal_grating_upwards_motion', ...
        'cartexp_vertical_grating_leftwards_motion','cartexp_horizontal_grating_downwards_motion', ...
        'cartexp_leftleaning_grating_upperrightwards_motion','cartexp_rightleaning_grating_upperleftwards_motion', ...
        'cartexp_leftleaning_grating_lowerleftwards_motion','cartexp_rightleaning_grating_lowerrightwards_motion', ...
        'cartexp_horizontal_stationary','cartexp_vertical_stationary','cartexp_rightleaning_grating_stationary', ...
        'cartexp_leftleaning_grating_stationary','cartexp_blank'};
    daStim = {'polexp_pinwheel_grating_clockwise_motion','polexp_annulus_grating_outwards_motion', ...
        'polexp_pinwheel_grating_cclockwise_motion','polexp_annulus_grating_inwards_motion', ...
        'polexp_ccspiral_grating_clockoutwards_motion','polexp_cspiral_grating_cclockoutwards_motion', ...
        'polexp_ccspiral_grating_cclockinwards_motion','polexp_cspiral_grating_clockinwards_motion', ...
        'polexp_annulus_grating_stationary','polexp_pinwheel_grating_stationary','polexp_cspiral_grating_stationary', ...
        'polexp_ccspiral_grating_stationary','polexp_blank'};

    % ROI names as createTables.m walks them, with its two renames. ORDER MATTERS: later
    % labels overwrite earlier ones where they overlap, so changing this list changes
    % which area a shared vertex reports.
    %
    % V2d IS DELIBERATELY OMITTED, and it should never be added back.
    %
    % left_V2d is not an analysis ROI. It was requested (JW) as a DIAGNOSTIC, to check
    % the polar-angle convention: V1/V2/V3 are bilateral and span dorsal and ventral, so
    % between them they cover the whole visual field and cannot discriminate one angle
    % convention from another. LEFT V2d is dorsal and one hemisphere, so it should
    % represent only the LOWER RIGHT visual field -- 270-360 deg if angles are
    % conventional (0 at the right horizontal, increasing counterclockwise), and 90-180
    % if they are still Benson. It came out at 96-99% in 270-360 with a circular mean
    % near 315, which is the confirmation recorded in ../STIMULUS_CONVENTIONS.md section
    % 3, point 3. Its job was to answer that question, and it did.
    %
    % It is also probably ATLAS-derived rather than hand-drawn, unlike every other label
    % here, so its boundaries need not agree with the hand-drawn V2. createTables.m lists
    % it last, which means applying it would OVERWRITE hand-drawn V2 with atlas
    % boundaries -- 1,190 V2 vertices and 66 V1 vertices for sub-0037 alone. The shipped
    % allsubjectsTable.csv contains zero left_V2d rows, which is the correct state, not
    % an omission to be reproduced for consistency's sake.
    %
    % opts.roinames can override this, but adding V2d to an analysis table is a mistake
    % rather than a preference.
    roinames = {'V1','V2','V3','hV4','V3a','V3b','hMTcomplex','pMT','pMST'};
    if isfield(opts,'roinames') && ~isempty(opts.roinames), roinames = opts.roinames; end

    R = load(fullfile(opts.root, sprintf('ret_%s.mat', subject)), ...
             'angle_adj','eccen','vexpl','sigma','nLH','nRH');
    nLH = double(R.nLH);
    N   = double(R.nLH) + double(R.nRH);

    T = table(repmat(string(subject), N, 1), strings(N,1), nan(N,1), nan(N,1), ...
              nan(N,1), nan(N,1), nan(N,1), false(N,1), ...
              'VariableNames', {'subject','visual_area','pRF_angle_bin','pRF_angle', ...
                                'pRF_ecc','pRF_r2','pRF_sigma','included'});

    % map_theta converts Benson to conventional degrees. createTables.m applies it to
    % angle_adj -- NOT to angle -- and everything downstream assumes conventional.
    T.pRF_angle = map_theta(double(R.angle_adj(:))')';
    T.pRF_ecc   = double(R.eccen(:));
    T.pRF_r2    = double(R.vexpl(:));
    T.pRF_sigma = double(R.sigma(:));

    binCenters = 0:45:315;
    circDist = abs(mod(T.pRF_angle - binCenters + 180, 360) - 180);
    [minDist, idx] = min(circDist, [], 2);
    T.pRF_angle_bin = nan(N,1);
    valid = minDist <= 22.5;
    T.pRF_angle_bin(valid) = binCenters(idx(valid));

    T.included = ~isnan(T.pRF_angle_bin) & ...
                 T.pRF_ecc >= 4 & T.pRF_ecc <= 8 & T.pRF_r2 >= 0.1;

    L = load(fullfile(opts.root, sprintf('labels_%s.mat', subject)));
    for ri = 1:numel(roinames)
        roiname = roinames{ri};
        fl = sprintf('lh_%s_REmanual', roiname);
        fr = sprintf('rh_%s_REmanual', roiname);
        if ~isfield(L, fl), continue; end            % absent label: leave those vertices as-is
        % DO NOT add 1 here. The .label files are 0-based, but COLLECT_LABELS already
        % stored them as raw(:,1) + 1, i.e. MATLAB 1-based -- unlike createTables.m,
        % which adds the 1 itself because it reads the raw files. Adding it twice
        % shifts every vertex by one, which does not change any label's SIZE and so
        % survives a count-based check; it shows up only as boundary vertices landing
        % in the neighbouring area. Only the right-hemisphere offset is applied here.
        if strcmp(roiname, 'V2d')
            label_idx = double(L.(fl)(:));            % left hemisphere only
            roiname   = 'left_V2d';
        else
            rh = [];
            if isfield(L, fr), rh = double(L.(fr)(:)) + nLH; end
            label_idx = [double(L.(fl)(:)); rh];
        end
        if strcmp(roiname,'pMT'),  roiname = 'MT';  end
        if strcmp(roiname,'pMST'), roiname = 'MST'; end
        T.visual_area(label_idx) = roiname;
    end

    for pj = {'dg','da'}
        pn = pj{1};
        stim = dgStim; if strcmp(pn,'da'), stim = daStim; end
        f = fullfile(opts.root, sprintf('betamaps_%s_%s.mat', subject, pn));
        if ~ismember(pn, opts.projects) || ~isfile(f)
            % No session for this observer in this experiment. NaN rather than zero:
            % zero is a response, NaN is "not measured", and every consumer filters on
            % finiteness. This is the case for the five dg-only observers.
            for k = 1:numel(stim), T.(stim{k}) = nan(N,1); end
            T.(sprintf('%s_beta_mean', pn)) = nan(N,1);
            T.(sprintf('%s_beta_std', pn))  = nan(N,1);
            continue
        end
        B = load(f, 'betamaps');
        if size(B.betamaps,1) ~= N
            error('build_table_rows:size', ...
                  '%s %s: betamaps has %d rows, retinotopy has %d vertices.', ...
                  subject, pn, size(B.betamaps,1), N);
        end
        for k = 1:numel(stim), T.(stim{k}) = B.betamaps(:, k); end
        T.(sprintf('%s_beta_mean', pn)) = mean(B.betamaps, 2);
        T.(sprintf('%s_beta_std',  pn)) = std(B.betamaps, 0, 2);
    end
end
