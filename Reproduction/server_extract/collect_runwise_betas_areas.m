function collect_runwise_betas_areas(varargin)
% COLLECT_RUNWISE_BETAS_AREAS  Per-RUN condition betas for EIGHT visual areas.
%
%   collect_runwise_betas_areas()
%   collect_runwise_betas_areas('subjects', {'sub-0037'}, 'experiments', {'dg'})
%   collect_runwise_betas_areas('dryrun', true)     % check labels, stream nothing
%
% The multi-area version of COLLECT_RUNWISE_BETAS. Same source, same arithmetic; it
% keeps V1 V2 V3 V3a V3b hV4 MT MST instead of V1 alone, and records which areas each
% retained vertex belongs to.
%
% WHY ALL AREAS AT ONCE. The expensive part is streaming modelmd -- 425 MB per
% subject-experiment, about 6.8 GB over the 16 files -- and that cost is IDENTICAL
% whether one area is kept or eight, because the whole variable has to come across the
% wire before anything can be subset (see the COST note in COLLECT_RUNWISE_BETAS: do
% not try to read scattered vertices through matfile, it is latency-bound and far
% slower). Extracting V1 now and V2/V3 later means paying 6.8 GB twice. The output
% grows from about 3 MB to roughly 10 MB per file, which is nothing.
%
% AREA NAMING. Output names follow Support/allsubjectsTable.csv, so anything keyed on
% its visual_area column joins without translation: V3a and V3b are lower-case, and
% MT / MST are the pMT / pMST labels. The label field actually read is in AREA_LABELS
% below. hMTcomplex is deliberately not included -- it overlaps pMT and pMST, and
% keeping all three would make the membership mask ambiguous for no gain.
%
% OVERLAPS ARE PRESERVED, NOT RESOLVED. Vertices are stored once, as a sorted union,
% with a logical nVertex x nArea membership mask. A vertex claimed by two labels stays
% in both columns rather than being assigned to whichever came first; the count of
% such vertices is reported per subject so it can be judged rather than discovered
% later. Select an area downstream with vertIndex(areaMask(:, strcmp(areaNames,'V2'))).
%
% PRE-FLIGHT. Every label is checked for every subject BEFORE any streaming starts, so
% a missing ROI fails in seconds rather than after half an hour of transfer. Use
% 'dryrun' to run only that check.
%
% TIMING A SINGLE FILE FIRST is worthwhile if the mount is slow or remote:
%   collect_runwise_betas_areas('subjects', {'sub-0037'}, 'experiments', {'dg'})
% prints the elapsed seconds and the achieved MB/s, which extrapolates to all 16.
%
% Requires the Vision volume mounted; run under caffeinate, the mount drops on sleep.
% Writes runbetas_areas_<subject>_<exp>.mat, leaving the V1-only runbetas_*.mat alone.

    AREA_LABELS = { ...
        'V1',  'V1'   ; ...
        'V2',  'V2'   ; ...
        'V3',  'V3'   ; ...
        'V3a', 'V3a'  ; ...
        'V3b', 'V3b'  ; ...
        'hV4', 'hV4'  ; ...
        'MT',  'pMT'  ; ...
        'MST', 'pMST' };

    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(thisDir), 'cleanroom'));

    p = inputParser;
    p.addParameter('deriv', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
                             'data_bids/derivatives'], @ischar);
    p.addParameter('root',   dg_collect_dir(), @ischar);
    p.addParameter('outdir', dg_collect_dir(), @ischar);
    p.addParameter('subjects', {}, @iscell);
    p.addParameter('experiments', {'dg','da'}, @iscell);
    p.addParameter('force', false, @(z) islogical(z) || isnumeric(z));
    p.addParameter('dryrun', false, @(z) islogical(z) || isnumeric(z));
    p.parse(varargin{:});
    opt = p.Results;

    subs = opt.subjects;
    if isempty(subs), cfg = config_repro(); subs = cfg.subjects; end
    areaNames = AREA_LABELS(:,1).';
    nA = numel(areaNames);

    % ---- pre-flight: labels for every subject, before any streaming ----------
    fprintf('pre-flight: labels for %d subjects x %d areas\n', numel(subs), nA);
    V = cell(numel(subs), 1);  A = cell(numel(subs), 1);  nLHs = nan(numel(subs),1);
    bad = {};
    for ii = 1:numel(subs)
        s = subs{ii};
        lf = fullfile(opt.root, sprintf('labels_%s.mat', s));
        rf = fullfile(opt.root, sprintf('ret_%s.mat', s));
        if ~isfile(lf) || ~isfile(rf)
            bad{end+1} = sprintf('%s: missing labels_/ret_ in %s', s, opt.root); continue %#ok<AGROW>
        end
        Lb = load(lf);  Rt = load(rf, 'nLH');
        nLH = double(Rt.nLH);  nLHs(ii) = nLH;

        sets = cell(1, nA);  missing = {};
        for a = 1:nA
            fl = ['lh_' AREA_LABELS{a,2} '_REmanual'];
            fr = ['rh_' AREA_LABELS{a,2} '_REmanual'];
            if ~isfield(Lb, fl) || ~isfield(Lb, fr)
                missing{end+1} = areaNames{a};  sets{a} = []; continue %#ok<AGROW>
            end
            sets{a} = [double(Lb.(fl)(:)); double(Lb.(fr)(:)) + nLH];
        end
        if ~isempty(missing)
            bad{end+1} = sprintf('%s: no label for %s', s, strjoin(missing, ' ')); %#ok<AGROW>
        end

        v = unique(vertcat(sets{:}));
        m = false(numel(v), nA);
        for a = 1:nA, m(:,a) = ismember(v, sets{a}); end
        V{ii} = v;  A{ii} = m;

        pairs = [areaNames; num2cell(sum(m,1))];
        fprintf('  %-14s %6d vertices  (', s, numel(v));
        fprintf('%s %d  ', pairs{:});
        fprintf(') overlap %d\n', nnz(sum(m,2) > 1));
    end
    if ~isempty(bad)
        error('collect_runwise_betas_areas:labels', ...
              'pre-flight failed:\n  %s', strjoin(bad, '\n  '));
    end
    fprintf('pre-flight OK\n\n');
    if opt.dryrun, fprintf('dryrun: stopping before any transfer.\n'); return; end

    expMap = struct('dg','dgGLM','da','daGLM');
    if ~isfolder(opt.deriv)
        error('collect_runwise_betas_areas:mount', 'Not mounted: %s', opt.deriv);
    end

    for ii = 1:numel(subs)
        s = subs{ii};  v = V{ii};  m = A{ii};  nLH = nLHs(ii);
        for e = opt.experiments(:).'
            en  = e{1};
            out = fullfile(opt.outdir, sprintf('runbetas_areas_%s_%s.mat', s, en));
            if isfile(out) && ~opt.force
                fprintf('  %-14s %s  skipped (exists)\n', s, en);  continue
            end
            t0 = tic;

            d = dir(fullfile(opt.deriv, expMap.(en), 'hRF_glmsingle', s, 'ses-*'));
            d = d([d.isdir]);
            if isempty(d), fprintf('  %-14s %s  NO ses-* -- skipped\n', s, en); continue; end
            src = fullfile(d(1).folder, d(1).name);

            DI = load(fullfile(src, 'DESIGNINFO.mat'), 'stimorder', 'numtrialrun');
            runOfTrial = repelem(1:numel(DI.numtrialrun), DI.numtrialrun);
            cond  = DI.stimorder(:).';
            nRun  = numel(DI.numtrialrun);
            nCond = max(cond);

            f  = fullfile(src, 'TYPED_FITHRF_GLMDENOISE_RR.mat');
            fi = dir(f);  mb = fi.bytes/1e6;
            M  = load(f, 'modelmd');
            tLoad = toc(t0);
            B = squeeze(M.modelmd);  clear M              % nVertex x nTrials

            if size(B,1) < max(v)
                error('collect_runwise_betas_areas:size', ...
                      '%s %s: modelmd has %d vertices, max area index is %d.', ...
                      s, en, size(B,1), max(v));
            end
            if size(B,2) ~= numel(cond)
                error('collect_runwise_betas_areas:trials', ...
                      '%s %s: modelmd has %d trials, DESIGNINFO has %d.', ...
                      s, en, size(B,2), numel(cond));
            end
            B = B(v, :);

            runBeta = nan(numel(v), nCond, nRun, 'single');
            nTrials = zeros(nCond, nRun);
            for c = 1:nCond
                for r = 1:nRun
                    sel = (cond == c) & (runOfTrial == r);
                    nTrials(c,r) = nnz(sel);
                    if any(sel), runBeta(:,c,r) = mean(B(:,sel), 2); end
                end
            end

            S = struct('subject', s, 'project', en, 'sourceFolder', src, ...
                       'vertIndex', v, 'areaMask', m, 'areaNames', {areaNames}, ...
                       'nLH', nLH, 'runBeta', runBeta, ...
                       'nTrialsPerCondRun', nTrials, 'nRun', nRun, 'nCond', nCond, ...
                       'note', ['runBeta is nVertex x nCond x nRun, mean single-trial ' ...
                                'beta, rows aligned with vertIndex. areaMask(:,k) says ' ...
                                'whether each vertex is in areaNames{k}; a vertex may be ' ...
                                'in more than one. Conditions 1-8 motion, 9=S_0 10=S_90 ' ...
                                '11=S_45 12=S_135, 13=blank.'], ...
                       'collected', datetime('now'));
            save(out, '-struct', 'S', '-v7.3');
            o = dir(out);
            fprintf(['  %-14s %s  %5d verts, %d cond x %d runs, ' ...
                     'load %.0f s (%.1f MB/s), total %.0f s, out %.1f MB\n'], ...
                    s, en, numel(v), nCond, nRun, tLoad, mb/tLoad, toc(t0), o.bytes/1e6);
        end
    end
    fprintf('collect_runwise_betas_areas: done.\n');
end
