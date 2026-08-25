function collect_runwise_betas(varargin)
% COLLECT_RUNWISE_BETAS  Per-RUN condition betas for V1, from GLMsingle single-trial fits.
%
%   collect_runwise_betas()
%   collect_runwise_betas('subjects', {'sub-0037'}, 'force', true)
%
% WHY. Every within-observer error estimate in this project needs to resample the
% MEASUREMENT, which means resampling runs. Resampling vertices does not do it: the
% GLM betas are held fixed and only reshuffled, so it characterises which patch of V1
% was sampled rather than the reliability of the measurement (and it ignores spatial
% autocorrelation). See ../supplement/SUPPLEMENT_harmonic_model.md and
% ../cleanroom/diagnose_context_asymmetry.m.
%
% WHAT IS AVAILABLE. GLMsingle's TYPED_FITHRF_GLMDENOISE_RR.mat holds
%   modelmd : [nVertex x 1 x 1 x nTrials] SINGLE-TRIAL betas
% (this is GLMsingle, not the older GLMdenoise -- there is no 100-sample bootstrap
% dimension; single-trial betas are strictly more useful, since run-level or
% bootstrapped estimates can be built from them.) DESIGNINFO.mat holds stimorder
% (condition of each trial) and numtrialrun (trials per run). The design is fully
% balanced: 8 runs x 52 trials, 13 conditions x 32 trials, exactly 4 trials per
% condition per run, so a 4-vs-4 split-half over runs is balanced by construction.
%
% WHAT THIS SAVES. For each subject and experiment, the mean beta per (V1 vertex,
% condition, run): a [nV1 x 13 x 8] single array, about 3 MB, versus 425 MB for the
% source modelmd. That granularity supports split-half over runs, bootstrap over runs,
% and per-run reliability, without ever touching the network again.
%
% COST. modelmd is 425 MB per subject-experiment and the volume streams at about
% 3.3 MB/s, so each file takes roughly 2 minutes and all 16 about 35 minutes. The cost
% is bandwidth, not latency -- do NOT try to read scattered vertices through matfile,
% which is latency-bound and far slower (and matfile requires equally-spaced indices
% anyway). Load the whole variable and subset in memory.
%
% Requires /Volumes/Vision mounted; run under caffeinate, the mount drops on sleep.

    p = inputParser;
    p.addParameter('deriv', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
                             'data_bids/derivatives'], @ischar);
    p.addParameter('root',   '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('outdir', '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('subjects', {}, @iscell);
    p.addParameter('force', false, @(z) islogical(z) || isnumeric(z));
    p.parse(varargin{:});
    opt = p.Results;

    if ~isfolder(opt.deriv)
        error('collect_runwise_betas:mount', 'Not mounted: %s', opt.deriv);
    end
    subs = opt.subjects;
    if isempty(subs)
        thisDir = fileparts(mfilename('fullpath'));
        addpath(fullfile(fileparts(thisDir), 'cleanroom'));
        cfg  = config_repro();
        subs = cfg.subjects;
    end
    expMap = struct('dg','dgGLM','da','daGLM');

    for ii = 1:numel(subs)
        s = subs{ii};

        % V1 vertices, from the local retinotopy + label files
        Rt = load(fullfile(opt.root, sprintf('ret_%s.mat', s)), 'nLH');
        Lb = load(fullfile(opt.root, sprintf('labels_%s.mat', s)), ...
                  'lh_V1_REmanual', 'rh_V1_REmanual');
        nLH = double(Rt.nLH);
        v1  = sort([double(Lb.lh_V1_REmanual(:)); double(Lb.rh_V1_REmanual(:)) + nLH]);

        for e = {'dg','da'}
            en  = e{1};
            out = fullfile(opt.outdir, sprintf('runbetas_%s_%s.mat', s, en));
            if isfile(out) && ~opt.force
                fprintf('  %-14s %s  skipped (exists)\n', s, en);
                continue
            end
            t0 = tic;

            d = dir(fullfile(opt.deriv, expMap.(en), 'hRF_glmsingle', s, 'ses-*'));
            d = d([d.isdir]);
            if isempty(d)
                fprintf('  %-14s %s  NO ses-* -- skipped\n', s, en);  continue
            end
            src = fullfile(d(1).folder, d(1).name);

            DI = load(fullfile(src, 'DESIGNINFO.mat'), 'stimorder', 'numtrialrun');
            runOfTrial = repelem(1:numel(DI.numtrialrun), DI.numtrialrun);
            cond       = DI.stimorder(:).';
            nRun       = numel(DI.numtrialrun);
            nCond      = max(cond);

            % Stream the whole variable: bandwidth-bound, ~2 min. Subset after.
            M = load(fullfile(src, 'TYPED_FITHRF_GLMDENOISE_RR.mat'), 'modelmd');
            B = squeeze(M.modelmd);                       % nVertex x nTrials
            clear M
            if size(B,1) < max(v1)
                error('collect_runwise_betas:size', ...
                      '%s %s: modelmd has %d vertices, V1 index max is %d.', ...
                      s, en, size(B,1), max(v1));
            end
            if size(B,2) ~= numel(cond)
                error('collect_runwise_betas:trials', ...
                      '%s %s: modelmd has %d trials, DESIGNINFO has %d.', ...
                      s, en, size(B,2), numel(cond));
            end
            B = B(v1, :);

            % mean beta per (vertex, condition, run)
            runBeta = nan(numel(v1), nCond, nRun, 'single');
            nTrials = zeros(nCond, nRun);
            for c = 1:nCond
                for r = 1:nRun
                    sel = (cond == c) & (runOfTrial == r);
                    nTrials(c,r) = nnz(sel);
                    if any(sel), runBeta(:,c,r) = mean(B(:,sel), 2); end
                end
            end

            S = struct('subject', s, 'project', en, 'sourceFolder', src, ...
                       'v1Index', v1, 'nLH', nLH, 'runBeta', runBeta, ...
                       'nTrialsPerCondRun', nTrials, 'nRun', nRun, 'nCond', nCond, ...
                       'note', ['runBeta is nV1 x nCond x nRun, mean single-trial beta. ' ...
                                'Conditions 1-8 motion, 9=S_0 10=S_90 11=S_45 12=S_135, 13=blank.'], ...
                       'collected', datetime('now'));
            save(out, '-struct', 'S', '-v7.3');
            fprintf('  %-14s %s  %d V1 verts, %d cond x %d runs, %.0f s -> %s\n', ...
                    s, en, numel(v1), nCond, nRun, toc(t0), out);
        end
    end
    fprintf('collect_runwise_betas: done.\n');
end
