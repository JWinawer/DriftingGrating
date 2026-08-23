function collect_timeseries(varargin)
% COLLECT_TIMESERIES  Preprocessed fsnative BOLD time series for the analysed vertices.
%
%   collect_timeseries('preflightOnly', true)          % seconds, writes nothing
%   collect_timeseries()                               % default: sub-0037, dg only
%   collect_timeseries('subjects', {'sub-0037'}, 'experiments', {'dg','da'})
%
% WHY. To quantify the Fig 4A run-mismatch control across all runs, and to calibrate
% ../cleanroom/run_mismatch_local.m, which estimates the mismatched-design R2 from local
% data but has to assume <e,yhat_m> = 0. That assumption is expected to be biased
% (event timing is shared across runs), so one session measured properly says whether
% the local estimate can be trusted for the rest.
%
% COST. The source files are fsnative .mgh -- UNCOMPRESSED -- at ~147 MB per hemisphere
% per run, so ~294 MB per run and ~2.35 GB per subject-experiment (8 runs); all 16
% subject-experiments are 37.6 GB. Measured mount throughput has ranged 0.8-3.3 MB/s
% (Abu Dhabi -> New York), i.e. roughly 50 min to 12 min per subject-experiment. Run
% PREFLIGHT first: it prints the byte count and an ETA at the currently measured rate.
%
% RESUMABLE, per (subject, experiment, RUN). Each run writes its own small output file
% and an existing one is skipped unless 'force' is set, so a dropped mount costs at most
% the run in flight (~6 min at 0.8 MB/s), never the whole job. Outputs are written with
% a .part suffix and renamed only on success, so a truncated file can never be mistaken
% for a complete one. A dropped mount returns EMPTY rather than erroring, so every read
% is size-checked against the expected vertex/TR count before it is accepted.
%
% OUTPUT, per run: ts_<subject>_<exp>_run<K>.mat holding
%   ts       nV1 x nTR single -- polynomial-projected % signal change, V1 vertices
%   v1Index  the vertex indices (same convention as runbetas_*.mat)
% about 7 MB per run, 56 MB per subject-experiment. Denoising mirrors
% AnalysisCode/02_ttave/createTTaveTable.m (computeConditionTTA_rawdata): GLMsingle's own
% per-run polynomial basis projected out, then % signal change against GLMsingle's meanvol.
%
% Requires /Volumes/Vision mounted; run under caffeinate (the mount drops on sleep).

    p = inputParser;
    p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids', @ischar);
    p.addParameter('root',    '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('outdir',  '/Users/jaw288/dg_collect/timeseries', @ischar);
    p.addParameter('subjects', {'sub-0037'}, @iscell);
    p.addParameter('experiments', {'dg'}, @iscell);
    p.addParameter('freesurfer', ['/Users/jaw288/repos/Code/Projects/' ...
                   'equivalent_input_noise_marc/external/freesurfer'], @ischar);
    p.addParameter('glmsingle', '/Users/jaw288/repos/Code/Toolboxes/GLMsingle', @ischar);
    p.addParameter('mbPerSec', 0.8, @isnumeric);      % for the ETA only
    p.addParameter('preflightOnly', false, @islogical);
    p.addParameter('force', false, @islogical);
    p.parse(varargin{:});
    opt = p.Results;

    addpath(opt.freesurfer);                      % MRIread: setup_user has no case for this Mac
    addpath(genpath(opt.glmsingle));              % constructpolynomialmatrix, projectionmatrix
    if isempty(which('MRIread')), error('collect_timeseries:mri', 'MRIread not on path'); end

    % ---------------- preflight ----------------
    fprintf('\n=============== PREFLIGHT ===============\n');
    if ~isfolder(opt.bidsDir)
        fprintf('[FAIL] bidsDir not reachable: %s\n', opt.bidsDir);
        fprintf('       Mount /Volumes/Vision and retry.\n=========================================\n'); return
    end
    fprintf('[ ok ] bidsDir reachable\n');

    jobs = {};  totalBytes = 0;  nSkip = 0;
    for si = 1:numel(opt.subjects)
        for ei = 1:numel(opt.experiments)
            s = opt.subjects{si};  en = opt.experiments{ei};
            fb = fullfile(opt.root, sprintf('runbetas_%s_%s.mat', s, en));
            if ~isfile(fb)
                fprintf('[FAIL] %s %s: no local runbetas (need v1Index) -- %s\n', s, en, fb); continue
            end
            B = load(fb, 'v1Index', 'nRun');
            ses = find_session(opt.bidsDir, s, en);
            if isempty(ses)
                fprintf('[FAIL] %s %s: no fmriprep session with task-%s\n', s, en, en); continue
            end
            for r = 1:B.nRun
                out = fullfile(opt.outdir, sprintf('ts_%s_%s_run%d.mat', s, en, r));
                if isfile(out) && ~opt.force, nSkip = nSkip + 1; continue; end
                f = hemi_files(opt.bidsDir, s, ses, en, r);
                if any(cellfun(@isempty, f))
                    fprintf('[FAIL] %s %s run %d: missing hemi file\n', s, en, r); continue
                end
                b = sum(cellfun(@(x) x.bytes, f));
                totalBytes = totalBytes + b;
                jobs(end+1,:) = {s, en, ses, r, f, out, double(B.v1Index)}; %#ok<AGROW>
            end
        end
    end

    gb = totalBytes / 1e9;
    fprintf('[ ok ] %d run(s) to fetch, %d already done (skipped)\n', size(jobs,1), nSkip);
    fprintf('[ ok ] %.2f GB to read; ETA %.1f h at %.1f MB/s\n', gb, gb*1000/opt.mbPerSec/3600, opt.mbPerSec);
    fprintf('-----------------------------------------\nPreflight OK.\n=========================================\n');
    if opt.preflightOnly || isempty(jobs), return; end
    if ~isfolder(opt.outdir), mkdir(opt.outdir); end

    % ---------------- fetch ----------------
    for k = 1:size(jobs,1)
        [s, en, ses, r, f, out, v1Index] = jobs{k,:};
        t0 = tic;
        fprintf('  %-14s %s run %d ... ', s, en, r);
        try
            [ts, nTR] = read_run(f, v1Index, opt, s, en, ses, r);
        catch e
            fprintf('FAILED (%s) -- rerun to resume\n', e.message);  continue
        end
        tmp = [out '.part'];
        S = struct('ts', single(ts), 'v1Index', v1Index, 'subject', s, ...
                   'experiment', en, 'run', r, 'nTR', nTR, 'collected', datetime('now'));
        save(tmp, '-struct', 'S', '-v7.3');
        movefile(tmp, out);
        fprintf('ok  %d vertices x %d TRs  (%.0f s)\n', size(ts,1), nTR, toc(t0));
    end
    fprintf('\nDone. Outputs in %s\n', opt.outdir);
end

% ------------------------------------------------------------------------
function [ts, nTR] = read_run(f, v1Index, opt, s, en, ses, r)
% Read both hemispheres, concatenate to the whole-surface vertex order, subset to V1,
% polynomial-denoise per run, convert to % signal change against GLMsingle's meanvol.
    D = cell(1,2);
    for h = 1:2
        m = MRIread(fullfile(f{h}.folder, f{h}.name));
        d = squeeze(m.vol);                                  % nVertexHemi x nTR
        if isempty(d) || size(d,2) < 10
            error('empty/truncated read (dropped mount?): %s', f{h}.name);
        end
        D{h} = d;
    end
    full = [D{1}; D{2}];                                     % lh then rh, GLM convention
    G = load(fullfile(opt.root, sprintf('glm_%s_%s.mat', s, en)), 'meanvol');
    if size(full,1) ~= numel(G.meanvol)
        error('vertex count %d ~= meanvol %d', size(full,1), numel(G.meanvol));
    end
    d   = double(full(v1Index, :));
    nTR = size(d, 2);

    maxpolydeg = round(((nTR * 1) / 60) / 2);                % TR = 1 s, as GLMsingle does
    pm = projectionmatrix(constructpolynomialmatrix(nTR, 0:maxpolydeg));
    d  = (pm * d.').';
    ts = ((d ./ double(G.meanvol(v1Index))) - 1) * 100;
end

% ------------------------------------------------------------------------
function ses = find_session(bidsDir, subj, en)
    ses = '';
    d = dir(fullfile(bidsDir, 'derivatives', 'fmriprep', subj, 'ses-*'));
    d = d([d.isdir]);
    for i = 1:numel(d)
        if ~isempty(dir(fullfile(d(i).folder, d(i).name, 'func', ...
                sprintf('*task-%s_*space-fsnative_hemi-L_bold.func.mgh', en))))
            ses = d(i).name;  return
        end
    end
end

% ------------------------------------------------------------------------
function f = hemi_files(bidsDir, subj, ses, en, r)
    f = cell(1,2);  hemis = {'L','R'};
    for h = 1:2
        d = dir(fullfile(bidsDir, 'derivatives', 'fmriprep', subj, ses, 'func', ...
            sprintf('*task-%s_*run-%d_space-fsnative_hemi-%s_bold.func.mgh', en, r, hemis{h})));
        if ~isempty(d), f{h} = d(1); end
    end
end
