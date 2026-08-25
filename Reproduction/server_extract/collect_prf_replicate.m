function collect_prf_replicate(varargin)
% COLLECT_PRF_REPLICATE  Copy the SECOND, independent pRF solution to ~/dg_collect.
%
%   collect_prf_replicate()
%   collect_prf_replicate('subjects', {'sub-0037'}, 'force', true)
%
% Every observer was fitted twice, independently, from two pRF runs with different
% stimuli. The derivatives tree holds them as two sibling folders:
%
%   prfvista_mov  the solution the whole analysis uses. Already local, as
%                 ret_<subject>.mat (verified identical: max|difference| = 0).
%   prfvista      the independent replicate. NOT local until this script runs.
%
% The disagreement between the two estimates the pRF measurement error, which is what
% DIAGNOSE_PRF_ANGLE_ERROR needs in order to bound the polar-angle error sigma. See
% ../supplement/SUPPLEMENT_harmonic_model.md for why sigma matters: angle error
% attenuates b1 in the polar
% experiment but not in the Cartesian one, so it inflates the cross-experiment gap.
%
% REQUIRES /Volumes/Vision mounted. That server is thousands of miles away and each
% .mgz costs a second or two, so this pulls 14 files per subject once (about 2-4
% minutes for all eight) and then never touches the network again. Run it under
% caffeinate -- the mount drops when the Mac sleeps.
%
% Writes ~/dg_collect/ret_prfvista_<subject>.mat with the same field names and the
% same lh-then-rh vertex order as ret_<subject>.mat, so the two drop into the same
% indexing without further care. Idempotent: existing outputs are skipped unless
% 'force' is set.

    p = inputParser;
    p.addParameter('deriv', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
                             'data_bids/derivatives'], @ischar);
    p.addParameter('source', 'prfvista', @ischar);
    p.addParameter('outdir', '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('subjects', {}, @iscell);
    p.addParameter('force', false, @(z) islogical(z) || isnumeric(z));
    p.parse(varargin{:});
    opt = p.Results;

    ensure_mriread();
    if ~isfolder(opt.deriv)
        error('collect_prf_replicate:mount', ...
              'Not mounted: %s. Mount /Volumes/Vision and retry.', opt.deriv);
    end
    if ~isfolder(opt.outdir), mkdir(opt.outdir); end

    subs = opt.subjects;
    if isempty(subs)
        thisDir = fileparts(mfilename('fullpath'));
        addpath(fullfile(fileparts(thisDir), 'cleanroom'));
        cfg  = config_repro();
        subs = cfg.subjects;
    end

    maps = {'x','y','eccen','angle','angle_adj','sigma','vexpl'};
    fprintf('collect_prf_replicate: %s -> %s\n', opt.source, opt.outdir);

    for ii = 1:numel(subs)
        s   = subs{ii};
        out = fullfile(opt.outdir, sprintf('ret_%s_%s.mat', opt.source, s));
        if isfile(out) && ~opt.force
            fprintf('  %-14s skipped (exists)\n', s);
            continue
        end

        t0 = tic;
        d  = dir(fullfile(opt.deriv, opt.source, s, 'ses-*'));
        d  = d([d.isdir]);
        if isempty(d)
            fprintf('  %-14s NO ses-* FOUND under %s/%s -- skipped\n', s, opt.source, s);
            continue
        end
        src = fullfile(d(1).folder, d(1).name);

        S = struct('subject', s, 'sourceFolder', src, 'scope', 'whole surface', ...
                   'source', opt.source, 'collected', datetime('now'));
        okAll = true;
        for mi = 1:numel(maps)
            m = maps{mi};
            fl = fullfile(src, sprintf('lh.%s.mgz', m));
            fr = fullfile(src, sprintf('rh.%s.mgz', m));
            if ~isfile(fl) || ~isfile(fr)
                fprintf('  %-14s missing %s -- skipped\n', s, m);
                okAll = false; break
            end
            vl = grab(fl);  vr = grab(fr);
            if mi == 1
                S.nLH = numel(vl);  S.nRH = numel(vr);
            elseif numel(vl) ~= S.nLH || numel(vr) ~= S.nRH
                error('collect_prf_replicate:vertexcount', ...
                      '%s: %s has %d/%d vertices, expected %d/%d.', ...
                      s, m, numel(vl), numel(vr), S.nLH, S.nRH);
            end
            S.(m) = [vl; vr];
        end
        if ~okAll, continue; end

        S.mapsSaved = maps;
        save(out, '-struct', 'S');
        fprintf('  %-14s %d vertices, %d maps, %.1f s -> %s\n', ...
                s, S.nLH + S.nRH, numel(maps), toc(t0), out);
    end
    fprintf('collect_prf_replicate: done.\n');
end

% ------------------------------------------------------------------------
function v = grab(f)
    m = MRIread(f);
    v = double(m.vol(:));
end

% ------------------------------------------------------------------------
function ensure_mriread()
    if exist('MRIread', 'file') == 2, return; end
    fsPath = ['/Users/jaw288/repos/Code/Projects/equivalent_input_noise_marc/' ...
              'external/freesurfer'];
    if isfolder(fsPath)
        addpath(fsPath);
    else
        error('collect_prf_replicate:mriread', ...
              'MRIread not on the path and %s does not exist.', fsPath);
    end
end
