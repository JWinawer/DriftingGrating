function collect_everything(bidsDir, expOutDir, outDir, opts)
% COLLECT_EVERYTHING  One liberal pass over the server. Ask for everything cheap.
%
%   collect_everything
%   collect_everything(bidsDir, expOutDir, outDir, opts)
%
% SELF-CONTAINED: no dependency on the DriftingGrating analysis code beyond the helpers
% already on the path for meanWithinLabel.m -- setup_user, read_label, MRIread.
%
% WHY THIS EXISTS, AND WHY IT FILTERS NOTHING
% -------------------------------------------
% An earlier script (extract_for_transfer.m) restricted to V1, 4-8 deg, pRF R2 > 0.1
% before saving. Every one of those choices later turned out to be a question we wanted
% to ask rather than an answer we could assume:
%
%   * 4-8 deg is a STIMULUS-MATCHING constraint (cartesian gratings hold spatial
%     frequency constant across the aperture, polar gratings' SF scales inversely with
%     eccentricity, so the two match only near 6 deg). It has nothing to do with where
%     the stimulus was -- the annulus is roughly 1-12 deg -- and nothing to do with
%     whether a GLM fit is sound.
%   * Thresholding pRF R2 selects vertices with good SNR, which raises GLM R2 whether or
%     not the vertex responded. Comparing a filtered patch against an unfiltered baseline
%     is biased upward by the selection alone.
%   * Restricting to V1 left no clean non-V1 baseline to compare V1 against.
%
% Each discovery cost a round trip to a server thousands of miles away. So this script
% takes the opposite approach: SAVE WHOLE-SURFACE, FILTER NOTHING, and let every
% restriction be a one-line choice made afterwards, locally, where it can be varied.
%
% This is affordable because the expensive part of results.mat is `modelmd` (nVertices x
% 416 single, ~450 MB, ~90% of each file) and we do not need it. Everything else is
% roughly 15 MB per subject x experiment.
%
% WHAT IT COLLECTS
%   1. GLMsingle quality, WHOLE SURFACE, per subject x experiment: R2, R2run, FRACvalue,
%      noisepool, HRFindex, meanvol, xvaltrend, pcnum, and any other non-huge field.
%   2. Every prfvista_mov retinotopy map, WHOLE SURFACE as full vectors (not medians).
%   3. Every FreeSurfer label found for the subject, as vertex index lists.
%   4. The experimental design and timing files per run (small ones only -- the ~36 MB
%      S09_const_file_* are skipped, see opts.includeConstFiles).
%   5. A manifest of everything SEEN, including what was skipped and why, so the next
%      question does not need another trip.
%
% Output is v7.3 (HDF5) so single variables can be read locally without loading whole
% files. Do not read these over a network mount -- HDF5 partial reads make many small
% dependent round trips and are pathologically slow over a high-latency link. Copy the
% folder across with rsync/ftp (bulk streaming, which long-haul links handle well) and
% read it from local disk.
%
% Expected size, measured against sub-0255's real results.mat: 71 MB per subject x
% experiment (of which ~30 MB is the 31 contrast maps), so ~1.15 GB for the GLM output,
% plus ~60 MB of retinotopy and a few MB of labels and design files. Call it 1.2 GB.
% With opts.includeBetas it is roughly 8 GB, because modelmd alone is 506 MB of each
% 528 MB source file -- 96%, not the 90% previously assumed.
%
% Nothing here modifies anything on the server. It only reads.

    if nargin < 1 || isempty(bidsDir)
        bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
    end
    if nargin < 2 || isempty(expOutDir)
        expOutDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/experimentalOutput';
    end
    if nargin < 3 || isempty(outDir)
        outDir = fullfile(pwd, 'dg_collect');
    end
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'includeBetas'),      opts.includeBetas      = false; end
    if ~isfield(opts,'includeConstFiles'), opts.includeConstFiles = false; end

    subjects = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123', ...
                'sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
    projects = {'dg','da'};
    hRF_setting = 'glmsingle';

    if ~isfolder(outDir), mkdir(outDir); end
    if ~isfolder(bidsDir)
        error('bidsDir not reachable:\n  %s\nEdit the default at the top of this file.', bidsDir);
    end
    try
        setup_user('rania', bidsDir);
    catch
        warning('setup_user failed -- MRIread/read_label may be missing. Continuing.');
    end

    fprintf('\ncollect_everything\n  bids  : %s\n  expOut: %s\n  out   : %s\n', ...
            bidsDir, expOutDir, outDir);
    fprintf('  betas : %d   constFiles: %d\n\n', opts.includeBetas, opts.includeConstFiles);

    manifest = {};   % {subject, project, item, status, detail, sizeMB}

    for si = 1:numel(subjects)
        subj = subjects{si};
        fprintf('%s\n', subj);

        manifest = [manifest; collect_retinotopy(bidsDir, subj, outDir)]; %#ok<AGROW>
        manifest = [manifest; collect_labels(bidsDir, subj, outDir)];     %#ok<AGROW>

        for pj = 1:numel(projects)
            proj = projects{pj};
            manifest = [manifest; collect_glm(bidsDir, subj, proj, hRF_setting, outDir, opts)]; %#ok<AGROW>
            manifest = [manifest; collect_design(expOutDir, subj, proj, outDir, opts)];         %#ok<AGROW>
        end
    end

    M = cell2table(manifest, 'VariableNames', ...
        {'subject','project','item','status','detail','sizeMB'});
    writetable(M, fullfile(outDir,'manifest.csv'));

    total = sum(M.sizeMB(isfinite(M.sizeMB)));
    fprintf('\n---------------------------------------------\n');
    fprintf('%d items, %d ok, %d missing/skipped\n', height(M), ...
            sum(strcmp(M.status,'ok')), sum(~strcmp(M.status,'ok')));
    fprintf('Total written: %.0f MB\n', total);
    fprintf('Zip or rsync and return: %s\n', outDir);
    fprintf('Check manifest.csv for any row whose status is not "ok".\n\n');
end

% =====================================================================
function rows = collect_glm(bidsDir, subj, proj, hRF_setting, outDir, opts)
% Whole-surface GLMsingle output, everything except the betas.
    rows = {};
    glmDir = fullfile(bidsDir,'derivatives',strcat(proj,'GLM'), ...
                      strcat('hRF_',hRF_setting), subj);
    f = dir(fullfile(glmDir,'**','results.mat'));
    if isempty(f)
        fprintf('    %-3s glm       MISSING under %s\n', proj, glmDir);
        rows(end+1,:) = {subj, proj, 'glm', 'missing', glmDir, NaN};
        return
    end

    S = load(fullfile(f(1).folder, f(1).name));
    if ~isfield(S,'results') || ~isfield(S.results,'allevents')
        fprintf('    %-3s glm       UNEXPECTED structure\n', proj);
        rows(end+1,:) = {subj, proj, 'glm', 'bad-structure', f(1).folder, NaN};
        return
    end
    a = S.results.allevents;

    g = struct('subject',subj,'project',proj,'sourceFolder',f(1).folder, ...
               'hRF_setting',hRF_setting,'scope','WHOLE SURFACE, unfiltered');
    fn = fieldnames(a);
    kept = {}; skipped = {};
    for k = 1:numel(fn)
        v = a.(fn{k});
        if strcmp(fn{k},'modelmd') && ~opts.includeBetas
            skipped{end+1} = fn{k}; %#ok<AGROW>
            continue
        end
        if isnumeric(v) && ~islogical(v), v = single(v); end
        g.(fn{k}) = v;
        kept{end+1} = fn{k}; %#ok<AGROW>
    end
    g.fieldsKept = kept; g.fieldsSkipped = skipped;

    % Anything else sitting alongside results.allevents is worth having too. This is
    % NOT a formality: results.contrasts is a struct of ~31 whole-surface contrast maps
    % (allmVblank, cardsVblank, s0Vb, ...), i.e. the analysis-relevant maps themselves,
    % and no previous extraction ever brought them off the server. Handle struct-valued
    % top-level fields explicitly -- an isnumeric() test silently drops them.
    otherTop = setdiff(fieldnames(S.results), {'allevents'});
    for k = 1:numel(otherTop)
        v = S.results.(otherTop{k});
        if isstruct(v) && isscalar(v)
            sub = struct(); sf = fieldnames(v);
            for j = 1:numel(sf)
                w = v.(sf{j});
                if isnumeric(w) && ~islogical(w), w = single(w); end
                sub.(sf{j}) = w;
            end
            g.(['top_' otherTop{k}]) = sub;
        elseif isnumeric(v) || islogical(v) || iscell(v)
            if isnumeric(v) && ~islogical(v), v = single(v); end
            g.(['top_' otherTop{k}]) = v;
        end
    end
    g.resultsTopFields = fieldnames(S.results);

    outFile = fullfile(outDir, sprintf('glm_%s_%s.mat', subj, proj));
    save(outFile, '-struct', 'g', '-v7.3');
    mb = fileMB(outFile);
    fprintf('    %-3s glm       ok  %5.0f MB  (%d fields, skipped: %s)\n', ...
            proj, mb, numel(kept), strjoin(skipped,','));
    rows(end+1,:) = {subj, proj, 'glm', 'ok', f(1).folder, mb};

    % Note whether the rawer modelOutput.mat exists, without reading it.
    mo = dir(fullfile(f(1).folder,'modelOutput.mat'));
    if isempty(mo)
        rows(end+1,:) = {subj, proj, 'modelOutput', 'absent', f(1).folder, NaN};
    else
        rows(end+1,:) = {subj, proj, 'modelOutput', 'present-not-copied', ...
                         fullfile(mo(1).folder,mo(1).name), mo(1).bytes/1e6};
    end
end

% ---------------------------------------------------------------------
function rows = collect_retinotopy(bidsDir, subj, outDir)
% EVERY prfvista_mov map, as full whole-surface vectors.
    rows = {};
    d = dir(fullfile(bidsDir,'derivatives','prfvista_mov',subj,'**','stimfiles.mat'));
    if isempty(d)
        d = dir(fullfile(bidsDir,'derivatives','prfvista_mov',subj,'**','lh.eccen.mgz'));
    end
    if isempty(d)
        fprintf('    ret       MISSING prfvista_mov\n');
        rows(end+1,:) = {subj, '', 'retinotopy', 'missing', '', NaN};
        return
    end
    retDir = d(1).folder;

    lhFiles = dir(fullfile(retDir,'lh.*.mgz'));
    R = struct('subject',subj,'sourceFolder',retDir, ...
               'mapsFound',{{lhFiles.name}},'scope','WHOLE SURFACE, unfiltered');
    ok = {};
    for k = 1:numel(lhFiles)
        nm = erase(erase(lhFiles(k).name,'lh.'),'.mgz');
        try
            lh = MRIread(fullfile(retDir,['lh.' nm '.mgz']));
            rh = MRIread(fullfile(retDir,['rh.' nm '.mgz']));
            R.(matlab.lang.makeValidName(nm)) = single([lh.vol(:); rh.vol(:)]);
            R.nLH = numel(lh.vol); R.nRH = numel(rh.vol);   % hemisphere sizes, no surfaces needed
            ok{end+1} = nm; %#ok<AGROW>
        catch ME
            rows(end+1,:) = {subj, '', ['ret:' nm], 'read-failed', ME.message, NaN}; %#ok<AGROW>
        end
    end
    R.mapsSaved = ok;

    outFile = fullfile(outDir, sprintf('ret_%s.mat', subj));
    save(outFile, '-struct', 'R', '-v7.3');
    mb = fileMB(outFile);
    fprintf('    ret       ok  %5.0f MB  (%s)\n', mb, strjoin(ok,','));
    rows(end+1,:) = {subj, '', 'retinotopy', 'ok', retDir, mb};
end

% ---------------------------------------------------------------------
function rows = collect_labels(bidsDir, subj, outDir)
% EVERY label under the subject's FreeSurfer label dir, as vertex index lists.
% Not just V1 -- other ROIs cost nothing and we have already wanted V2 once.
    rows = {};
    labelDir = fullfile(bidsDir,'derivatives','freesurfer',subj,'label');
    if ~isfolder(labelDir)
        fprintf('    labels    MISSING %s\n', labelDir);
        rows(end+1,:) = {subj, '', 'labels', 'missing', labelDir, NaN};
        return
    end
    f = dir(fullfile(labelDir,'**','*.label'));
    L = struct('subject',subj,'sourceFolder',labelDir,'labelsFound',{{f.name}});
    names = {};
    for k = 1:numel(f)
        full = fullfile(f(k).folder, f(k).name);
        try
            raw = readLabelFile(full);
            key = matlab.lang.makeValidName(erase(f(k).name,'.label'));
            % keep the sub-path too, since lh.V1 can exist in several subfolders
            rel = erase(f(k).folder, labelDir);
            L.(key) = uint32(raw(:,1) + 1);          % label files are 0-based
            L.([key '_from']) = rel;
            names{end+1} = f(k).name; %#ok<AGROW>
        catch ME
            rows(end+1,:) = {subj, '', ['label:' f(k).name], 'read-failed', ME.message, NaN}; %#ok<AGROW>
        end
    end
    L.labelsSaved = names;

    outFile = fullfile(outDir, sprintf('labels_%s.mat', subj));
    save(outFile, '-struct', 'L', '-v7.3');
    mb = fileMB(outFile);
    fprintf('    labels    ok  %5.0f MB  (%d labels)\n', mb, numel(names));
    rows(end+1,:) = {subj, '', 'labels', 'ok', labelDir, mb};
end

% ---------------------------------------------------------------------
function rows = collect_design(expOutDir, subj, proj, outDir, opts)
% Per-run design / timing / response files. These live in a SEPARATE tree from
% data_bids (setup.json calls it expoutputdir) and answer whether a flat session is a
% processing error -- wrong design matrix, misaligned timing -- rather than a real
% non-response. The ~36 MB S09_const_file_* are skipped by default.
    rows = {};
    if ~isfolder(expOutDir)
        rows(end+1,:) = {subj, proj, 'design', 'missing-tree', expOutDir, NaN};
        return
    end
    % Session folders are named per subject in a way we cannot assume, so glob widely
    % and keep whatever design files turn up for this subject.
    pat = {'*design*.mat','*timestamps*.mat','*responses*.mat','runlog.txt'};
    if opts.includeConstFiles, pat{end+1} = '*const_file*.mat'; end

    found = [];
    for k = 1:numel(pat)
        found = [found; dir(fullfile(expOutDir, proj, '**', pat{k}))]; %#ok<AGROW>
    end
    if isempty(found)
        fprintf('    %-3s design    none found under %s\n', proj, fullfile(expOutDir,proj));
        rows(end+1,:) = {subj, proj, 'design', 'none-found', fullfile(expOutDir,proj), NaN};
        return
    end

    dst = fullfile(outDir, 'design', proj);
    if ~isfolder(dst), mkdir(dst); end
    n = 0; mb = 0;
    for k = 1:numel(found)
        src = fullfile(found(k).folder, found(k).name);
        rel = strrep(erase(found(k).folder, expOutDir), filesep, '_');
        try
            copyfile(src, fullfile(dst, sprintf('%s_%s', rel, found(k).name)));
            n = n + 1; mb = mb + found(k).bytes/1e6;
        catch ME
            rows(end+1,:) = {subj, proj, ['design:' found(k).name], 'copy-failed', ME.message, NaN}; %#ok<AGROW>
        end
    end
    fprintf('    %-3s design    ok  %5.0f MB  (%d files)\n', proj, mb, n);
    rows(end+1,:) = {subj, proj, 'design', 'ok', fullfile(expOutDir,proj), mb};
end

% ---------------------------------------------------------------------
function v = readLabelFile(f)
% Minimal FreeSurfer .label reader, so this does not depend on read_label being
% on the path (and does not need the subject-relative naming read_label expects).
    fid = fopen(f,'r');
    if fid < 0, error('cannot open %s', f); end
    c = onCleanup(@() fclose(fid));
    fgetl(fid);                      % comment line
    n = str2double(fgetl(fid));      % vertex count
    v = fscanf(fid, '%d %f %f %f %f', [5 n])';
end

% ---------------------------------------------------------------------
function mb = fileMB(f)
    d = dir(f);
    if isempty(d), mb = NaN; else, mb = d(1).bytes/1e6; end
end
