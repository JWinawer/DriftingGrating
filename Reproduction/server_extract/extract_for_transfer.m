function extract_for_transfer(bidsDir, outDir)
% EXTRACT_FOR_TRANSFER  Pull everything the Fig 7 follow-up needs off the server.
%
%   extract_for_transfer                       % uses the default bidsDir below
%   extract_for_transfer(bidsDir, outDir)
%
% SELF-CONTAINED: this file has no dependencies on the DriftingGrating analysis code.
% It needs only the helpers already on the path for meanWithinLabel.m -- setup_user,
% read_label, get_surfsize, MRIread.
%
% It writes a small folder (a few MB total) that can be zipped and sent back. It does
% NOT need the ~500 MB results.mat files to be copied anywhere: it reads them in place
% and keeps only the quality fields, discarding `modelmd`, which is ~90% of each file.
%
% WHAT IT COLLECTS, per subject x experiment:
%   (a) GLMsingle fit quality -- R2, R2run, FRACvalue, noisepool, HRFindex, meanvol.
%       Used to decide whether sub-0201 and sub-0037 should be in the polar analysis.
%       Two specific things to look for: a bad run in sub-0201's R2run, and a da-only
%       quality drop for sub-0037 (which responds normally in dg).
%   (b) An inventory of the prfvista_mov retinotopy outputs, plus the value of every
%       .mgz map found there in the analysed V1 patch. We need a per-observer BOLD gain
%       estimated INDEPENDENTLY of the 13 GLM conditions, and the retinotopy scan is a
%       separate session, so it is the natural candidate -- but we do not know from here
%       which maps that pipeline saves. The inventory answers that in the same trip.
%
% Everything is restricted to the V1 label so the FreeSurfer labels and retinotopy
% maps never have to be copied either. Each vertex's eccentricity comes along in
% `patchEccen`, so the published 4-8 deg band -- or any other -- is a local filter rather
% than a reason to re-run this. Whole-surface percentile summaries are kept alongside so
% no context is lost.
%
% Nothing here modifies anything on the server. It only reads.

    if nargin < 1 || isempty(bidsDir)
        bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
    end
    if nargin < 2 || isempty(outDir)
        outDir = fullfile(pwd, 'glm_qc_for_transfer');
    end
    if ~isfolder(outDir), mkdir(outDir); end

    % --- settings ------------------------------------------------------------------
    subjects = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123', ...
                'sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
    projects    = {'dg','da'};
    hRF_setting = 'glmsingle';
    % NO eccentricity restriction: keep the whole V1 label and save each vertex's
    % eccentricity alongside (`patchEccen`).
    %
    % The published analysis restricts to 4-8 deg, but that restriction is a
    % STIMULUS-MATCHING constraint, not the stimulus extent. The stimulus is a much
    % larger annulus (roughly 1-12 deg). Cartesian gratings have uniform spatial
    % frequency across the aperture while polar gratings' SF scales inversely with
    % eccentricity, so the two stimulus types are SF-matched only near 6 deg; 4-8 is the
    % band where the mismatch stays acceptable.
    %
    % That constraint governs the ASYMMETRY analysis. It has no bearing on a GLM
    % FIT-QUALITY question, which only asks whether the session produced a usable
    % response anywhere the stimulus drove V1. Restricting quality assessment to 4-8
    % throws away most of the stimulated cortex and most of the statistical power.
    %
    % So: extract everything, cut locally. Any band -- 4-8, 2-8, 1-12 -- is then a
    % one-line filter on patchEccen at the receiving end, with no further server run.
    eccRange    = [-inf inf];
    % NO pRF-R2 restriction either, for the same reason plus a subtler one. vexpl is the
    % pRF fit from the RETINOTOPY scan, so thresholding it is not literally thresholding
    % the GLM R2 under test -- but it selects vertices with good SNR (clean timeseries,
    % no dropout, good registration), and those have higher GLM R2 whether or not they
    % responded to this stimulus. Any comparison of a vexpl-filtered patch against an
    % unfiltered baseline is therefore biased upward by the selection alone.
    %
    % Keeping every V1 vertex and saving vexpl per vertex (`patchVexpl`) makes the
    % threshold a local choice, and makes sensitivity to it testable rather than baked in.
    r2min       = -inf;
    qcFields    = {'R2','R2run','FRACvalue','noisepool','HRFindex','meanvol','xvaltrend','pcnum'};

    if ~isfolder(bidsDir)
        error('bidsDir not reachable:\n  %s\nEdit the default at the top of this file.', bidsDir);
    end
    try
        setup_user('rania', bidsDir);
    catch
        warning(['setup_user failed. If read_label / get_surfsize / MRIread are not on ' ...
                 'the path, the V1 restriction will be skipped (metrics are still saved).']);
    end

    fprintf('\nextract_for_transfer\n  bidsDir: %s\n  outDir : %s\n\n', bidsDir, outDir);
    rows = {};

    for si = 1:numel(subjects)
        subjectname = subjects{si};

        % ---- V1 patch, shared by both experiments --------------------------------
        [v1idx, patchNote, retInventory, retPatch, patchEccen, patchVexpl] = ...
            build_patch(bidsDir, subjectname, eccRange, r2min);
        fprintf('%-14s %s\n', subjectname, patchNote);

        for pj = 1:numel(projects)
            proj = projects{pj};
            tag  = sprintf('%s_%s', subjectname, proj);

            glmDir = fullfile(bidsDir,'derivatives',strcat(proj,'GLM'), ...
                              strcat('hRF_',hRF_setting), subjectname);
            f = dir(fullfile(glmDir,'**','results.mat'));
            if isempty(f)
                fprintf('    %-6s MISSING results.mat under %s\n', proj, glmDir);
                rows(end+1,:) = {subjectname, proj, 'missing', '', NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
                continue
            end

            S = load(fullfile(f(1).folder, f(1).name));
            if ~isfield(S,'results') || ~isfield(S.results,'allevents')
                fprintf('    %-6s UNEXPECTED structure -- different hRF_setting?\n', proj);
                rows(end+1,:) = {subjectname, proj, 'bad-structure', f(1).folder, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
                continue
            end
            a = S.results.allevents;

            qc = struct('subject',subjectname,'project',proj,'sourceFolder',f(1).folder, ...
                        'nVertices',numel(a.R2),'v1idx',uint32(v1idx),'patchNote',patchNote, ...
                        'retInventory',{retInventory},'retPatchMedians',retPatch, ...
                        'eccRange',eccRange,'vexplMin',r2min, ...
                        'patchEccen',patchEccen,'patchVexpl',patchVexpl);
            qc.fieldsMissing = setdiff(qcFields, qcFields(isfield(a,qcFields)));
            qc.wholeSurface  = summarise(a, qcFields);

            % A BASELINE THE PATCH CAN HONESTLY BE COMPARED AGAINST.
            % wholeSurface holds six pooled percentiles and includes the V1 vertices
            % themselves, so it cannot support a matched comparison. Keep instead a fine
            % percentile grid of R2 over NON-V1 cortex, whole-experiment and per run.
            % 401 doubles is nothing, and it makes both directions available by
            % interpolation: an arbitrary percentile of the baseline, or the fraction of
            % baseline vertices below any given patch value.
            qc.baselinePct = 0:0.25:100;
            R2all = double(a.R2(:));
            outside = true(numel(R2all),1);
            if ~isempty(v1idx), outside(v1idx) = false; end
            qc.baselineR2 = prctile(R2all(outside & isfinite(R2all)), qc.baselinePct);
            qc.baselineN  = sum(outside & isfinite(R2all));
            if isfield(a,'R2run')
                rr = double(reshape(a.R2run, numel(R2all), []));
                qc.baselineR2run = nan(size(rr,2), numel(qc.baselinePct));
                for rIdx = 1:size(rr,2)
                    col = rr(:,rIdx);
                    qc.baselineR2run(rIdx,:) = prctile(col(outside & isfinite(col)), qc.baselinePct);
                end
            end

            for k = 1:numel(qcFields)
                fn = qcFields{k};
                if ~isfield(a,fn), continue; end
                v = a.(fn);
                if ~islogical(v), v = single(v); end
                if ~isempty(v1idx) && size(v,1) == numel(a.R2)
                    v = v(v1idx,:,:,:);              % keep only the analysed patch
                end
                qc.(fn) = v;
            end

            save(fullfile(outDir, sprintf('glmqc_%s.mat', tag)), 'qc', '-v7');

            R2 = double(qc.R2(:));
            rr = NaN; rs = NaN;
            if isfield(qc,'R2run')
                perRun = median(double(squeeze(qc.R2run)), 1, 'omitnan');
                rr = min(perRun); rs = max(perRun) - min(perRun);
            end
            rows(end+1,:) = {subjectname, proj, 'ok', f(1).folder, numel(R2), ...
                             median(R2,'omitnan'), mean(R2>5), rr, rs}; %#ok<AGROW>
            fprintf('    %-6s ok  nPatch=%-6d medR2=%.2f%%  worstRun=%.2f\n', ...
                    proj, numel(R2), median(R2,'omitnan'), rr);
        end
    end

    M = cell2table(rows, 'VariableNames', {'subject','project','status','sourceFolder', ...
        'nPatchVertices','medianR2','fracR2gt5','worstRunR2','runSpread'});
    writetable(M, fullfile(outDir,'summary.csv'));

    fprintf('\nDone. %d of %d ok.\n', sum(strcmp(M.status,'ok')), height(M));
    fprintf('Zip and return: %s\n', outDir);
    fprintf('If it is too big to email, summary.csv alone answers most of the question.\n\n');
end

% =========================================================================
function [v1idx, note, inventory, retPatch, patchEccen, patchVexpl] = ...
        build_patch(bidsDir, subjectname, eccRange, r2min)
% V1 label intersected with the eccentricity and pRF-R2 filters. Deliberately uses only
% eccen and vexpl, never angle -- polar angle has a Benson-vs-conventional convention
% trap that is irrelevant here (see AUDIT.md), so this sidesteps it entirely.
%
% Returns each retained vertex's eccentricity and vexpl as well, so the receiving end can
% narrow the range (e.g. back to the published 4-8 deg) without re-reading the server.
    v1idx = []; inventory = {}; retPatch = struct();
    patchEccen = single([]); patchVexpl = single([]);
    try
        hSize = get_surfsize(subjectname);
        lh = read_label(subjectname, 'retinotopy_RE/lh.V1_REmanual');
        rh = read_label(subjectname, 'retinotopy_RE/rh.V1_REmanual');
        labelIdx = [lh(:,1)+1 ; rh(:,1)+hSize(1)+1];        % label files are 0-based

        d = dir(fullfile(bidsDir,'derivatives','prfvista_mov',subjectname,'**','stimfiles.mat'));
        retDir = d(1).folder;

        ecc  = [getfield(MRIread(fullfile(retDir,'lh.eccen.mgz')),'vol'), ...
                getfield(MRIread(fullfile(retDir,'rh.eccen.mgz')),'vol')];
        vexp = [getfield(MRIread(fullfile(retDir,'lh.vexpl.mgz')),'vol'), ...
                getfield(MRIread(fullfile(retDir,'rh.vexpl.mgz')),'vol')];

        nVert = numel(ecc);
        inLabel = false(nVert,1); inLabel(labelIdx) = true;
        % Apply each filter only if it is actually bounded. A blanket comparison would
        % silently drop NaN-eccen / NaN-vexpl vertices, which is itself a selection on
        % pRF fit quality -- exactly what this extraction is trying not to bake in.
        % Vertices whose pRF failed are kept, with NaN in patchEccen/patchVexpl, so any
        % local filter drops them explicitly rather than by accident.
        keep = inLabel;
        if all(isfinite(eccRange))
            keep = keep & ecc(:) >= eccRange(1) & ecc(:) <= eccRange(2);
        end
        if isfinite(r2min)
            keep = keep & vexp(:) > r2min;
        end
        v1idx = find(keep);
        patchEccen = single(ecc(v1idx));
        patchVexpl = single(vexp(v1idx));
        if all(isfinite(eccRange)) || isfinite(r2min)
            note = sprintf('V1 patch: %d vertices, ecc [%g %g], vexpl > %g', ...
                           numel(v1idx), eccRange(1), eccRange(2), r2min);
        else
            note = sprintf(['V1 patch: %d vertices, unfiltered ' ...
                            '(cut locally on patchEccen / patchVexpl)'], numel(v1idx));
        end

        % Inventory every retinotopy map, and its median in the patch. We are looking for
        % a response-amplitude / gain map to use as an effect-independent normaliser, but
        % do not know from off-server which maps this pipeline writes.
        lhFiles = dir(fullfile(retDir,'lh.*.mgz'));
        inventory = {lhFiles.name};
        for k = 1:numel(lhFiles)
            nm = erase(erase(lhFiles(k).name,'lh.'),'.mgz');
            try
                v = [getfield(MRIread(fullfile(retDir,['lh.' nm '.mgz'])),'vol'), ...
                     getfield(MRIread(fullfile(retDir,['rh.' nm '.mgz'])),'vol')];
                retPatch.(matlab.lang.makeValidName(nm)) = median(double(v(v1idx)),'omitnan');
            catch
            end
        end
    catch ME
        note = sprintf('V1 patch UNAVAILABLE (%s) -- whole-surface metrics still saved', ME.message);
    end
end

% -------------------------------------------------------------------------
function s = summarise(a, qcFields)
    s = struct();
    for k = 1:numel(qcFields)
        f = qcFields{k};
        if ~isfield(a,f), continue; end
        v = double(a.(f));
        if numel(v) < 100, s.(f) = v; continue; end
        v = v(isfinite(v));
        s.(f) = struct('p5',prctile(v,5),'p25',prctile(v,25),'median',median(v), ...
                       'p75',prctile(v,75),'p95',prctile(v,95),'mean',mean(v));
    end
end
