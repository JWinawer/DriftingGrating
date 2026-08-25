function extract_glm_qc(bidsDir, outDir, saveFullSurface)
% EXTRACT_GLM_QC  Pull GLMsingle fit-quality metrics out of each subject's results.mat.
%
%   extract_glm_qc()                       % uses the default bidsDir below
%   extract_glm_qc(bidsDir, outDir)
%   extract_glm_qc(bidsDir, outDir, true)  % also keep whole-surface vectors
%
% RUN THIS ON THE MACHINE WHERE /Volumes/Vision IS MOUNTED. It writes one small .mat
% per subject x experiment into outDir (default Reproduction/_glmqc/); those small files
% are what you copy back, NOT the results.mat files themselves.
%
% Why: each results.mat is ~500 MB, almost all of it `modelmd` (nVertices x 416 single),
% and the 16 files come to ~8 GB. The quality fields are a rounded ~1% of that, and
% restricting to the analysed V1 patch here means the FreeSurfer labels and pRF maps do
% not need copying either. The files are v7 (not v7.3), so partial loading is impossible
% -- each one is read into memory whole, which is fine locally and wasteful over a wire.
%
% Session IDs are discovered, not assumed: main_singlesub.m sets `ses` by hand per
% subject, so this globs for results.mat the same way meanWithinLabel.m:100 does.
%
% STATUS: the results.mat half is tested against Support/sub-0255. The V1-restriction
% half follows meanWithinLabel.m:96-183 but is UNTESTED -- no labels or pRF maps are on
% the machine this was written on. If it errors, the whole-surface metrics are still
% written and the V1 fields are left empty; check the warning it prints.

    if nargin < 1 || isempty(bidsDir)
        bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
    end
    cfg = config_repro();
    if nargin < 2 || isempty(outDir), outDir = fullfile(cfg.reproDir, '_glmqc'); end
    if nargin < 3 || isempty(saveFullSurface), saveFullSurface = false; end
    if ~isfolder(outDir), mkdir(outDir); end

    if ~isfolder(bidsDir)
        error('extract_glm_qc:noBidsDir', ...
              'bidsDir not reachable: %s\nMount the volume, or pass the path explicitly.', bidsDir);
    end

    % Same path setup meanWithinLabel.m relies on (read_label / get_surfsize / MRIread).
    try, setup_user('rania', bidsDir); catch, warning('setup_user failed; V1 restriction may not work.'); end

    hRF_setting = 'glmsingle';
    qcFields = {'R2','R2run','FRACvalue','noisepool','HRFindex','meanvol','xvaltrend','pcnum'};
    projects = {'dg','da'};

    fprintf('\nextract_glm_qc: bidsDir = %s\n            outDir  = %s\n\n', bidsDir, outDir);
    manifest = {};

    for pi = 1:numel(projects)
        proj = projects{pi};
        for si = 1:numel(cfg.subjects)
            subjectname = cfg.subjects{si};
            tag = sprintf('%s_%s', subjectname, proj);

            glmResultsfolder = fullfile(bidsDir, 'derivatives', ...
                strcat(proj,'GLM'), strcat('hRF_', hRF_setting), subjectname);
            f = dir(fullfile(glmResultsfolder, '**', 'results.mat'));

            if isempty(f)
                fprintf('  %-24s MISSING results.mat under %s\n', tag, glmResultsfolder);
                manifest(end+1,:) = {tag, 'missing', '', NaN}; %#ok<AGROW>
                continue
            end
            if numel(f) > 1
                fprintf('  %-24s NOTE: %d results.mat found, using %s\n', tag, numel(f), f(1).folder);
            end

            S = load(fullfile(f(1).folder, f(1).name));
            if ~isfield(S,'results') || ~isfield(S.results,'allevents')
                fprintf('  %-24s UNEXPECTED structure (no results.allevents) -- different hRF_setting?\n', tag);
                manifest(end+1,:) = {tag, 'bad-structure', f(1).folder, NaN}; %#ok<AGROW>
                continue
            end
            a = S.results.allevents;

            have = qcFields(isfield(a, qcFields));
            miss = setdiff(qcFields, have);
            if ~isempty(miss)
                fprintf('  %-24s WARNING missing fields: %s\n', tag, strjoin(miss, ', '));
            end

            qc = struct('subject', subjectname, 'project', proj, ...
                        'sourceFolder', f(1).folder, 'fieldsPresent', {have}, ...
                        'fieldsMissing', {miss});
            for k = 1:numel(have)
                v = a.(have{k});
                if islogical(v), qc.(have{k}) = v; else, qc.(have{k}) = single(v); end
            end
            nVert = numel(a.R2);
            qc.nVertices = nVert;

            % --- restrict to the analysed V1 patch (mirrors meanWithinLabel.m) -------
            qc.v1idx = []; qc.v1note = '';
            try
                hSize = get_surfsize(subjectname);
                lh = read_label(subjectname, 'retinotopy_RE/lh.V1_REmanual');
                rh = read_label(subjectname, 'retinotopy_RE/rh.V1_REmanual');
                label_idx = [lh(:,1)+1 ; rh(:,1)+hSize(1)+1];   % 0-based -> 1-based

                retDir = dir(fullfile(bidsDir,'derivatives','prfvista_mov',subjectname,'**','stimfiles.mat'));
                retDir = retDir(1).folder;
                ecc  = [MRIread(fullfile(retDir,'lh.eccen.mgz')).vol, ...
                        MRIread(fullfile(retDir,'rh.eccen.mgz')).vol];
                vexp = [MRIread(fullfile(retDir,'lh.vexpl.mgz')).vol, ...
                        MRIread(fullfile(retDir,'rh.vexpl.mgz')).vol];
                % NB: only ecc and vexpl are used, never angle -- this deliberately avoids
                % the Benson-vs-conventional polar-angle convention
                % (see ../STIMULUS_CONVENTIONS.md).

                inPatch = false(nVert,1);
                inPatch(label_idx) = true;
                keep = inPatch(:) & ecc(:) >= cfg.eccRange(1) & ecc(:) <= cfg.eccRange(2) ...
                                  & vexp(:) > cfg.r2min;
                qc.v1idx = uint32(find(keep));
                qc.v1note = sprintf('V1 label + ecc in [%g %g] + vexpl > %g', ...
                                    cfg.eccRange(1), cfg.eccRange(2), cfg.r2min);
                fprintf('  %-24s ok, %d vertices, V1 patch %d\n', tag, nVert, numel(qc.v1idx));
            catch ME
                qc.v1note = sprintf('V1 restriction FAILED: %s', ME.message);
                warning('%s: %s', tag, qc.v1note);
                fprintf('  %-24s ok, %d vertices, V1 patch UNAVAILABLE\n', tag, nVert);
            end

            % Trim to keep the transfer small: V1-patch vectors in full, plus
            % whole-surface percentile summaries. modelmd is never carried.
            if ~saveFullSurface && ~isempty(qc.v1idx)
                qc.wholeSurfaceSummary = summarise(a, qcFields);
                for k = 1:numel(have)
                    if ismember(have{k}, {'xvaltrend','pcnum'}), continue; end
                    v = qc.(have{k});
                    if size(v,1) == nVert, qc.(have{k}) = v(qc.v1idx, :, :, :); end
                end
                qc.scope = 'V1 patch only (+ wholeSurfaceSummary)';
            else
                qc.scope = 'whole surface';
            end

            save(fullfile(outDir, sprintf('glmqc_%s.mat', tag)), 'qc', '-v7');
            manifest(end+1,:) = {tag, 'ok', f(1).folder, numel(qc.v1idx)}; %#ok<AGROW>
        end
    end

    M = cell2table(manifest, 'VariableNames', {'tag','status','sourceFolder','nV1'});
    writetable(M, fullfile(outDir, 'manifest.csv'));
    fprintf('\n%d/%d extracted. Manifest: %s\n', sum(strcmp(M.status,'ok')), height(M), ...
            fullfile(outDir, 'manifest.csv'));
    fprintf('Copy back the whole %s folder (small), then run audit_glm_quality.\n\n', outDir);
end

% -------------------------------------------------------------------------
function s = summarise(a, qcFields)
% Percentile summary of the whole surface, so trimming to V1 loses no context.
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
