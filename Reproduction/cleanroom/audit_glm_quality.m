function A = audit_glm_quality(qcDir)
% AUDIT_GLM_QUALITY  Per-subject GLMsingle fit quality, from EXTRACT_GLM_QC output.
%
%   A = audit_glm_quality()          % reads Reproduction/_glmqc/
%   A = audit_glm_quality(qcDir)
%
% Consumes the small per-subject files written by EXTRACT_GLM_QC and answers the
% question in ../NEXT_STEPS.md: are the observers that per-vertex z-scoring up-weights
% most also the ones with the worst GLM fits?
%
% Two specific predictions are worth checking, from ZSCORE_FIG7.md section 3a:
%   sub-0201  its blank beta exceeds all 12 stimulus betas in BOTH experiments, which is
%             what a bad run or a motion artefact produces -> expect low R2, and expect
%             R2run to single out one or more bad runs.
%   sub-0037  no condition differentiation at all in `da` but a strong response in `dg`
%             -> expect a da-specific quality drop, dg comparatively normal.
%
% Reports metrics within the analysed V1 patch where available, and flags any subject
% more than 1.5 IQR outside the group on each metric. Nothing here decides exclusion --
% it produces the table that decision should be made from.

    cfg = config_repro();
    if nargin < 1 || isempty(qcDir), qcDir = fullfile(cfg.reproDir, '_glmqc'); end
    files = dir(fullfile(qcDir, 'glmqc_*.mat'));
    if isempty(files)
        error('audit_glm_quality:noData', ...
            ['No glmqc_*.mat in %s.\nRun EXTRACT_GLM_QC on the machine where the ' ...
             'data volume is mounted, then copy that folder here.'], qcDir);
    end

    rows = struct('subject',{},'project',{},'nVertex',{},'medR2',{},'fracR2gt5',{}, ...
                  'medFRAC',{},'fracNoisepool',{},'medMeanvol',{},'worstRunR2',{}, ...
                  'runSpread',{},'scope',{});
    for i = 1:numel(files)
        S = load(fullfile(files(i).folder, files(i).name), 'qc'); q = S.qc;
        R2 = double(q.R2(:)); R2 = R2(isfinite(R2));
        r.subject = q.subject; r.project = q.project; r.nVertex = numel(R2);
        r.medR2 = median(R2); r.fracR2gt5 = mean(R2 > 5);
        r.medFRAC = pick(q,'FRACvalue', @(v) median(double(v(isfinite(v)))));
        r.fracNoisepool = pick(q,'noisepool', @(v) mean(logical(v)));
        r.medMeanvol = pick(q,'meanvol', @(v) median(double(v(isfinite(v)))));
        if isfield(q,'R2run')
            rr = double(squeeze(q.R2run)); if size(rr,1) ~= numel(R2) && size(rr,2) == numel(R2), rr = rr.'; end
            perRun = median(rr, 1, 'omitnan');
            r.worstRunR2 = min(perRun); r.runSpread = max(perRun) - min(perRun);
        else
            r.worstRunR2 = NaN; r.runSpread = NaN;
        end
        % EXTRACT_GLM_QC calls this field `scope`; the standalone server script
        % ../server_extract/extract_for_transfer.m calls it `patchNote`. Accept either,
        % so output from whichever extractor was actually run can be audited.
        if isfield(q,'scope'),   r.scope = q.scope;
        elseif isfield(q,'patchNote'), r.scope = q.patchNote;
        else,                    r.scope = 'unknown';
        end
        rows(end+1) = r; %#ok<AGROW>
    end
    A = struct2table(rows);

    for proj = {'dg','da'}
        P = A(strcmp(A.project, proj{1}), :);
        if isempty(P), continue; end
        P = sortrows(P, 'medR2');
        fprintf('\n================ %s : GLM fit quality ================\n', proj{1});
        % scope is per-subject when it came from extract_for_transfer (it embeds that
        % subject's patch size), so print it only if every row agrees. Per-row vertex
        % counts are in the nVert column either way.
        if numel(unique(P.scope)) == 1
            fprintf('scope: %s\n', P.scope{1});
        else
            fprintf('scope: V1 patch, per-subject size (see nVert)\n');
        end
        fprintf('%-14s %8s %8s %9s %9s %9s %9s %9s\n', 'subject','nVert', ...
                'medR2%','R2>5%','medFRAC','noisepool','worstRun','runSpread');
        for i = 1:height(P)
            fprintf('%-14s %8d %8.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n', P.subject{i}, ...
                P.nVertex(i), P.medR2(i), P.fracR2gt5(i), P.medFRAC(i), ...
                P.fracNoisepool(i), P.worstRunR2(i), P.runSpread(i));
        end
        flagOutliers(P, {'medR2','fracR2gt5','medFRAC','fracNoisepool','worstRunR2'});
    end

    % --- dg-vs-da within subject: isolates session-specific failures -----------
    if all(ismember({'dg','da'}, A.project))
        fprintf('\n================ within-subject dg vs da ================\n');
        fprintf('%-14s %9s %9s %9s\n','subject','dg medR2','da medR2','da - dg');
        subs = unique(A.subject);
        for i = 1:numel(subs)
            d = A.medR2(strcmp(A.subject,subs{i}) & strcmp(A.project,'dg'));
            p = A.medR2(strcmp(A.subject,subs{i}) & strcmp(A.project,'da'));
            if isscalar(d) && isscalar(p)
                fprintf('%-14s %9.2f %9.2f %9.2f\n', subs{i}, d, p, p-d);
            end
        end
        fprintf(['(A large negative "da - dg" for sub-0037 would support a session-specific\n' ...
                 ' failure in the polar experiment; see ZSCORE_FIG7.md section 3a.)\n']);
    end

    fprintf('\nReminder: this table informs an inclusion decision, it does not make one.\n');
    fprintf('State whatever is decided in the Methods -- see ZSCORE_FIG7.md recommendation 7.\n\n');
end

% -------------------------------------------------------------------------
function v = pick(q, f, fn)
    if isfield(q, f) && ~isempty(q.(f)), v = fn(q.(f)); else, v = NaN; end
end

function flagOutliers(P, metrics)
    any_ = false;
    for m = 1:numel(metrics)
        x = P.(metrics{m});
        if all(isnan(x)), continue; end
        lo = prctile(x,25) - 1.5*iqr(x); hi = prctile(x,75) + 1.5*iqr(x);
        bad = find(x < lo | x > hi);
        for b = bad(:).'
            fprintf('  FLAG %-14s %s = %.2f (group median %.2f)\n', ...
                    P.subject{b}, metrics{m}, x(b), median(x,'omitnan'));
            any_ = true;
        end
    end
    if ~any_, fprintf('  (no subject outside 1.5 IQR on any metric)\n'); end
end
