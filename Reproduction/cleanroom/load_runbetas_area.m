function [A, ok] = load_runbetas_area(S, cfg, en, root, area)
% LOAD_RUNBETAS_AREA  Restrict a run-wise beta extraction to the analysed vertices.
%
%   [A, ok] = load_runbetas_area(S, cfg, en, root, area)
%
% S is a loaded ~/dg_collect/runbetas[_areas]_<subject>_<exp>.mat. Applies the
% analysis inclusion filter (cfg.eccRange, pRF R2 > cfg.r2min) within `area` and
% attaches both polar-angle representations the settled specification needs:
% A.thetaV (continuous, conventional degrees) and A.wedge (nearest of cfg.paBins).
%
% ok is false when no vertex survives, which happens in the sparser maps and is a
% reportable fact rather than an error -- see ../SPECIFICATION.md section 7.
%
% EXTRACTED from DIAGNOSE_WITHIN_OBSERVER_ERROR so SPEC_PROFILES fits the same
% vertices from the same files. One definition of the analysed vertex set, not two.
%
% Accepts either extraction layout:
%   runbetas_areas_*  vertIndex + areaMask + areaNames, eight visual areas
%   runbetas_*        v1Index, V1 only (the original COLLECT_RUNWISE_BETAS)
% so the V1 results are reproducible from whichever files are present.
%
% Output A
%   .runBeta  nVertex x nCond x nRun
%   .thetaV   nVertex x 1  conventional polar angle, continuous
%   .wedge    nVertex x 1  index into cfg.paBins
%   .expn     the experiment name
%
% ON THE ANGLE. ret_*.mat stores Benson angle_adj, which is what meanWithinLabel.m
% bins; conv = mod(90 - ang, 360) puts it in the conventional frame cfg.paBins is
% written in. Getting this backwards reflects the wedges about 45 deg and swaps the
% four cardinal meridians -- see ../STIMULUS_CONVENTIONS.md sections 3 and 5.

    ok = false;  A = struct();
    R = load(fullfile(root, sprintf('ret_%s.mat', S.subject)), 'eccen','vexpl','angle_adj');

    if isfield(S, 'vertIndex')
        ai = find(strcmp(S.areaNames, area), 1);
        if isempty(ai)
            error('load_runbetas_area:area', ...
                  '%s holds %s, not %s.', S.subject, strjoin(S.areaNames, '/'), area);
        end
        sel = S.areaMask(:, ai);
        v   = S.vertIndex(sel);
        rb  = S.runBeta(sel, :, :);
    else
        if ~strcmp(area, 'V1')
            error('load_runbetas_area:v1only', ...
                  ['%s has only the V1-only extraction, so area ''%s'' is unavailable. ' ...
                   'Run collect_runwise_betas_areas.m.'], S.subject, area);
        end
        v  = S.v1Index;
        rb = S.runBeta;
    end

    good = double(R.eccen(v)) >= cfg.eccRange(1) & double(R.eccen(v)) <= cfg.eccRange(2) ...
         & double(R.vexpl(v)) > cfg.r2min;
    if ~any(good), return; end
    A.runBeta = rb(good, :, :);
    ang  = double(R.angle_adj(v(good)));            % Benson deg, as meanWithinLabel bins
    conv = mod(90 - ang, 360);                      % conventional, matches cfg.paBins
    [~, A.wedge] = min(abs(mod(conv - cfg.paBins(:).' + 180, 360) - 180), [], 2);
    A.thetaV = conv(:);                             % continuous, for the harmonic route
    A.expn = en;
    ok = true;
end
