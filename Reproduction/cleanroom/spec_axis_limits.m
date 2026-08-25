function L = spec_axis_limits(S)
% SPEC_AXIS_LIMITS  One set of axis limits for every panel that invites comparison.
%
%   L = spec_axis_limits(S)          % S from SPEC_PROFILES
%   L = spec_axis_limits({S1,S2,S3}) % span several variants, so they are comparable
%
% Panels drawn in the same units, showing the same kind of comparison, must share a
% scale. Otherwise a small effect drawn on a small axis looks like a large one, and
% the reader has to check the tick labels to avoid being misled -- which is exactly
% what a figure is supposed to save them from.
%
% The limits span BOTH experiments, so Figure 5 and Figure 6 are comparable to each
% other and not merely internally consistent. dg-versus-da IS the paper's claim, so
% putting the two figures on different scales would undercut the one comparison the
% reader most needs to make by eye.
%
% Returns
%   L.rmax   scalar, polar-plot radius; every polar panel gets [-rmax rmax]
%   L.dot    1x2, the y limits for every per-observer difference panel
%   L.prof   1x2, the y limits for every panel of the polar-angle profile figure
%
% Panels that are NOT the same kind of comparison keep their own scale: in the
% hierarchy figure the asymmetries (top row) and the context effects (bottom row) are
% different quantities, so each row is linked internally and the two rows are not.
% That grouping lives in PLOT_SPEC_HIERARCHY, which has the values; this function
% covers everything derivable from one SPEC_PROFILES result.

    if ~iscell(S), S = {S}; end
    Ss = S;
    expn = {'dg','da'};

    % --- polar panels: model traces and observed markers, both experiments -----
    v = [];
    for q = 1:numel(Ss)
        for ei = 1:2
            E = Ss{q}.(expn{ei});
            for f = {'mPro','mCon','oPro','oCon'}
                m = squeeze(mean(E.(f{1}), 1, 'omitnan'));
                v = [v; m(:)]; %#ok<AGROW>
            end
        end
    end
    L.rmax = pad1(max(abs(v), [], 'omitnan'), 1.15);

    % --- per-observer difference panels: observers, and the group interval -----
    lo = []; hi = [];
    for q = 1:numel(Ss)
        for ei = 1:2
            E = Ss{q}.(expn{ei});
            for j = 1:size(E.asym,2)
                x = E.asym(isfinite(E.asym(:,j)), j);
                if numel(x) < 3, continue; end
                for w = {'equal','precision'}
                    G = spec_group(E.asym(:,j), E.sigma(:,j), w{1});
                    if ~isfinite(G.mean), continue; end
                    lo(end+1) = min([x; G.lo; 0]); %#ok<AGROW>
                    hi(end+1) = max([x; G.hi; 0]); %#ok<AGROW>
                end
            end
        end
    end
    L.dot = pad2([min(lo) max(hi)], 0.08);

    % --- profile panels: observed bin means with their SEM band, plus the model -
    lo = []; hi = [];
    for q = 1:numel(Ss)
    for ei = 1:2
        E = Ss{q}.(expn{ei});
        for k = 1:3
            o = E.fine.obs(:,:,k);
            n = sum(isfinite(o), 1);
            mu = mean(o, 1, 'omitnan');  se = std(o, 0, 1, 'omitnan') ./ max(sqrt(n), 1);
            mu(n < 3) = NaN;  se(n < 3) = NaN;
            mm = mean(E.fine.mdlDense(:,:,k), 1, 'omitnan');
            wv = mean(E.fine.wedge(:,:,k), 1, 'omitnan');
            lo(end+1) = min([mu - se, mm, wv], [], 'omitnan'); %#ok<AGROW>
            hi(end+1) = max([mu + se, mm, wv], [], 'omitnan'); %#ok<AGROW>
        end
    end
    end
    L.prof = pad2([min(lo) max(hi)], 0.05);
end

function r = pad1(r, f)
    r = r * f;
    if ~isfinite(r) || r <= 0, r = 1; end
end

function y = pad2(y, f)
    if any(~isfinite(y)), y = [-1 1]; return; end
    p = max(f * diff(y), eps);
    y = y + [-p p];
end
