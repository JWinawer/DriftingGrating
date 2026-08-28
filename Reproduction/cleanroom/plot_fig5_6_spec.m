function out = plot_fig5_6_spec(S, en, figNum, figDir, L, V)
% PLOT_FIG5_6_SPEC  Figure 5 (dg) / Figure 6 (da) under the settled specification.
%
%   out = plot_fig5_6_spec(S, 'dg', 5, figDir)          % S from SPEC_PROFILES
%   out = plot_fig5_6_spec(S, 'dg', 5, figDir, L)       % L from SPEC_AXIS_LIMITS
%   out = plot_fig5_6_spec(S, 'dg', 5, figDir, L, V)    % V from SPEC_VARIANTS
%
% SCALES ARE SHARED, and shared across BOTH experiments, not just within one figure.
% Every polar panel gets the same radius and every difference panel the same y limits,
% here and in the companion figure, so Figure 5 and Figure 6 can be compared by eye.
% dg-versus-da is the paper's claim; drawing the two on different scales would undercut
% the one comparison the reader most needs to make. Pass L to fix the limits explicitly;
% omitted, it is computed from S, which already holds both experiments.
%
% Row 1, polar plots: the pro and con class responses. LINES are the fitted four-term
%   harmonic model (continuous thetaV, equal coverage at 45 deg, per observer x map
%   gain); MARKERS are the observed wedge means over the same vertices, at the eight
%   polar-angle ROI centres. Where they separate, the difference is the within-wedge
%   local-orientation term that binning thetaV conflates with context -- visible here
%   rather than argued (../supplement/SUPPLEMENT_harmonic_model.md section S5.1).
%
%   THE MODEL IS DRAWN AT 0.5 DEG, the data at 45. The fit is continuous in thetaV,
%   so sampling it at the eight marker positions and joining the dots would draw an
%   octagon that is an artefact of the display grid, not of the model -- and the
%   second and fourth harmonics it is made of have their extrema BETWEEN the wedge
%   centres, so the octagon understates them. SPEC_PROFILES supplies the curve as
%   .dense, and asserts that it passes exactly through the eight centres, so the
%   smooth line and the markers cannot come apart. The data are still shown only
%   where they were measured. Display sampling only: nothing fitted or tested
%   changes, and the ROI route (which fits no model) still joins its wedge means.
%
%   RADIUS IS DEMEANED. The model is fitted to each vertex's four orientation
%   responses with that vertex's mean over the four removed, so it carries no overall
%   response level and the radial axis is a deviation from the local mean, not a raw
%   percent signal change. That is deliberate: demeaning removes the blank, which is
%   full-field pink noise rather than a baseline (../local_qc/DATA_QUALITY.md section 1),
%   so only orientation DIFFERENCES enter. Zero is drawn as a grey circle.
%
% Row 2: each observer's pro-minus-con as one point, with the group mean and a t
%   interval on n-1 df. NOT the paired pro-vs-con lines of the original figure: under
%   this parameterisation the two classes are demeaned and partition the four
%   orientations, so con = -pro exactly for the two second-harmonic asymmetries and
%   the paired plot would show the same information mirrored. The difference is the
%   estimand, so the difference is what is plotted.
%
%   t, NOT the percentile bootstrap: at n = 8 the percentile method has poor coverage
%   and disagrees with t on exactly the two polar cells this figure would be read for
%   (../METHOD_DECISIONS.md section 5). SPEC_TABLES reports both, so the choice is visible.

    if nargin < 5 || isempty(L), L = spec_axis_limits(S); end
    if nargin < 6 || isempty(V)
        V = spec_variants('spec');
        if ~S.hasModel, V = spec_variants('roi'); end
    end
    E   = S.(en);
    nS  = numel(E.subjects);      % this experiment's observers, not a shared list
    nP  = numel(S.paBins);
    titles = {'horizontal vs vertical','cardinal vs oblique', ...
              'radial vs tangential','polar-card vs polar-obl'};
    proName = {'horizontal','cardinal','radial','polar-card'};
    conName = {'vertical','oblique','tangential','polar-obl'};
    proC = {[0.16 0.20 0.52],[0.13 0.47 0.24],[0.20 0.55 0.75],[0.55 0.33 0.12]};
    conC = {[0.93 0.60 0.11],[0.55 0.40 0.72],[0.78 0.16 0.16],[0.18 0.58 0.58]};

    th = deg2rad([S.paBins, S.paBins(1)]);
    cl = @(v) [v(:); v(1)];

    fig = figure('Color','w','Position',[60 60 1560 780],'Visible','off');
    tl  = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
    expLbl = 'Cartesian gratings (dg)';
    if strcmp(en,'da'), expLbl = 'Polar gratings (da)'; end
    hdr = 'settled specification';
    if ~strcmp(V.tag,'spec'), hdr = 'ALTERNATIVE ROUTE'; end
    title(tl, {sprintf('Figure %d, %s  \\bf%s\\rm  |  %s, %g-%g%s  |  n = %d', ...
               figNum, hdr, expLbl, S.area, S.eccRange(1), S.eccRange(2), char(176), nS), ...
               sprintf(['\\rm\\fontsize{9}\\bf%s\\rm  |  equal coverage at 45%s, gain per ' ...
               'observer \\times map, fit per observer then combined'], V.label, char(176)), ...
               [tern(S.hasModel, ...
                 '\rm\fontsize{9}lines = fitted model (0.5\circ), markers = observed wedge means (45\circ); ', ...
                 '\rm\fontsize{9}lines join the observed wedge means (this route fits no model); ') ...
               'radius is the deviation from each vertex''s mean over the four orientations']}, ...
               'Interpreter','tex','FontWeight','bold');

    mPro = squeeze(mean(E.mPro,1,'omitnan'));   mCon = squeeze(mean(E.mCon,1,'omitnan'));
    oPro = squeeze(mean(E.oPro,1,'omitnan'));   oCon = squeeze(mean(E.oCon,1,'omitnan'));
    rmax = L.rmax;

    % The model curve, on its own fine grid. Averaged across observers exactly as
    % mPro/mCon are -- the mean of eight fitted curves, not a curve through the mean.
    hasDense = S.hasModel && isfield(E,'dense');
    if hasDense
        thD  = deg2rad(E.dense.centres(:));      % 0:0.5:360 already closes the loop
        dPro = squeeze(mean(E.dense.mPro,1,'omitnan'));
        dCon = squeeze(mean(E.dense.mCon,1,'omitnan'));
    end

    % The zero reference is a CIRCLE, so draw it on a dense grid. Drawn at the eight
    % wedge angles it renders as a grey OCTAGON sitting among the model curves, which
    % in a figure whose whole point is a smooth fit against 45-deg samples reads as if
    % the model itself had been sampled every 45 deg and joined by straight lines.
    thZero = deg2rad(0:0.5:360);

    for i = 1:4
        nexttile(i);
        polarplot(thZero, zeros(size(thZero)), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.75); hold on;
        if hasDense
            % Lines are the fitted model at 0.5 deg, markers the observed wedge means
            % at 45. Where they separate, the gap is the within-wedge
            % local-orientation term.
            hP = polarplot(thD, dPro(:,i), '-', 'Color', proC{i}, 'LineWidth', 2.2);
            hC = polarplot(thD, dCon(:,i), '-', 'Color', conC{i}, 'LineWidth', 2.2);
        elseif S.hasModel
            % No dense grid in this result (an older SPEC_PROFILES): fall back to the
            % eight centres, which is what this figure drew before.
            hP = polarplot(th, cl(mPro(:,i)), '-', 'Color', proC{i}, 'LineWidth', 2.2);
            hC = polarplot(th, cl(mCon(:,i)), '-', 'Color', conC{i}, 'LineWidth', 2.2);
        else
            % The ROI route fits no per-vertex model, so there is nothing to overlay:
            % the lines simply join the wedge means, which is all this route computes.
            hP = polarplot(th, cl(oPro(:,i)), '-', 'Color', proC{i}, 'LineWidth', 2.2);
            hC = polarplot(th, cl(oCon(:,i)), '-', 'Color', conC{i}, 'LineWidth', 2.2);
        end
        polarplot(th(1:nP), oPro(:,i), 'o', 'Color', proC{i}, 'MarkerFaceColor', proC{i}, 'MarkerSize', 5);
        polarplot(th(1:nP), oCon(:,i), 'o', 'Color', conC{i}, 'MarkerFaceColor', conC{i}, 'MarkerSize', 5);
        ax = gca; ax.RLim = [-rmax rmax]; ax.ThetaZeroLocation = 'right';
        ax.ThetaDir = 'counterclockwise'; ax.RTickLabel = []; ax.FontSize = 8;
        ax.ThetaTick = 0:45:315;
        title(titles{i}, 'FontWeight','bold', 'FontSize', 10);
        lg = legend([hP hC], {proName{i}, conName{i}}, 'FontSize', 8.5, ...
                    'Location','southoutside', 'Orientation','horizontal', 'Box','off');
        lg.ItemTokenSize = [16 8];
    end

    for i = 1:4
        d  = E.asym(:,i);
        ok = isfinite(d);
        G  = spec_group(d, E.sigma(:,i), V.weighting);
        n  = G.n;  g = G.mean;  ci = [G.lo G.hi];  pv = G.p;

        nexttile(4+i); hold on;
        yline(0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.9);
        rng(0);                                   % reproducible jitter
        jx = 1 + (rand(nS,1)-0.5)*0.30;
        plot(jx(ok), d(ok), 'o', 'MarkerSize', 6, 'MarkerFaceColor', [1 1 1], ...
             'MarkerEdgeColor', [0.45 0.45 0.45], 'LineWidth', 0.9);
        cc = proC{i};
        plot([1.55 1.55], ci, '-', 'Color', cc, 'LineWidth', 2.4);
        plot([1.48 1.62], [ci(1) ci(1)], '-', 'Color', cc, 'LineWidth', 1.6);
        plot([1.48 1.62], [ci(2) ci(2)], '-', 'Color', cc, 'LineWidth', 1.6);
        plot(1.55, g, 'o', 'MarkerSize', 10, 'MarkerFaceColor', cc, 'MarkerEdgeColor','k', 'LineWidth', 1);
        xlim([0.65 1.85]);
        ylim(L.dot);
        set(gca,'XTick',[1 1.55],'XTickLabel',{'observers','group'},'FontSize',9);
        grid on; box off; ax = gca; ax.YGrid = 'on'; ax.XGrid = 'off';
        star = ''; if ci(1) > 0 || ci(2) < 0, star = ' *'; end
        nAgree = G.obsAgree;
        title(sprintf('%s %s %s\\rm\\fontsize{8}\\newline\\Delta = %.3f [%.3f, %.3f]%s,  \\itp\\rm %s,  %d/%d obs', ...
              proName{i}, char(8722), conName{i}, g, ci(1), ci(2), star, pstr(pv), nAgree, n), ...
              'Interpreter','tex', 'FontSize', 9.5, 'FontWeight','bold');
        if i == 1
            ylabel({'\Delta demeaned response','(% signal change)'}, 'FontSize', 9, 'Interpreter','tex');
        end
    end

    if ~isfolder(figDir), mkdir(figDir); end
    base = fullfile(figDir, sprintf('Figure_%d_%s_%s%s', figNum, V.tag, en, areaSfx(S)));
    exportgraphics(fig, [base '.png'], 'Resolution', 200);
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector');
    close(fig);
    out = [base '.png'];
    fprintf('plot_fig5_6_spec: wrote %s\n', out);
end

function s = areaSfx(S)
% The supplement figures are V1; anything else gets the map in its name so the two
% cannot be confused for one another.
    if strcmpi(S.area,'V1') && isequal(S.eccRange(:).', [4 8]), s = '';
    else, s = sprintf('_%s_%g-%g', lower(S.area), S.eccRange(1), S.eccRange(2)); end
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

function s = pstr(p)
    if p < 0.001, s = '< .001'; else, s = sprintf('= %.3f', p); end
end
