function out = plot_spec_hierarchy(A, figDir)
% PLOT_SPEC_HIERARCHY  The extrastriate supplement figure.
%
%   out = plot_spec_hierarchy(A, figDir)       % A from SPEC_AREAS_SUMMARY
%
% Top row: each asymmetry in each experiment, across the maps that pass the coverage
% criterion. Bottom row: the context effect (dg - da) across the same maps.
% Filled markers are reportable maps at 4-8 deg; the open marker is V3a, which
% qualifies only at 2-10 deg and so is NOT on the same footing as the others.
%
% Maps that fail the criterion are absent by design, not by oversight: the criterion
% is a statement about what this design can resolve BY POLAR ANGLE (../SPECIFICATION.md
% section 7), and their numbers are still in spec_areas_*.csv with a reportable flag.
%
% SCALES ARE SHARED WITHIN EACH ROW, and deliberately not between them. All four top
% panels are asymmetries and all four bottom panels are context effects, so within a
% row the panels are the same units and the same kind of comparison and must be drawn
% against one axis. Across rows they are different quantities -- an asymmetry and a
% difference of asymmetries -- so a common scale there would assert a comparison that
% is not being made.

    nm = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    ttl = {'horizontal - vertical','cardinal - oblique','radial - tangential', ...
           'polar-card - polar-obl'};
    order = {'V1','V2','V3','V3a'};
    band  = {'4-8','4-8','4-8','2-10'};
    cDg = [0.16 0.20 0.52]; cDa = [0.85 0.42 0.10]; cCx = [0.30 0.30 0.30];

    % --- one scale per row, spanning every panel that will be drawn in it -------
    yTop = rowLimits(A.asym,    order, band, {'dg','da'}, nm);
    yBot = rowLimits(A.context, order, band, {'dg-da'},   nm);

    fig = figure('Color','w','Position',[60 60 1400 720],'Visible','off');
    tl = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
    vl = 'settled specification';
    if isfield(A,'variant') && ~strcmp(A.variant.tag,'spec')
        vl = ['ALTERNATIVE ROUTE -- ' A.variant.label];
    end
    title(tl, {['Orientation asymmetries across the visual hierarchy, ' vl], ...
        ['\rm\fontsize{9}' wt(A) ', t intervals on n-1 df.  V1/V2/V3 at 4-8' char(176) ...
         ';  V3a at 2-10' char(176) ' (open), the only band where it passes the coverage ' ...
         'criterion and so not on the same footing.']}, ...
        'Interpreter','tex','FontWeight','bold');

    for j = 1:4
        nexttile(j); hold on; yline(0,'-','Color',[0.75 0.75 0.75]);
        [xg, lbl] = deal([], {});
        for k = 1:numel(order)
            r1 = pick(A.asym, order{k}, band{k}, 'dg', nm{j});
            r2 = pick(A.asym, order{k}, band{k}, 'da', nm{j});
            if isempty(r1), continue; end
            xg(end+1) = k; lbl{end+1} = order{k}; %#ok<AGROW>
            open = strcmp(band{k},'2-10');
            draw(k-0.12, r1, cDg, open); draw(k+0.12, r2, cDa, open);
        end
        finish(xg, lbl, ttl{j}, j==1, 'asymmetry (% signal change)', yTop);
        if j == 1
            h(1) = plot(nan,nan,'o','MarkerFaceColor',cDg,'MarkerEdgeColor','k','MarkerSize',7);
            h(2) = plot(nan,nan,'o','MarkerFaceColor',cDa,'MarkerEdgeColor','k','MarkerSize',7);
            lg = legend(h, {'Cartesian (dg)','Polar (da)'}, 'FontSize',8,'Box','off','Location','best');
            lg.ItemTokenSize = [10 8];
        end
    end

    for j = 1:4
        nexttile(4+j); hold on; yline(0,'-','Color',[0.75 0.75 0.75]);
        [xg, lbl] = deal([], {});
        for k = 1:numel(order)
            r = pick(A.context, order{k}, band{k}, 'dg-da', nm{j});
            if isempty(r), continue; end
            xg(end+1) = k; lbl{end+1} = order{k}; %#ok<AGROW>
            draw(k, r, cCx, strcmp(band{k},'2-10'));
        end
        finish(xg, lbl, ['context: ' ttl{j}], j==1, 'dg - da (% signal change)', yBot);
    end

    if ~isfolder(figDir), mkdir(figDir); end
    tag = 'spec';
    if isfield(A,'variant'), tag = A.variant.tag; end
    base = fullfile(figDir, ['Figure_S5_' tag '_hierarchy']);
    exportgraphics(fig, [base '.png'], 'Resolution', 200);
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector');
    close(fig);
    out = [base '.png'];
    fprintf('plot_spec_hierarchy: wrote %s\n', out);
end

function s = wt(A)
    if isfield(A,'variant') && strcmp(A.variant.weighting,'precision')
        s = 'precision weighting';
    else
        s = 'equal weighting';
    end
end

function r = pick(T, area, band, exp, asym)
    m = strcmp(T.area,area) & strcmp(T.band,band) & strcmp(T.experiment,exp) & ...
        strcmp(T.asymmetry,asym) & T.reportable;
    if ~any(m), r = []; else, r = T(find(m,1), :); end
end

function draw(x, r, c, open)
    plot([x x], [r.t_lo r.t_hi], '-', 'Color', c, 'LineWidth', 1.8);
    plot([x-0.06 x+0.06], [r.t_lo r.t_lo], '-', 'Color', c, 'LineWidth', 1.2);
    plot([x-0.06 x+0.06], [r.t_hi r.t_hi], '-', 'Color', c, 'LineWidth', 1.2);
    fc = c; if open, fc = [1 1 1]; end
    plot(x, r.mean, 'o', 'MarkerSize', 8, 'MarkerFaceColor', fc, 'MarkerEdgeColor', c, 'LineWidth', 1.4);
end

function y = rowLimits(T, order, band, exps, nm)
% Span every interval the row will draw, over all maps, experiments and asymmetries.
    lo = []; hi = [];
    for k = 1:numel(order)
        for e = 1:numel(exps)
            for j = 1:numel(nm)
                m = strcmp(T.area,order{k}) & strcmp(T.band,band{k}) & ...
                    strcmp(T.experiment,exps{e}) & strcmp(T.asymmetry,nm{j}) & T.reportable;
                if ~any(m), continue; end
                r = T(find(m,1), :);
                if ~isfinite(r.t_lo) || ~isfinite(r.t_hi), continue; end
                lo(end+1) = min(r.t_lo, 0); %#ok<AGROW>
                hi(end+1) = max(r.t_hi, 0); %#ok<AGROW>
            end
        end
    end
    if isempty(lo), y = [-1 1]; return; end
    y = [min(lo) max(hi)];
    p = max(0.06 * diff(y), eps);
    y = y + [-p p];
end

function finish(xg, lbl, ttl, first, yl, ylims)
    if isempty(xg), axis off; title(ttl,'FontSize',10); return; end
    xlim([min(xg)-0.5 max(xg)+0.5]); ylim(ylims);
    set(gca,'XTick',xg,'XTickLabel',lbl,'FontSize',9);
    grid on; box off; ax = gca; ax.YGrid = 'on'; ax.XGrid = 'off';
    title(ttl,'FontSize',10,'FontWeight','bold');
    if first, ylabel(yl,'FontSize',9); end
end
