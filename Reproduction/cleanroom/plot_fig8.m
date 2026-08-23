function outFile = plot_fig8(M, res, cfg, expCfg, variantName, figLabel)
% PLOT_FIG8  Per-location data + model polar plots (Figure 8, one experiment).
%
%   outFile = plot_fig8(M, res, cfg, expCfg, variantName, figLabel)
%
% One small polar plot per visual-field location (8 polar-angle wedges). Within
% each, the angular axis is the local orientation of each stimulus and the radius
% is the V1 response: black = observer data (mean across subjects), red = LME
% fixed-effect prediction (res.pred). The 4 orientations are duplicated (+180 deg)
% to close each plot. figLabel is e.g. 'A' (dg) or 'B' (da).

    paBins = cfg.paBins;
    nP = numel(paBins);
    dataMean = mean(M, 3, 'omitnan');           % nOri x nPA
    rmax = max([dataMean(:); res.pred(:)]) * 1.1;
    rmin = min([0; dataMean(:); res.pred(:)]);

    fig = figure('Color','w','Position',[100 100 1400 720],'Visible','off');
    tl  = tiledlayout(fig, 2, 4, 'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Figure 8%s reproduction (%s, %s)  \\bf clean-room Path A', ...
          figLabel, expCfg.name, variantName), 'Interpreter','tex');

    for p = 1:nP
        loc = mod(expCfg.oriAngle(:).' + expCfg.isPolar*(paBins(p) - 90), 180);  % 1 x nOri
        dV  = dataMean(:, p).';
        mV  = res.pred(:, p).';
        th  = [loc, loc + 180];
        dd  = [dV, dV]; mm = [mV, mV];
        [ths, ord] = sort(th);
        thc = deg2rad([ths, ths(1)]);

        nexttile(p);
        polarplot(thc, [dd(ord), dd(ord(1))], 'k-o', 'LineWidth', 1.5, ...
                  'MarkerFaceColor','k','MarkerSize',4); hold on;
        polarplot(thc, [mm(ord), mm(ord(1))], 'r-o', 'LineWidth', 1.5, ...
                  'MarkerFaceColor','r','MarkerSize',4);
        ax = gca; ax.RLim = [min(rmin,0) rmax];
        ax.ThetaTick = 0:45:315; ax.RTickLabel = [];
        title(sprintf('location %d\\circ', paBins(p)), 'Interpreter','tex');
        if p == 1
            legend({'data','model'}, 'Location','southoutside', 'Orientation','horizontal');
        end
    end

    if ~isfolder(cfg.figDir), mkdir(cfg.figDir); end
    base = fullfile(cfg.figDir, sprintf('fig8%s_%s_%s', figLabel, expCfg.name, variantName));
    exportgraphics(fig, [base '.png'], 'Resolution', 150);
    exportgraphics(fig, [base '.pdf']);
    close(fig);
    outFile = [base '.png'];
    fprintf('plot_fig8: wrote %s\n', outFile);
end
