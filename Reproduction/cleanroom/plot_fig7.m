function outFile = plot_fig7(resDG, resDA, cfg, variantName)
% PLOT_FIG7  Reproduce Figure 7: joint-LME asymmetry weights, dg (A) and da (B).
%
% For each asymmetry the model contributes +beta to its "pro" condition and -beta
% to its "con" condition, so we draw a pro bar at +beta and a con bar at -beta
% (beta = delta/2). Error bars are the 68% bootstrap CI (subject resampling); a '*'
% marks asymmetries whose delta CI excludes 0.

    labels = {'horiz vs vert','card vs obl','rad vs tang','polC vs polO'};
    proC = [0.15 0.15 0.45; 0.20 0.55 0.25; 0.45 0.70 0.90; 0.55 0.35 0.15];
    conC = [0.95 0.60 0.10; 0.60 0.45 0.75; 0.80 0.15 0.15; 0.20 0.60 0.60];

    fig = figure('Color','w','Position',[100 100 1100 480],'Visible','off');
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Figure 7 reproduction (%s)  \\bf clean-room Path A', variantName), ...
          'Interpreter','tex');

    panels = {resDG, 'A: Cartesian gratings (dg)'; resDA, 'B: Polar gratings (da)'};
    ymax = 0;
    for pi = 1:2, ymax = max(ymax, max(abs(panels{pi,1}.delta))/2*1.6); end
    ymax = max(ymax, 0.3);

    for pi = 1:2
        res = panels{pi,1};
        nexttile(pi); hold on;
        for j = 1:4
            beta = res.delta(j)/2;
            ciB  = res.ci68(j,:)/2;                 % CI on beta scale
            xc = j;
            bar(xc-0.19,  beta, 0.36, 'FaceColor', proC(j,:), 'EdgeColor','none');
            bar(xc+0.19, -beta, 0.36, 'FaceColor', conC(j,:), 'EdgeColor','none');
            errorbar(xc-0.19, beta, beta-ciB(1), ciB(2)-beta, 'k', 'LineWidth',1, 'CapSize',6);
            sig = ~(res.ci68(j,1) <= 0 && 0 <= res.ci68(j,2));
            if sig, text(xc, ymax*0.9, '*', 'HorizontalAlignment','center','FontSize',18); end
        end
        yline(0,'k-');
        xlim([0.5 4.5]); ylim([-ymax ymax]);
        set(gca,'XTick',1:4,'XTickLabel',labels); xtickangle(20); grid on;
        ylabel(sprintf('BOLD (%s), pro at +\\beta', variantName),'Interpreter','tex');
        title(panels{pi,2});
    end

    if ~isfolder(cfg.figDir), mkdir(cfg.figDir); end
    base = fullfile(cfg.figDir, sprintf('fig7_LME_%s', variantName));
    exportgraphics(fig, [base '.png'], 'Resolution', 150);
    exportgraphics(fig, [base '.pdf']);
    close(fig);
    outFile = [base '.png'];
    fprintf('plot_fig7: wrote %s\n', outFile);
end
