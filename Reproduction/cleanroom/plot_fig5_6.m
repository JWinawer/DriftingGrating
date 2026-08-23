function outFile = plot_fig5_6(A, cfg, expCfg, variantName, figNum)
% PLOT_FIG5_6  Reproduce a Figure 5 (dg) or Figure 6 (da) panel set.
%
%   outFile = plot_fig5_6(A, cfg, expCfg, variantName, figNum)
%
% Top row : one polar plot per asymmetry (theta = 8 polar-angle wedges, radius =
%           across-subject mean of the pro and con category responses).
% Bottom  : per-subject pairwise plot (pro vs con, each = mean over the 8 wedges),
%           with the group mean and a bootstrap 95% CI on the difference; a '*'
%           marks CIs excluding 0.
%
% Saves a PNG + PDF to cfg.figDir and returns the PNG path.

    asyms  = A.order;
    titles = {'horizontal vs vertical','cardinal vs oblique', ...
              'radial vs tangential','polar-card vs polar-obl'};
    % pro/con colors per asymmetry (loosely matching the manuscript scheme)
    proC = {[0.15 0.15 0.45],[0.20 0.55 0.25],[0.45 0.70 0.90],[0.55 0.35 0.15]};
    conC = {[0.95 0.60 0.10],[0.60 0.45 0.75],[0.80 0.15 0.15],[0.20 0.60 0.60]};

    th   = deg2rad([cfg.paBins, cfg.paBins(1)]);   % close the loop
    nS   = numel(cfg.subjects);

    fig = figure('Color','w','Position',[100 100 1500 700], 'Visible','off');
    tl  = tiledlayout(fig, 2, 4, 'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Figure %d reproduction (%s, %s)  \\bf clean-room Path A', ...
          figNum, expCfg.name, variantName), 'Interpreter','tex');

    % rho limits: symmetric, common across the polar plots of this figure
    allvals = [];
    for i = 1:numel(asyms)
        allvals = [allvals; mean(A.(asyms{i}).pro,2); mean(A.(asyms{i}).con,2)]; %#ok<AGROW>
    end
    rmax = max(abs(allvals)) * 1.1;

    % --- top row: polar plots ---
    for i = 1:numel(asyms)
        s = A.(asyms{i});
        pro = mean(s.pro, 2); con = mean(s.con, 2);
        nexttile(i);
        polarplot(th, [pro; pro(1)], '-o', 'Color', proC{i}, 'LineWidth', 2, ...
                  'MarkerFaceColor', proC{i}, 'MarkerSize', 4); hold on;
        polarplot(th, [con; con(1)], '-o', 'Color', conC{i}, 'LineWidth', 2, ...
                  'MarkerFaceColor', conC{i}, 'MarkerSize', 4);
        ax = gca; ax.RLim = [-rmax rmax]; ax.ThetaZeroLocation = 'right';
        ax.ThetaDir = 'counterclockwise'; ax.RTickLabel = [];
        title(titles{i}, 'FontWeight','bold');
    end

    % --- bottom row: pairwise ---
    for i = 1:numel(asyms)
        s = A.(asyms{i});
        proSubj = mean(s.pro, 1);          % 1 x nSubj (avg over PA)
        conSubj = mean(s.con, 1);
        diffSubj = proSubj - conSubj;
        grp = mean(diffSubj);
        ci  = bootstrap_ci(diffSubj, cfg.nBoot, cfg.ciLevel, 0);
        sig = ~(ci(1) <= 0 && 0 <= ci(2));

        nexttile(4 + i); hold on;
        for si = 1:nS
            plot([1 2], [proSubj(si) conSubj(si)], '-', 'Color', [0.6 0.6 0.6]);
        end
        plot(1, mean(proSubj), 'o', 'MarkerFaceColor', proC{i}, 'MarkerEdgeColor','k', 'MarkerSize', 9);
        plot(2, mean(conSubj), 'o', 'MarkerFaceColor', conC{i}, 'MarkerEdgeColor','k', 'MarkerSize', 9);
        xlim([0.5 2.5]); set(gca,'XTick',[1 2],'XTickLabel',{s.proName, s.conName});
        xtickangle(20); grid on;
        ttl = sprintf('\\Delta=%.3f [%.2f, %.2f]%s', grp, ci(1), ci(2), repmat(' *',1,sig));
        title(ttl, 'Interpreter','tex');
        if i == 1, ylabel(sprintf('V1 %% change (%s)', variantName)); end
    end

    if ~isfolder(cfg.figDir), mkdir(cfg.figDir); end
    base = fullfile(cfg.figDir, sprintf('fig%d_%s_%s', figNum, expCfg.name, variantName));
    exportgraphics(fig, [base '.png'], 'Resolution', 150);
    exportgraphics(fig, [base '.pdf']);
    close(fig);
    outFile = [base '.png'];
    fprintf('plot_fig5_6: wrote %s\n', outFile);
end
