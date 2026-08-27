function plotAsymmetryAcrossROIs(projectSettings)
% PLOTASYMMETRYACROSSROIS  Successor to lme1_fit.m's "plot across ROIs
% per figure" visualization: one figure PER ASYMMETRY, with cortical
% areas along the x-axis (the inverse layout of plotROISummary.m, which
% is one figure per cortical area with asymmetries on the x-axis).
% Sourced from fitAsymmetryRegression.m's cached joint regression fit
% instead of the LME. Fully independent of lme1_fit.m -- does not read or
% write anything under LME_results/; lme1_fit.m and its output remain
% untouched and independently runnable.
%
%   plotAsymmetryAcrossROIs(projectSettings)
%
% projectSettings - same struct plot1/plot2_experimentalCond.m and
%                    plotROISummary.m use (must include .rois, .roi_idx,
%                    .comparisonName, .projectName, .colors_data,
%                    .figureDir, .gainSummaryFile).
%
% Pro/con points are MODEL-DERIVED, zero-centered (pro = +beta,
% con = -beta), and error bars are the 68% CI of the cached fit's
% bootstrapped difference, full width, applied symmetrically to both --
% identical conventions to plotROISummary.m (see that file for the full
% derivation of why: matches lme1_fit.m's ACTUAL plotted behavior given
% its meanRelative=1 setting, and the full-width CI approximates the
% conventional ~95% "non-overlap" significance heuristic). Significance
% asterisks at a fixed y=0.35 above each cortical area's point: ** if the
% 95% CI of the difference excludes 0, * if only the 68% CI does.

projectName = projectSettings.projectName;
comparisonName = projectSettings.comparisonName;
colors_data = projectSettings.colors_data;
rois = projectSettings.rois;
figureDir = projectSettings.figureDir;

[summaryTablesDir,~,~] = fileparts(projectSettings.gainSummaryFile);
[derivativesDir,~,~] = fileparts(summaryTablesDir);
[bidsDir,~,~] = fileparts(derivativesDir);

% Concept label and COLORS.json style key per raw term (termIdx 1..4 =
% mainCardinal, derivedCardinal, mainSubset, derivedSubset) -- same
% mapping as plotROISummary.m's colorKeys, plus the human-readable
% concept label used there for xLabels but here indexed directly by
% termIdx (no cross-asymmetry x-axis reordering needed here, since each
% asymmetry gets its own figure).
if strcmp(projectName, 'dg')
    conceptLabels = {'Cardinal vs Oblique','Polar Cardinal vs Polar Oblique','Horizontal vs Vertical','Radial vs Tangential'};
    colorKeys = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','verticalVsHorizontal','radialVsTangential'};
elseif strcmp(projectName, 'da')
    conceptLabels = {'Polar Cardinal vs Polar Oblique','Cardinal vs Oblique','Radial vs Tangential','Horizontal vs Vertical'};
    colorKeys = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','radialVsTangential','verticalVsHorizontal'};
else
    error('plotAsymmetryAcrossROIs:project', 'projectName must be ''dg'' or ''da''.');
end

nROIs = length(rois);
statsRows = {};

% Load all cortical areas' cached fits once, up front (each figure below
% needs all of them, not just one).
% fitLabel: which regressionResults/<label>/ subfolder to read --
% defaults to projectName, but for dg 'matched' mode (7 subjects, same
% ones run in da) set projectSettings.fitLabel = 'dgMatched7' explicitly.
if isfield(projectSettings, 'fitLabel') && ~isempty(projectSettings.fitLabel)
    fitLabel = projectSettings.fitLabel;
else
    fitLabel = projectName;
end
F = cell(nROIs,1);
for ri = 1:nROIs
    fitFile = fullfile(bidsDir,'derivatives','summaryTables','regressionResults',fitLabel,sprintf('%s.mat',rois{ri}));
    if ~isfile(fitFile)
        error('plotAsymmetryAcrossROIs:missingFit', ...
            'Cached fit not found for %s / %s at %s -- run fitAsymmetryRegression(''%s'') first.', ...
            fitLabel, rois{ri}, fitFile, projectName);
    end
    F{ri} = load(fitFile);
end

for termIdx = 1:4
    figure; hold on

    styleInfo = colors_data.conditions.(projectName).(colorKeys{termIdx});
    colorPro = styleInfo.color_pro';
    colorCon = 0.5*styleInfo.color_con' + 0.5*[1 1 1]; % same 50%-white blend used throughout plot1/plot2/plotROISummary

    for ri = 1:nROIs
        beta = F{ri}.estimates(termIdx) / 2;
        dotPro = beta;
        dotCon = -beta;

        bootDraws = F{ri}.coeffs(termIdx, :)';
        estDiff = F{ri}.estimates(termIdx);
        ci68 = prctile(bootDraws, [16 84]);
        ci95 = prctile(bootDraws, [2.5 97.5]);
        ci68_halfwidth = (ci68(2) - ci68(1)) / 2;

        statsRows(end+1,:) = {rois{ri}, conceptLabels{termIdx}, dotPro, dotCon, estDiff, ci68(1), ci68(2), ci95(1), ci95(2)}; %#ok<AGROW>

        x = ri;
        % boxchart, matching this figure's original style. No separate
        % edge color: BoxEdgeColor = BoxFaceColor always. Whichever of
        % pro/con is bigger gets the full-strength color; the smaller one
        % gets the 50%-white-blended color (colorCon, reused regardless
        % of whether it's pro or con that ends up smaller, since
        % color_pro==color_con by design). Since dotCon=-dotPro, "con
        % bigger" reduces to dotPro<0.
        if dotPro >= dotCon
            proColor = colorPro; conColor = colorCon;
        else
            proColor = colorCon; conColor = colorPro;
        end
        boxchart(x, dotPro, 'BoxFaceColor', proColor, 'BoxEdgeColor', proColor, 'LineWidth', 4, 'BoxWidth', 0.6);
        errorbar(x, dotPro, ci68_halfwidth, ci68_halfwidth, 'LineStyle','none', 'LineWidth', 2, 'Color', proColor);

        boxchart(x, dotCon, 'BoxFaceColor', conColor, 'BoxEdgeColor', conColor, 'LineWidth', 4, 'BoxWidth', 0.6);
        errorbar(x, dotCon, ci68_halfwidth, ci68_halfwidth, 'LineStyle','none', 'LineWidth', 2, 'Color', conColor);

        sig95 = ci95(1) > 0 || ci95(2) < 0;
        sig68 = ci68(1) > 0 || ci68(2) < 0;
        if sig95
            sigStr = '**';
        elseif sig68
            sigStr = '*';
        else
            sigStr = '';
        end
        if ~isempty(sigStr)
            text(x, 0.35, sigStr, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 14, 'Color', [0 0 0]);
        end
    end

    yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)
    ylim([-0.4 0.4])
    xlim([0.5 nROIs+0.5])
    set(gca,'XTick',[])
    box off
    set(gca,'linewidth',2, 'YColor',[0 0 0], 'XColor',[0 0 0]);
    set(gca,'FontName','Arial','FontSize',20);
    xticks(1:nROIs)
    xticklabels(rois)
    xtickangle(25)
    ylabel('zscored BOLD psc', 'FontSize', 20);
    title(conceptLabels{termIdx}, 'FontSize', 20, 'Interpreter', 'none');

    f = gcf;
    f.Position = [298 843 651 494];
    drawnow;

    print(fullfile(figureDir, sprintf('AsymmetryAcrossROIs_%s_%s_%s', comparisonName, projectName, colorKeys{termIdx})), '-dpdf', '-bestfit');
    close all;
end

statsTable = cell2table(statsRows, 'VariableNames', ...
    {'roi','asymmetry','pro_mean','con_mean','diff_estimate','ci68_lower','ci68_upper','ci95_lower','ci95_upper'});
writetable(statsTable, fullfile(figureDir, sprintf('AsymmetryAcrossROIs_%s_%s_stats.csv', comparisonName, projectName)));

end
