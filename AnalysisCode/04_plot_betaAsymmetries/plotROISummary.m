function plotROISummary(projectSettings)
% PLOTROISUMMARY  Successor to lme1_fit.m's per-cortical-area "master
% figure" visualization: one figure per cortical area, all 4 asymmetries
% along the x-axis, pro vs con each shown as a point + error bar. Sourced
% from fitAsymmetryRegression.m's cached joint regression fit instead of
% the LME. Fully independent of lme1_fit.m -- does not read or write
% anything under LME_results/; lme1_fit.m and its output remain untouched
% and independently runnable.
%
%   plotROISummary(projectSettings)
%
% projectSettings - same struct plot1/plot2_experimentalCond.m use (must
%                    include .rois, .roi_idx, .comparisonName,
%                    .projectName, .colors_data, .figureDir,
%                    .gainSummaryFile -- the last only to locate bidsDir,
%                    same convention as plot2_experimentalCond.m).
%
% Pro/con points are MODEL-DERIVED, zero-centered: pro = +beta,
% con = -beta (beta = estimates(termIdx)/2, since estimates is already
% the beta*2 pro-minus-con scale) -- the direct readout of the +-1/0
% design coding, with no intercept added back in. This matches what
% lme1_fit.m's master figure ACTUALLY plots: it locally computes
% Gintercept+beta/Gintercept-beta, but then (since meanRelative=1 is its
% configured setting) subtracts baselineSub=Gintercept back out before
% storing/plotting (see its mean_pro(ai) = y1(x)-baselineSub, line ~543),
% so the intercept never actually appears in its figure -- only beta and
% -beta do. This is a DELIBERATE difference from
% plot1_experimentalCond.m/plot2_experimentalCond.m, whose dots are
% data-derived (absolute BOLD levels) rather than model-derived
% (asymmetry-only, zero-centered) -- this script's whole purpose is to
% reproduce lme1_fit.m's own design (all 4 asymmetries sharing a common
% center at exactly 0, representing the asymmetry effect alone, decoupled
% from the cortical area's overall response magnitude), not plot2's
% "dots stay data-derived" convention.
%
% Error bars: 68% CI of the cached fit's bootstrapped difference
% (coeffs), FULL half-width applied symmetrically to both pro and con --
% the SAME convention plot1_experimentalCond.m/plot2_experimentalCond.m
% already use (confirmed deliberately, not lme1_fit.m's own convention of
% mirroring the raw, undoubled coefficient's own CI onto each side, which
% corresponds to a more liberal ~68%-of-the-difference significance
% boundary rather than the ~95%-ish one this convention approximates).

projectName = projectSettings.projectName;
comparisonName = projectSettings.comparisonName;
colors_data = projectSettings.colors_data;
rois = projectSettings.rois;
figureDir = projectSettings.figureDir;

% Locate fitAsymmetryRegression.m's cached output (same convention as
% plot2_experimentalCond.m).
[summaryTablesDir,~,~] = fileparts(projectSettings.gainSummaryFile);
[derivativesDir,~,~] = fileparts(summaryTablesDir);
[bidsDir,~,~] = fileparts(derivativesDir);

% x-axis order and display labels: matches lme1_fit.m's plotOrder/
% asymLabel exactly (desired order: vertical-vs-horizontal,
% cardinal-vs-oblique, radial-vs-tangential, polar-cardinal-vs-oblique,
% relabeled per project since which raw term is which concept swaps
% between dg/da -- see fitAsymmetryRegression_dgVsDa.m's
% dgTermForConcept/daTermForConcept for the same mapping). termIdx here
% indexes termNames = {mainCardinal, derivedCardinal, mainSubset,
% derivedSubset}, the order fitAsymmetryRegression.m saves estimates/
% coeffs in.
if strcmp(projectName, 'dg')
    xOrderTermIdx = [3, 1, 4, 2];
    xLabels = {'Horizontal vs Vertical','Cardinal vs Oblique','Radial vs Tangential','Polar Cardinal vs Polar Oblique'};
    colorKeys = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','verticalVsHorizontal','radialVsTangential'}; % indexed by termIdx (1..4), not x-position
elseif strcmp(projectName, 'da')
    xOrderTermIdx = [4, 2, 3, 1];
    xLabels = {'Horizontal vs Vertical','Cardinal vs Oblique','Radial vs Tangential','Polar Cardinal vs Polar Oblique'};
    colorKeys = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','radialVsTangential','verticalVsHorizontal'};
else
    error('plotROISummary:project', 'projectName must be ''dg'' or ''da''.');
end
nA = 4;

% fitLabel: which regressionResults/<label>/ subfolder to read --
% defaults to projectName, but for dg 'matched' mode (7 subjects, same
% ones run in da) set projectSettings.fitLabel = 'dgMatched7' explicitly.
if isfield(projectSettings, 'fitLabel') && ~isempty(projectSettings.fitLabel)
    fitLabel = projectSettings.fitLabel;
else
    fitLabel = projectName;
end

statsRows = {};

for ri = 1:length(rois)
    roiname = rois{ri};

    fitFile = fullfile(bidsDir,'derivatives','summaryTables','regressionResults',fitLabel,sprintf('%s.mat',roiname));
    if ~isfile(fitFile)
        error('plotROISummary:missingFit', ...
            'Cached fit not found for %s / %s at %s -- run fitAsymmetryRegression(''%s'') first.', ...
            fitLabel, roiname, fitFile, projectName);
    end
    F = load(fitFile);

    figure; hold on

    for ai = 1:nA
        termIdx = xOrderTermIdx(ai);
        beta = F.estimates(termIdx) / 2;
        dotPro = beta;
        dotCon = -beta;

        bootDraws = F.coeffs(termIdx, :)';
        estDiff = F.estimates(termIdx);
        ci68 = prctile(bootDraws, [16 84]);
        ci95 = prctile(bootDraws, [2.5 97.5]);
        ci68_halfwidth = (ci68(2) - ci68(1)) / 2;

        statsRows(end+1,:) = {roiname, xLabels{ai}, dotPro, dotCon, estDiff, ci68(1), ci68(2), ci95(1), ci95(2)}; %#ok<AGROW>

        styleInfo = colors_data.conditions.(projectName).(colorKeys{termIdx});
        colorPro = styleInfo.color_pro';
        % con: same 50%-white blend used throughout plot1/plot2 for the
        % con condition (color_pro == color_con by design in COLORS.json
        % -- pro/con are meant to be distinguished by this blend, not by
        % hue). This is con's own "identity" color, used whenever con is
        % filled below -- not pro's full-strength color.
        colorCon = 0.5*styleInfo.color_con' + 0.5*[1 1 1];

        x = ai;
        % boxchart, matching this figure's original style. No separate
        % edge color: BoxEdgeColor = BoxFaceColor always. Whichever of
        % pro/con is bigger gets the full-strength color; the smaller one
        % gets the 50%-white-blended color (colorCon, reused regardless
        % of whether it's pro or con that ends up smaller, since
        % color_pro==color_con by design -- same hue either way). Since
        % dotCon=-dotPro, "con bigger" reduces to dotPro<0.
        if dotPro >= dotCon
            proColor = colorPro; conColor = colorCon;
        else
            proColor = colorCon; conColor = colorPro;
        end
        boxchart(x, dotPro, 'BoxFaceColor', proColor, 'BoxEdgeColor', proColor, 'LineWidth', 4, 'BoxWidth', 0.6);
        errorbar(x, dotPro, ci68_halfwidth, ci68_halfwidth, 'LineStyle','none', 'LineWidth', 2, 'Color', proColor);

        boxchart(x, dotCon, 'BoxFaceColor', conColor, 'BoxEdgeColor', conColor, 'LineWidth', 4, 'BoxWidth', 0.6);
        errorbar(x, dotCon, ci68_halfwidth, ci68_halfwidth, 'LineStyle','none', 'LineWidth', 2, 'Color', conColor);

        % Significance asterisk at a fixed y=0.35 for this asymmetry: two
        % asterisks if the 95% CI of the difference excludes 0, one if
        % only the 68% CI does (95% always implies 68%, since it's the
        % wider, nested interval) -- same convention as
        % plot2_experimentalCond.m's asterisks.
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
    xlim([0.5 nA+0.5])
    set(gca,'XTick',[])
    box off
    set(gca,'linewidth',2, 'YColor',[0 0 0], 'XColor',[0 0 0]);
    set(gca,'FontName','Arial','FontSize',20);
    xticks(1:nA)
    xticklabels(xLabels)
    xtickangle(25)
    ylabel('zscored BOLD psc', 'FontSize', 20);
    title(roiname, 'FontSize', 20, 'Interpreter', 'none');

    f = gcf;
    f.Position = [298 843 651 494];
    drawnow;

    print(fullfile(figureDir, sprintf('ROIsummary_%s_%s_%s', comparisonName, projectName, roiname)), '-dpdf', '-bestfit');
    close all;
end

statsTable = cell2table(statsRows, 'VariableNames', ...
    {'roi','asymmetry','pro_mean','con_mean','diff_estimate','ci68_lower','ci68_upper','ci95_lower','ci95_upper'});
writetable(statsTable, fullfile(figureDir, sprintf('ROIsummary_%s_%s_stats.csv', comparisonName, projectName)));

end
