function plot2_experimentalCond(medianBOLDpa, asymmetryName, projectSettings, varargin)

    rng('default')
    %rng(0)
    % Plot mean across polar angles
    % keep in mind that this equally weighs each PA, whereas there could be
    % differential # of voxels representing the PAs

    % for now, ensure a third argument is always provided
    if nargin < 3 || isempty(varargin{1})
        error('A third input is required when projectName is "%s" and asymmetryName is "%s".', projectSettings.projectName, asymmetryName);
    else
        subset = varargin{1};
    end

    % Check project name & request
    if strcmp(projectSettings.projectName, 'da') || strcmp(projectSettings.projectName, 'dots')
        if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
            derivedVals = 0; % this is needed because radialVsTang occurs for both conditions
            % Ensure a third argument is provided
            if subset == 1
                asymmetryName = 'radialVsTangential';
                %radialvstang = 1;
            end
        elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
             % Third argument is ignored
             derivedVals = 1;
             if subset == 1
                asymmetryName = 'verticalVsHorizontal';
                %radialvstang = 0;
             end
        end
    elseif strcmp(projectSettings.projectName, 'dg')
        if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
            % Third argument is ignored
            derivedVals = 0; 
            if subset == 1
                asymmetryName = 'verticalVsHorizontal';
                %radialvstang = 0;
            end
        elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
            derivedVals = 1;
            if subset == 1
                asymmetryName = 'radialVsTangential';
                %radialvstang = 1;
            end
        end
    else
        error('Unknown projectName. Expected "da" or "dg".');
    end
    
    projectName = projectSettings.projectName;
    comparisonName = projectSettings.comparisonName;
    colors_data = projectSettings.colors_data;
    rois = projectSettings.rois;
    %pairaxes_limits = projectSettings.pairaxes_limits;
    pairaxes_PAew_limits = projectSettings.pairaxes_PAew_limits;
    figureDir = projectSettings.figureDir;

    styleInfo = colors_data.conditions.(projectName).(asymmetryName);
    colors = {styleInfo.color_pro', styleInfo.color_con'};

    % pro = filled marker, full color; con = unfilled marker, same
    % 50%-white-blended color used in plot1_experimentalCond.m's polar
    % plots (color_pro == color_con by design -- pro/con distinguished by
    % fill state and this alpha-equivalent blend, not by hue). Blending
    % toward white rather than true alpha for the same reason as the
    % polar plots: painters-rendered vector PDF export doesn't support
    % line/marker transparency.
    proFaceColor = colors{1};
    conFaceColor = [1 1 1]; % white fill (not unfilled)
    proEdgeColor = colors{1};
    conEdgeColor = 0.5*colors{2} + 0.5*[1 1 1];
    subjectLineColor = [0.8196, 0.8275, 0.8314]; % RGB(209,211,212)
    % Match plot1_experimentalCond.m's polar-plot marker size exactly.
    % plot()'s MarkerSize is a linear (roughly diameter, in points)
    % measure; scatter()'s size argument is marker AREA (points^2). To
    % render the same apparent diameter d, a circular scatter marker
    % needs area = pi/4 * d^2.
    polarMarkerSize = 6 * 0.8; % must match plot1_experimentalCond.m's markerSize
    meanDotSize = (pi/4) * polarMarkerSize^2;
    axisLineWidth = 1; % matches plot1_experimentalCond.m's polar-plot axis lines
    axisLineColor = [0.25 0.25 0.25]; % matches plot1_experimentalCond.m's polar-plot axis/grid color
    errorbarLineWidth = styleInfo.errorbar_lineWidth; % same as plot1_experimentalCond.m's polar-plot error bars
    proLineWidth = styleInfo.pro_lineWidth; % marker outline thickness, matches polar plot's markers
    conLineWidth = styleInfo.con_lineWidth;
    showTitleLegend = false;

    % Narrow xlim (a shift, cropping in on the existing x-range) rather
    % than shrinking the axis Position width directly, so the physical
    % cm-per-data-unit scale stays exactly what it was -- i.e. the plotted
    % points/spacing are not compressed, only how much surrounding margin
    % is shown changes. previousXlim/previousPairwisePlotWidth_cm
    % describe the last (already-tuned) width=3.525cm-at-xlim=[0.5,2.5]
    % configuration; solving for the xlim that yields exactly 3cm at that
    % same scale, centered on the same midpoint.
    previousXlim = [0.5, 2.5];
    previousPairwisePlotWidth_cm = 3 * 1.25 * 0.94; % ~3.525cm
    targetPairwisePlotWidth_cm = 3;
    xScale_cm_per_unit = previousPairwisePlotWidth_cm / diff(previousXlim);
    newXlimRange = targetPairwisePlotWidth_cm / xScale_cm_per_unit;
    xlimCenter = mean(previousXlim);
    pairwiseXlim = xlimCenter + [-newXlimRange/2, newXlimRange/2];

    % retrieve the indices for specific asymmetries (e.g., motion -
    % orientation for main cardinal v main oblique)
    if ~derivedVals
        [proConditions, conConditions, ~] = retrieveProConIdx(projectName, comparisonName, subset);
    else
        proConditions = 1; conConditions = 2;
    end

    % Accumulates the exact numbers underlying this figure (one row per
    % ROI) across the loop below, so a full top-to-bottom run leaves a
    % persistent, inspectable record of the plotted values -- not just
    % console output that scrolls away, and not something you have to
    % re-run (with its own fresh bootstrap draws) to see again.
    statsRows = {};

    figure

    for ii=1:length(rois)
        rois{ii}
    
        % Specify the region index
        regionIndex = projectSettings.roi_idx{ii};
    
        % Extract the relevant conditions for the specified region
        conditions1 = medianBOLDpa(proConditions, :, regionIndex, :);
        conditions2 = medianBOLDpa(conConditions, :, regionIndex, :);
        conditions1 = squeeze(conditions1);
        conditions2 = squeeze(conditions2);
        
        % Average across the conditions within subjects
        avgConditions1 = nanmean(conditions1, 1);
        avgConditions2 = nanmean(conditions2, 1);
        avgConditions1 = squeeze(avgConditions1);
        avgConditions2 = squeeze(avgConditions2);

        % Gain-weight each observer: divide their value (irrespective of
        % polar angle) by their own mean pRF gain (prfvista_mov/prfvista
        % average), BEFORE any averaging across observers. This down-weights
        % high-gain observers and up-weights low-gain observers; the
        % displayed statistics, plotted per-subject lines/scatter points,
        % and bootstrapped error bars below all inherit the adjustment
        % automatically since they are computed from these arrays. See
        % projectSettings.observerGain / retrieveObserverGainWeights.m
        %
        % The across-observer average gain is then multiplied back in, so
        % the plotted scale/units resemble the original (unweighted) data.
        % Dividing by gain_i and then multiplying by groupGain is the same
        % as dividing by gain_i normalized to the group mean (gain_i /
        % groupGain), so the relative weighting -- and therefore the
        % relative pattern across per-subject lines, group markers, and
        % error bars -- is unchanged; only the overall scale shifts.
        %
        % groupGain uses the geometric, not arithmetic, mean: the applied
        % factor (groupGain ./ gainWeights) is itself a multiplicative
        % scale factor, and geometric mean is the choice under which the
        % *geometric* mean of that factor across observers is exactly 1
        % (magnitude-neutral in the multiplicative sense, matching what
        % the factor actually is) -- arithmetic mean does not have this
        % property. Implemented via exp(mean(log(.))) to avoid a
        % dependency on the Statistics and Machine Learning Toolbox.
        gainWeights = projectSettings.observerGain;
        if any(gainWeights <= 0)
            error('gainWeights must be strictly positive to take log() for the geometric mean (found %d non-positive value(s))', ...
                sum(gainWeights <= 0));
        end
        groupGain = exp(mean(log(gainWeights)));
        avgConditions1 = avgConditions1 .* (groupGain ./ gainWeights);
        avgConditions2 = avgConditions2 .* (groupGain ./ gainWeights);

        % Already averaged across the conditions within subjects < -- already did this in
        % the loop above
        
        % Extract polar angles
        %anglevals = [90, 135, 180, 225, 270, 315, 0, 45];
        anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
        
        vals_1 = nanmean(avgConditions1,1)';
        vals_2 = nanmean(avgConditions2,1)';
    
        vals_1_overall = nanmean(avgConditions1,'all')';
        vals_2_overall = nanmean(avgConditions2,'all')';
        %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
        %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
        
        % Plot the data on a polar plot
        if length(rois) == 1
            subplot(1,1,ii)
        else
            subplot(2,5,ii)
        end
        
        grandMean = mean([vals_1, vals_2], 'all');


        for subjectIndex = 1:size(medianBOLDpa, 4)
            temp = [vals_1(subjectIndex) vals_2(subjectIndex)];
            tempMean = mean(temp);
            %scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); % [127/255, 191/255, 123/255]
%             scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none'); % 
%             hold on
%             %scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %  [175/255, 141/255, 195/255]
%             scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', 'none'); 
            
            % adding grandMean back, subtracting subject mean, otherwise:
            % [vals_1 vals_2]
            plot([1 2], grandMean + [vals_1(subjectIndex)-tempMean vals_2(subjectIndex)-tempMean], 'Color', subjectLineColor); %colors_rgb(subjectIndex,:))
            xlim(pairwiseXlim)
    %         ylim([-0.15 0.25])
        hold on
        end
        
        % if i need to plot asymmetry itself

        
%         differences = vals_1-vals_2;
%         for subjectIndex = 1:size(medianBOLDpa, 4)
%             scatter(1.5, differences(subjectIndex), 30, 'MarkerFaceColor', [.85 .85 .85], 'MarkerEdgeColor', 'none');
%         end
        hold on

%         errDiff = std(vals_1-vals_2)/(sqrt((size(medianBOLDpa, 4))));
%         errorbar(1.5, mean(vals_1-vals_2), errDiff, 'k', 'LineWidth', 3);
%         [h, p] = ttest(vals_1-vals_2, 0)

        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        % Assume: vals_1 and vals_2 are column vectors [nSubjects x 1]
        differences = vals_1 - vals_2;
        nBoot = 10000;
        nSubs = length(differences);
        
        %% Bootstrap Mean
        bootMeans = zeros(nBoot, 1);
        for i = 1:nBoot
            resample = datasample(differences, nSubs);
            bootMeans(i) = mean(resample);
        end
        % Observed mean difference
        meanDiff = mean(differences);

        % 95% bootstrap confidence interval
        ci_mean = prctile(bootMeans, [2.5 97.5]);

        % 68% bootstrap confidence interval, plotted below as the error
        % bar (half-width applied symmetrically to both vals_1 and
        % vals_2, same simplification used elsewhere in this file: this
        % CI describes the (pro-con) difference, not either condition's
        % own uncertainty separately, so there's no unique way to split
        % it between the two markers).
        ci_mean_68 = prctile(bootMeans, [16 84]);
        ci68_halfwidth = (ci_mean_68(2) - ci_mean_68(1)) / 2;

        % record this ROI's aggregate stats for the export below
        statsRows(end+1,:) = {rois{ii}, meanDiff, ci_mean_68(1), ci_mean_68(2), ci_mean(1), ci_mean(2)}; %#ok<SAGROW>

%         errLower_mean = meanDiff - ci_mean(1);
%         errUpper_mean = ci_mean(2) - meanDiff;
        
%         %% Bootstrap Median
%         bootMedians = zeros(nBoot, 1);
%         for i = 1:nBoot
%             resample = datasample(differences, nSubs);
%             bootMedians(i) = median(resample);
%         end
%         medianDiff = median(differences);
%         ci_median = prctile(bootMedians, [2.5 97.5]);
%         errLower_median = medianDiff - ci_median(1);
%         errUpper_median = ci_median(2) - medianDiff;


        
        %% Print to console
        fprintf('Mean difference: %.4f, 95%% CI: [%.4f, %.4f]\n', meanDiff, ci_mean(1), ci_mean(2));
%         fprintf('Median difference: %.4f, 95%% CI: [%.4f, %.4f]\n', medianDiff, ci_median(1), ci_median(2));
        differences
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        

        %plot([1 2], [vals_1_overall vals_2_overall], 'k', 'LineWidth', 3)
        %hold on
        % Mean dots drawn first; marker outline thickness matches
        % plot1_experimentalCond.m's polar-plot markers (proLineWidth/
        % conLineWidth); size is unaffected -- meanDotSize (SizeData,
        % i.e. marker area) is set independently of LineWidth (the edge
        % stroke width).
        scatter(1, vals_1_overall, meanDotSize, 'MarkerFaceColor', proFaceColor, 'MarkerEdgeColor', proEdgeColor, 'LineWidth',proLineWidth); %, 'MarkerFaceAlpha', 0.5);
        hold on
        scatter(2, vals_2_overall,  meanDotSize, 'MarkerFaceColor', conFaceColor, 'MarkerEdgeColor', conEdgeColor, 'LineWidth',conLineWidth); %, 'MarkerFaceAlpha', 0.5);
        hold on

        % Error bars -- top layer, drawn last -- matches
        % plot1_experimentalCond.m's polar-plot draw order.
        errorbar(1, mean(vals_1), ci68_halfwidth, 'Color', proEdgeColor, 'LineWidth', errorbarLineWidth, 'CapSize', 0);
        hold on
        errorbar(2, mean(vals_2), ci68_halfwidth, 'Color', conEdgeColor, 'LineWidth', errorbarLineWidth, 'CapSize', 0);

        if showTitleLegend
            title(rois{ii});
        end
        %ylabel('zscored PSC')
        set(gca, 'XTick', []);
        ax = gca;
        ax.LineWidth = axisLineWidth;
        ax.XColor = axisLineColor; % matches plot1_experimentalCond.m's polar-plot ThetaColor/RColor
        ax.YColor = axisLineColor;
        % ax.Layer left at its default ('bottom') -- matches
        % plot1_experimentalCond.m's polar plots, where 'top' drew grid
        % lines over solid data points.
        box off

        if ismember(rois{ii}, {'pMT', 'pMST', 'hMTcomplex'})
            ROI_category = 'ROIs_motion';
        else
            ROI_category = 'ROIs_early';
        end

        ylim([pairaxes_PAew_limits.(projectName).(comparisonName).(ROI_category).min ...
                     pairaxes_PAew_limits.(projectName).(comparisonName).(ROI_category).max])
    
        if ii==1 && showTitleLegend
            if strcmp(asymmetryName, 'radialVsTangential')
                lg1 = legend('Radial', 'Tangential', 'Location', 'northeast');
            elseif strcmp(asymmetryName, 'verticalVsHorizontal')
                lg1 = legend('Horizontal', 'Vertical', 'Location', 'northeast');
            else
                lg1 = legend('Card', 'Obl', 'Location', 'northeast');
            end
            lg1.Position = [0.7946 0.2614 0.1089 0.0631];
        end
    
        hold off;
    end
    
    if showTitleLegend
        sgtitle(sprintf('%s: %s', projectName, strrep(comparisonName, '_', ' ')), 'FontSize', 40)
    end

%     fa = gcf;
%     fa.Position = [1000 555 1514 782];

    filename = fullfile(figureDir,sprintf('pairwise_PAequalweight_%s_%s_%s', comparisonName, projectName, asymmetryName));

    % Save the exact values underlying this figure alongside the PDF, so
    % a full top-to-bottom run leaves a persistent, reproducible record
    % (diffable across runs) rather than requiring a re-run -- with its
    % own fresh bootstrap draws -- to see the numbers again.
    statsTable = cell2table(statsRows, 'VariableNames', ...
        {'roi','meanDiff','ci68_lower','ci68_upper','ci95_lower','ci95_upper'});
    writetable(statsTable, [filename, '_stats.csv']);

    % The axis (plot) box itself is the sized element (3.75 x 3.6 cm) --
    % the PDF page is padded larger around it (padding_cm each side) so
    % the axis box is guaranteed exactly this size regardless of how much
    % room MATLAB's auto layout would otherwise reserve for tick labels.
    % NOT fitFig2Page, which scales the figure to fill a full
    % letter-landscape page and would override any custom size set here.
    pairwisePlotWidth_cm = targetPairwisePlotWidth_cm; % 3cm, achieved via pairwiseXlim above rather than by scaling the axis box directly
    pairwisePlotHeight_cm = 3.6 * 0.94; % 6% smaller
    basePadding_cm = 0.5; % previous padding, on each side
    areaScale = 1.25; % total figure area, not just padding, should grow by this factor
    % Solve for the (still-uniform, all 4 sides) padding p that gives
    % (axisW+2p)*(axisH+2p) = areaScale * (axisW+2*basePadding)*(axisH+2*basePadding),
    % i.e. a quadratic in p: 4p^2 + 2p(axisW+axisH) + axisW*axisH - targetArea = 0
    targetArea = areaScale * (pairwisePlotWidth_cm + 2*basePadding_cm) * (pairwisePlotHeight_cm + 2*basePadding_cm);
    qa = 4; qb = 2*(pairwisePlotWidth_cm + pairwisePlotHeight_cm); qc = pairwisePlotWidth_cm*pairwisePlotHeight_cm - targetArea;
    padding_cm = (-qb + sqrt(qb^2 - 4*qa*qc)) / (2*qa);
    figWidth_cm = pairwisePlotWidth_cm + 2*padding_cm;
    figHeight_cm = pairwisePlotHeight_cm + 2*padding_cm;

    gcf_edit = gcf;
    gcf_edit.Units = 'centimeters';
    gcf_edit.Position(3:4) = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperUnits = 'centimeters';
    gcf_edit.PaperPositionMode = 'manual'; % otherwise 'auto' ignores PaperPosition below
    gcf_edit.PaperSize = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperPosition = [0, 0, figWidth_cm, figHeight_cm];

    ax.Units = 'centimeters';
    ax.Position = [padding_cm, padding_cm, pairwisePlotWidth_cm, pairwisePlotHeight_cm];

    ylim([-.25 0.75]) % if zero-meaning the data %%%%
    % Save as PDF
    set(gcf_edit,'Renderer','painters'); % new
    print(gcf_edit, filename, '-dpdf', '-painters');
    close all;

end
