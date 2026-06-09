%ttave_combinefits.m

clear all; close all; clc;

% must be in DriftingGrating directory to run

% setup path
addpath(genpath(pwd));
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
figDir = strrep(bidsDir, 'data_bids', 'figures');
%bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.2.0';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir)

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois;
clear projectSettings;
comparisonName = 'motion_minus_baseline';

projects = {'dg', 'da'};
hRF_setting = 'glmsingle'; %
asymmetries = {'cartesianCardinalvsCartesianOblique', 'polarCardinalvsPolarOblique', ...
    'verticalVsHorizontal','radialVsTangential'};

% for bookkeeping
fitted_allStaticVals_array = []; raw_allStaticVals_array = [];

% if strcmp(projectName, 'dg')
%     % to leave out one subject
%     subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-0397', ...
%         'sub-0442', 'sub-wlsubj121', ...
%         'sub-wlsubj123', 'sub-wlsubj124', 'sub-wlsubj127', 'sub-0395', 'sub-0426', ...
%         'sub-0250'};
% elseif strcmp(projectName, 'da')
    % load subjects
subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250'}; % , 
% end

% make a three 3D arrays (1 for each asymmetry)
% rows (4x) 
% cols (x25) timepoints

nCurves = 4; % conditions (motion pro, motion con, orientation pro, orientation con)
nSubjects = numel(subjects);
nAsymmetries = numel(asymmetries);
nTimepoints = 25;

%% Aggregate

for ri=1:numel(rois)  % specific to the ROI

    roiName = rois{ri};

    for si=1:numel(subjects)
    
        subj = subjects{si};
    
        ttave_combinedPath = fullfile(bidsDir, 'derivatives', 'ttaveCombinedfits', subj);
    
        % this overrides 'projectSettings'
        % save the variables: 'ttave_arrays', 'projectSettings', 'subj'
        load(fullfile(ttave_combinedPath,sprintf('ttaveCombined_%s',comparisonName)));

        for pi=1:numel(projects)    % will take the mean over projects
        
            projectName = projects{pi};
    
            for as=1:numel(asymmetries)

                asymmetryName = asymmetries{as};
    
                % initialize group array if does not exist
                try
                 retrived_fits = group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits'));
                 retrived_rawvals = group_ttave_arrays.(roiName).(asymmetryName).(projectName);
                catch
                 group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits')) = nan(nCurves, nTimepoints, nSubjects);
                 
                 group_ttave_arrays.(roiName).(asymmetryName).(projectName) = nan(nCurves, nTimepoints, nSubjects);

                 % intercept, beta (scale)
                 group_ttave_metrics.(roiName).(asymmetryName).(projectName) = nan(nCurves, 2, nSubjects);

                end

                metrics = ttave_metrics.(roiName).(asymmetryName).(projectName);
                group_ttave_metrics.(roiName).(asymmetryName).(projectName)(:,:,si) = metrics;

                fits = ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits'));
                group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits'))(:,:,si) = fits;

                rawvals = ttave_arrays.(roiName).(asymmetryName).(projectName);
                group_ttave_arrays.(roiName).(asymmetryName).(projectName)(:,:,si) = rawvals;

                group_ttave_arrays.(roiName).(asymmetryName).colors = ttave_arrays.(roiName).(asymmetryName).colors;
                group_ttave_arrays.(roiName).(asymmetryName).legendLabels = ttave_arrays.(roiName).(asymmetryName).legendLabels;

            end
        
        end
    end
end

% example of querying stationary condition (3:4), beta fits (scale) 
mean(group_ttave_metrics.(roiName).('cartesianCardinalvsCartesianOblique').('dg')(3:4,2,:), 'all')
mean(group_ttave_metrics.(roiName).('cartesianCardinalvsCartesianOblique').('da')(3:4,2,:), 'all')

%% Now average and plot

% asymmetries

ttave_combinedGroupPath = fullfile(figDir, 'combined');

if ~isfolder(ttave_combinedGroupPath)
    mkdir(ttave_combinedGroupPath)
end

pointer=1;
for ri=1:numel(rois)  % specific to the ROI

    roiName = rois{ri};

    iter=0;

    for ri=1:numel(asymmetries)

        figure
        asymmetryName = asymmetries{ri};

        colors = group_ttave_arrays.(roiName).(asymmetryName).colors;
        legendLabels = group_ttave_arrays.(roiName).(asymmetryName).legendLabels;

        for pi=1:numel(projects)    % will take the mean over projects

            projectName = projects{pi};
        
            iter=iter+1;
            figure(1)
            subplot(4,2,iter)

            % copied over for now; condense with ttave_fit LATER
            rawVals = mean(group_ttave_arrays.(roiName).(asymmetryName).(projectName),3);
            fittedVals = mean(group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits')),3);

            advMotionVals = rawVals(1,:);
            disadvMotionVals = rawVals(2,:);
            advStaticVals = rawVals(3,:);
            disadvStaticVals = rawVals(4,:);
    
            fitted_advMotionVals = fittedVals(1,:);
            fitted_disadvMotionVals = fittedVals(2,:);
            fitted_advStaticVals = fittedVals(3,:);
            fitted_disadvStaticVals = fittedVals(4,:);
    
            %figure(1)
            hold on
            %p1=plot(fitted_advMotionVals, 'b-', 'Linewidth',1, 'Color', colors{1});
            %p2=plot(fitted_disadvMotionVals, 'b-', 'Linewidth',1, 'Color', colors{2});
            p3=plot(fitted_advStaticVals, 'b-', 'Linewidth',1, 'Color', colors{1});
            p4=plot(fitted_disadvStaticVals, 'b-', 'Linewidth',1, 'Color', colors{2});

            if ri==1
                overallFits = [fitted_advStaticVals + fitted_disadvStaticVals]/2;
                fitted_allStaticVals_array = [fitted_allStaticVals_array; overallFits];
                overallVals = [advStaticVals + disadvStaticVals]/2;
                raw_allStaticVals_array = [raw_allStaticVals_array; overallVals];
            end

            %scatter(1:length(advMotionVals), advMotionVals, 25, 'filled', 'Marker', 'd', 'MarkerFaceColor', colors{1});
            %scatter(1:length(disadvMotionVals), disadvMotionVals, 25, 'filled', 'Marker', 'd', 'MarkerFaceColor', colors{2});
            scatter(1:length(advStaticVals), advStaticVals, 25, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{1}, 'MarkerFaceAlpha', 0.4, 'LineWidth',1);
            scatter(1:length(disadvStaticVals), disadvStaticVals, 25, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{2}, 'MarkerFaceAlpha', 0.4, 'LineWidth',1);
            hold off

            if ismember(iter, [1,3,5,7])
                %legend([p1,p2,p3,p4],legendLabels, 'Location', 'northeast', 'Fontsize', 6)
                legend([p3,p4],legendLabels(3:4), 'Location', 'northeast', 'Fontsize', 6)
                ylabel('BOLD (psc)')
            end

            if strcmp(projectName, 'dg')
                projectDescription = 'Cartesian grating stimulus';
            elseif strcmp(projectName, 'da')
                projectDescription = 'Polar grating stimulus';
            end

            if ismember(iter, [1,2])
                title(sprintf('%s', projectDescription))
            end

            if ismember(iter, [7,8])
                xlabel('TRs')
                current_ticks = xticks();
                new_labels = current_ticks - projectSettings.eventTRs_prior;
                xticklabels(string(new_labels));
            else
                xticklabels([]);
            end

            %plot(grandMean, '-', 'Linewidth',3, 'Color', 'k')
            ylim([-.5 1])
            set(gca, 'FontName', 'Arial', 'FontSize', 12)

            % for peak datapoints (4 x 8 array - 4 conditions; 8 subjects)
            fittedpeakVals = squeeze(max(group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits')),[],2));

            x = 1:size(fittedpeakVals, 1);
            %jitterAmount = .2; % Maximum jitter range
            %jitter = x + (rand(size(x)) - 0.5) * 2 * jitterAmount; % Add jitter
            
%             figure(2)
%             subplot(3,2,iter)
%             hold on;
%             for i = 1:size(fittedpeakVals, 1)
%                 if i<=2
%                     currentColor = colors{i};
%                     hold on
%                     scatter(repmat(x(i), 1, size(fittedpeakVals, 2)), fittedpeakVals(i, :), 100, 'filled', 'MarkerFaceColor', currentColor); % Plot each row
%                 else
%                     currentColor = colors{i-2};
%                     scatter(repmat(x(i), 1, size(fittedpeakVals, 2)), fittedpeakVals(i, :), 100, 'MarkerEdgeColor', currentColor); % Plot each row
%                 end
%             end
%             for i = 1:size(mot_minus_ori_peaks, 2) % Loop through each column (pair of values)
%                 plot([1,2], fittedpeakVals(1:2, i), 'k-', 'LineWidth', 1.5); % Connect the two points with a line
%                 hold on
%                 plot([3,4], fittedpeakVals(3:4, i), 'k-', 'LineWidth', 1.5); % Connect the two points with a line
%                 hold on
%             end
%             hold off;
%             
%             % Add labels
%             xlim([0 5])
%             xlabel('condition');
%             ylabel('peak');
%             xticks([])

%             figure(3)
%             if ~strcmp(asymmetryName, 'radialVsTangential') && ~strcmp(asymmetryName, 'verticalVsHorizontal') % RE added verticalVHor
%                 mot_minus_ori_peaks = fittedpeakVals(1:2,:)-fittedpeakVals(3:4,:);
%             else
%                 mot_minus_ori_peaks = fittedpeakVals(1:2,:)-flipud(fittedpeakVals(3:4,:));
%             end
%             x = 1:size(mot_minus_ori_peaks, 1);
%             %x_jittered = x + (rand(size(x)) - 0.5) * 2 * jitterAmount; % Add jitter
%             subplot(3,2,iter)
%             hold on;
%             for i = 1:size(mot_minus_ori_peaks, 1)
%                 scatter(repmat(x(i), 1, size(mot_minus_ori_peaks, 2)), mot_minus_ori_peaks(i, :), 100, 'filled', 'MarkerFaceColor', colors{i}); % Plot each row
%             end
%             for i = 1:size(mot_minus_ori_peaks, 2) % Loop through each column (pair of values)
%                 plot([1,2], mot_minus_ori_peaks(:, i), 'k-', 'LineWidth', 1.5); % Connect the two points with a line
%             end
%             hold off;
%             
%             % Add labels
%             xlim([0 3])
%             xlabel('condition');
%             ylabel('peak');
%             %title('Drifting Grating (left) vs Annulus (Right): Motion - Ori')
%             xticks([])
        end

    end


    sgtitle(roiName)
    fig = gcf;
    fig.Position = [1 306 1211 1031];
    gcf_edit = fitFig2Page(gcf);

    % Save as PDF
    print(gcf_edit, fullfile(ttave_combinedGroupPath,sprintf('ttaveCombined_%s_%s',roiName,projectSettings.comparisonName)), '-dpdf');

    %close all;

    figure
    subplot(1,2,1)
    p1=plot(fitted_allStaticVals_array(1,:), 'k-', 'Linewidth',1);
    hold on
    scatter(1:length(raw_allStaticVals_array), raw_allStaticVals_array(1,:), 25, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.4, 'LineWidth',1);
    xlabel('TRs')
    current_ticks = xticks();
    new_labels = current_ticks - projectSettings.eventTRs_prior;
    xticklabels(string(new_labels)); 
    ylim([-.5 1])
    ax = gca; 
    set(ax, 'FontName', 'Arial', 'FontSize', 12)  
    subplot(1,2,2)
    p2=plot(fitted_allStaticVals_array(2,:), 'k-', 'Linewidth',1);
    hold on
    scatter(1:length(raw_allStaticVals_array), raw_allStaticVals_array(2,:), 25, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.4, 'LineWidth',1);
    xlabel('TRs')
    current_ticks = xticks();
    new_labels = current_ticks - projectSettings.eventTRs_prior;
    xticklabels(string(new_labels)); 
    ylim([-.5 1])
    ax = gca; 
    set(ax, 'FontName', 'Arial', 'FontSize', 12)  
    %fig = gcf;
    %fig.Position = [2 1139 459 197];

    fig = gcf;
    fig.Position = [28 1119 1212 213];
    gcf_edit = fitFig2Page(gcf);
    print(gcf_edit, fullfile(ttave_combinedGroupPath,sprintf('ALLttaveCombined_%s_%s',roiName,projectSettings.comparisonName)), '-dpdf');

end
