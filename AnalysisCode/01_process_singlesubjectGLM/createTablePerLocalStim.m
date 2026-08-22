clc; clear all; close all

%% create new table with columns per local direction / orientation
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
savedir = fullfile(bidsDir, 'derivatives', 'summaryTables');
load(fullfile(savedir, 'allsubjectsTable.mat'));

% new table, reorganized
% lookup

[stimNames, stimMap] = createStimMap;

% print out mapping

angleBins = 0:45:315;

for iStim = 1:numel(stimNames)

    stim = stimNames{iStim};
    info = stimMap.(stim);

    fprintf('\n%s\n', stim);
    fprintf('-----------------------------\n');

    for angle = angleBins

        if isfield(info,'blank') && info.blank

            fprintf('pRF %3d -> blank\n', angle);

        elseif isfield(info,'motdir') && isempty(info.motdir)

            % stationary
            if strcmp(info.type,'cart')
                ori = info.ori;
            else
                ori = mod(angle + info.oriOffset,180);
            end

            fprintf('pRF %3d -> ori %3d, mot none\n', ...
                angle, ori);

        else

            % moving
            if strcmp(info.type,'cart')

                ori = info.ori;
                mot = info.motdir;

            else

                ori = mod(angle + info.oriOffset,180);
                mot = mod(angle + info.motOffset,360);

            end

            fprintf('pRF %3d -> ori %3d, mot %3d\n', ...
                angle, ori, mot);

        end
    end
end


%% Create all possible output condition names

oriBins = [0 45 90 135];

conditions = {};

% Define orthogonal motion directions for each orientation
orthogonalMotions = struct();
orthogonalMotions.ori0   = [90 270];
orthogonalMotions.ori45  = [135 315];
orthogonalMotions.ori90  = [0 180];
orthogonalMotions.ori135 = [45 225];

% Moving
for ori = oriBins

    motBins = orthogonalMotions.(sprintf('ori%d',ori));

    for mot = motBins

        conditions{end+1} = sprintf(...
            'cart_localori_%d_localmotdir_%d', ori, mot);

        conditions{end+1} = sprintf(...
            'pol_localori_%d_localmotdir_%d', ori, mot);

    end
end

% Stationary
for ori = oriBins

    conditions{end+1} = sprintf(...
        'cart_localori_%d_localmotdir_none', ori);

    conditions{end+1} = sprintf(...
        'pol_localori_%d_localmotdir_none', ori);

end

% Blank
conditions{end+1} = 'cart_blank';
conditions{end+1} = 'pol_blank';

conditions = string(conditions);
conditions = matlab.lang.makeValidName(conditions);

T = allsubjectsTable;

localData = nan(height(T), numel(conditions));

betaData = T{:,stimNames};   % rows = voxels, cols = stimuli
angles = T.pRF_angle_bin;

%% Unique pRF angle bins

angleBins = 0:45:315;

angleLookup = struct();

for iAngle = 1:numel(angleBins)

    angle = angleBins(iAngle);

    angleLookup(iAngle).angle = angle;

    for iStim = 1:numel(stimNames)

        stim = stimNames{iStim};
        info = stimMap.(stim);

        if isfield(info,'blank')

            colName = sprintf('%s_blank',info.type);

        elseif isfield(info,'motdir') && isempty(info.motdir)

            % stationary

            if strcmp(info.type,'cart')
                ori = info.ori;
            else
                ori = mod(angle + info.oriOffset,180);
            end

            colName = sprintf('%s_localori_%d_localmotdir_none',...
                info.type,ori);

        else

            % moving

            if strcmp(info.type,'cart')

                ori = info.ori;
                mot = info.motdir;

            else

                ori = mod(angle + info.oriOffset,180);
                mot = mod(angle + info.motOffset,360);

            end

            colName = sprintf('%s_localori_%d_localmotdir_%d',...
                info.type,ori,mot);

        end

        angleLookup(iAngle).stim(iStim).column = ...
            find(strcmp(conditions,matlab.lang.makeValidName(colName)));

    end
end

%% Apply mapping

for iAngle = 1:numel(angleBins)

    rows = angles == angleBins(iAngle);

    for iStim = 1:numel(stimNames)

        col = angleLookup(iAngle).stim(iStim).column;

        localData(rows,col) = betaData(rows,iStim);

    end
end

%% Make final table, organized by local absolute motiondir/orientation (instead of stimulus)

localTable = array2table(localData,...
    'VariableNames',conditions);

finalTable = [T(:,1:8), T(:,{'dg_beta_mean','dg_beta_std','da_beta_mean','da_beta_std'}), localTable];


%%

% compute each asymmetry:

finalTable.cart_oriasym_HvV = finalTable.cart_localori_0_localmotdir_none - finalTable.cart_localori_90_localmotdir_none;
finalTable.pol_oriasym_HvV  = finalTable.pol_localori_0_localmotdir_none  - finalTable.pol_localori_90_localmotdir_none;

finalTable.cart_oriasym_CvO = mean([finalTable.cart_localori_0_localmotdir_none finalTable.cart_localori_90_localmotdir_none],2) - ...
                              mean([finalTable.cart_localori_45_localmotdir_none finalTable.cart_localori_135_localmotdir_none],2);

finalTable.pol_oriasym_CvO  = mean([finalTable.pol_localori_0_localmotdir_none finalTable.pol_localori_90_localmotdir_none],2) - ...
                              mean([finalTable.pol_localori_45_localmotdir_none finalTable.pol_localori_135_localmotdir_none],2);


a = finalTable.pRF_angle_bin;

% Map pRF angle bin -> radial orientation, tangential orientation
radOri = nan(size(a));
tanOri = nan(size(a));

radOri(ismember(a,[0 180]))   = 0;
tanOri(ismember(a,[0 180]))   = 90;

radOri(ismember(a,[90 270]))  = 90;
tanOri(ismember(a,[90 270]))  = 0;

radOri(ismember(a,[45 225]))  = 45;
tanOri(ismember(a,[45 225]))  = 135;

radOri(ismember(a,[135 315])) = 135;
tanOri(ismember(a,[135 315])) = 45;


% Initialize output
cartRad = nan(height(finalTable),1);
cartTan = nan(height(finalTable),1);

polRad = nan(height(finalTable),1);
polTan = nan(height(finalTable),1);

% Fill according to orientation
for ori = [0 45 90 135]

    idxRad = radOri == ori;
    idxTan = tanOri == ori;

    cartRad(idxRad) = finalTable.(sprintf('cart_localori_%d_localmotdir_none',ori))(idxRad);
    cartTan(idxTan) = finalTable.(sprintf('cart_localori_%d_localmotdir_none',ori))(idxTan);

    polRad(idxRad) = finalTable.(sprintf('pol_localori_%d_localmotdir_none',ori))(idxRad);
    polTan(idxTan) = finalTable.(sprintf('pol_localori_%d_localmotdir_none',ori))(idxTan);

end

% Radial minus tangential
finalTable.cart_oriasym_RvT = cartRad - cartTan;
finalTable.pol_oriasym_RvT = polRad - polTan;


isCardinal = ismember(finalTable.pRF_angle_bin,[0 90 180 270]);
isOblique  = ismember(finalTable.pRF_angle_bin,[45 135 225 315]);

cartPC = nan(height(finalTable),1);
cartPO = nan(height(finalTable),1);
polPC  = nan(height(finalTable),1);
polPO  = nan(height(finalTable),1);

% Cardinal pRF locations
cartPC(isCardinal) = mean([ ...
    finalTable.cart_localori_0_localmotdir_none(isCardinal), ...
    finalTable.cart_localori_90_localmotdir_none(isCardinal)],2);

cartPO(isCardinal) = mean([ ...
    finalTable.cart_localori_45_localmotdir_none(isCardinal), ...
    finalTable.cart_localori_135_localmotdir_none(isCardinal)],2);

polPC(isCardinal) = mean([ ...
    finalTable.pol_localori_0_localmotdir_none(isCardinal), ...
    finalTable.pol_localori_90_localmotdir_none(isCardinal)],2);

polPO(isCardinal) = mean([ ...
    finalTable.pol_localori_45_localmotdir_none(isCardinal), ...
    finalTable.pol_localori_135_localmotdir_none(isCardinal)],2);


% Oblique pRF locations
cartPC(isOblique) = mean([ ...
    finalTable.cart_localori_45_localmotdir_none(isOblique), ...
    finalTable.cart_localori_135_localmotdir_none(isOblique)],2);

cartPO(isOblique) = mean([ ...
    finalTable.cart_localori_0_localmotdir_none(isOblique), ...
    finalTable.cart_localori_90_localmotdir_none(isOblique)],2);

polPC(isOblique) = mean([ ...
    finalTable.pol_localori_45_localmotdir_none(isOblique), ...
    finalTable.pol_localori_135_localmotdir_none(isOblique)],2);

polPO(isOblique) = mean([ ...
    finalTable.pol_localori_0_localmotdir_none(isOblique), ...
    finalTable.pol_localori_90_localmotdir_none(isOblique)],2);


finalTable.cart_oriasym_PCvPO = cartPC - cartPO;
finalTable.pol_oriasym_PCvPO = polPC - polPO;

%% Compute motion vs static per experiment

% Get column names
vars = finalTable.Properties.VariableNames;


% Cartesian motion vs static
% 8 motion-direction columns
cartMotCols = vars(startsWith(vars,'cart_') & contains(vars,'localmotdir_') & ~contains(vars,'localmotdir_none'));

% 4 static columns
cartStatCols = vars(startsWith(vars,'cart_') & contains(vars,'localmotdir_none'));

finalTable.cart_MotvStat = ...
    (mean(finalTable{:,cartMotCols},2,'omitnan') - ...
    mean(finalTable{:,cartStatCols},2,'omitnan')) ./ ...
    (mean(finalTable{:,cartMotCols},2,'omitnan') + ...
    mean(finalTable{:,cartStatCols},2,'omitnan'));

% Polar motion vs static
% 8 motion-direction columns
polMotCols = vars(startsWith(vars,'pol_') & contains(vars,'localmotdir_') & ~contains(vars,'localmotdir_none'));

% 4 static columns
polStatCols = vars(startsWith(vars,'pol_') & contains(vars,'localmotdir_none'));

finalTable.pol_MotvStat = ...
    (mean(finalTable{:,polMotCols},2,'omitnan') - ...
    mean(finalTable{:,polStatCols},2,'omitnan')) ./ ...
    (mean(finalTable{:,polMotCols},2,'omitnan') + ...
    mean(finalTable{:,polStatCols},2,'omitnan'));


%% Replot polar plot with horizontal vs vertical

% Included V1 voxels
idx = strcmp(finalTable.visual_area,'V1') & finalTable.included==1;
T = finalTable(idx,:);

subjects = unique(T.subject);
angleBins = 0:45:315;

% Variables to compute
vars = { ...
    'cart_localori_0_localmotdir_none', ...
    'cart_localori_90_localmotdir_none', ...
    'pol_localori_0_localmotdir_none', ...
    'pol_localori_90_localmotdir_none'};

groupMean = cell(size(vars));

for v = 1:numel(vars)

    subjMedian = nan(numel(subjects), numel(angleBins));

    for s = 1:numel(subjects)

        subjIdx = strcmp(T.subject, subjects{s});

        for b = 1:numel(angleBins)

            binIdx = subjIdx & T.pRF_angle_bin == angleBins(b);

            subjMedian(s,b) = median( ...
                T.(vars{v})(binIdx), ...
                'omitnan');

        end
    end

    groupMean{v} = mean(subjMedian,1,'omitnan');

end

theta = deg2rad([angleBins angleBins(1)]);

figure

% Cartesian
subplot(1,2,1)

polarplot(theta,[groupMean{1} groupMean{1}(1)],'-o','LineWidth',2, 'Color', [0, 0, 0.5])
hold on
polarplot(theta,[groupMean{2} groupMean{2}(1)],'-o','LineWidth',2, 'Color', [1, 0.5, 0])
rlim([-1 1])
pax = gca;
pax.ThetaZeroLocation = 'right';
pax.ThetaDir = 'counterclockwise';

title('Cartesian')
legend({'Local ori 0','Local ori 90'},'Location','best')

% Polar
subplot(1,2,2)

polarplot(theta,[groupMean{3} groupMean{3}(1)],'-o','LineWidth',2, 'Color', [0, 0, 0.5])
hold on
polarplot(theta,[groupMean{4} groupMean{4}(1)],'-o','LineWidth',2, 'Color', [1, 0.5, 0])
rlim([-1 1])
pax = gca;
pax.ThetaZeroLocation = 'right';
pax.ThetaDir = 'counterclockwise';

title('Polar')
legend({'Local ori 0','Local ori 90'},'Location','best')

%% Replot polar plot with radial vs tangential

% Compute radial and tangential responses per voxel

n = height(finalTable);

finalTable.cart_radial = nan(n,1);
finalTable.cart_tangential = nan(n,1);
finalTable.pol_radial = nan(n,1);
finalTable.pol_tangential = nan(n,1);

a = finalTable.pRF_angle_bin;

% Horizontal meridian
idx = ismember(a,[0 180]);
finalTable.cart_radial(idx)      = finalTable.cart_localori_0_localmotdir_none(idx);
finalTable.cart_tangential(idx)  = finalTable.cart_localori_90_localmotdir_none(idx);
finalTable.pol_radial(idx)       = finalTable.pol_localori_0_localmotdir_none(idx);
finalTable.pol_tangential(idx)   = finalTable.pol_localori_90_localmotdir_none(idx);

% Vertical meridian
idx = ismember(a,[90 270]);
finalTable.cart_radial(idx)      = finalTable.cart_localori_90_localmotdir_none(idx);
finalTable.cart_tangential(idx)  = finalTable.cart_localori_0_localmotdir_none(idx);
finalTable.pol_radial(idx)       = finalTable.pol_localori_90_localmotdir_none(idx);
finalTable.pol_tangential(idx)   = finalTable.pol_localori_0_localmotdir_none(idx);

% 45° diagonal
idx = ismember(a,[45 225]);
finalTable.cart_radial(idx)      = finalTable.cart_localori_45_localmotdir_none(idx);
finalTable.cart_tangential(idx)  = finalTable.cart_localori_135_localmotdir_none(idx);
finalTable.pol_radial(idx)       = finalTable.pol_localori_45_localmotdir_none(idx);
finalTable.pol_tangential(idx)   = finalTable.pol_localori_135_localmotdir_none(idx);

% 135° diagonal
idx = ismember(a,[135 315]);
finalTable.cart_radial(idx)      = finalTable.cart_localori_135_localmotdir_none(idx);
finalTable.cart_tangential(idx)  = finalTable.cart_localori_45_localmotdir_none(idx);
finalTable.pol_radial(idx)       = finalTable.pol_localori_135_localmotdir_none(idx);
finalTable.pol_tangential(idx)   = finalTable.pol_localori_45_localmotdir_none(idx);

% Subject-averaged radial and tangential responses

idx = strcmp(finalTable.visual_area,'V1') & finalTable.included==1;
T = finalTable(idx,:);

subjects = unique(T.subject);
angleBins = 0:45:315;

vars = { ...
    'cart_radial', ...
    'cart_tangential', ...
    'pol_radial', ...
    'pol_tangential'};

groupMean = cell(size(vars));

for v = 1:numel(vars)

    subjMedian = nan(numel(subjects), numel(angleBins));

    for s = 1:numel(subjects)

        subjIdx = strcmp(T.subject,subjects{s});

        for b = 1:numel(angleBins)

            binIdx = subjIdx & T.pRF_angle_bin==angleBins(b);

            subjMedian(s,b) = median(T.(vars{v})(binIdx),'omitnan');

        end
    end

    groupMean{v} = mean(subjMedian,1,'omitnan');

end

theta = deg2rad([angleBins angleBins(1)]);

figure

% Cartesian
subplot(1,2,1)
polarplot(theta,[groupMean{1} groupMean{1}(1)],'-o','LineWidth',2,'Color',[0.5725, 0.7725, 0.8706])
hold on
polarplot(theta,[groupMean{2} groupMean{2}(1)],'-o','LineWidth',2,'Color',[0.7922, 0, 0.1255])
rlim([-1 1])
pax = gca;
pax.ThetaZeroLocation = 'right';
pax.ThetaDir = 'counterclockwise';
title('Cartesian')
legend({'Radial','Tangential'},'Location','best')

% Polar
subplot(1,2,2)
polarplot(theta,[groupMean{3} groupMean{3}(1)],'-o','LineWidth',2,'Color',[0.5725, 0.7725, 0.8706])
hold on
polarplot(theta,[groupMean{4} groupMean{4}(1)],'-o','LineWidth',2,'Color',[0.7922, 0, 0.1255])
rlim([-1 1])
pax = gca;
pax.ThetaZeroLocation = 'right';
pax.ThetaDir = 'counterclockwise';
title('Polar')
legend({'Radial','Tangential'},'Location','best')


%%

plotAsymmetry(finalTable,...
    'cart_oriasym_HvV',...
    'pol_oriasym_HvV',...
    'H - V orientation asymmetry')


plotAsymmetry(finalTable,...
    'cart_oriasym_CvO',...
    'pol_oriasym_CvO',...
    'C - O orientation asymmetry')


plotAsymmetry(finalTable,...
    'cart_oriasym_RvT',...
    'pol_oriasym_RvT',...
    'R - T orientation asymmetry')


plotAsymmetry(finalTable,...
    'cart_oriasym_PCvPO',...
    'pol_oriasym_PCvPO',...
    'PC - PO orientation asymmetry')



%% Plot by polar angle the mean and std

% Select included V1 voxels
idx = strcmp(finalTable.visual_area,'V1') & finalTable.included == 1;

angleBins = 0:45:315;
thetaMed = deg2rad(angleBins);

figure

% DG mean
subplot(2,2,1)

theta = deg2rad(finalTable.pRF_angle_bin(idx));
rho = finalTable.dg_beta_mean(idx);

% polarscatter(theta, rho, 10, 'filled', ...
%     'MarkerFaceAlpha',0.3,...
%     'MarkerEdgeAlpha',0.3)
% rlim([-2 2])
% hold on

rhoMed = arrayfun(@(a) ...
    median(rho(finalTable.pRF_angle_bin(idx)==a),'omitnan'), ...
    angleBins);

theta = linspace(0, 2*pi, 361);
polarplot(theta, zeros(size(theta)), 'r-', 'LineWidth', 2)
hold on
polarscatter(thetaMed, rhoMed, 120, 'k+', 'LineWidth',2)
rlim([-.5 .5])
hold off
title('Cartesian experiment: beta mean')


% DA mean
subplot(2,2,2)

rho = finalTable.da_beta_mean(idx);

% polarscatter(theta, rho, 10, 'filled', ...
%     'MarkerFaceAlpha',0.3,...
%     'MarkerEdgeAlpha',0.3)
% rlim([-2 2])
% hold on

rhoMed = arrayfun(@(a) ...
    median(rho(finalTable.pRF_angle_bin(idx)==a),'omitnan'), ...
    angleBins);


polarplot(theta, zeros(size(theta)), 'r-', 'LineWidth', 2)
hold on
polarscatter(thetaMed, rhoMed, 120, 'k+', 'LineWidth',2)
rlim([-.5 .5])
hold off
title('Polar experiment: beta mean')


% DG std
subplot(2,2,3)

rho = finalTable.dg_beta_std(idx);

% polarscatter(theta, rho, 10, 'filled', ...
%     'MarkerFaceAlpha',0.3,...
%     'MarkerEdgeAlpha',0.3)
% 
% hold on

rhoMed = arrayfun(@(a) ...
    median(rho(finalTable.pRF_angle_bin(idx)==a),'omitnan'), ...
    angleBins);

polarscatter(thetaMed, rhoMed, 120, 'k+', 'LineWidth',2)
rlim([0 1])
hold off
title('Cartesian experiment: beta std')


% DA std
subplot(2,2,4)

rho = finalTable.da_beta_std(idx);

% polarscatter(theta, rho, 10, 'filled', ...
%     'MarkerFaceAlpha',0.3,...
%     'MarkerEdgeAlpha',0.3)
% 
% hold on

rhoMed = arrayfun(@(a) ...
    median(rho(finalTable.pRF_angle_bin(idx)==a),'omitnan'), ...
    angleBins);

polarscatter(thetaMed, rhoMed, 120, 'k+', 'LineWidth',2)
rlim([0 1])
hold off
title('Polar experiment: beta std')
f1 = gcf;
f1.Position = [1072 289 1207 941];

%% Save table
filename = fullfile(savedir, 'allsubjectsTable_local_AbsOriDir.mat');
save(filename, 'finalTable');

filename = fullfile(savedir, 'allsubjectsTable_local_AbsOriDir.csv');
writetable(finalTable, filename);


%% For V1, included voxels:
% average each subjects' voxels per unique pRF_bin, for horizontal orientations only
% then do vertical only


%%

function plotAsymmetryCombinedSubjects(finalTable, asymColCart, asymColPol, titleText)

idx = strcmp(finalTable.visual_area,'V1') & finalTable.included == 1;

figure;

plots = {
    finalTable.dg_beta_mean(idx), finalTable.(asymColCart)(idx), 'beta mean', 'Cartesian'
    finalTable.da_beta_mean(idx), finalTable.(asymColPol)(idx),  'beta mean', 'Polar'
    finalTable.dg_beta_std(idx),  finalTable.(asymColCart)(idx), 'beta std',  'Cartesian'
    finalTable.da_beta_std(idx),  finalTable.(asymColPol)(idx),  'beta std',  'Polar'
    };

for i = 1:4

    subplot(2,2,i)

    x = plots{i,1};
    y = plots{i,2};


    if i <= 2
        xlim([-3 3])
        ylim([-3 3])
        xline(0,'w-','LineWidth',1)
        xedges = linspace(-3, 3, 81);
        yedges = linspace(-3, 3, 81);
        hold on
    else
        xlim([0 3])
        ylim([-3 3])
        xedges = linspace(-3, 3, 81);
        yedges = linspace(-3, 3, 81);
    end

    yline(0,'w-','LineWidth',1)
    hold on

    %histogram2(x, y, xedges, yedges, 'DisplayStyle', 'tile', 'ShowEmptyBins','on'); 
    %set(gca, "ZScale", 'linear')

    h=scatter(x,y,10,'filled');
    h.MarkerFaceAlpha = 0.3;
    h.MarkerEdgeAlpha = 0.3;

%     plot(median(x,'omitnan'), median(y,'omitnan'), ...
%         'k+', 'MarkerSize',16,'LineWidth',2)

    hold off

    xlabel(plots{i,3})
    ylabel(plots{i,4} + " " + titleText)


end

sgtitle(strcat(titleText, ' per V1 voxel (ecc 4-8, R^2>.1)'))

end


%%


function plotAsymmetryPersubject(finalTable, asymColCart, asymColPol, titleText)

betaType = 'std';

% Included V1 voxels
idx = strcmp(finalTable.visual_area,'V1') & finalTable.included == 1;

T = finalTable(idx,:);

subjects = unique(T.subject,'stable');
nSubjects = numel(subjects);

figure;

% Choose beta variable
switch betaType
    case 'mean'
        dgBeta = 'dg_beta_mean';
        daBeta = 'da_beta_mean';
        xlimVals = [-3 3];
    case 'std'
        dgBeta = 'dg_beta_std';
        daBeta = 'da_beta_std';
        xlimVals = [0 3];
end

% Histogram bins
xedges = linspace(xlimVals(1), xlimVals(2), 41);
yedges = linspace(-3,3,41);

for s = 1:nSubjects

    subjIdx = strcmp(T.subject, subjects{s});

    % DG row
    subplot(2,nSubjects,s)

    x = T.(dgBeta)(subjIdx);
    y = T.(asymColCart)(subjIdx);

    histogram2(x,y,xedges,yedges,...
        'DisplayStyle','tile',...
        'ShowEmptyBins','on');

    hold on
    xline(0,'w-','LineWidth',1,'HandleVisibility','off')
    yline(0,'w-','LineWidth',1,'HandleVisibility','off')
    hold off

    xlim(xlimVals)
    ylim([-3 3])

    title(subjects{s})
    
    if s == 1
        ylabel('DG asymmetry')
    end


    % DA row
    subplot(2,nSubjects,nSubjects+s)

    x = T.(daBeta)(subjIdx);
    y = T.(asymColPol)(subjIdx);

    histogram2(x,y,xedges,yedges,...
        'DisplayStyle','tile',...
        'ShowEmptyBins','on');

    hold on
    xline(0,'w-','LineWidth',1,'HandleVisibility','off')
    yline(0,'w-','LineWidth',1,'HandleVisibility','off')
    hold off

    xlim(xlimVals)
    ylim([-3 3])

    if s == 1
        ylabel('DA asymmetry')
    end

end

% Same color scale across all plots
clim([0 20])   % adjust depending on max counts

sgtitle(sprintf('%s %s per subject', titleText, betaType))

end

%%

% function plotAsymmetry(finalTable, asymColCart, asymColPol, titleText)
% 
% betaType = 'std';
% 
% % Included V1 voxels
% idx = strcmp(finalTable.visual_area,'V1') & finalTable.included == 1;
% 
% T = finalTable(idx,:);
% 
% subjects = unique(T.subject,'stable');
% nSubjects = numel(subjects);
% 
% % Choose beta variable
% switch betaType
%     case 'mean'
%         dgBeta = 'dg_beta_mean';
%         daBeta = 'da_beta_mean';
%         xlimVals = [-3 3];
%     case 'std'
%         dgBeta = 'dg_beta_std';
%         daBeta = 'da_beta_std';
%         xlimVals = [0 3];
% end
% 
% % pRF angle colors
% angleBins = 0:45:315;
% colors = hsv(numel(angleBins));
% 
% figure;
% 
% for s = 1:nSubjects
% 
%     subjIdx = strcmp(T.subject, subjects{s});
% 
%     % DG row
%     subplot(2,nSubjects,s)
% 
%     x = T.(dgBeta)(subjIdx);
%     y = T.(asymColCart)(subjIdx);
%     angle = T.pRF_angle_bin(subjIdx);
% 
%     hold on
% 
%     for a = 1:numel(angleBins)
% 
%         aIdx = angle == angleBins(a);
% 
%         scatter(x(aIdx), y(aIdx), ...
%             10, colors(a,:), 'filled', ...
%             'MarkerFaceAlpha',0.3,...
%             'MarkerEdgeAlpha',0.3);
% 
%     end
% 
%     xline(0,'k-','LineWidth',1,'HandleVisibility','off')
%     yline(0,'k-','LineWidth',1,'HandleVisibility','off')
% 
%     hold off
% 
%     xlim(xlimVals)
%     ylim([-3 3])
% 
%     title(subjects{s})
% 
%     if s == 1
%         ylabel('DG asymmetry')
%     end
% 
% 
%     % DA row
%     subplot(2,nSubjects,nSubjects+s)
% 
%     x = T.(daBeta)(subjIdx);
%     y = T.(asymColPol)(subjIdx);
%     angle = T.pRF_angle_bin(subjIdx);
% 
%     hold on
% 
%     for a = 1:numel(angleBins)
% 
%         aIdx = angle == angleBins(a);
% 
%         scatter(x(aIdx), y(aIdx), ...
%             10, colors(a,:), 'filled', ...
%             'MarkerFaceAlpha',0.3,...
%             'MarkerEdgeAlpha',0.3);
% 
%     end
% 
%     xline(0,'k-','LineWidth',1,'HandleVisibility','off')
%     yline(0,'k-','LineWidth',1,'HandleVisibility','off')
% 
%     hold off
% 
%     xlim(xlimVals)
%     ylim([-3 3])
% 
%     if s == 1
%         ylabel('DA asymmetry')
%     end
% 
% end
% 
% 
% % Add one shared legend
% figure(gcf)
% hold on
% h = gobjects(numel(angleBins),1);
% 
% for a = 1:numel(angleBins)
%     h(a) = scatter(nan,nan,30,colors(a,:),'filled');
% end
% 
% legend(h,string(angleBins)+"°",...
%     'Location','eastoutside')
% 
% sgtitle(sprintf('%s %s per subject (colored by pRF angle)', ...
%     titleText,betaType))
% 
% end

%%

function plotAsymmetry(finalTable, asymColCart, asymColPol, titleText)

betaType = 'std';

% Included V1 voxels
idx = strcmp(finalTable.visual_area,'V1') & finalTable.included == 1;

T = finalTable(idx,:);

subjects = unique(T.subject,'stable');
nSubjects = numel(subjects);

% Choose beta variable
switch betaType
    case 'mean'
        dgBeta = 'dg_beta_mean';
        daBeta = 'da_beta_mean';
        xlimVals = [-3 3];
    case 'std'
        dgBeta = 'dg_beta_std';
        daBeta = 'da_beta_std';
        xlimVals = [0 4];
end

figure;

for s = 1:nSubjects

    subjIdx = strcmp(T.subject, subjects{s});


    % DG row (Cartesian)
    subplot(2,nSubjects,s)

    x = T.(dgBeta)(subjIdx);
    y = T.(asymColCart)(subjIdx);
    motStat = T.cart_MotvStat(subjIdx);

    h = plotMotHighlight(x,y,motStat);
    legend({'Remaining voxels','',''})
    

    xline(0,'k-','LineWidth',1,'HandleVisibility','off')
    yline(0,'k-','LineWidth',1,'HandleVisibility','off')

    xlim(xlimVals)
    ylim([-3 3])

    title(subjects{s})

    if s == 1
        legend(h, {'Remaining voxels','Bottom 20% motion-selective','Top 20% motion-selective'}, ...
            'Location','best');
        ylabel('DG asymmetry')
    else
        legend('off');
    end


    % DA row (Polar)
    subplot(2,nSubjects,nSubjects+s)

    x = T.(daBeta)(subjIdx);
    y = T.(asymColPol)(subjIdx);
    motStat = T.pol_MotvStat(subjIdx);

    plotMotHighlight(x,y,motStat);

    xline(0,'k-','LineWidth',1,'HandleVisibility','off')
    yline(0,'k-','LineWidth',1,'HandleVisibility','off')

    xlim(xlimVals)
    ylim([-3 3])

    if s == 1
        ylabel('DA asymmetry')
    end

end

sgtitle(sprintf('%s %s per subject (highlighting motion sensitivity)', ...
    titleText,betaType))
f = gcf;

f.Position = [32 274 2153 719];
end


function h = plotMotHighlight(x,y,motStat)

% percentile thresholds
lowThresh = prctile(motStat,20);
highThresh = prctile(motStat,80);

lowIdx = motStat <= lowThresh;
highIdx = motStat >= highThresh;
midIdx = ~lowIdx & ~highIdx;

hold on

% Middle 60%
h(1) = scatter(x(midIdx),y(midIdx),10,[0.3 0.3 0.3],'filled',...
    'MarkerFaceAlpha',0.3,...
    'MarkerEdgeAlpha',0.3);

% Lowest 20%
h(2) = scatter(x(lowIdx),y(lowIdx),10,[0 0.45 1],'filled',...
    'MarkerFaceAlpha',0.5,...
    'MarkerEdgeAlpha',0.5);

% Highest 20%
h(3) = scatter(x(highIdx),y(highIdx),10,[1 0 0],'filled',...
    'MarkerFaceAlpha',0.5,...
    'MarkerEdgeAlpha',0.5);

hold off

end

