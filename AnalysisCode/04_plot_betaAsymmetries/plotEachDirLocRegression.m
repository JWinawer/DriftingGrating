function plotEachDirLocRegression(projectName, roiname, varargin)
% PLOTEACHDIRLOCREGRESSION  Offshoot of lme2_ploteachDirLoc.m: the same
% 8-location "compass" grid of polar subplots (model prediction overlaid
% on empirical data, one subplot per location, spatially arranged so each
% subplot sits in its own direction from the center), but sourced from
% fitAsymmetryRegression.m's cached joint regression fit instead of
% LME_bold.mat/modeldata.mat, and with the empirical data gain- and
% precision-weighted across subjects (the original script's "meta-subject"
% averaging was an unweighted mean across all subjects' rows -- this
% version applies the same per-subject gain correction and per-(subject,
% cortical area) precision weighting used everywhere else in this
% pipeline). Fully independent of lme1_fit.m/lme2_ploteachDirLoc.m and
% their LME_results/ output.
%
%   plotEachDirLocRegression('dg', 'V1')
%   plotEachDirLocRegression('da', 'V1', 'precisionWeights', T)
%
% Model prediction (red line/dots): M * [grandInterceptFE; beta1..beta4],
% where beta_k = estimates(k)/2 (estimates is already the beta*2
% pro-minus-con scale) and grandInterceptFE is
% fitAsymmetryRegression.m's fixed-subject-intercept refit -- the analog
% of lme1_fit.m's fitlme-derived intercept, reconstructing the SAME kind
% of design-matrix-times-coefficients prediction the original script
% computed from the LME's estimates.
%
% Empirical data (black dots): gain-corrected then precision-weighted
% mean across subjects at each (direction, location) cell -- same
% weightedNanMean logic used in plot1/plot2/plotROISummary -- replacing
% the original script's unweighted varfun(@mean,...) across subjects.
%
% Orientation stimuli are axial (a 0-degree grating looks the same
% rotated 180 degrees), so each of the 4 real stimulus directions
% (0/90/45/135) is duplicated at direction+180 for display across all 8
% compass positions, exactly matching the original script's "if
% orientation_minus_baseline, duplicate direction+180" step -- valid
% because every one of the 4 asymmetry predictors is invariant to a
% direction+180 shift by construction (verified: {0,90,180,270} and
% {45,135,225,315} are each closed under +180, so mainCardinal etc. don't
% change), so the duplicated rows correctly reuse the same term values.

p = inputParser;
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.addParameter('githubDir', '~/Documents/GitHub', @ischar);
p.addParameter('precisionWeights', [], @(x) isempty(x) || istable(x));
p.addParameter('figureDir', '', @ischar); % '' = default production location (strrep(bidsDir,'data_bids','figures')/<projectName>)
p.addParameter('dgSubjectMode', 'all', @(x) ismember(x, {'all','matched'})); % only used when projectName='dg'; 'matched' = the same 7 subjects also run in da (sub-0395 excluded), matching fitAsymmetryRegression.m's dgSubjectMode
p.parse(varargin{:});
opt = p.Results;

githubDir = opt.githubDir;
bidsDir = opt.bidsDir;

addpath(genpath(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')));
cd(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode'));
setup_user('rania', bidsDir);

comparisonName = 'orientation_minus_baseline';
projectSettings = loadConfig(githubDir);
roi_idx = projectSettings.roi_idx;
rois = projectSettings.rois;
riSel = find(strcmp(rois, roiname));
if isempty(riSel)
    error('plotEachDirLocRegression:roi', 'roiname ''%s'' not found in projectSettings.rois.', roiname);
end
regionIndex = roi_idx{riSel};

if isempty(opt.figureDir)
    figureDirSuffix = '';
    if strcmp(projectName, 'dg')
        figureDirSuffix = ['_', opt.dgSubjectMode]; % matches plot_NeuralAsymmetries.m's figureDir convention (e.g. figures/dg_all or figures/dg_matched)
    end
    figureDir = [strrep(bidsDir, 'data_bids', 'figures'), projectName, figureDirSuffix];
else
    figureDir = opt.figureDir;
end
if ~isfolder(figureDir), mkdir(figureDir); end

contrasts_dict = projectSettings.contrasts_dict;
contrastnames = {contrasts_dict.contrasts.('dg_contrast_name')};
s0_idx = find(strcmp(contrastnames,'s0_v_b'));
s90_idx = find(strcmp(contrastnames,'s90_v_b'));
s45_idx = find(strcmp(contrastnames,'s45_v_b'));
s135_idx = find(strcmp(contrastnames,'s135_v_b'));
mdirvals_dg = [0, 90, 45, 135];
anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % index order matching meanBOLDpa's 2nd dim
maincardinalmDir = [0,90,180,270];
primaryMeridians = [90,0,270,180];

if strcmp(projectName, 'dg')
    if strcmp(opt.dgSubjectMode, 'all')
        subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
            'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
            'sub-0397', 'sub-0427'};
        fitLabel = 'dg';
    else % 'matched'
        subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
            'sub-0426', 'sub-0250'};
        fitLabel = 'dgMatched7';
    end
elseif strcmp(projectName, 'da')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0426', 'sub-0250'};
    fitLabel = 'da';
else
    error('plotEachDirLocRegression:project', 'projectName must be ''dg'' or ''da''.');
end

gainWeightsFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummaryByROI.mat');
Ggain = load(gainWeightsFile, 'gainTable');
gainWeights = retrieveObserverGainWeights2(subjects, roiname, Ggain.gainTable);
groupGain = exp(mean(log(gainWeights), 'omitnan')); % omitnan: see retrieveObserverGainWeights2.m
subjectScale = groupGain ./ gainWeights;
% Row vector, matching vals' orientation below -- MATLAB vector indexing
% preserves the INDEXED array's own orientation (not the index's), so a
% column here would silently broadcast into a matrix instead of an
% elementwise product against vals(valid) (same bug found and fixed
% earlier in plot1_experimentalCond.m).
precisionW = retrieveObserverPrecisionWeights(subjects, roiname, opt.precisionWeights);

fitFile = fullfile(bidsDir,'derivatives','summaryTables','regressionResults',fitLabel,sprintf('%s.mat',roiname));
if ~isfile(fitFile)
    error('plotEachDirLocRegression:missingFit', ...
        'Cached fit not found for %s / %s at %s -- run fitAsymmetryRegression(''%s'') first.', ...
        fitLabel, roiname, fitFile, projectName);
end
F = load(fitFile);
betaVec = F.estimates(:) / 2; % raw (undoubled) coefficients, same scale as M below expects
estimatesVec = [F.grandInterceptFE; betaVec]; % [intercept; mainCardinal; derivedCardinal; mainSubset; derivedSubset]

glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), 'hRF_glmsingle');
S1 = load(fullfile(glmResultsfolder, 'meanBOLDpa'));
meanBOLDpa_full = S1.meanBOLDpa;
if strcmp(projectName,'dg')
    dg_full13 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
        'sub-0397', 'sub-0427'};
    subjIdx = cellfun(@(s) find(strcmp(dg_full13,s)), subjects);
else
    da_full8 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
    subjIdx = cellfun(@(s) find(strcmp(da_full8,s)), subjects);
end
data = squeeze(meanBOLDpa_full([s0_idx,s90_idx,s45_idx,s135_idx], :, regionIndex, subjIdx)); % 4(direction) x 8(location) x nSubj

% weighted empirical mean across subjects, and the 4 asymmetry term
% values, at every (direction, location) cell -- including the
% direction+180 duplicate, since orientation stimuli are axial.
N = [0 45 90 135 180 225 270 315]; % the 8 compass "direction" slots
weightedBold = nan(numel(N), 8); % rows=direction(N), cols=location(pa, in anglevals order)
mainCardinalGrid = nan(numel(N), 8);
derivedCardinalGrid = nan(numel(N), 8);
mainSubsetGrid = nan(numel(N), 8);
derivedSubsetGrid = nan(numel(N), 8);

for li = 1:8
    pa = anglevals(li);
    for mi = 1:4
        md = mdirvals_dg(mi);
        vals = squeeze(data(mi, li, :))' .* subjectScale; % 1 x nSubj, gain-corrected
        valid = ~isnan(vals);
        wBold = sum(precisionW(valid) .* vals(valid)) / sum(precisionW(valid));

        mainCardinal = 2*ismember(md, maincardinalmDir) - 1;
        derivedCardinal = 2*((ismember(md,maincardinalmDir) & ismember(pa,primaryMeridians)) | ...
            (~ismember(md,maincardinalmDir) & ~ismember(pa,primaryMeridians))) - 1;
        mainSubset = 0; derivedSubset = 0;
        if strcmp(projectName, 'dg')
            if ismember(md,[0,180]); mainSubset=1; elseif ismember(md,[90,270]); mainSubset=-1; end
            if (abs(md-pa)==0 || abs(md-pa)==180); derivedSubset=1;
            elseif (abs(md-pa)==90 || abs(md-pa)==270); derivedSubset=-1; end
        else
            if ismember(md,[90,270]); mainSubset=1; elseif ismember(md,[0,180]); mainSubset=-1; end
            isCardMd = ismember(md,[0,90]); isOblMd = ismember(md,[45,135]);
            diffv = abs(md-pa);
            proH = (isCardMd && (diffv==90||diffv==270)) || (isOblMd && (diffv==0||diffv==180));
            conV = (isCardMd && (diffv==0||diffv==180)) || (isOblMd && (diffv==90||diffv==270));
            if proH; derivedSubset=1; elseif conV; derivedSubset=-1; end
        end

        % same (direction, direction+180) duplication as the original
        % script -- both compass slots get the identical term values and
        % weighted BOLD value, since all 4 terms are +-180-invariant.
        for d = [md, md+180]
            rowIdx = find(N == d);
            weightedBold(rowIdx, li) = wBold;
            mainCardinalGrid(rowIdx, li) = mainCardinal;
            derivedCardinalGrid(rowIdx, li) = derivedCardinal;
            mainSubsetGrid(rowIdx, li) = mainSubset;
            derivedSubsetGrid(rowIdx, li) = derivedSubset;
        end
    end
end

%% Plot: same "compass" spatial layout as lme2_ploteachDirLoc.m

% Square figure (Position width=height), so the axes' own normalized
% width/height fractions (p(3), p(4) below) already represent the SAME
% physical scale once 'axis square' is applied -- this is what makes the
% trigonometric shifts below give genuinely equal on-page distances, and
% is also what's needed for the saved PDF (below) to come out square.
figSize = 1250;
figure
set(gcf, 'Position', [1118 87 figSize figSize])
plot(0,0,'+k','MarkerSize',12, 'linewidth',3)
axis square
hold on
xlim([-4 4])
ylim([-4 4])
xticks([])
yticks([])
xticklabels({})
yticklabels({})
p = get(gca, 'Position');

box on
set(gca,'linewidth',1, 'YColor', [0 0 0]);
set(gca,'linewidth',1, 'XColor', [0 0 0]);

w = p(3)-p(3)*.835; h = p(4)-p(4)*.8;
left = p(1)+p(3)/2-(w/2); bottom = p(2)+p(4)/2-(h/2);
origin = [left bottom w h];

% Equal-radius placement (by construction, via trigonometry) for all 8
% locations, replacing the original script's separate "straight"/"diag"
% magnitudes (1.6 vs 1.1*sqrt(2)~=1.556), which were close but not
% exactly equal -- radius chosen to roughly match that original visual
% spacing (average of the two).
% Each subplot is placed at its own location's actual polar angle
% (anglevals(li)) as the compass direction from center -- same design
% intent as the original script (location value N(pos) doubled as the
% spatial placement angle), just computed directly via trigonometry
% instead of 8 hand-picked shift vectors.
radius = mean([w, h]) * mean([1.6, 1.1*sqrt(2)]);
shifts = cell(1,8);
for k = 1:8
    shifts{k} = [radius*cosd(anglevals(k)), radius*sind(anglevals(k)), 0, 0];
end

dotsize = 8;
globalMin = -0.5;
globalMax = 1;

modelPlotVals = nan(8,8); % rows=location(anglevals order), cols=direction(N order)

for li = 1:8 % location, matches anglevals(li) -- same subplot spatial position convention as the original (positioned by location index, not by the location's own polar angle)
    shift = shifts{li};

    pax = polaraxes(gcf);
    set(pax, 'Position', origin+shift)

    M = [ones(8,1), mainCardinalGrid(:,li), derivedCardinalGrid(:,li), mainSubsetGrid(:,li), derivedSubsetGrid(:,li)];
    predicted = M * estimatesVec;
    modelPlotVals(li,:) = predicted';

    rhoModel = [predicted; predicted(1)];
    polarplot(deg2rad([N N(1)]), rhoModel, 'r', 'linewidth', 2)
    hold on
    q = polarplot(deg2rad([N N(1)]), rhoModel, 'o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', dotsize);
    q.MarkerEdgeColor = 'w';
    hold on

    dataDot = weightedBold(:,li);
    pdot = polarplot(deg2rad(N), dataDot, 'ok', 'LineWidth', 1, 'MarkerFaceColor', 'black', 'MarkerSize', dotsize/1.5);
    pdot.MarkerEdgeColor = 'w';
    pax.FontSize = 6;
    pax.RTickLabel = {''};
    pax.ThetaTickLabel = {''};
    thetaticks(0:45:315);
    rlim([globalMin globalMax])
    hold on
end

sgtitle(sprintf('Regression-based est %s %s %s (gain+precision-weighted)', projectName, strrep(comparisonName,'_','-'), roiname), 'Interpreter', 'none')

% Explicit square PaperSize/PaperPosition, matching the figure's own
% square Position -- NOT '-bestfit', which scales/fits into a standard
% (non-square) page size and would undo the squareness set up above.
gcf_edit = gcf;
gcf_edit.Units = 'inches';
figDims = gcf_edit.Position(3:4); % square by construction (figSize x figSize)
gcf_edit.PaperUnits = 'inches';
gcf_edit.PaperPositionMode = 'manual';
gcf_edit.PaperSize = figDims;
gcf_edit.PaperPosition = [0, 0, figDims(1), figDims(2)];

print(gcf_edit, fullfile(figureDir, sprintf('EachDirLocRegression_%s_%s_%s', projectName, comparisonName, roiname)), '-dpdf');
close all;

end
