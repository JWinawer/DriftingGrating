clc; clear all; 

% load in the nonNormalized 13 betas and characterize it within V1:
projectName = 'dg'; % 'dg', or 'dgl' or 'da'

subjects = {'sub-0037', 'sub-0426', 'sub-0250', 'sub-0255', 'sub-0201', 'sub-0395', 'sub-wlsubj123', 'sub-wlsubj124'};
averageOri_array = []; averageSigma_array = []; meanTime_array = [];
averageCard_array = []; averageObl_array = []; medianSigma_array = [];
minSigma_array = []; maxSigma_array = [];

for si=1:length(subjects)

    subj = subjects{si}
    %ses = 'ses-01'; %'ses-01'; %'ses-nyu3t02'; %'ses-01'; <-- can I make this not manual?
    space = 'fsnative';
    hRF_setting = 'glmsingle';
    
    % checks/sets working directory and adds dependencies
    bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
    setup_user('rania', bidsDir)
    jsonparamfile = 'setup.json';
    jsonParams = jsondecode(fileread(jsonparamfile));
    
    % load predefined paths
    bidsDir = jsonParams.bidsdir.Path; 

    % find session
    findSessionInDir = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
        sprintf('hRF_%s', hRF_setting), subj); 
    sesDirs = dir(fullfile(findSessionInDir, 'ses-*'));
    sesNames = {sesDirs([sesDirs.isdir]).name};
    if length(sesNames) > 1
        error("more than one session found")
    else
        ses = sesNames{1};
    end

    derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
        sprintf('hRF_%s', hRF_setting), subj, ses);
    
    
    load(fullfile(derivativesFolder,'modelOutput.mat'));
    
    %
    results.(char('allevents')) = modelOut{1,4};
    
    %load(fullfile(derivativesFolder, 'nonNormalized13betas.mat'))
    betamaps = [];
    for ci=1:13
        condSelect = designSINGLE.stimorder==ci;
        betas = results.allevents.modelmd(:,:,:,condSelect);
        newbetas = squeeze(betas);
        
        % trial mean
        betamaps(:,ci) = mean(newbetas,2);
    end
    
    betamaps_condOnly = betamaps(:,1:13);
    
    % also load the ROI
    fileName = sprintf('ttaveSignal_%s_V1_mainCardinalVsMainOblique_motion_minus_baseline.mat', subj);
    load(fullfile(derivativesFolder, 'ttaveData','motion_minus_baseline',fileName))
    idx = find(~isnan(projectSettings.filteredPrfBins)); % ROI that meets ecc and r^2 criteria
    
    % undo PSC from GLMsingle that /mean and * 100
    %fixMean = modelOut{1,4}.meanvol;
    %betamaps_condOnly = betamaps_condOnly.*fixMean;

    % this should be computed within V1
    mu1 = mean(betamaps_condOnly(idx,1), 'omitnan');
    mu2 = mean(betamaps_condOnly(idx,2), 'omitnan');
    mu3 = mean(betamaps_condOnly(idx,3), 'omitnan');
    mu4 = mean(betamaps_condOnly(idx,4), 'omitnan');
    mu5 = mean(betamaps_condOnly(idx,5), 'omitnan');
    mu6 = mean(betamaps_condOnly(idx,6), 'omitnan');
    mu7 = mean(betamaps_condOnly(idx,7), 'omitnan');
    mu8 = mean(betamaps_condOnly(idx,8), 'omitnan');
    mu9 = mean(betamaps_condOnly(idx,9), 'omitnan');
    mu10 = mean(betamaps_condOnly(idx,10), 'omitnan');
    mu11 = mean(betamaps_condOnly(idx,11), 'omitnan');
    mu12 = mean(betamaps_condOnly(idx,12), 'omitnan');
    mu13 = mean(betamaps_condOnly(idx,13), 'omitnan');
    
    % these do not have a mean of 0 (b/c were computed as % relative to mean
    % timeseries including padding)
    mu = [mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10, mu11, mu12, mu13];
    stim_minus_blank = mu - mu(13);
    averageOri = mean(stim_minus_blank(9:12))
    averageCard = mean(stim_minus_blank(9:10));
    averageObl = mean(stim_minus_blank(11:12));
    averageOri_array = [averageOri_array, averageOri];
    averageCard_array = [averageCard_array, averageCard];
    averageObl_array = [averageObl_array, averageObl];

    meanTime = mean(modelOut{1,4}.meanvol);
    meanTime_array = [meanTime_array meanTime];
    
    % sigma
    sigma = nanmean(std(betamaps_condOnly, 0, 2))
    mediansigma = nanmedian(std(betamaps_condOnly, 0, 2))
    minsigma = nanmin(std(betamaps_condOnly, 0, 2));
    maxsigma = nanmax(std(betamaps_condOnly, 0, 2));
    averageSigma_array = [averageSigma_array, sigma];
    medianSigma_array = [medianSigma_array, mediansigma];
    minSigma_array = [minSigma_array, minsigma];
    maxSigma_array = [maxSigma_array, maxsigma];
    
    % pseudo z-score (but ok for now)
    zscore = (mu - mean(mu)) / std(mu);
    stim_minus_blank_zscore = zscore - zscore(13); % this makes blank 0, and most of the negative z-scores positive
    averageOri_zscore = mean(stim_minus_blank_zscore(9:12)); % sigma values
end