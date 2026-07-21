clc; clear all; close all;

% create table
% ROW = vertex
% COLUMN = vertex idx, subject, ROI, 13 betas DG (stimuli), 13 betas DA (polar stimuli), pRF ecc, pRF
% PA, pRF R^2

bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
allTables = cell(length(subjects),1);
projects = {'dg','da'};
savedir = fullfile(bidsDir, 'derivatives', 'summaryTables');

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

hemis = {'lh'; 'rh'};
roinames = {'V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b', 'hMTcomplex', 'pMT', 'pMST'};

dg_stimNames = {'cartexp_vertical_grating_rightwards_motion', 'cartexp_horizontal_grating_upwards_motion', ...
    'cartexp_vertical_grating_leftwards_motion', 'cartexp_horizontal_grating_downwards_motion', ...
    'cartexp_leftleaning_grating_upperrightwards_motion', 'cartexp_rightleaning_grating_upperleftwards_motion', ...
    'cartexp_leftleaning_grating_lowerleftwards_motion', 'cartexp_rightleaning_grating_lowerrightwards_motion', ...
    'cartexp_horizontal_stationary', 'cartexp_vertical_stationary', 'cartexp_rightleaning_grating_stationary', ...
    'cartexp_leftleaning_grating_stationary', 'cartexp_blank'};

da_stimNames = {'polexp_pinwheel_grating_clockwise_motion', 'polexp_annulus_grating_outwards_motion', ...
    'polexp_pinwheel_grating_cclockwise_motion', 'polexp_annulus_grating_inwards_motion', ...
    'polexp_ccspiral_grating_clockoutwards_motion', 'polexp_cspiral_grating_cclockoutwards_motion', ...
    'polexp_ccspiral_grating_cclockinwards_motion', 'polexp_cspiral_grating_clockinwards_motion', ...
    'polexp_annulus_grating_stationary', 'polexp_pinwheel_grating_stationary', 'polexp_cspiral_grating_stationary', ...
    'polexp_ccspiral_grating_stationary', 'polexp_blank'};

for si=1:numel(subjects)
    subjectname = subjects{si};

    retDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subjectname, '**/stimfiles.mat'));
    retDir = retDir.folder;
    
    for hi=1:numel(hemis)
        hemi = hemis{hi};
        ret_moving.(sprintf('%s_pa', hemi)) = MRIread(fullfile(retDir, sprintf('%s.angle_adj.mgz', hemi)));
        ret_moving.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(retDir, sprintf('%s.eccen.mgz', hemi)));
        ret_moving.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(retDir, sprintf('%s.vexpl.mgz', hemi)));
        ret_moving.(sprintf('%s_sigma', hemi)) = MRIread(fullfile(retDir, sprintf('%s.sigma.mgz', hemi)));
    end
    
    retData = [ret_moving.lh_pa.vol, ret_moving.rh_pa.vol ; ...
            ret_moving.lh_ecc.vol, ret_moving.rh_ecc.vol ; ...
            ret_moving.lh_vexp.vol, ret_moving.rh_vexp.vol; ...
            ret_moving.lh_sigma.vol, ret_moving.rh_sigma.vol]';
    
    hSize = [size(ret_moving.lh_pa.vol,2), size(ret_moving.rh_pa.vol,2)];
    size_combinedSurfs = sum(hSize);
    
    % create table with data
    N = size_combinedSurfs;
    
    T = table( ...
        strings(N,1), ...    % subject
        strings(N,1), ...    % visual_area
        nan(N,1), ...        % pRF_angle_bin
        nan(N,1), ...        % pRF_angle
        nan(N,1), ...        % pRF_ecc
        nan(N,1), ...        % pRF_r2
        nan(N,1), ...        % pRF_sigma
        false(N,1), ...      % included
        'VariableNames', {'subject','visual_area', 'pRF_angle_bin', 'pRF_angle','pRF_ecc', ...
        'pRF_r2','pRF_sigma', 'included'});
    
    T.subject = repmat(string(subjectname), N, 1);
    % angle converted from Benson coordinates to 0 to 360 (cc from right horizontal)
    T.pRF_angle = map_theta([ret_moving.lh_pa.vol, ret_moving.rh_pa.vol])'; 
    T.pRF_ecc = [ret_moving.lh_ecc.vol, ret_moving.rh_ecc.vol]';
    T.pRF_r2 = [ret_moving.lh_vexp.vol, ret_moving.rh_vexp.vol]';
    T.pRF_sigma = [ret_moving.lh_sigma.vol, ret_moving.rh_sigma.vol]';
    
    
    % visual field bin
    binCenters = 0:45:315;    % Bin labels (degrees)
    angles = T.pRF_angle;
    
    % Circular distance from each angle to each bin center
    circDist = abs(mod(angles - binCenters + 180, 360) - 180);
    
    % Closest bin and its distance
    [minDist, idx] = min(circDist, [], 2);
    
    % Initialize as NaN
    T.pRF_angle_bin = nan(height(T),1);
    
    % Assign only if within +/-22.5°
    valid = minDist <= 22.5;
    T.pRF_angle_bin(valid) = binCenters(idx(valid));
    
    % included in analysis?
    T.included = ...
        ~isnan(T.pRF_angle_bin) & ...
        T.pRF_ecc >= 4 & T.pRF_ecc <= 8 & ...
        T.pRF_r2 >= 0.1;
    

    % read ROI labels    
    for ri=1:length(roinames)
        roiname = roinames{ri};
        
        lh_label = read_label(subjectname, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
        rh_label = read_label(subjectname, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));
        
        % plus one because matlab starts from 1 not 0
        label_idx = [lh_label(:,1)+1 ; rh_label(:,1)+hSize(1)+1];
    
        if strcmp(roiname, 'pMT')
            roiname = 'MT';
        elseif strcmp(roiname, 'pMST')
            roiname = 'MST';
        end
        
        T.visual_area(label_idx) = roiname;
    end
    
    % betas
    
    % (1) M_0,      cartexp_vertical_grating_rightwards_motion
    %               polexp_pinwheel_grating_clockwise_motion
    % (2) M_90,     cartexp_horizontal_grating_upwards_motion 
    %               polexp_annulus_grating_outwards_motion
    % (3) M_180,    cartexp_vertical_grating_leftwards_motion  
    %               polexp_pinwheel_grating_cclockwise_motion
    % (4) M_270,    cartexp_horizontal_grating_downwards_motion  
    %               polexp_annulus_grating_inwards_motion
    % (5) M_45,     cartexp_leftleaning_grating_upperrightwards_motion  
    %               polexp_ccspiral_grating_clockoutwards_motion
    % (6) M_135,    cartexp_rightleaning_grating_upperleftwards_motion  
    %               polexp_cspiral_grating_cclockoutwards_motion
    % (7) M_225,    cartexp_leftleaning_grating_lowerleftwards_motion  
    %               polexp_ccspiral_grating_cclockinwards_motion
    % (8) M_315,    cartexp_rightleaning_grating_lowerrightwards_motion  
    %               polexp_cspiral_grating_clockinwards_motion
    % (9) S_0,      cartexp_horizontal_stationary             
    %               polexp_annulus_grating_stationary
    % (10) S_90,    cartexp_vertical_stationary            
    %               polexp_pinwheel_grating_stationary
    % (11) S_45,    cartexp_rightleaning_grating_stationary          
    %               polexp_cspiral_grating_stationary
    % (12) S_135,   cartexp_leftleaning_grating_stationary        
    %               polexp_ccspiral_grating_stationary
    % (13) B,       cartexp_blank  
    %               polexp_blank
    
    for iProj = 1:numel(projects)
    
        projectName = projects{iProj};
    
        if strcmp(projectName, 'dg')
            stimNames = dg_stimNames;
        elseif strcmp(projectName, 'da')
            stimNames = da_stimNames;
        end
    
        % Path to betafile
        subjdir = fullfile(bidsDir, 'derivatives', ...
            [projectName 'GLM'], 'hRF_glmsingle', subjectname);
    
        sesFolder = findses(subjdir);
    
        betafile = load(fullfile(subjdir, sesFolder, 'betas_nonzscored.mat'));
    
        for iStim = 1:numel(stimNames)
            T.(stimNames{iStim}) = betafile.betamaps(:, iStim);
        end
    
        T.(strcat(projectName, '_beta_mean')) = mean(betafile.betamaps, 2);
        T.(strcat(projectName, '_beta_std')) = std(betafile.betamaps, 0, 2);
    
    end
    
    % save this table with columns per beta / stimulus
    allTables{si} = T;

end

% Concatenate all tables vertically
allsubjectsTable = vertcat(allTables{:});

filename = fullfile(savedir, 'allsubjectsTable.csv');
writetable(allsubjectsTable, filename);

filename = fullfile(savedir, 'allsubjectsTable.mat');
save(filename, 'allsubjectsTable');

%%

function sesFolder = findses(subjdir)
    % find session
    d = dir(fullfile(subjdir, 'ses-*'));
    d = d([d.isdir]);    % Keep only directories
    if isempty(d)
        error('No ses-* directory found in %s', subjdir);
    elseif numel(d) > 1
        error('Multiple ses-* directories found in %s', subjdir);
    end
    sesFolder = d.name;
end

