clc; clear all; close all;

% computeMeanMTplus2PerSubject.m
%
% For each of the 8 subjects, computes a single scalar: the mean mt+2.mgz
% response (lh+rh combined) across valid V1 vertices. "Valid" uses the same
% voxel-filtering method as createTTaveTable.m/createTables.m -- pRF
% eccentricity in [4,8] deg, pRF R^2 >= 0.1, and a defined pRF polar angle
% -- restricted here to the V1 ROI only (not all 10 ROIs). Prints one
% number per subject.

%% SET UP

bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir)

subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250'};

hemis = {'lh'; 'rh'};
roiname = 'V1';

minECC = 4;
maxECC = 8;
minVAREXP = 0.1;

mtDir = fullfile(bidsDir, 'derivatives', 'dotsGLM', 'hRF_canonical');

meanResponse = nan(numel(subjects), 1);
nValidVertices = nan(numel(subjects), 1);

%% PER SUBJECT

for si = 1:numel(subjects)

    subjectname = subjects{si};
    fprintf('=== %s ===\n', subjectname);

    %% pRF metadata (same source/logic as createTTaveTable.m)

    retDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subjectname, '**/stimfiles.mat'));
    retDir = retDir.folder;

    for hi = 1:numel(hemis)
        hemi = hemis{hi};
        ret_moving.(sprintf('%s_pa', hemi)) = MRIread(fullfile(retDir, sprintf('%s.angle_adj.mgz', hemi)));
        ret_moving.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(retDir, sprintf('%s.eccen.mgz', hemi)));
        ret_moving.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(retDir, sprintf('%s.vexpl.mgz', hemi)));
    end

    hSize = [size(ret_moving.lh_pa.vol,2), size(ret_moving.rh_pa.vol,2)];
    N = sum(hSize);

    pRF_angle = map_theta([ret_moving.lh_pa.vol, ret_moving.rh_pa.vol])';
    pRF_ecc = [ret_moving.lh_ecc.vol, ret_moving.rh_ecc.vol]';
    pRF_r2 = [ret_moving.lh_vexp.vol, ret_moving.rh_vexp.vol]';

    % same "valid angle bin" check used elsewhere (excludes only vertices
    % with an undefined pRF angle, since the 8 bins otherwise tile the
    % full circle with no gaps)
    binCenters = 0:45:315;
    circDist = abs(mod(pRF_angle - binCenters + 180, 360) - 180);
    minDist = min(circDist, [], 2);
    validAngleBin = minDist <= 22.5;

    included = validAngleBin & pRF_ecc >= minECC & pRF_ecc <= maxECC & pRF_r2 >= minVAREXP;

    %% V1 ROI label

    lh_label = read_label(subjectname, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
    rh_label = read_label(subjectname, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));
    label_idx = [lh_label(:,1)+1 ; rh_label(:,1)+hSize(1)+1];

    isV1 = false(N, 1);
    isV1(label_idx) = true;

    validV1 = isV1 & included;
    nValidVertices(si) = sum(validV1);

    %% mt+2 response

    lh_mt = MRIread(fullfile(mtDir, subjectname, 'lh.mt+2.mgz'));
    rh_mt = MRIread(fullfile(mtDir, subjectname, 'rh.mt+2.mgz'));

    if numel(lh_mt.vol) ~= hSize(1) || numel(rh_mt.vol) ~= hSize(2)
        warning('%s: mt+2.mgz vertex count (lh=%d, rh=%d) does not match pRF surface size (lh=%d, rh=%d) -- check both are in the same (fsnative) space.', ...
            subjectname, numel(lh_mt.vol), numel(rh_mt.vol), hSize(1), hSize(2));
    end

    mtResponse = [lh_mt.vol', rh_mt.vol']';

    meanResponse(si) = mean(mtResponse(validV1), 'omitnan');

    fprintf('  %d valid V1 vertices, mean mt+2 response = %.4f\n', nValidVertices(si), meanResponse(si));
end

%% SUMMARY

fprintf('\nMean mt+2 response per subject (V1, valid voxels: ecc 4-8, R^2>=0.1):\n');
for si = 1:numel(subjects)
    fprintf('  %s: %.4f\n', subjects{si}, meanResponse(si));
end
