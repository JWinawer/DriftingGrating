function results = dg_compareGainROI(varargin)
% dg_compareGainROI - mean/median pRF gain in an ROI, moving vs stationary
%
%   results = dg_compareGainROI()
%   results = dg_compareGainROI('subjects', {'sub-0037', 'sub-0201'})
%
% For each subject, builds ONE vertex mask -- the union of the manually
% drawn lh/rh V1 labels, restricted to vertices with
%   variance explained > vexplMin
%   eccenMin <= eccentricity <= eccenMax
% -- and reports the mean/median gain within that mask for the moving
% ('movRoot', prfvista_mov) and stationary ('statRoot', prfvista) carrier.
%
% The pRF criteria (variance explained, eccentricity) always come from
% movRoot's lh/rh .vexpl.mgz and .eccen.mgz maps, regardless of which
% experiment gain is computed from. That keeps the vertex set identical
% between experiments, so gain differences reflect the stimulus, not a
% difference in which vertices were selected. See README.md.
%
% Gain itself is computed straight from each experiment's results.mat via
% RMMODELGAIN (same as DG_COMPUTEGAIN), not read from any precomputed file.
%
% OUTPUT
%   results - table with one row per subject: nVertices (the ROI+criteria
%             mask size), meanGain_mov, medianGain_mov, meanGain_stat,
%             medianGain_stat.
%
% REQUIREMENTS
%   Base MATLAB only. Reads FreeSurfer .label files and .mgz surface
%   overlays directly (no FreeSurfer matlab toolbox needed).
%
% See also DG_COMPUTEGAIN, RMMODELGAIN, DG_INSPECTPRFRESULTS

p = inputParser;
p.addParameter('movRoot', ['/Volumes/EXTERNAL_US/Project_dg/data_bids/' ...
    'derivatives/prfvista_mov'], @ischar);
p.addParameter('statRoot', ['/Volumes/EXTERNAL_US/Project_dg/data_bids/' ...
    'derivatives/prfvista'], @ischar);
p.addParameter('fsRoot', ['/Volumes/EXTERNAL_US/Project_dg/data_bids/' ...
    'derivatives/freesurfer'], @ischar);
p.addParameter('roiDir',  'label/retinotopy_RE', @ischar);
p.addParameter('roiName', 'V1_REmanual', @ischar);
p.addParameter('subjects', {}, @iscell);
p.addParameter('vexplMin', 0.1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('eccenMin', 4,   @(x) isnumeric(x) && isscalar(x));
p.addParameter('eccenMax', 8,   @(x) isnumeric(x) && isscalar(x));
p.addParameter('metric',   'maxabs', @ischar);
p.parse(varargin{:});
opt = p.Results;

if isempty(which('rmModelGain'))
    error('rmModelGain.m is not on the path. Add this handoff folder.');
end

if isempty(opt.subjects)
    opt.subjects = findCommonSubjects(opt);
end
fprintf('%d subject(s)\n\n', numel(opt.subjects));

rows = struct('subject', {}, 'nVertices', {}, ...
    'meanGain_mov', {}, 'medianGain_mov', {}, ...
    'meanGain_stat', {}, 'medianGain_stat', {});

for ii = 1:numel(opt.subjects)
    sub = opt.subjects{ii};
    fprintf('--- %s\n', sub);
    try
        rows(end+1) = computeOne(sub, opt); %#ok<AGROW>
    catch ME
        fprintf('    FAILED: %s\n', ME.message);
    end
end

fprintf('\n%d of %d subjects done.\n', numel(rows), numel(opt.subjects));
results = struct2table(rows);
return


% ------------------------------------------------------------------------
function subs = findCommonSubjects(opt)
d1 = dir(fullfile(opt.movRoot,  'sub-*')); d1 = d1([d1.isdir]);
d2 = dir(fullfile(opt.statRoot, 'sub-*')); d2 = d2([d2.isdir]);
d3 = dir(fullfile(opt.fsRoot,   'sub-*')); d3 = d3([d3.isdir]);
subs = intersect(intersect({d1.name}, {d2.name}), {d3.name});
return


% ------------------------------------------------------------------------
function row = computeOne(sub, opt)

% Identify the mask BEFORE computing any gain: it only takes label + mgz
% reads, so this is cheap. rmModelGain then only has to build Gaussians and
% predicted time series for the few hundred/thousand selected vertices
% instead of the whole ~150-450k-vertex surface -- a large speedup, and the
% same trick works regardless of how many experiments are compared.
sesDir = findResultsSesDir(opt.movRoot, sub);

vexplLH = readMgz(fullfile(sesDir, 'lh.vexpl.mgz'));
vexplRH = readMgz(fullfile(sesDir, 'rh.vexpl.mgz'));
eccenLH = readMgz(fullfile(sesDir, 'lh.eccen.mgz'));
eccenRH = readMgz(fullfile(sesDir, 'rh.eccen.mgz'));

nL = numel(vexplLH);
nR = numel(vexplRH);

% ROI: union of manually drawn lh/rh V1
roiLH = readLabelMask(fullfile(opt.fsRoot, sub, opt.roiDir, ...
    sprintf('lh.%s.label', opt.roiName)), nL);
roiRH = readLabelMask(fullfile(opt.fsRoot, sub, opt.roiDir, ...
    sprintf('rh.%s.label', opt.roiName)), nR);

critLH = vexplLH > opt.vexplMin & eccenLH >= opt.eccenMin & eccenLH <= opt.eccenMax;
critRH = vexplRH > opt.vexplMin & eccenRH >= opt.eccenMin & eccenRH <= opt.eccenMax;

% Concatenated lh-then-rh, matching the gain vector's vertex order
mask = [roiLH & critLH; roiRH & critRH];
nROI = sum(roiLH) + sum(roiRH);
voxels = find(mask);
fprintf('    %d vertices in V1, %d survive vexpl>%.2g & %g<=eccen<=%g\n', ...
    nROI, numel(voxels), opt.vexplMin, opt.eccenMin, opt.eccenMax);
if isempty(voxels)
    error('no vertices survive the ROI + pRF criteria');
end

nTotal = nL + nR;
gainMov  = gainFromRoot(opt.movRoot,  sub, opt.metric, voxels, nTotal);
gainStat = gainFromRoot(opt.statRoot, sub, opt.metric, voxels, nTotal);

row.subject         = sub;
row.nVertices       = numel(voxels);
row.meanGain_mov    = mean(gainMov);
row.medianGain_mov  = median(gainMov);
row.meanGain_stat   = mean(gainStat);
row.medianGain_stat = median(gainStat);
return


% ------------------------------------------------------------------------
function sesDir = findResultsSesDir(root, sub)
hits = dir(fullfile(root, sub, 'ses-*', '*results.mat'));
if isempty(hits)
    error('no *results.mat under %s', fullfile(root, sub));
end
sesDir = hits(1).folder;
return


% ------------------------------------------------------------------------
function gain = gainFromRoot(root, sub, metric, voxels, expectedNVerts)
hits = dir(fullfile(root, sub, 'ses-*', '*results.mat'));
if isempty(hits)
    error('no *results.mat under %s', fullfile(root, sub));
end
if numel(hits) > 1
    fprintf('    %d results files under %s, using %s\n', ...
        numel(hits), root, hits(1).name);
end
srcFile = fullfile(hits(1).folder, hits(1).name);

V = load(srcFile);
[model, params] = findModelParams(V);

if ~isfield(params, 'analysis') || ~isfield(params.analysis, 'allstimimages')
    error('params.analysis.allstimimages is missing in %s', srcFile);
end

m = model; if iscell(m), m = m{1}; end
nV = numel(m.x0);
if nV ~= expectedNVerts
    error('%s has %d vertices, expected %d (from movRoot''s mgz maps)', ...
        srcFile, nV, expectedNVerts);
end

rm.model  = model;
rm.params = params;
gain = rmModelGain(rm, 'voxels', voxels, 'metric', metric);
gain = double(gain(:));
return


% ------------------------------------------------------------------------
function [model, params] = findModelParams(V)
% The results file may store model/params at the top level, or nested
% inside a variable such as 'results'. Handle both.
if isfield(V, 'model') && isfield(V, 'params')
    model = V.model; params = V.params; return
end
top = fieldnames(V);
for k = 1:numel(top)
    n = V.(top{k});
    if isstruct(n) && isscalar(n) && isfield(n, 'model') && isfield(n, 'params')
        model = n.model; params = n.params; return
    end
end
error(['could not find model and params in the results file. ' ...
    'Top-level variables were: %s'], strjoin(top', ', '));


% ------------------------------------------------------------------------
function vol = readMgz(mgzPath)
% Minimal MGH/MGZ reader for the per-hemisphere surface overlays written by
% vistaPRF2MAP (angle/eccen/sigma/vexpl/x/y .mgz files). Base MATLAB only:
% gunzip plus fread against the standard MGH header layout. These files are
% always [nVertices x 1 x 1 x 1], so no reshape beyond the column read below
% is needed.
if ~exist(mgzPath, 'file')
    error('mgz file not found: %s', mgzPath);
end

tmpDir = tempname;
mkdir(tmpDir);
cleanupDir = onCleanup(@() rmdir(tmpDir, 's'));

files = gunzip(mgzPath, tmpDir);
fid = fopen(files{1}, 'rb', 'b');
cleanupFid = onCleanup(@() fclose(fid));

fread(fid, 1, 'int32'); % version
ndim1   = fread(fid, 1, 'int32');
ndim2   = fread(fid, 1, 'int32');
ndim3   = fread(fid, 1, 'int32');
nframes = fread(fid, 1, 'int32');
type    = fread(fid, 1, 'int32');
fread(fid, 1, 'int32'); % dof

UNUSED_SPACE_SIZE = 256;
USED_SPACE_SIZE   = (3*4 + 4*3*4);
unused_space_size = UNUSED_SPACE_SIZE - 2;

ras_good_flag = fread(fid, 1, 'short');
if ras_good_flag == 1
    fread(fid, 3, 'float32');
    fread(fid, 9, 'float32');
    fread(fid, 3, 'float32');
    unused_space_size = unused_space_size - USED_SPACE_SIZE;
end
fseek(fid, unused_space_size, 'cof');

nv = ndim1 * ndim2 * ndim3 * nframes;
switch type
    case 0, vol = fread(fid, nv, 'uchar');
    case 1, vol = fread(fid, nv, 'int32');
    case 3, vol = fread(fid, nv, 'float32');
    case 4, vol = fread(fid, nv, 'short');
    otherwise, error('unsupported mgz data type %d', type);
end
vol = double(vol);
return


% ------------------------------------------------------------------------
function mask = readLabelMask(labelPath, nVerts)
% FreeSurfer ASCII label: comment line, vertex count, then one row per
% vertex "index x y z value" with 0-based indices.
if ~exist(labelPath, 'file')
    error('label file not found: %s', labelPath);
end
fid = fopen(labelPath, 'r');
fgetl(fid);
n = str2double(fgetl(fid));
data = textscan(fid, '%f %f %f %f %f', n);
fclose(fid);

mask = false(nVerts, 1);
mask(data{1} + 1) = true;
return
