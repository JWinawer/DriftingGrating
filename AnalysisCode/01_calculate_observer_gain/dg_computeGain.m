function summary = dg_computeGain(varargin)
% dg_computeGain - STEP 2: compute pRF gain at the server, save small files
%
%   summary = dg_computeGain()
%   summary = dg_computeGain('subjects', {'sub-0037'})
%
% For each subject, runs rmModelGain on both the moving (prfvista_mov) and
% stationary (prfvista) results files and writes one small gain_<subject>.mat
% per subject per protocol. Each output is a vector over the whole
% concatenated surface (left hemisphere then right), so it lines up
% one-to-one with the ret_<subject>.mat maps already in New York -- but for
% speed, rmModelGain is only asked to compute the V1 ROI + pRF-criteria
% vertices (see below), not the full ~150-450k-vertex surface. All other
% entries in the saved vector are 0 (not computed, not necessarily unfit --
% see info.note in the saved file).
%
% It also builds and saves a per-subject SUMMARY TABLE of mean/median gain
% within a V1 ROI:
%   subject, nVertices,
%   meanGain_mov, medianGain_mov,     - prfvista_mov (moving carrier)
%   meanGain_stat, medianGain_stat,   - prfvista (stationary carrier)
%   meanGain_avg, medianGain_avg      - average of the mov/stat columns
%
% The ROI is the union of the manually drawn lh/rh V1 labels, restricted to
% vertices with variance explained > vexplMin and eccenMin <= eccentricity
% <= eccenMax. Those criteria are always read from prfvista_mov's lh/rh
% .vexpl.mgz and .eccen.mgz maps, regardless of which protocol gain is
% computed from -- so the same vertex set is used for both, since that
% vertex set is what all subsequent analyses were selected on. See README.md
% and DG_COMPAREGAINROI, which uses the same rule.
%
% The point of doing this at the server: the results files are large, the
% gain vectors are about 1 MB each. Send back only the gain_*.mat files and
% the summary table.
%
% REQUIREMENTS
%   Base MATLAB only. No vistasoft, no toolboxes. Everything needed is in
%   this folder. Reads FreeSurfer .label files and .mgz surface overlays
%   directly (no FreeSurfer matlab toolbox needed).
%
% USAGE
%   addpath('/path/to/this/handoff/folder');
%   dg_computeGain('outdir', '~/dg_gain');
%
% Then send back everything in outdir.
%
% See also DG_INSPECTPRFRESULTS, DG_COMPAREGAINROI, RMMODELGAIN

p = inputParser;
p.addParameter('movRoot', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
    'data_bids/derivatives/prfvista_mov'], @ischar);
p.addParameter('statRoot', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
    'data_bids/derivatives/prfvista'], @ischar);
p.addParameter('fsRoot', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
    'data_bids/derivatives/freesurfer'], @ischar);
p.addParameter('roiDir',  'label/retinotopy_RE', @ischar);
p.addParameter('roiName', 'V1_REmanual', @ischar);
p.addParameter('outdir',   ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
    'data_bids/derivatives/summaryTables'], @ischar);
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

movOutdir  = fullfile(opt.outdir, 'prfvista_mov');
statOutdir = fullfile(opt.outdir, 'prfvista');
if ~exist(movOutdir, 'dir'),  mkdir(movOutdir);  end
if ~exist(statOutdir, 'dir'), mkdir(statOutdir); end

if isempty(opt.subjects)
    opt.subjects = findCommonSubjects(opt);
end

fprintf('movRoot : %s\n', opt.movRoot);
fprintf('statRoot: %s\n', opt.statRoot);
fprintf('fsRoot  : %s\n', opt.fsRoot);
fprintf('outdir  : %s\n', opt.outdir);
fprintf('%d subject(s)\n\n', numel(opt.subjects));

rows = struct('subject', {}, 'nVertices', {}, ...
    'meanGain_mov', {}, 'medianGain_mov', {}, ...
    'meanGain_stat', {}, 'medianGain_stat', {}, ...
    'meanGain_avg', {}, 'medianGain_avg', {});

nOK = 0;
for ii = 1:numel(opt.subjects)
    sub = opt.subjects{ii};
    fprintf('--- %s\n', sub);
    try
        rows(end+1) = computeOne(sub, opt, movOutdir, statOutdir); %#ok<AGROW>
        nOK = nOK + 1;
    catch ME
        fprintf('    FAILED: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('    at %s line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
end

fprintf('\n%d of %d subjects done.\n', nOK, numel(opt.subjects));

summary = struct2table(rows);
summaryMatFile = fullfile(opt.outdir, 'gainSummary.mat');
summaryCsvFile = fullfile(opt.outdir, 'gainSummary.csv');
save(summaryMatFile, 'summary');
writetable(summary, summaryCsvFile);
fprintf('wrote summary table to:\n  %s\n  %s\n', summaryMatFile, summaryCsvFile);
return


% ------------------------------------------------------------------------
function subs = findCommonSubjects(opt)
d1 = dir(fullfile(opt.movRoot,  'sub-*')); d1 = d1([d1.isdir]);
d2 = dir(fullfile(opt.statRoot, 'sub-*')); d2 = d2([d2.isdir]);
d3 = dir(fullfile(opt.fsRoot,   'sub-*')); d3 = d3([d3.isdir]);
subs = intersect(intersect({d1.name}, {d2.name}), {d3.name});
return


% ------------------------------------------------------------------------
function row = computeOne(sub, opt, movOutdir, statOutdir)

% ROI + pRF-criteria mask, ALWAYS from prfvista_mov's maps -- see header
% comment and README.md. Identifying the mask first is cheap (label + mgz
% reads only); it also lets us report the ROI-restricted mean/median without
% a second pass over the results files.
sesDir = findResultsSesDir(opt.movRoot, sub);

vexplLH = readMgz(fullfile(sesDir, 'lh.vexpl.mgz'));
vexplRH = readMgz(fullfile(sesDir, 'rh.vexpl.mgz'));
eccenLH = readMgz(fullfile(sesDir, 'lh.eccen.mgz'));
eccenRH = readMgz(fullfile(sesDir, 'rh.eccen.mgz'));

nL = numel(vexplLH);
nR = numel(vexplRH);

roiLH = readLabelMask(fullfile(opt.fsRoot, sub, opt.roiDir, ...
    sprintf('lh.%s.label', opt.roiName)), nL);
roiRH = readLabelMask(fullfile(opt.fsRoot, sub, opt.roiDir, ...
    sprintf('rh.%s.label', opt.roiName)), nR);

critLH = vexplLH > opt.vexplMin & eccenLH >= opt.eccenMin & eccenLH <= opt.eccenMax;
critRH = vexplRH > opt.vexplMin & eccenRH >= opt.eccenMin & eccenRH <= opt.eccenMax;

mask = [roiLH & critLH; roiRH & critRH];
nROI = sum(roiLH) + sum(roiRH);
voxels = find(mask);
fprintf('    %d vertices in V1, %d survive vexpl>%.2g & %g<=eccen<=%g (from prfvista_mov)\n', ...
    nROI, numel(voxels), opt.vexplMin, opt.eccenMin, opt.eccenMax);
if isempty(voxels)
    error('no vertices survive the ROI + pRF criteria');
end

nTotal = nL + nR;

gainMov  = computeAndSaveGain(opt.movRoot,  sub, opt.metric, nTotal, voxels, movOutdir);
gainStat = computeAndSaveGain(opt.statRoot, sub, opt.metric, nTotal, voxels, statOutdir);

gMov  = gainMov(voxels);
gStat = gainStat(voxels);

row.subject         = sub;
row.nVertices        = numel(voxels);
row.meanGain_mov     = mean(gMov);
row.medianGain_mov   = median(gMov);
row.meanGain_stat    = mean(gStat);
row.medianGain_stat  = median(gStat);
row.meanGain_avg     = mean([row.meanGain_mov,   row.meanGain_stat]);
row.medianGain_avg   = mean([row.medianGain_mov, row.medianGain_stat]);
return


% ------------------------------------------------------------------------
function gain = computeAndSaveGain(root, sub, metric, expectedNVerts, voxels, outdir)
hits = dir(fullfile(root, sub, 'ses-*', '*results.mat'));
if isempty(hits)
    error('no *results.mat under %s', fullfile(root, sub));
end
if numel(hits) > 1
    fprintf('    %d results files under %s, using %s\n', ...
        numel(hits), root, hits(1).name);
end
srcFile = fullfile(hits(1).folder, hits(1).name);

%% Load and normalise into the {model, params} shape rmModelGain expects
V = load(srcFile);
[model, params] = findModelParams(V);
clear V

if ~isfield(params, 'analysis') || ~isfield(params.analysis, 'allstimimages')
    error(['params.analysis.allstimimages is missing, so the predicted ' ...
        'time series cannot be built. Run dg_inspectPrfResults and send ' ...
        'the output back before going further.']);
end

rm.model  = model;
rm.params = params;

%% Provenance checks
m = model; if iscell(m), m = m{1}; end
nV = numel(m.x0);
if nV ~= expectedNVerts
    error('%s has %d vertices, expected %d (from movRoot''s mgz maps)', ...
        srcFile, nV, expectedNVerts);
end

%% Gain, restricted to the ROI + pRF-criteria vertices only (for speed --
% the summary table only needs these; see dg_compareGainROI for the same
% pattern). Embedded back into a full-surface-length vector so indexing
% still matches ret_<subject>.mat; everything outside voxels is left at 0.
gROI = rmModelGain(rm, 'voxels', voxels, 'metric', metric);
gROI = double(gROI(:));
gain = zeros(nV, 1);
gain(voxels) = gROI;

%% Sanity report, over the computed (ROI) vertices only
g = sort(gROI);
pct = @(q) g(max(1, min(numel(g), round(q/100*numel(g)))));
fprintf('    %s: %d vertices computed (ROI), gain median %.3f, 5th-95th pct %.3f - %.3f %%BOLD\n', ...
    root, numel(voxels), median(g), pct(5), pct(95));

info.subject    = sub;
info.sourceFile = srcFile;
info.metric     = metric;
info.created    = datestr(now, 'yyyy-mm-dd HH:MM:SS');
info.nVertices  = nV;
info.units      = 'percent BOLD, peak excursion of the predicted time series from baseline';
info.computedVoxels = voxels;
info.note       = ['Vertex order matches ret_<subject>.mat: left hemisphere ' ...
    'first, then right. Index it with the same indices. Only the V1 ROI + ' ...
    'pRF-criteria vertices in info.computedVoxels were computed (for speed); ' ...
    'all other entries are 0, which does NOT mean "unfit" here -- it means ' ...
    '"not requested".'];

outFile = fullfile(outdir, sprintf('gain_%s.mat', sub));
save(outFile, 'gain', 'info', '-v7');
dd = dir(outFile);
fprintf('    wrote %s (%.1f MB)\n', outFile, dd.bytes/1e6);
return


% ------------------------------------------------------------------------
function [model, params] = findModelParams(V)
% The results file may store model/params at the top level, or nested inside
% a variable such as 'results'. Handle both.
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
function sesDir = findResultsSesDir(root, sub)
hits = dir(fullfile(root, sub, 'ses-*', '*results.mat'));
if isempty(hits)
    error('no *results.mat under %s', fullfile(root, sub));
end
sesDir = hits(1).folder;
return


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
