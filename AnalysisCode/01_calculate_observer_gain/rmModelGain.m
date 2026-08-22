function [gain, prediction, voxels] = rmModelGain(rm, varargin)
% rmModelGain - pRF gain defined from the predicted time series
%
%   [gain, prediction, voxels] = rmModelGain(rm, 'param', value, ...)
%
% Returns one gain value per voxel (or vertex), in units of percent BOLD,
% derived from the model's own predicted time series rather than from the
% raw beta weight.
%
% WHY NOT JUST USE BETA
%   vistasoft builds pRFs with rfGaussian2d, which returns a Gaussian of
%   unit HEIGHT (the area normalization in that function is commented out
%   on purpose). The stimulus is scaled by sampleRate^2 and convolved with
%   a unit-AREA HRF, so the fitted beta has units of percent-BOLD/deg^2:
%
%       prediction(t) = beta(1) * (stim .* sampleRate^2 (*) hrf * RF)(t)
%                       + DC
%
%   Because RF has unit height, its "volume" grows as 2*pi*sigma^2, so
%   beta(1) shrinks as roughly 1/sigma^2 for a fixed BOLD amplitude. beta
%   is therefore not comparable across voxels with different pRF sizes.
%
%   The predicted time series does not have this problem: it is already in
%   percent BOLD, and already convolved with the HRF. Since the fMRI data
%   are converted to percent BOLD (mean 0) before fitting, the DC term of
%   the model IS the baseline, and the stimulus-driven part of the
%   prediction is exactly the displacement from that baseline. So we define
%
%       gain = max_t | beta(1) * pred(t) |
%
%   i.e. the largest excursion the model predicts away from baseline. The
%   trends/DC columns are deliberately excluded.
%
%   Note that we do NOT multiply by params.analysis.HrfMaxResponse here.
%   That factor is only needed when displaying the RF amplitude as an
%   instantaneous quantity (as rmPredictedTSeries and rmPlotGUI do), to
%   undo the unit-area HRF normalization. The predicted time series is
%   already convolved with the HRF, so its peak is already percent BOLD.
%
% INPUTS
%   rm      - path to a saved retinotopy model .mat file (as written by
%             rmSave, containing 'model' and 'params'), or a struct with
%             fields .model and .params already loaded.
%
%   Optional parameter/value pairs:
%   'modelnum'   - which entry of the model cell array. [1]
%   'voxels'     - linear indices of voxels to compute. [] = all voxels.
%   'metric'     - how to reduce the prediction to one number:
%                    'maxabs' largest absolute excursion from baseline
%                             (default; equals 'max' for positive-only fits)
%                    'max'    largest positive excursion
%                    'range'  max minus min. Use this if the stimulus never
%                             leaves the pRF unstimulated, so that there is
%                             no time point at which the prediction
%                             actually returns to baseline.
%   'chunksize'  - voxels processed per pass, to bound memory. [2000]
%   'prediction' - true to also return the per-voxel predicted time series
%                  (nTimePoints x nVoxels). Only do this for a subset;
%                  it is large. [false]
%   'plot'       - true to plot the predicted time series of the requested
%                  voxels (capped at 10 traces). [false]
%
% OUTPUTS
%   gain       - percent BOLD. Same size as rmGet(model,'x0') when all
%                voxels are computed (so it can be dropped straight into a
%                parameter map), otherwise a column vector matching
%                'voxels'.
%   prediction - nTimePoints x nVoxels stimulus-driven prediction, with the
%                DC/trend terms excluded. Empty unless 'prediction' is true.
%   voxels     - the linear voxel indices that were computed.
%
% Voxels with zero variance explained (unfit) or a zero/NaN sigma get a
% gain of 0.
%
% EXAMPLE
%   gain = rmModelGain('Gray/Averages/retModel-fFit.mat');
%
%   % one voxel, with its predicted time series plotted
%   [g, p] = rmModelGain(rmFile, 'voxels', 1234, 'prediction', true, ...
%                        'plot', true);
%
% DEPENDENCIES
%   None. This function is deliberately self-contained: it reads the model
%   fields directly and builds its own Gaussians, so it can be handed to
%   someone who does not have vistasoft installed. The local modelGet and
%   gauss2d subfunctions reproduce rmGet and rfGaussian2d exactly for the
%   parameters used here; test_rmModelGain checks that they agree.
%
% See also RMPREDICTEDTSERIES, RMGET, RFGAUSSIAN2D, RMMAKESTIMULUS.

%% Parse inputs
if nargin < 1 || isempty(rm), error('Need a model file or struct'); end

p = inputParser;
p.addParameter('modelnum',   1,        @(x) isnumeric(x) && isscalar(x));
p.addParameter('voxels',     [],       @isnumeric);
p.addParameter('metric',     'maxabs', @ischar);
p.addParameter('chunksize',  2000,     @(x) isnumeric(x) && isscalar(x));
p.addParameter('prediction', false,    @(x) islogical(x) || isnumeric(x));
p.addParameter('plot',       false,    @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

metric = lower(opt.metric);
if ~ismember(metric, {'maxabs' 'max' 'range'})
    error('Unknown metric "%s". Use maxabs, max, or range.', opt.metric);
end

%% Load model and params
if ischar(rm) || isstring(rm)
    if ~exist(rm, 'file'), error('Cannot find model file %s', rm); end
    rm = load(rm, 'model', 'params');
end
if ~isfield(rm, 'model') || ~isfield(rm, 'params')
    error('Expected a struct with fields "model" and "params".');
end

model  = rm.model;
if iscell(model), model = model{opt.modelnum}; end
params = rm.params;

if ~isfield(params, 'analysis') || ~isfield(params.analysis, 'X')
    error(['params.analysis is missing the sampling grid. This model file ' ...
        'was saved without its stimulus; regenerate it with rmMakeStimulus.']);
end

%% pRF parameters, flattened to one row per voxel
% All of these are reshaped by rmSave into the same view-shaped array, so
% linear indexing is consistent across them.
x0    = modelGet(model, 'x0');
sz    = size(x0);
x0    = x0(:);
y0    = modelGet(model, 'y0');         y0    = y0(:);
sMaj  = modelGet(model, 'sigmamajor'); sMaj  = sMaj(:);
sMin  = modelGet(model, 'sigmaminor'); sMin  = sMin(:);
theta = modelGet(model, 'sigmatheta'); theta = theta(:);
beta1 = modelGet(model, 'bcomp1');     beta1 = beta1(:);
ve    = modelGet(model, 'varexp');     ve    = ve(:);

nVoxTotal = numel(x0);

isCSS = isfield(model, 'exponent') && ~isempty(model.exponent);
if isCSS
    expt = modelGet(model, 'exponent'); expt = expt(:);
end

%% Which voxels
if isempty(opt.voxels)
    voxels  = (1:nVoxTotal)';
    allVox  = true;
else
    voxels  = opt.voxels(:);
    allVox  = false;
    if any(voxels < 1 | voxels > nVoxTotal)
        error('voxels must be linear indices in 1:%d', nVoxTotal);
    end
end
nVox = numel(voxels);

% Unfit voxels: all parameters are zero, which would make the Gaussian
% divide by zero. Skip them and report a gain of zero.
valid = ve(voxels) > 0 & sMaj(voxels) > 0 & isfinite(beta1(voxels));

%% Stimulus
X = params.analysis.X(:);
Y = params.analysis.Y(:);

if isCSS
    A = double(params.analysis.allstimimages_unconvolved);
else
    A = double(params.analysis.allstimimages);
end
nT = size(A, 1);

%% Loop over voxels in chunks
gainVec = zeros(nVox, 1);
if opt.prediction
    prediction = zeros(nT, nVox);
else
    prediction = [];
end

chunkStarts = 1:opt.chunksize:nVox;
for c = 1:numel(chunkStarts)
    lo   = chunkStarts(c);
    hi   = min(lo + opt.chunksize - 1, nVox);
    here = lo:hi;
    here = here(valid(here));
    if isempty(here), continue; end

    v = voxels(here);

    % Unit-height Gaussians, one column per voxel
    RFs = gauss2d(X, Y, sMaj(v), sMin(v), theta(v), x0(v), y0(v));

    if isCSS
        % Static nonlinearity acts on the neural (unconvolved) response,
        % then we convolve, scan by scan.
        pred = (A * RFs) .^ repmat(expt(v)', nT, 1);
        pred = convolveWithHrf(pred, params);
    else
        % allstimimages is already HRF-convolved
        pred = A * RFs;
    end

    % Scale by the fitted pRF beta. The DC/trend terms are the baseline
    % and are deliberately left out.
    pred = pred .* repmat(beta1(v)', nT, 1);

    switch metric
        case 'maxabs', g = max(abs(pred), [], 1);
        case 'max',    g = max(pred, [], 1);
        case 'range',  g = max(pred, [], 1) - min(pred, [], 1);
    end
    gainVec(here) = g(:);

    if opt.prediction
        prediction(:, here) = pred;
    end
end

%% Shape the output
if allVox
    gain = reshape(gainVec, sz);
else
    gain = gainVec;
end

%% Optional plot
if opt.plot
    if isempty(prediction)
        [~, prediction] = rmModelGain(rm, 'modelnum', opt.modelnum, ...
            'voxels', voxels, 'metric', metric, 'prediction', true);
    end
    nShow = min(10, nVox);
    figure('Color', 'w'); hold on;
    plot(prediction(:, 1:nShow), 'LineWidth', 1);
    plot([1 nT], [0 0], 'k-');
    xlabel('Time (frames)');
    ylabel('Predicted BOLD (% signal, baseline removed)');
    title(sprintf('pRF predictions, %d voxel(s); gain = %s', nShow, metric));
    for ii = 1:nShow
        plot([1 nT], gainVec(ii)*[1 1], '--');
    end
end

return


% ------------------------------------------------------------------------
function val = modelGet(model, param)
% Local stand-in for rmGet, covering only the parameters this function uses.
% Semantics deliberately match rmGet case for case, so that rmModelGain has
% no vistasoft dependency and can be shipped on its own.
switch lower(param)
    case 'x0'
        val = model.x0;
    case 'y0'
        val = model.y0;
    case 'sigmamajor'
        val = model.sigma.major;
    case 'sigmaminor'
        val = model.sigma.minor;
    case 'sigmatheta'
        val = model.sigma.theta;
    case 'bcomp1'
        % Betas live in the last dimension; component 1 is the pRF gain.
        switch ndims(model.beta)
            case 2, val = model.beta(1, :);
            case 3, val = model.beta(:, :, 1);
            case 4, val = model.beta(:, :, :, 1);
            otherwise
                error('Unexpected beta with %d dimensions', ndims(model.beta));
        end
    case 'varexp'
        % rmGet clamps to [0 1] and maps non-finite to 0
        val = 1 - (model.rss ./ model.rawrss);
        val(~isfinite(val)) = 0;
        val = min(max(val, 0), 1);
    case 'exponent'
        if isfield(model, 'exponent') && ~isempty(model.exponent)
            val = model.exponent;
        else
            val = ones(size(model.sigma.major));
        end
    otherwise
        error('modelGet: unknown parameter %s', param);
end
val = double(val);
return


% ------------------------------------------------------------------------
function RF = gauss2d(X, Y, sigmaMajor, sigmaMinor, theta, x0, y0)
% Local stand-in for rfGaussian2d. Unit HEIGHT, not unit area, matching
% vistasoft: the area normalization is intentionally absent there, and the
% fitted beta absorbs the pRF volume.
if numel(sigmaMajor) ~= 1
    sz1 = numel(X);
    sz2 = numel(sigmaMajor);
    X   = repmat(X(:), 1, sz2);
    Y   = repmat(Y(:), 1, sz2);
    sigmaMajor = repmat(sigmaMajor(:)', sz1, 1);
    sigmaMinor = repmat(sigmaMinor(:)', sz1, 1);
    if any(theta(:)), theta = repmat(theta(:)', sz1, 1); end
    x0 = repmat(x0(:)', sz1, 1);
    y0 = repmat(y0(:)', sz1, 1);
end

X = X - x0;
Y = Y - y0;

if any(theta(:))
    Xold = X;
    Yold = Y;
    X = Xold .* cos(theta) - Yold .* sin(theta);
    Y = Xold .* sin(theta) + Yold .* cos(theta);
end

RF = exp(-.5 * ((Y ./ sigmaMajor).^2 + (X ./ sigmaMinor).^2));
return


% ------------------------------------------------------------------------
function pred = convolveWithHrf(pred, params)
% Convolve each scan's segment with that scan's HRF, matching the way
% rmPredictedTSeries reconvolves the CSS prediction.
if isfield(params.analysis, 'scan_number') && ~isempty(params.analysis.scan_number)
    scanNumber = params.analysis.scan_number;
else
    scanNumber = ones(size(pred, 1), 1);
end
for scan = 1:numel(params.stim)
    these = scanNumber == scan;
    if ~any(these), continue; end
    hrf = params.analysis.Hrf{scan};
    pred(these, :) = filter(hrf, 1, pred(these, :));
end
return
