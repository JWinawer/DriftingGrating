function [X, info] = harmonic_predictors(thetaV, expCfg, opts)
% HARMONIC_PREDICTORS  Per-vertex harmonic design matrix for the orientation model.
%
%   [X, info] = harmonic_predictors(thetaV, expCfg, opts)
%
% Builds the design for the per-vertex model of the four stationary-orientation
% betas (see ../supplement/SUPPLEMENT_harmonic_model.md):
%
%   y_vk = b1*cos(2*theta) + b2*cos(4*theta)
%        + b3*cos(2*(theta-thetaV)) + b4*cos(4*(theta-thetaV))
%
% ANGLE CONVENTION (settled; see ../STIMULUS_CONVENTIONS.md section 2).
%   theta  = orientation of the grating BARS, in conventional visual-field degrees
%            (0 = rightward horizontal meridian, CCW positive), mod 180.
%            So theta=0 is a HORIZONTAL grating and theta=90 a VERTICAL one.
%            NOTE this is NOT the "direction of luminance variation" convention,
%            which is rotated 90 deg and would flip the sign of BOTH first-harmonic
%            terms (b1 would mean vertical>horizontal, b3 tangential>radial).
%   thetaV = vertex pRF polar angle, same frame. The CSV's pRF_angle has already
%            been through map_theta (createTables.m:75), so it is conventional.
%   A stimulus is RADIAL when its bars lie along the radius, i.e. theta == thetaV,
%   giving cos(2*(theta-thetaV)) = +1.
%
% Each stimulus's local orientation follows the repo rule
%   theta = mod(oriAngle + isPolar*(thetaV - 90), 180)
% (config_repro.m:59-65, compute_asymmetries.m:36, lme_codes.m:21), so the corrected
% spiral identities are inherited from expCfg rather than restated here.
%
% Inputs
%   thetaV : nVertex x 1 polar angles in degrees (conventional).
%   expCfg : cfg.dg or cfg.da (uses .oriAngle and .isPolar).
%   opts   : optional struct.
%            .expanded (default false) -- also return the four sin columns, giving
%                      the complete harmonic basis at harmonics 2 and 4. The core
%                      model imposes two testable constraints relative to this
%                      (see ../supplement/SUPPLEMENT_harmonic_model.md); the sin
%                      columns should be ~0
%                      under left-right visual-field symmetry.
%
% Outputs
%   X    : (nVertex*nOri) x nPred design. Rows are the column-major unfolding of an
%          (nVertex x nOri) array, i.e. row = v + (k-1)*nVertex, matching Y(:) for
%          a Y that is nVertex x nOri. nPred = 4, or 8 if opts.expanded.
%   info : .names   1 x nPred predictor names. The core four are named for the
%                   asymmetry they code, matching LME_CODES / FIT_LME_FIG7:
%                   hVv, cVo, rVt, pcVpo.
%          .theta   nVertex x nOri local orientation in degrees.
%          .d       nVertex x nOri local orientation relative to the radius.
%          .nOri, .nVertex
%          .isZeroCol  1 x nPred, true where a column is identically zero (this
%                   happens for sin(4*theta) in the Cartesian experiment, because
%                   sin(4*theta)==0 at all four 45-deg-spaced orientations).
%
% At the eight wedge centres the four core columns reduce exactly to the +1/0/-1
% and +/-1 codes of LME_CODES; TEST_HARMONIC_MODEL asserts this.

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'expanded'), opts.expanded = false; end

    tv  = double(thetaV(:));                        % nVertex x 1
    ori = double(expCfg.oriAngle(:)).';             % 1 x nOri
    nV  = numel(tv);
    nO  = numel(ori);

    % Local orientation of each stimulus at each vertex, and its offset from the
    % radial axis. mod 180 is cosmetic here (cos(2x), cos(4x) both have period 180)
    % but keeps info.theta directly comparable to expCfg.oriAngle.
    theta = mod(ori + double(expCfg.isPolar) * (tv - 90), 180);   % nVertex x nOri
    d     = theta - tv;                                            % nVertex x nOri

    cols  = {cosd(2*theta), cosd(4*theta), cosd(2*d), cosd(4*d)};
    names = {'hVv', 'cVo', 'rVt', 'pcVpo'};

    if opts.expanded
        cols  = [cols,  {sind(2*theta), sind(4*theta), sind(2*d), sind(4*d)}];
        names = [names, {'sin2ori', 'sin4ori', 'sin2rad', 'sin4rad'}];
    end

    X = zeros(nV*nO, numel(cols));
    for j = 1:numel(cols)
        c = cols{j};
        X(:,j) = c(:);
    end

    info.names     = names;
    info.theta     = theta;
    info.d         = d;
    info.nOri      = nO;
    info.nVertex   = nV;
    info.isZeroCol = all(abs(X) < 1e-12, 1);
end
