function Yhat = predict_harmonic(b, thetaV, expCfg, opts)
% PREDICT_HARMONIC  Apply harmonic coefficients to an experiment's geometry.
%
%   Yhat = predict_harmonic(b, thetaV, expCfg, opts)
%
% b      : 1 x nPred coefficients from FIT_HARMONIC_VERTEX (NaN entries, i.e. columns
%          that were dropped as identically zero, are treated as 0).
% thetaV : nVertex x 1 polar angles, degrees (conventional).
% expCfg : cfg.dg or cfg.da -- the geometry to predict INTO. Passing coefficients
%          fitted on one experiment together with the other experiment's expCfg is
%          the cross-experiment prediction: the same four mechanisms acting on
%          locally-rotated stimuli, with no context effect.
%
% Returns Yhat, nVertex x nOri predicted demeaned responses, columns in
% expCfg.oriCols order.

    if nargin < 4, opts = struct(); end

    [X, info] = harmonic_predictors(thetaV, expCfg, opts);
    bb = b(:);
    if numel(bb) ~= size(X, 2)
        error('predict_harmonic:size', ...
              'b has %d entries but the design has %d columns.', numel(bb), size(X,2));
    end
    bb(~isfinite(bb)) = 0;

    Yhat = reshape(X * bb, [info.nVertex, info.nOri]);
end
