function [scale, gain] = observer_gain_weights(cfg)
% OBSERVER_GAIN_WEIGHTS  Per-observer pRF-gain rescaling factors.
%
%   [scale, gain] = observer_gain_weights(cfg)
%
% Reproduces AnalysisCode/04_plot_betaAsymmetries/retrieveObserverGainWeights.m plus the
% rescaling applied in lme1_fit.m / plot1_ / plot2_experimentalCond.m: each observer's
% BOLD is divided by their own mean pRF gain (gainSummary.csv column meanGain_avg,
% written by dg_computeGain.m from the independent prfvista/prfvista_mov protocols) and
% the group mean gain is multiplied back in, so
%
%       scale_i = groupGain / gain_i
%
% cfg.gainMean selects the group summary. The repo code uses the ARITHMETIC mean; the
% manuscript Methods say GEOMETRIC, which is the intended version (JW, 2026-08-19) --
% so 'geometric' is the default here and the repo should be updated to match.
% The choice is a single scalar common to all observers: it multiplies every rescaled
% value by the same constant (geometric/arithmetic = 0.9895 for these 8 observers), so
% it shifts reported effect sizes by ~1% and leaves every correlation, variance RATIO,
% t statistic and p value exactly unchanged.
%
% Returns scale as nSubj x 1 (ones, with a warning, if the gain file is unavailable).

    nS    = numel(cfg.subjects);
    scale = ones(nS, 1);
    gain  = nan(nS, 1);

    if isempty(cfg.gainFile) || ~isfile(cfg.gainFile)
        warning('observer_gain_weights:missing', ...
                'gain file not found (%s); returning unit weights', cfg.gainFile);
        return
    end

    G = readtable(cfg.gainFile);
    for si = 1:nS
        row = strcmp(G.subject, cfg.subjects{si});
        if any(row), gain(si) = G.meanGain_avg(row); end
    end
    if ~all(isfinite(gain))
        warning('observer_gain_weights:subject', ...
                'gain missing for some subjects; returning unit weights');
        gain(:) = NaN;  return
    end

    switch lower(cfg.gainMean)
        case 'geometric', groupGain = exp(mean(log(gain)));
        case 'arithmetic', groupGain = mean(gain);
        otherwise, error('observer_gain_weights:mean', 'cfg.gainMean must be geometric|arithmetic');
    end
    scale = groupGain ./ gain;
end
