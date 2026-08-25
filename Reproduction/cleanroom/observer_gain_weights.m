function [scale, gain] = observer_gain_weights(cfg, area, eccBand)
% OBSERVER_GAIN_WEIGHTS  Per-observer pRF-gain rescaling factors.
%
%   [scale, gain] = observer_gain_weights(cfg)                  % V1, legacy
%   [scale, gain] = observer_gain_weights(cfg, 'V2')            % per map
%   [scale, gain] = observer_gain_weights(cfg, 'hV4', '2-10')   % per map and band
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
% PER-MAP GAIN. With no `area`, this reads gainSummary.csv, whose single gain per
% observer is computed over V1 alone (DG_GAININV1) -- applying that to V2/V3/hV4
% special-cases V1 as the source. With an `area`, it reads gain_areas_summary.csv
% instead (written by ../server_extract/COLLECT_GAIN_AREAS), which holds one gain per
% observer per map per eccentricity band, and combines the two pRF protocols by
% GEOMETRIC mean. eccBand defaults to '4-8', matching cfg.eccRange; pass '2-10' or
% '1-12' where the analysis relaxes the criteria.
%
% How much this matters, measured on the current spec: per-map gain differs from the
% V1-derived scalar by at most 0.0073 in the group estimate (hV4; 0.0045 V2, 0.0051 V3),
% against 0.0127 for gain-versus-no-gain. So the per-map refinement is real but small,
% and smaller than every other choice in SPECIFICATION.md's review.
%
% Returns scale as nSubj x 1 (ones, with a warning, if the gain file is unavailable).

    if nargin < 2, area = ''; end
    if nargin < 3 || isempty(eccBand), eccBand = '4-8'; end

    nS    = numel(cfg.subjects);
    scale = ones(nS, 1);
    gain  = nan(nS, 1);

    if ~isempty(area)
        f = fullfile(fileparts(cfg.gainFile), 'gain_areas_summary.csv');
        if ~isfile(f)
            warning('observer_gain_weights:noAreaFile', ...
                    ['per-map gain not found (%s); run ' ...
                     '../server_extract/collect_gain_areas.m'], f);
            return
        end
        A = readtable(f);
        for si = 1:nS
            r = strcmp(A.subject, cfg.subjects{si}) & strcmp(A.area, area) & ...
                strcmp(A.eccBand, eccBand);
            if any(r), gain(si) = A.meanGain_avg(find(r, 1)); end
        end
        if ~all(isfinite(gain))
            warning('observer_gain_weights:areaSubject', ...
                    'per-map gain missing for some subjects in %s %s', area, eccBand);
            gain(:) = NaN;  return
        end
        scale = geo_or_arith(gain, cfg.gainMean) ./ gain;
        return
    end

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

    scale = geo_or_arith(gain, cfg.gainMean) ./ gain;
end

% ------------------------------------------------------------------------
function g = geo_or_arith(gain, how)
    switch lower(how)
        case 'geometric',  g = exp(mean(log(gain)));
        case 'arithmetic', g = mean(gain);
        otherwise, error('observer_gain_weights:mean', 'cfg.gainMean must be geometric|arithmetic');
    end
end
