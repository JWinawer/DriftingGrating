function [M, counts] = bin_and_aggregate(T, cfg, expCfg, doZscore)
% BIN_AND_AGGREGATE  Aggregate orientation contrast within each V1 polar-angle wedge.
%
%   [M, counts] = bin_and_aggregate(T, cfg, expCfg, doZscore)
%
% Applies the analysis inclusion filter (ecc in cfg.eccRange, pRF_r2 > cfg.r2min;
% V1 already enforced by LOAD_AND_FILTER), computes per-vertex contrasts, then takes
% cfg.aggregator ('mean', the manuscript choice, or 'median') across vertices within each
% polar-angle wedge, separately per subject. Each subject's wedge values are then scaled
% by their pRF-gain factor (OBSERVER_GAIN_WEIGHTS), matching the manuscript route.
%
% M      : nOri x nPA x nSubj, ordering follows expCfg.oriCols, cfg.paBins, and
%          SUBJECTS_FOR(cfg, expCfg) -- the observer list for THIS experiment, which
%          since 2026-08-27 is not the same for dg and da. The subject dimension of a
%          dg result and a da result therefore need not be the same length or the same
%          people: pass both through ASSERT_SAME_OBSERVERS before pairing them.
% counts : nPA x nSubj, number of vertices contributing to each wedge/subject.

    % Recorded BEFORE the inclusion filter, so the check below can tell "this observer
    % is not in the table at all" apart from "this observer is in the table but no
    % vertex survived ecc/R2". They need different fixes and must not report the same.
    inTable = unique(cellstr(T.subject));

    keep = T.pRF_ecc >= cfg.eccRange(1) & T.pRF_ecc <= cfg.eccRange(2) & T.pRF_r2 > cfg.r2min;
    T = T(keep, :);

    C    = compute_vertex_contrasts(T, expCfg, doZscore);
    subj = string(T.subject);
    pab  = T.pRF_angle_bin;

    % Resolve the observer list for THIS experiment, then rebind cfg.subjects on the
    % local copy so everything downstream that still reads cfg.subjects -- notably
    % OBSERVER_GAIN_WEIGHTS, whose scale vector is positional -- stays aligned with the
    % subject dimension built here. Rebinding rather than threading an extra argument
    % is deliberate: a gain vector silently ordered by a different list than M is
    % exactly the failure this split exists to remove.
    subs = subjects_for(cfg, expCfg);
    cfg.subjects = subs;

    absent = subs(~ismember(subs, inTable));
    if ~isempty(absent)
        error('bin_and_aggregate:subjectsNotInTable', ...
              ['%s observers requested but not present in the data table at all: %s\n' ...
               'Support/allsubjectsTable.csv holds only the original 8 observers, so ' ...
               'the 13-observer dg set cannot be aggregated until that table is ' ...
               'rebuilt (AnalysisCode/01_process_singlesubjectGLM/createTables.m). ' ...
               'Without this check they would aggregate to silent NaN rows.'], ...
              upper(expCfg.name), strjoin(absent, ', '));
    end

    emptied = subs(ismember(subs, inTable) & ~ismember(subs, unique(cellstr(subj))));
    if ~isempty(emptied)
        error('bin_and_aggregate:subjectsFilteredOut', ...
              ['%s observers are in the data table but no vertex survived the ' ...
               'inclusion filter (ecc %g-%g deg, pRF R2 > %g): %s\n' ...
               'That is a data or filter problem, not a subject-list problem.'], ...
              upper(expCfg.name), cfg.eccRange(1), cfg.eccRange(2), cfg.r2min, ...
              strjoin(emptied, ', '));
    end

    nO = numel(expCfg.oriCols);
    nP = numel(cfg.paBins);
    nS = numel(subs);

    M      = nan(nO, nP, nS);
    counts = zeros(nP, nS);
    for si = 1:nS
        inSubj = subj == subs{si};
        for pi = 1:nP
            idx = inSubj & (pab == cfg.paBins(pi));
            counts(pi, si) = nnz(idx);
            if any(idx)
                if strcmpi(cfg.aggregator, 'median')
                    M(:, pi, si) = median(C(idx, :), 1);
                else
                    M(:, pi, si) = mean(C(idx, :), 1);
                end
            end
        end
    end

    % Observer gain rescaling: divide by each observer's own pRF gain, multiply the
    % group gain back in. A per-observer scalar, so it leaves every within-observer
    % quantity (split-half correlations, variance ratios) untouched.
    scale = observer_gain_weights(cfg);
    for si = 1:nS
        M(:, :, si) = M(:, :, si) * scale(si);
    end
end
