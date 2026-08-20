function R = asymmetry_tables(varargin)
% ASYMMETRY_TABLES  The four asymmetries per experiment, and the four context effects,
% on the published route (cfg.aggregator = 'mean' + observer pRF-gain rescaling).
%
%   R = asymmetry_tables()
%   R = asymmetry_tables('csv', '../supplement/asymmetry_tables')
%
% Reports, per experiment and per asymmetry: the group mean, a bootstrap percentile CI
% over observers (1000 resamples -- the manuscript's Figs 5/6 method) and a t interval
% on 7 df, plus the per-observer and per-ROI consistency counts. Then the four context
% effects (dg minus da) with the same intervals and a paired t test.

    p = inputParser;
    p.addParameter('csv', '', @ischar);
    p.addParameter('nBoot', 1000, @isnumeric);
    p.parse(varargin{:});
    opt = p.Results;

    cfg  = config_repro();
    T    = load_and_filter(cfg);
    expn = {'dg','da'};
    nm   = {'horiz-vert','card-obl','rad-tang','polc-polo'};

    D = struct();     % per-observer wedge-averaged asymmetry, nSubj x 4 x 2
    P = struct();     % per-ROI, observer-averaged,           nPA   x 4 x 2
    for ei = 1:2
        M = bin_and_aggregate(T, cfg, cfg.(expn{ei}), false);   % nOri x nPA x nSubj
        A = compute_asymmetries(M, cfg, cfg.(expn{ei}));
        for j = 1:4
            d = A.(A.order{j}).diff;              % nPA x nSubj
            D.v(:,j,ei) = mean(d, 1, 'omitnan').';
            P.v(:,j,ei) = mean(d, 2, 'omitnan');
        end
    end
    R.perObserver = D.v;  R.perROI = P.v;  R.names = nm;  R.experiments = expn;
    R.aggregator  = cfg.aggregator;  R.gainMean = cfg.gainMean;

    bar = repmat('=',1,100);
    fprintf('\n%s\nFOUR ASYMMETRIES PER EXPERIMENT  (aggregator=%s, gain rescaling=%s mean)\n%s\n', ...
            bar, cfg.aggregator, cfg.gainMean, bar);
    fprintf('%-4s %-11s %8s %22s %22s %8s %9s %9s\n', 'exp','asymmetry','mean', ...
            'bootstrap 95% CI','t 95% CI (7 df)','t','p','obs/ROI');
    rows = {};
    for ei = 1:2
        for j = 1:4
            y  = D.v(:,j,ei);
            [ci, ~] = bootstrap_ci(y.', opt.nBoot, 95, 0);
            [~, pv, tci, st] = ttest(y);
            nObs = sum(sign(y) == sign(mean(y)));
            nROI = sum(sign(P.v(:,j,ei)) == sign(mean(y)));
            fprintf('%-4s %-11s %8.3f   [%7.3f, %7.3f]   [%7.3f, %7.3f] %8.2f %9.4f %6d/%d\n', ...
                    expn{ei}, nm{j}, mean(y), ci(1), ci(2), tci(1), tci(2), st.tstat, pv, nObs, nROI);
            rows(end+1,:) = {expn{ei}, nm{j}, mean(y), ci(1), ci(2), tci(1), tci(2), ...
                             st.tstat, pv, nObs, nROI}; %#ok<AGROW>
        end
    end
    R.asym = cell2table(rows, 'VariableNames', {'experiment','asymmetry','mean', ...
        'bootLo','bootHi','tLo','tHi','t','p','nObsSameSign','nROIsameSign'});

    fprintf('\n%s\nCONTEXT EFFECTS  (dg minus da), paired across observers\n%s\n', bar, bar);
    fprintf('%-11s %8s %22s %22s %8s %9s %8s\n', 'asymmetry','mean', ...
            'bootstrap 95% CI','t 95% CI (7 df)','t','p','n same');
    rows = {};
    for j = 1:4
        y  = D.v(:,j,1) - D.v(:,j,2);
        [ci, ~] = bootstrap_ci(y.', opt.nBoot, 95, 0);
        [~, pv, tci, st] = ttest(y);
        nObs = sum(sign(y) == sign(mean(y)));
        fprintf('%-11s %8.3f   [%7.3f, %7.3f]   [%7.3f, %7.3f] %8.2f %9.4f %8d\n', ...
                nm{j}, mean(y), ci(1), ci(2), tci(1), tci(2), st.tstat, pv, nObs);
        rows(end+1,:) = {nm{j}, mean(y), ci(1), ci(2), tci(1), tci(2), st.tstat, pv, nObs}; %#ok<AGROW>
    end
    R.context = cell2table(rows, 'VariableNames', {'asymmetry','mean','bootLo','bootHi', ...
        'tLo','tHi','t','p','nObsSameSign'});

    fprintf('\n%s\nPER-OBSERVER VALUES (wedge-averaged)\n%s\n', bar, bar);
    for ei = 1:2
        fprintf('\n%s:\n%-15s %10s %10s %10s %10s\n', expn{ei}, 'observer', nm{:});
        for si = 1:numel(cfg.subjects)
            fprintf('%-15s %10.3f %10.3f %10.3f %10.3f\n', cfg.subjects{si}, D.v(si,:,ei));
        end
    end
    fprintf('\n%s\nPER-ROI VALUES (observer-averaged), polar angle in cfg.paBins order\n%s\n', bar, bar);
    for ei = 1:2
        fprintf('\n%s:\n%-15s %10s %10s %10s %10s\n', expn{ei}, 'polar angle', nm{:});
        for pi = 1:numel(cfg.paBins)
            fprintf('%-15d %10.3f %10.3f %10.3f %10.3f\n', cfg.paBins(pi), P.v(pi,:,ei));
        end
    end

    if ~isempty(opt.csv)
        writetable(R.asym,    [opt.csv '_asymmetries.csv']);
        writetable(R.context, [opt.csv '_context.csv']);
        fprintf('\nwrote %s_asymmetries.csv and %s_context.csv\n', opt.csv, opt.csv);
    end
end
