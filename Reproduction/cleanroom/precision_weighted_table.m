function P = precision_weighted_table(varargin)
% PRECISION_WEIGHTED_TABLE  The four asymmetries in both experiments, and the four
% context effects, estimated with and without precision (inverse-variance) weighting.
%
%   P = precision_weighted_table()
%   P = precision_weighted_table('csv', 'out.csv')
%
% This IS the mixed-effects answer to "shouldn't noisier observers count less?" -- a
% random-effects model in which the within-observer variances are MEASURED rather than
% inferred from a design that has no replication:
%
%   y_i        observer i's asymmetry (or dg-minus-da context difference), from the
%              wedge-average route (BIN_AND_AGGREGATE -> COMPUTE_ASYMMETRIES)
%   sigma_i^2  observer i's within-observer measurement variance, from resampling RUNS
%              (DIAGNOSE_WITHIN_OBSERVER_ERROR). Runs, not vertices: vertex responses are
%              not independent, and an earlier vertex-resampling estimate was withdrawn.
%   tau^2      between-observer variance of the TRUE effects -- the heterogeneity, i.e.
%              the random-effect variance component, one number for the sample. By method
%              of moments, tau^2 = max(0, var(y_i) - mean(sigma_i^2)).
%
%   w_i = 1/(tau^2 + sigma_i^2), normalised. Equal weighting is the w_i = 1/n special
%   case, exact in the limit where tau^2 dominates every sigma_i^2.
%
% Both estimates get a 95% CI on n-1 df; the equal-weighted one is the ordinary
% across-observer t interval, the weighted one uses SE = sqrt(sum(w_i^2 (tau^2+sigma_i^2))).
%
% The result of the table is that the two agree. Precision weighting is available and
% principled and changes nothing, because tau^2 is common to every observer and here
% exceeds the mean sigma_i^2 -- so a large spread in RELIABILITY becomes a small spread in
% WEIGHT. See ../LME.md section 5.
%
% Returns P.asym and P.context, each a table. Optional 'csv' writes them out.

    p = inputParser;
    p.addParameter('csv', '', @ischar);
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('area', 'V1', @ischar);
    p.addParameter('route', 'harmonic', @ischar);
    p.addParameter('thetaV', 'continuous', @ischar);
    p.addParameter('gain', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('weighting', 'equalcoverage', @ischar);
    p.addParameter('eccRange', [], @(x) isempty(x) || numel(x)==2);
    p.addParameter('quiet', false, @islogical);
    p.parse(varargin{:});
    opt = p.Results;

    W    = diagnose_within_observer_error('root', opt.root, 'area', opt.area, ...
                                          'eccRange', opt.eccRange, 'route', opt.route, ...
                                          'thetaV', opt.thetaV, 'gain', opt.gain, ...
                                          'weighting', opt.weighting, 'quiet', true);
    expn = {'dg','da'};
    nS   = numel(W.subjects);

    rows = {};
    for e = 1:2
        for a = 1:4
            rows(end+1,:) = est(W.full(:,a,e), W.seBoot(:,a,e).^2, ...
                                expn{e}, W.names{a}, nS);       %#ok<AGROW>
        end
    end
    P.asym = mk(rows);

    rows = {};
    for a = 1:4
        d  = W.full(:,a,1) - W.full(:,a,2);
        s2 = W.seBoot(:,a,1).^2 + W.seBoot(:,a,2).^2;
        rows(end+1,:) = est(d, s2, 'dg-da', W.names{a}, nS);     %#ok<AGROW>
    end
    P.context = mk(rows);

    if ~opt.quiet
        show('ASYMMETRIES, per experiment', P.asym);
        show('CONTEXT EFFECTS (dg - da)',   P.context);
        fprintf(['\nThe two weightings agree on every context effect and on all six Cartesian-frame\n' ...
                 'entries. They differ in status for da rad-tang, which is marginal by any route.\n' ...
                 'tau^2 is common to every observer, so a large spread in reliability compresses\n' ...
                 'into a small spread in weight; and tau^2 is estimated on %d df, so the weights\n' ...
                 'are themselves imprecise.\n'], nS-1);
    end

    if ~isempty(opt.csv)
        T = [P.asym; P.context];
        writetable(T, opt.csv);
        fprintf('wrote %s\n', opt.csv);
    end
end

% ------------------------------------------------------------------------
function r = est(y, s2, expLabel, name, nS)
    tau2 = max(0, var(y,0) - mean(s2));
    w    = 1./(tau2 + s2);  w = w/sum(w);
    tcrit = tinv(0.975, nS-1);

    mEq  = mean(y);          seEq = std(y)/sqrt(nS);
    mW   = sum(w.*y);        seW  = sqrt(sum(w.^2 .* (tau2 + s2)));
    [~, pEq] = ttest(y);
    pW  = 2*(1 - tcdf(abs(mW/seW), nS-1));

    r = {expLabel, name, mEq, mEq-tcrit*seEq, mEq+tcrit*seEq, pEq, ...
         mW,  mW-tcrit*seW,  mW+tcrit*seW,  pW, ...
         sqrt(tau2), sqrt(mean(s2)), max(w)/min(w)};
end

function T = mk(rows)
    T = cell2table(rows, 'VariableNames', {'experiment','asymmetry', ...
        'equal','equal_lo','equal_hi','equal_p', ...
        'weighted','weighted_lo','weighted_hi','weighted_p', ...
        'tau','mean_sigma','weight_ratio'});
end

function show(ttl, T)
    fprintf('\n%s\n%s\n', ttl, repmat('=',1,86));
    fprintf('%-6s %-11s %22s %22s %7s %7s %6s\n', ...
            'exp','asymmetry','equal-weighted [95% CI]','precision-w [95% CI]', ...
            'tau','sigma','w max/min');
    for i = 1:height(T)
        fprintf('%-6s %-11s %7.3f [%6.3f %6.3f] %7.3f [%6.3f %6.3f] %7.3f %7.3f %6.2fx\n', ...
            T.experiment{i}, T.asymmetry{i}, ...
            T.equal(i), T.equal_lo(i), T.equal_hi(i), ...
            T.weighted(i), T.weighted_lo(i), T.weighted_hi(i), ...
            T.tau(i), T.mean_sigma(i), T.weight_ratio(i));
    end
end
