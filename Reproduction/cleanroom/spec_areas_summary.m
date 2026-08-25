function A = spec_areas_summary(varargin)
% SPEC_AREAS_SUMMARY  The extrastriate supplement, under the settled specification.
%
%   A = spec_areas_summary()                      % all eight maps, both bands
%   A = spec_areas_summary('maps',{'V1','V2','V3'})
%
% One analysis, applied unchanged to every map (../SPECIFICATION.md section 1). Runs
% SPEC_PROFILES + SPEC_TABLES per map x band, records coverage, applies the section 6
% reportability criterion, and tests the hierarchy trend.
%
% COVERAGE CRITERION (section 6), evaluated rather than assumed: at most 2 empty
% (observer x ROI) cells of 64, a median of at least 20 vertices per cell, and a
% maximum precision-weight ratio below 25. Maps that fail are still computed and
% still written out, flagged rather than dropped, so the excluded numbers are
% inspectable. Failing it is a statement about what THIS DESIGN can resolve by polar
% angle, not a claim that the map has no asymmetries.
%
% THE HIERARCHY TREND IS ALWAYS EQUAL-WEIGHTED, whatever `variant` asks for. Precision
% weighting it would need the within-observer covariance of the SAME observer's estimate
% in two different maps, and those come from the same runs, so their errors are
% correlated by an unmeasured amount. Assuming independence would understate the SE.
% Equal weighting needs no such assumption.
%
% THE HIERARCHY TREND. ../RESULTS.md section 5 records that the monotonic
% V1 -> V2 -> V3 decline rested on six individually significant cells falling in the
% same order, and that the trend itself had never been computed. It is computed here,
% WITHIN observer: for each observer, the V1-minus-V3 difference of the context
% effect, then a t test across the 8. Same observers in both maps, so the pairing is
% exact and the between-observer variance that limits everything else cancels.
%
% Returns A.cov, A.asym, A.context, A.trend, A.perObs.

    p = inputParser;
    p.addParameter('maps', {'V1','V2','V3','V3a','V3b','hV4','MT','MST'}, @iscell);
    p.addParameter('bands', {[4 8],[2 10]}, @iscell);
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('write', true, @islogical);
    p.addParameter('variant', 'spec', @ischar);
    p.parse(varargin{:});
    opt = p.Results;

    cfg = config_repro();
    Vr  = spec_variants(opt.variant);
    outDir = fullfile(cfg.reproDir, 'supplement');
    covRows = {};  asym = [];  ctx = [];  perObs = [];
    store = struct();

    for mi = 1:numel(opt.maps)
        for bi = 1:numel(opt.bands)
            m = opt.maps{mi};  b = opt.bands{bi};
            bandLbl = sprintf('%g-%g', b(1), b(2));
            fprintf('\n--- %s, %s deg ---\n', m, bandLbl);
            lastwarn('');
            try
                S = spec_profiles('area', m, 'eccRange', b, 'root', opt.root, ...
                                  'route', Vr.route, 'verify', false);
            catch ME
                fprintf('  skipped: %s\n', ME.message);
                continue
            end
            % A map with no per-map gain would be scored with a V1-derived scalar,
            % which is the V1 special-casing the specification set out to remove
            % (../SPECIFICATION.md section 4). Record it, and do not report the map.
            [~, wid] = lastwarn();
            gainFellBack = strcmp(wid, 'spec_profiles:gainFallback') || ...
                           strcmp(wid, 'observer_gain_weights:missing');

            % --- coverage -----------------------------------------------------
            % Counted over ONE experiment's 64 (observer x ROI) cells, which is what
            % the section 6 criterion is written against. dg and da give identical
            % counts by construction -- the inclusion mask is built from eccentricity
            % and pRF vexpl, neither of which depends on the experiment -- and that
            % identity is asserted rather than assumed.
            assert(isequal(S.dg.nVert, S.da.nVert), 'spec_areas_summary:nVert', ...
                   '%s %s: dg and da select different vertices.', m, bandLbl);
            nV = S.dg.nVert;                          % nSubj x nPA
            nEmpty = nnz(nV == 0);
            medV   = median(nV(:));
            allNaN = all(~isfinite(S.dg.asym(:))) && all(~isfinite(S.da.asym(:)));

            T = [];
            if ~allNaN
                % Hand over the profiles just computed: recomputing them here would
                % repeat the run bootstrap, which is most of this pipeline's cost.
                T = spec_tables('area', m, 'eccRange', b, 'root', opt.root, ...
                                'write', false, 'quiet', true, 'variant', opt.variant, ...
                                'profiles', S);
                wr = max([T.asym.weight_ratio; T.context.weight_ratio]);
            else
                wr = Inf;
            end
            reportable = ~allNaN && nEmpty <= 2 && medV >= 20 && wr < 25 && ~gainFellBack;

            covRows(end+1,:) = {m, bandLbl, nEmpty, medV, wr, gainFellBack, reportable}; %#ok<AGROW>
            fprintf(['  empty cells %d/%d, median vertices/cell %g, max weight ratio %.1f, ' ...
                     'per-map gain %s -> %s\n'], nEmpty, numel(nV), medV, wr, ...
                     tern(gainFellBack,'MISSING','ok'), ...
                     tern(reportable,'REPORTABLE','not reportable'));

            if allNaN, continue; end
            T.asym.band    = repmat({bandLbl}, height(T.asym), 1);
            T.context.band = repmat({bandLbl}, height(T.context), 1);
            T.perObs.band  = repmat({bandLbl}, height(T.perObs), 1);
            T.asym.reportable    = repmat(reportable, height(T.asym), 1);
            T.context.reportable = repmat(reportable, height(T.context), 1);
            asym   = [asym;   T.asym];    %#ok<AGROW>
            ctx    = [ctx;    T.context]; %#ok<AGROW>
            perObs = [perObs; T.perObs];  %#ok<AGROW>
            store.(sprintf('%s_%d', matlab.lang.makeValidName(m), bi)) = S;
        end
    end

    A.cov = cell2table(covRows, 'VariableNames', ...
        {'map','band','empty_cells','median_vertices_per_cell','max_weight_ratio', ...
         'gain_fell_back','reportable'});
    A.asym = asym;  A.context = ctx;  A.perObs = perObs;  A.store = store;
    A.trend = hierarchy_trend(store, Vr);   % Vr ignored: see the function
    A.variant = Vr;

    report(A);
    if opt.write
        if ~isfolder(outDir), mkdir(outDir); end
        sfx = Vr.tag;
        writetable(A.cov,     fullfile(outDir, sprintf('spec_areas_coverage_%s.csv', sfx)));
        writetable(A.asym,    fullfile(outDir, sprintf('spec_areas_asymmetries_%s.csv', sfx)));
        writetable(A.context, fullfile(outDir, sprintf('spec_areas_context_%s.csv', sfx)));
        writetable(A.perObs,  fullfile(outDir, sprintf('spec_areas_perobserver_%s.csv', sfx)));
        writetable(A.trend,   fullfile(outDir, sprintf('spec_areas_trend_%s.csv', sfx)));
        fprintf('\nspec_areas_summary: wrote 5 CSVs to %s\n', outDir);
    end
end

% ------------------------------------------------------------------------
function Tr = hierarchy_trend(store, ~)
% The variant is accepted and ignored on purpose: the trend is ALWAYS equal-weighted,
% for the reason given in the header. Taking the argument keeps the call site honest
% about what was requested; ignoring it keeps the answer defensible.
% Within-observer V1 - V3 difference of the context effect, at 4-8 deg. Paired: the
% same 8 observers contribute to both maps, so this is a one-sample t on the
% per-observer difference of differences, not a comparison of two group means.
    nm = {'horiz-vert','card-obl','rad-tang','polc-polo'};
    rows = {};
    have = @(f) isfield(store, f);
    if ~(have('V1_1') && have('V3_1')), Tr = cell2table(cell(0,8), 'VariableNames', ...
        {'comparison','asymmetry','n','mean','t_lo','t_hi','t_p','obs_agree'}); return; end

    pairs = {'V1','V2'; 'V2','V3'; 'V1','V3'};
    for pi = 1:size(pairs,1)
        fa = [pairs{pi,1} '_1'];  fb = [pairs{pi,2} '_1'];
        if ~(have(fa) && have(fb)), continue; end
        Sa = store.(fa);  Sb = store.(fb);
        Ca = Sa.dg.asym - Sa.da.asym;         % context effect, per observer
        Cb = Sb.dg.asym - Sb.da.asym;
        for j = 1:4
            d = Ca(:,j) - Cb(:,j);
            if nnz(isfinite(d)) < 3, continue; end
            G = spec_group(d, zeros(numel(d),1), 'equal');   % see the note in the header
            if ~isfinite(G.mean), continue; end
            rows(end+1,:) = {sprintf('%s - %s', pairs{pi,1}, pairs{pi,2}), nm{j}, G.n, ...
                             G.mean, G.lo, G.hi, G.p, G.obsAgree}; %#ok<AGROW>
        end
    end
    Tr = cell2table(rows, 'VariableNames', ...
        {'comparison','asymmetry','n','mean','t_lo','t_hi','t_p','obs_agree'});
end

% ------------------------------------------------------------------------
function report(A)
    fprintf('\n%s\nCOVERAGE AND REPORTABILITY\n%s\n', repmat('=',1,78), repmat('=',1,78));
    fprintf('%-5s %-7s %11s %10s %12s %10s  %s\n','map','band','empty/64', ...
            'med vert','max w ratio','map gain','reportable');
    for r = 1:height(A.cov)
        fprintf('%-5s %-7s %11d %10g %12.1f %10s  %s\n', A.cov.map{r}, A.cov.band{r}, ...
            A.cov.empty_cells(r), A.cov.median_vertices_per_cell(r), ...
            A.cov.max_weight_ratio(r), tern(A.cov.gain_fell_back(r),'MISSING','ok'), ...
            tern(A.cov.reportable(r),'yes','no'));
    end
    fprintf(['\nCriterion (../SPECIFICATION.md section 7): at most 2 empty (observer x ROI)\n' ...
             'cells of 64, median at least 20 vertices per cell, max precision-weight ratio\n' ...
             'below 25, and a per-map gain that exists. Failing it says this DESIGN cannot\n' ...
             'resolve that map by polar angle -- not that the map has no asymmetries.\n']);

    fprintf('\n%s\nCONTEXT EFFECTS (dg - da), reportable maps only\n%s\n', ...
            repmat('=',1,78), repmat('=',1,78));
    C = A.context(A.context.reportable, :);
    fprintf('%-5s %-7s %-11s %8s %19s %8s %5s\n','map','band','asymmetry','mean','t 95% CI','t p','obs');
    for r = 1:height(C)
        fprintf('%-5s %-7s %-11s %8.3f  [%7.3f %7.3f] %8s %3d/%d\n', C.area{r}, C.band{r}, ...
            C.asymmetry{r}, C.mean(r), C.t_lo(r), C.t_hi(r), pstr(C.t_p(r)), C.obs_agree(r), C.n(r));
    end

    if ~isempty(A.trend)
        fprintf('\n%s\nHIERARCHY TREND, within observer, 4-8 deg\n%s\n', ...
                repmat('=',1,78), repmat('=',1,78));
        fprintf(['Per-observer difference of the CONTEXT EFFECT between two maps, reported as\n' ...
                 '(earlier map) - (later map), then a t test across the 8 observers. Same\n' ...
                 'observers in both maps, so the pairing is exact and the between-observer\n' ...
                 'variance that limits everything else cancels.\n\n' ...
                 'The two Cartesian context effects are NEGATIVE, so a negative difference here\n' ...
                 'means the effect is larger IN MAGNITUDE in the earlier map -- i.e. it\n' ...
                 'attenuates up the hierarchy. 12 tests, uncorrected.\n\n']);
        fprintf('%-10s %-11s %8s %19s %8s %5s\n','comparison','asymmetry','mean','t 95% CI','t p','obs');
        for r = 1:height(A.trend)
            fprintf('%-10s %-11s %8.3f  [%7.3f %7.3f] %8s %3d/%d\n', A.trend.comparison{r}, ...
                A.trend.asymmetry{r}, A.trend.mean(r), A.trend.t_lo(r), A.trend.t_hi(r), ...
                pstr(A.trend.t_p(r)), A.trend.obs_agree(r), A.trend.n(r));
        end
    end
end

function s = pstr(p)
    if p < 0.001, s = '<.001'; else, s = sprintf('%.3f', p); end
end
function s = tern(c,a,b)
    if c, s = a; else, s = b; end
end
