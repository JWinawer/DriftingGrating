function R = spec_tables(varargin)
% SPEC_TABLES  The Figure 5/6 tables, under the settled specification.
%
%   R = spec_tables()                       % V1, 4-8 deg, writes CSV + markdown
%   R = spec_tables('area','V2','write',false)
%
% Produces, for one map:
%   R.asym      the four asymmetries in each experiment
%   R.context   the four context effects (dg - da), paired across observers
%   R.perObs    every observer's value, so nothing in the tables is unreproducible
%
% Each estimate is reported four ways, deliberately:
%   equal-weighted with a t interval on n-1 df   <- PRIMARY, and what the figures show
%   equal-weighted with a percentile bootstrap   <- what the supplement currently quotes
%   precision-weighted with a t interval         <- the robustness check
%   obs, the count of observers sharing the group sign
%
% WHY BOTH INTERVALS. They disagree for two polar cells, and only those two
% (../METHOD_DECISIONS.md section 5). Reporting one silently would hide a real open decision;
% reporting both makes it a choice the reader can see. At n = 8 the percentile
% bootstrap has poor coverage and does not account for uncertainty in the spread, so
% t is the primary interval here.

    p = inputParser;
    p.addParameter('area', 'V1', @ischar);
    p.addParameter('eccRange', [], @(x) isempty(x) || numel(x)==2);
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('nBoot', 10000, @isscalar);
    p.addParameter('write', true, @islogical);
    p.addParameter('outDir', '', @ischar);
    p.addParameter('quiet', false, @islogical);
    p.addParameter('variant', 'spec', @ischar);
    % Accept an already-computed SPEC_PROFILES result. The fit and the run bootstrap
    % behind it are the expensive part of this pipeline (~75% of it is the 500-resample
    % bootstrap over runs), and a caller that already holds S for this map and route
    % would otherwise pay for it twice. Numbers are identical either way -- the guard
    % below refuses an S that does not match what was asked for.
    p.addParameter('profiles', [], @(x) isempty(x) || isstruct(x));
    p.parse(varargin{:});
    opt = p.Results;

    cfg = config_repro();
    if ~isempty(opt.eccRange), cfg.eccRange = opt.eccRange; end
    if isempty(opt.outDir), opt.outDir = fullfile(cfg.reproDir, 'supplement'); end

    V = spec_variants(opt.variant);
    if isempty(opt.profiles)
        S = spec_profiles('area', opt.area, 'eccRange', opt.eccRange, 'root', opt.root, ...
                          'route', V.route);
    else
        S = opt.profiles;
        assert(strcmp(S.area, opt.area) && strcmp(S.route, V.route) && ...
               isequal(S.eccRange(:).', cfg.eccRange(:).'), 'spec_tables:profiles', ...
               ['supplied profiles are %s/%s/%g-%g but %s/%s/%g-%g was asked for. ' ...
                'Reusing the wrong ones would silently mislabel the table.'], ...
               S.area, S.route, S.eccRange(1), S.eccRange(2), ...
               opt.area, V.route, cfg.eccRange(1), cfg.eccRange(2));
    end
    R.variant = V.tag;  R.route = V.route;  R.weighting = V.weighting;

    nm   = S.names;
    expn = {'dg','da'};
    nS   = numel(S.subjects);

    % --- per-observer long table ------------------------------------------
    rows = {};
    for ei = 1:2
        A = S.(expn{ei}).asym;
        for si = 1:nS
            for j = 1:4
                rows(end+1,:) = {opt.area, expn{ei}, nm{j}, S.subjects{si}, ...
                                 A(si,j), S.gainScale(si), sum(S.(expn{ei}).nVert(si,:))}; %#ok<AGROW>
            end
        end
    end
    R.perObs = cell2table(rows, 'VariableNames', ...
        {'area','experiment','asymmetry','observer','value','gainScale','nVertex'});

    % --- asymmetries -------------------------------------------------------
    rows = {};
    for ei = 1:2
        A = S.(expn{ei}).asym;
        for j = 1:4
            rows(end+1,:) = summarise(A(:,j), S.(expn{ei}).sigma(:,j), opt.area, ...
                                      expn{ei}, nm{j}, V, opt.nBoot); %#ok<AGROW>
        end
    end
    R.asym = mkTable(rows);

    % --- context effects, paired within observer ---------------------------
    rows = {};
    D = S.dg.asym - S.da.asym;                 % the per-observer difference FIRST
    for j = 1:4
        sg = sqrt(S.dg.sigma(:,j).^2 + S.da.sigma(:,j).^2);   % independent sessions
        rows(end+1,:) = summarise(D(:,j), sg, opt.area, 'dg-da', nm{j}, V, opt.nBoot); %#ok<AGROW>
    end
    R.context = mkTable(rows);

    R.area = opt.area;  R.eccRange = cfg.eccRange;  R.subjects = S.subjects;

    if ~opt.quiet, show(R); end
    if opt.write, writeOut(R, opt); end
end

% ------------------------------------------------------------------------
function row = summarise(d, sigma, area, en, name, V, nBoot)
% Both weightings are always computed and both are reported; V.weighting only decides
% which one lands in the primary columns (est/lo/hi/p). Reporting one silently would
% hide that the choice is a choice.
%
% NAN CONVENTION, matching PRECISION_WEIGHTED_TABLE: if any observer is missing, every
% summary is NaN rather than an estimate over the observers that remain. An estimate
% silently computed over 7 of 8 is not comparable with one computed over 8, and in the
% sparse maps the missing observer is missing because the map has no vertices there --
% exactly the fact worth surfacing (../SPECIFICATION.md section 7).
    Ge = spec_group(d, sigma, 'equal');
    Gp = spec_group(d, sigma, 'precision');
    G  = Ge;  if strcmp(V.weighting,'precision'), G = Gp; end

    nFinite = nnz(isfinite(d));
    if nFinite < numel(d) || nFinite < 3
        row = {area, en, name, V.tag, V.route, V.weighting, nFinite, ...
               NaN, NaN, NaN, NaN, NaN, NaN, ...
               Ge.mean, Ge.lo, Ge.hi, Ge.p, Gp.mean, Gp.lo, Gp.hi, Gp.p, ...
               NaN, NaN, NaN, 0, false};
        return
    end

    % Percentile bootstrap over observers, the Figs 5/6 method, for comparison only.
    n = numel(d);
    rng(0);
    bs  = mean(d(randi(n, n, nBoot)), 1);
    bci = prctile(bs, [2.5 97.5]);

    row = {area, en, name, V.tag, V.route, V.weighting, n, ...
           G.mean, G.lo, G.hi, G.p, bci(1), bci(2), ...
           Ge.mean, Ge.lo, Ge.hi, Ge.p, Gp.mean, Gp.lo, Gp.hi, Gp.p, ...
           Ge.tau, Ge.meanSigma, Gp.weightRatio, G.obsAgree, ...
           (bci(1) > 0 || bci(2) < 0) ~= (G.lo > 0 || G.hi < 0)};
end

function T = mkTable(rows)
    T = cell2table(rows, 'VariableNames', {'area','experiment','asymmetry', ...
        'variant','route','weighting','n', ...
        'mean','t_lo','t_hi','t_p','boot_lo','boot_hi', ...
        'eq_mean','eq_lo','eq_hi','eq_p','pw_mean','pw_lo','pw_hi','pw_p', ...
        'tau','mean_sigma','weight_ratio','obs_agree','ci_methods_disagree'});
end

% ------------------------------------------------------------------------
function show(R)
    for which = {'asym','context'}
        T = R.(which{1});
        if strcmp(which{1},'asym')
            hdr = sprintf('ASYMMETRIES, per experiment  --  %s, %g-%g deg  [%s: %s, %s wt]', ...
                          R.area, R.eccRange(1), R.eccRange(2), R.variant, R.route, R.weighting);
        else
            hdr = sprintf('CONTEXT EFFECTS (dg - da), paired  --  %s, %g-%g deg  [%s: %s, %s wt]', ...
                          R.area, R.eccRange(1), R.eccRange(2), R.variant, R.route, R.weighting);
        end
        fprintf('\n%s\n%s\n%s\n', repmat('=',1,96), hdr, repmat('=',1,96));
        fprintf('%-6s %-11s %3s %8s %19s %9s %19s %6s %s\n', ...
                'exp','asymmetry','n','estimate','95% CI','p','the other weighting','p','obs');
        for r = 1:height(T)
            if ~isfinite(T.mean(r))
                fprintf('%-6s %-11s %3d %8s   %-19s %9s\n', T.experiment{r}, T.asymmetry{r}, ...
                        T.n(r), 'NaN', '(observer with no vertices)', '');
                continue
            end
            fprintf('%-6s %-11s %3d %8.3f  [%7.3f %7.3f] %9s  [%7.3f %7.3f] %6s  %d/%d%s\n', ...
                T.experiment{r}, T.asymmetry{r}, T.n(r), T.mean(r), T.t_lo(r), T.t_hi(r), ...
                pstr(T.t_p(r)), oth(T,r,'lo'), oth(T,r,'hi'), pstr(oth(T,r,'p')), ...
                T.obs_agree(r), T.n(r), tern(T.ci_methods_disagree(r), '  <- CI methods disagree', ''));
        end
    end
    d = R.asym.ci_methods_disagree | false;
    if any(d) || any(R.context.ci_methods_disagree)
        % n comes from the table rather than a literal: the observer count is no longer
        % fixed at 8, and a note that states the wrong n is worse than no note.
        nObs = max([R.asym.n(:); R.context.n(:)]);
        fprintf(['\nThe percentile bootstrap and the t interval disagree where marked. t is ' ...
                 'primary\nhere; at n = %d the percentile method has poor coverage. See ' ...
                 '../METHOD_DECISIONS.md section 5.\n'], nObs);
    end
end

function v = oth(T, r, f)
% The weighting NOT used for the primary columns, so both are always visible.
    if strcmp(T.weighting{r}, 'precision'), pre = 'eq_'; else, pre = 'pw_'; end
    v = T.([pre f])(r);
end

function s = pstr(p)
    if p < 0.001, s = '<.001'; else, s = sprintf('%.3f', p); end
end
function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

% ------------------------------------------------------------------------
function writeOut(R, opt)
    if ~isfolder(opt.outDir), mkdir(opt.outDir); end
    tag = sprintf('%s_%s_%g-%g', R.variant, lower(R.area), R.eccRange(1), R.eccRange(2));
    f1 = fullfile(opt.outDir, sprintf('spec_asymmetries_%s.csv', tag));
    f2 = fullfile(opt.outDir, sprintf('spec_context_%s.csv', tag));
    f3 = fullfile(opt.outDir, sprintf('spec_perobserver_%s.csv', tag));
    writetable(R.asym, f1);  writetable(R.context, f2);  writetable(R.perObs, f3);
    fprintf('\nspec_tables: wrote\n  %s\n  %s\n  %s\n', f1, f2, f3);
end
