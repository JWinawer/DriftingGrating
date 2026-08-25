function R = diagnose_cell_loss(varargin)
% DIAGNOSE_CELL_LOSS  What an extrastriate map's missing cells would do to V1.
%
%   R = diagnose_cell_loss()                                   % MT's 4-8 deg pattern
%   R = diagnose_cell_loss('donor','MT','eccRange',[2 10])
%   R = diagnose_cell_loss('donor','hV4','nRand',50)
%
% THE QUESTION. An extrastriate map has 64 (observer x polar-angle ROI) cells and
% some hold no vertex that passes the inclusion filter. Take those same cells in V1,
% where every one of them is populated, delete them, and re-run the specification.
% How much do the results move? V1 is the right test bed because its answer is known
% on complete data, so any change is attributable to the loss and nothing else.
%
% ../SPECIFICATION.md section 3 argues that cell loss damages the ROI-average route
% and not the harmonic one, on the grounds that the harmonic route never bins. This
% measures that rather than asserting it: both routes are run on identical deletions.
%
% THREE CONDITIONS, all on the same V1 vertices:
%   baseline   nothing deleted
%   donor      the donor map's own empty cells, deleted from V1
%   null       nRand random deletions of the SAME SIZE, preserving each observer's
%              number of lost cells but randomising WHICH ROIs are lost
% The null is what separates the two things the donor pattern confounds: how much
% is lost, and where. Only the second is what section 3 warns about.
%
% MATCHED OBSERVER SET. If deletion empties an observer entirely, SPEC_GROUP's
% convention makes the whole summary NaN -- correctly, since an estimate over 7 of 8
% is not comparable with one over 8. That convention would also make baseline and
% donor incomparable here, so every summary below is computed over the observers
% finite in BOTH, with n reported. The count is the finding, not a nuisance: it is
% reported in R.nCommon and printed.
%
% NBOOT. The run bootstrap only feeds sigma, which only the precision weighting uses,
% and the specification's primary is equal weighting with a t interval. Nothing
% reported here depends on it, so this defaults to 50 rather than the specification's
% 500 and the sweep costs minutes instead of half an hour. Raise it if you want the
% precision-weighted columns of R.detail to be trustworthy.
%
% Returns R with .summary (one row per quantity per route), .cells, .perObs,
% .nCommon, .drop and .donorN.

    p = inputParser;
    p.addParameter('area', 'V1', @ischar);          % the map the deletion is applied TO
    p.addParameter('donor', 'MT', @ischar);         % the map the pattern comes FROM
    p.addParameter('eccRange', [4 8], @(x) numel(x)==2);
    p.addParameter('root', dg_collect_dir(), @ischar);
    p.addParameter('nRand', 20, @isscalar);
    p.addParameter('nBoot', 50, @isscalar);
    p.addParameter('routes', {'harmonic','roi'}, @iscell);
    p.addParameter('seed', 7, @isscalar);
    % The mechanism check: delete a whole POLAR-ANGLE CLASS from every observer and
    % see which asymmetry moves. It is what turns "the results changed" into "they
    % changed because the wedges no longer balance", and it costs six extra fits.
    p.addParameter('classTest', true, @islogical);
    p.parse(varargin{:});
    opt = p.Results;

    cfg  = config_repro();
    nS   = numel(cfg.subjects);  nP = numel(cfg.paBins);
    expn = {'dg','da'};
    nm   = {'horiz-vert','card-obl','rad-tang','polc-polo'};

    % --- the deletion pattern -------------------------------------------------
    donorN = cell_occupancy('area', opt.donor, 'eccRange', opt.eccRange, 'root', opt.root);
    drop   = donorN == 0;
    nDrop  = nnz(drop);
    baseN  = cell_occupancy('area', opt.area,  'eccRange', opt.eccRange, 'root', opt.root);

    fprintf('\n%s\n', repmat('=',1,86));
    fprintf('CELL LOSS: %s''s empty cells, deleted from %s   (%g-%g deg)\n', ...
            opt.donor, opt.area, opt.eccRange(1), opt.eccRange(2));
    fprintf('%s\n', repmat('=',1,86));
    fprintf('%s has %d empty cells of %d; %s has %d.\n', ...
            opt.donor, nDrop, numel(drop), opt.area, nnz(baseN == 0));
    fprintf('\ndeleted cells (x), and the %s vertices each one destroys:\n', opt.area);
    fprintf('%-15s', 'observer');  fprintf('%7g', cfg.paBins);
    fprintf('%9s %9s\n', 'lost', 'of');
    for si = 1:nS
        fprintf('%-15s', cfg.subjects{si});
        for pIdx = 1:nP
            if drop(si,pIdx), fprintf('%7s', sprintf('x%d', baseN(si,pIdx)));
            else,             fprintf('%7s', '.');
            end
        end
        fprintf('%9d %9d\n', sum(baseN(si,drop(si,:))), sum(baseN(si,:)));
    end
    fprintf('%-15s', 'cells lost');  fprintf('%7d', sum(drop,1));
    fprintf('%9d %9d\n', sum(baseN(drop)), sum(baseN(:)));
    fprintf(['\nThe column totals are the point: loss concentrated in a few ROIs is what\n' ...
             'section 3 warns about, loss spread evenly is not.\n']);

    % --- the runs -------------------------------------------------------------
    rows = {};  R.perObs = struct();  R.nCommon = struct();
    for ri = 1:numel(opt.routes)
        route = opt.routes{ri};
        fprintf('\n--- route %s: baseline ---\n', route);
        Sb = runOne(opt, route, false(nS,nP));
        fprintf('--- route %s: %s pattern ---\n', route, opt.donor);
        Sd = runOne(opt, route, drop);

        % Null: same number of lost cells per observer, random ROIs.
        % Patterns are generated UP FRONT, before any fitting call. The fitting path
        % reseeds the global stream (DIAGNOSE_WITHIN_OBSERVER_ERROR does rng(si) for
        % its run bootstrap), so drawing them inside the loop would silently hand
        % every iteration the same pattern -- the same trap ../SPECIFICATION.md
        % section 9 records for the bootstrap draws.
        rng(opt.seed);
        kPer = sum(drop, 2);
        nulls = cell(opt.nRand, 1);
        for b = 1:opt.nRand
            dr = false(nS,nP);
            for si = 1:nS
                if kPer(si) > 0, dr(si, randperm(nP, kPer(si))) = true; end
            end
            nulls{b} = dr;
        end

        Dn = nan(opt.nRand, 12);                      % 8 asymmetries + 4 context
        for b = 1:opt.nRand
            fprintf('--- route %s: null %d/%d ---\n', route, b, opt.nRand);
            Sr = runOne(opt, route, nulls{b});
            Dn(b,:) = deltaVec(Sb, Sr);
        end

        dDon = deltaVec(Sb, Sd);
        q = quantities(Sb, Sd, expn, nm);
        for j = 1:numel(q)
            rows(end+1,:) = {opt.area, opt.donor, route, q(j).exp, q(j).name, ...
                             q(j).nCommon, q(j).base, q(j).baseP, ...
                             q(j).drop, q(j).dropP, dDon(j), ...
                             mean(abs(Dn(:,j)),'omitnan'), ...
                             prctile(abs(Dn(:,j)), 95), q(j).maxObs}; %#ok<AGROW>
        end
        R.perObs.(route)  = struct('base', {{Sb.dg.asym, Sb.da.asym}}, ...
                                   'drop', {{Sd.dg.asym, Sd.da.asym}});
        R.nCommon.(route) = [q.nCommon];
    end

    R.summary = cell2table(rows, 'VariableNames', ...
        {'area','donor','route','experiment','asymmetry','n', ...
         'base_mean','base_p','drop_mean','drop_p','delta', ...
         'null_mean_abs_delta','null_p95_abs_delta','max_observer_delta'});
    R.drop = drop;  R.donorN = donorN;  R.baseN = baseN;  R.nDrop = nDrop;
    R.eccRange = opt.eccRange;  R.subjects = cfg.subjects;  R.area = opt.area;

    report(R, opt);
    if opt.classTest, R.classTest = class_test(opt, cfg, nS, nm); end
end

% ------------------------------------------------------------------------
function T = class_test(opt, cfg, nS, nm)
% WHY the donor pattern does what it does. Delete a whole polar-angle class from
% every observer and read which asymmetry moves.
%
% The ROI route's asymmetries are a mean over the eight wedges, and two of the four
% are measured in OPPOSITE Cartesian directions at different wedges: a radial grating
% is horizontal on the horizontal meridian and vertical on the vertical meridian, and
% polar-cardinal is exactly plus-or-minus cardinal with the sign alternating
% cardinal/oblique (../SPECIFICATION.md section 6). The mean is only the intended
% quantity while the classes balance. Delete one class and the large Cartesian
% asymmetry leaks straight into the polar-frame estimate.
%
% The harmonic route fits continuous thetaV and never forms that mean, so it has
% nothing to unbalance. This measures the difference rather than asserting it.
    cls = struct('vertical',   ismember(cfg.paBins, [90 270]), ...
                 'horizontal', ismember(cfg.paBins, [0 180]), ...
                 'oblique',    ismember(cfg.paBins, [45 135 225 315]));
    fn = fieldnames(cls);
    rows = {};
    fprintf('\n%s\nMECHANISM: delete a whole polar-angle class from EVERY observer\n%s\n', ...
            repmat('=',1,86), repmat('=',1,86));
    for ri = 1:numel(opt.routes)
        route = opt.routes{ri};
        fprintf('\n--- route: %s ---\n%-12s %-6s', route, 'deleted', 'exp');
        fprintf('%11s', nm{:});  fprintf('\n');
        S0 = runOne(opt, route, false(nS, numel(cfg.paBins)));
        printRow('(none)', S0, route, cfg.paBins);
        rows = [rows; mkRows('(none)', S0, route, nm)]; %#ok<AGROW>
        for k = 1:numel(fn)
            S = runOne(opt, route, repmat(cls.(fn{k}), nS, 1));
            printRow(fn{k}, S, route, cfg.paBins);
            rows = [rows; mkRows(fn{k}, S, route, nm)]; %#ok<AGROW>
        end
    end
    fprintf(['\nRead the two rad-tang columns and the two polc-polo columns down the ROI\n' ...
             'block: deleting the vertical meridian and deleting the horizontal meridian\n' ...
             'push them in OPPOSITE directions, which is the signature of an unbalanced\n' ...
             'mean rather than of lost precision.\n']);
    T = cell2table(rows, 'VariableNames', ...
        {'route','deleted','experiment','asymmetry','mean'});
end

function printRow(lbl, S, ~, ~)
    fprintf('%-12s %-6s', lbl, 'dg');  fprintf('%11.3f', mean(S.dg.asym,1,'omitnan'));
    fprintf('\n%-12s %-6s', '', 'da'); fprintf('%11.3f', mean(S.da.asym,1,'omitnan'));
    fprintf('\n');
end

function rows = mkRows(lbl, S, route, nm)
    rows = {};
    for en = {'dg','da'}
        m = mean(S.(en{1}).asym, 1, 'omitnan');
        for j = 1:4, rows(end+1,:) = {route, lbl, en{1}, nm{j}, m(j)}; end %#ok<AGROW>
    end
end

% ------------------------------------------------------------------------
function S = runOne(opt, route, dr)
    S = spec_profiles('area', opt.area, 'eccRange', opt.eccRange, 'root', opt.root, ...
                      'route', route, 'dropCells', dr, 'nBoot', opt.nBoot, ...
                      'verify', false);
end

% ------------------------------------------------------------------------
function q = quantities(Sb, Sd, expn, nm)
% The eight per-experiment asymmetries then the four context effects, each summarised
% over the observers finite in BOTH conditions. SPEC_GROUP does the combining, so
% there is still only one definition of the group estimator.
    q = struct('exp',{},'name',{},'nCommon',{},'base',{},'baseP',{}, ...
               'drop',{},'dropP',{},'maxObs',{});
    for ei = 1:2
        for j = 1:4
            q(end+1) = one(Sb.(expn{ei}).asym(:,j), Sd.(expn{ei}).asym(:,j), ...
                           Sb.(expn{ei}).sigma(:,j), Sd.(expn{ei}).sigma(:,j), ...
                           expn{ei}, nm{j}); %#ok<AGROW>
        end
    end
    Db = Sb.dg.asym - Sb.da.asym;   Dd = Sd.dg.asym - Sd.da.asym;
    sb = sqrt(Sb.dg.sigma.^2 + Sb.da.sigma.^2);
    sd = sqrt(Sd.dg.sigma.^2 + Sd.da.sigma.^2);
    for j = 1:4
        q(end+1) = one(Db(:,j), Dd(:,j), sb(:,j), sd(:,j), 'dg-da', nm{j}); %#ok<AGROW>
    end
end

function s = one(b, d, sgb, sgd, en, name)
    m  = isfinite(b) & isfinite(d);
    Gb = spec_group(b(m), sgb(m), 'equal');
    Gd = spec_group(d(m), sgd(m), 'equal');
    s = struct('exp', en, 'name', name, 'nCommon', nnz(m), ...
               'base', Gb.mean, 'baseP', Gb.p, 'drop', Gd.mean, 'dropP', Gd.p, ...
               'maxObs', max(abs(d(m) - b(m))));
end

function v = deltaVec(Sb, Sx)
    q = quantities(Sb, Sx, {'dg','da'}, {'1','2','3','4'});
    v = [q.drop] - [q.base];
end

% ------------------------------------------------------------------------
function report(R, opt)
    fprintf('\n%s\nWHAT THE DELETION DOES  (equal weighting, t test, %s %g-%g deg)\n%s\n', ...
            repmat('=',1,102), R.area, R.eccRange(1), R.eccRange(2), repmat('=',1,102));
    for ri = 1:numel(opt.routes)
        T = R.summary(strcmp(R.summary.route, opt.routes{ri}), :);
        fprintf('\n--- route: %s ---\n', opt.routes{ri});
        fprintf('%-6s %-11s %3s %9s %8s %9s %8s %9s %9s %9s\n', ...
                'exp','asymmetry','n','baseline','p','deleted','p','delta', ...
                'null|d|','null p95');
        for r = 1:height(T)
            fprintf('%-6s %-11s %3d %9.3f %8s %9.3f %8s %9.3f %9.3f %9.3f%s\n', ...
                T.experiment{r}, T.asymmetry{r}, T.n(r), T.base_mean(r), pstr(T.base_p(r)), ...
                T.drop_mean(r), pstr(T.drop_p(r)), T.delta(r), ...
                T.null_mean_abs_delta(r), T.null_p95_abs_delta(r), ...
                flag(T.base_p(r), T.drop_p(r)));
        end
    end
    fprintf(['\n<- marks a quantity that crosses p = .05 in one direction or the other.\n' ...
             '"null |d|" is the mean absolute shift under random deletions of the same\n' ...
             'size; "null p95" its 95th percentile. A delta inside the null band is the\n' ...
             'price of losing that many cells; one outside it is the price of losing\n' ...
             'THOSE cells.\n']);
end

function s = pstr(p)
    if ~isfinite(p),   s = 'NaN';
    elseif p < 0.001,  s = '<.001';
    else,              s = sprintf('%.3f', p);
    end
end
function s = flag(p1, p2)
    if isfinite(p1) && isfinite(p2) && (p1 < 0.05) ~= (p2 < 0.05), s = '  <-';
    else, s = '';
    end
end
