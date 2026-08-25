function R = fit_cell_meta(d, obs, wedge, S, varargin)
% FIT_CELL_META  Random-effects meta-analysis at the (observer x wedge) cell level.
%
%   R = fit_cell_meta(d, obs, wedge, S)
%   R = fit_cell_meta(..., 'tau2', [tau2obs tau2cell], 'wedgeEffects', false)
%
% Fits, for ONE asymmetry in ONE map and experiment,
%
%       d_ip = mu + alpha_p + u_i + v_ip + e_ip
%
%   alpha_p   wedge (polar-angle) effects, FIXED, sum-to-zero coded. The polar-angle
%             profile is systematic, not a draw from a population, so it is not a
%             random effect; and sum-to-zero coding is what makes mu the average over
%             the FULL set of wedge levels rather than over whichever cells survived.
%   u_i       ~ N(0, tau2obs), observer offset, shared across that observer's wedges
%   v_ip      ~ N(0, tau2cell), observer-by-wedge heterogeneity of the TRUE effect
%   e_ip      ~ N(0, S), sampling error with KNOWN covariance, from resampling runs
%
% mu is the estimand: the asymmetry averaged over the full wedge set, for the average
% observer. With complete data it is the ordinary grand mean; with empty cells it is
% the least-squares mean, recovered from the incomplete design instead of silently
% redefined as "the mean over whichever wedges this observer happened to have". That
% is the entire reason for fitting a model here -- see ../METHOD_DECISIONS.md section 4 and the
% note there that the wedge-average and pooled routes diverge as soon as balance goes.
%
% S IS KNOWN AND NOT DIAGONAL. The per-wedge asymmetries of one observer are computed
% from the SAME runs, so their sampling errors are correlated; the run bootstrap in
% DIAGNOSE_WITHIN_OBSERVER_ERROR returns the whole nP x nP covariance and all of it is
% used. Passing only the diagonal would understate the correlated part and overstate
% the precision of mu. This known-covariance structure is exactly what fitlme cannot
% express (METHOD_DECISIONS.md section 4), which is why the fit is done here
% rather than with it.
%
% Inputs
%   d      nCell x 1 cell-level asymmetries.
%   obs    nCell x 1 observer index (any integers; only grouping matters).
%   wedge  nCell x 1 wedge index (any integers; the distinct values PRESENT ANYWHERE
%          in the input define the level set that mu averages over -- see the
%          structural-gap warning below).
%   S      nCell x nCell known sampling covariance (block diagonal by observer), or
%          nCell x 1 of variances for the diagonal-only case.
%
% Options
%   'tau2'         [tau2obs tau2cell] to FIX the variance components instead of
%                  estimating them. Used by VALIDATE_CELL_META to make the reduction
%                  to the observer-level estimator an exact algebraic identity.
%   'wedgeEffects' true (default). false drops alpha_p, which is the wrong model
%                  whenever cells are missing; provided only for comparison.
%   'df'           degrees of freedom for the interval on mu. Default nObs-1.
%                  NOT the cell count: the row count is 8x the observer count and
%                  using it repeats the error diagnosed in METHOD_DECISIONS.md section 3.
%
% STRUCTURAL GAPS. A wedge absent for EVERY observer cannot contribute to mu -- alpha_p
% is not estimable there and no model recovers it. Such levels are dropped and named in
% R.droppedWedges, and mu is then the average over the REMAINING wedges. That is a
% different estimand from V1's, so any cross-map comparison must be recomputed on a
% matched wedge set. Sampling gaps (a wedge empty for some observers) are the case this
% function is for; structural gaps (a map that does not represent that part of the
% field) are not repairable here.
%
% Returns R with mu, se, ci, p, tau2obs, tau2cell, alpha, df, nObs, nCell,
% wedgeLevels, droppedWedges, and R.wObs, the effective weight each observer carries
% in mu (the summed GLS weights of that observer's rows, normalised).

    p = inputParser;
    p.addParameter('tau2', [], @(x) isempty(x) || numel(x)==2);
    p.addParameter('wedgeEffects', true, @islogical);
    p.addParameter('df', [], @(x) isempty(x) || isscalar(x));
    p.parse(varargin{:});
    opt = p.Results;

    d = d(:);  obs = obs(:);  wedge = wedge(:);
    ok = isfinite(d);
    if isvector(S) && numel(S) == numel(d), ok = ok & isfinite(S(:)); end
    d = d(ok);  obs = obs(ok);  wedge = wedge(ok);
    if isvector(S) && numel(ok) == numel(S), S = diag(S(ok));
    else,                                    S = S(ok, ok);
    end
    n = numel(d);
    assert(n > 0, 'fit_cell_meta:empty', 'no finite cells');

    [uo, ~, oi] = unique(obs);
    [uw, ~, wi] = unique(wedge);
    nO = numel(uo);  nW = numel(uw);

    % --- fixed effects: intercept + sum-to-zero wedge coding ---------------------
    X = ones(n, 1);
    if opt.wedgeEffects && nW > 1
        C = zeros(n, nW-1);
        for j = 1:nW-1, C(:,j) = (wi == j) - (wi == nW); end
        X = [X, C];
    end
    Z = sparse(1:n, oi, 1, n, nO);
    ZZt = full(Z*Z.');

    % --- variance components -----------------------------------------------------
    if isempty(opt.tau2)
        % Start from a method-of-moments split of the excess over the known sampling
        % variance, the same subtraction the observer-level route uses.
        excess = max(0, var(d, 0) - mean(diag(S)));
        th0    = [excess/2, excess/2];
        lo     = [0 0];
        hi     = [Inf Inf];
        f      = @(t) reml_deviance(t, d, X, ZZt, S);
        o      = optimoptions('fmincon', 'Display', 'off', ...
                              'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-10);
        th     = fmincon(f, th0, [], [], [], [], lo, hi, [], o);
        % Restart from the boundary: the REML surface for two components on this little
        % data is often maximised at tau2 = 0, which an interior start can miss.
        for st = {[0 excess], [excess 0]}
            th2 = fmincon(f, st{1}, [], [], [], [], lo, hi, [], o);
            if f(th2) < f(th), th = th2; end
        end
    else
        th = opt.tau2;
    end

    [~, beta, XtViX, Vi] = reml_deviance(th, d, X, ZZt, S);   % full return, once
    covB = inv(XtViX);

    R.mu       = beta(1);
    R.se       = sqrt(covB(1,1));
    R.tau2obs  = th(1);
    R.tau2cell = th(2);
    R.alpha    = alpha_full(beta, nW, opt.wedgeEffects);
    R.nObs     = nO;   R.nCell = n;
    R.wedgeLevels   = uw(:).';
    R.droppedWedges = setdiff(unique(wedge(:).'), uw(:).');
    R.df = opt.df;  if isempty(R.df), R.df = nO - 1; end

    tcrit = tinv(0.975, R.df);
    R.ci  = [R.mu - tcrit*R.se, R.mu + tcrit*R.se];
    R.p   = 2 * tcdf(-abs(R.mu/R.se), R.df);

    % Effective observer weights: mu = a'd with a = e1' (X'V^-1 X)^-1 X'V^-1.
    a = (covB(1,:) * X.' * Vi).';
    R.wObs = accumarray(oi, a, [nO 1]);
    R.wObs = R.wObs / sum(R.wObs);
    R.obsLevels = uo(:).';
end

% ------------------------------------------------------------------------
function [dev, beta, XtViX, Vi] = reml_deviance(th, d, X, ZZt, S)
% -2 log REML, up to an additive constant. th = [tau2obs tau2cell], both >= 0.
% The explicit V inverse is only formed when the caller asks for it (once, at the end);
% inside the optimiser this is called for `dev` alone and stays at triangular solves.
    n = numel(d);
    V = max(th(1),0)*ZZt + max(th(2),0)*eye(n) + S;
    V = (V + V.')/2;
    [L, flag] = chol(V, 'lower');
    if flag, dev = 1e12; beta = nan(size(X,2),1); XtViX = []; Vi = []; return; end

    ViX   = L.'\(L\X);
    XtViX = X.'*ViX;  XtViX = (XtViX + XtViX.')/2;
    [Lx, fx] = chol(XtViX, 'lower');
    if fx, dev = 1e12; beta = nan(size(X,2),1); Vi = []; return; end

    beta = XtViX \ (ViX.'*d);
    r    = d - X*beta;
    dev  = 2*sum(log(diag(L))) + 2*sum(log(diag(Lx))) + r.'*(L.'\(L\r));
    if nargout > 3, Vi = L.'\(L\eye(n)); end
end

% ------------------------------------------------------------------------
function a = alpha_full(beta, nW, useWedge)
% Recover all nW sum-to-zero wedge effects from the nW-1 estimated contrasts.
    if ~useWedge || nW <= 1, a = zeros(1, nW); return; end
    a = beta(2:nW).';
    a = [a, -sum(a)];
end
