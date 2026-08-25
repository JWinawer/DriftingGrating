function E = diagnose_prf_angle_error(varargin)
% DIAGNOSE_PRF_ANGLE_ERROR  Measure pRF polar-angle precision from two independent fits.
%
%   E = diagnose_prf_angle_error()
%
% REQUIRES the NYUAD volume mounted at /Volumes/Vision (thousands of miles away, so
% each .mgz read costs a second or two). Reads 4 files per subject.
%
% WHY THIS EXISTS. The per-vertex harmonic model (../HARMONIC_MODEL.md) regresses on
% each vertex's pRF polar angle thetaV. Error in thetaV attenuates whichever terms
% depend on it, and the two experiments depend on it through DIFFERENT terms:
%
%            dg (Cartesian)                  da (polar)
%   b1, b2   regressor is thetaV-free  ->    depends on thetaV -> ATTENUATED
%   b3, b4   depends on thetaV -> ATTEN.     regressor is thetaV-free
%
% So angle error deflates b1(da) but NOT b1(dg), which inflates the very gap the
% context claim rests on. A forward simulation of a context-free world reproduces the
% observed b1 ratio at about 40 deg of angle error. This function measures the error
% instead of assuming it.
%
% HOW. Every observer was fitted twice, independently, from two pRF runs with
% different stimuli:
%   prfvista_mov  - the solution the analysis actually uses. ret_<subject>.mat in
%                   ~/dg_collect IS this solution (verified: max|difference| = 0),
%                   so it is read locally and costs nothing.
%   prfvista      - the independent replicate, read from the server.
% The disagreement between them estimates the measurement error. Treating the two as
% independent with equal variance, sigma_single = circSD(difference)/sqrt(2).
%
% This is an UPPER BOUND on sigma. The two stimuli differ (moving/flickering vs not),
% so any genuine stimulus-driven difference in the retinotopic estimate inflates the
% disagreement beyond pure noise. An upper bound is what the argument needs: if even
% the bound is far below 40 deg, angle error cannot manufacture the result.
%
% Vertices: V1_REmanual, 4-8 deg eccentricity, vexpl > 0.1 -- the manuscript inclusion
% criteria -- and by default vexpl > 0.1 in the REPLICATE too, since a vertex the
% replicate failed to fit tells us nothing about precision.
%
% Options: 'root' (~/dg_collect), 'deriv' (server derivatives), 'veThresh' (0.1),
%          'eccRange' ([4 8]), 'requireBothVE' (true), 'verbose' (true).
%
% Returns E with per-subject and pooled circular statistics, and the implied
% second- and fourth-harmonic attenuation factors lambda2 = exp(-2 sigma^2),
% lambda4 = exp(-8 sigma^2) (sigma in radians).

    p = inputParser;
    p.addParameter('root',  '/Users/jaw288/dg_collect', @ischar);
    p.addParameter('deriv', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
                             'data_bids/derivatives'], @ischar);
    p.addParameter('veThresh', 0.10, @isnumeric);
    p.addParameter('eccRange', [4 8], @isnumeric);
    p.addParameter('requireBothVE', true, @(z) islogical(z) || isnumeric(z));
    p.addParameter('verbose', true, @(z) islogical(z) || isnumeric(z));
    p.parse(varargin{:});
    opt = p.Results;

    ensure_mriread();
    cfg  = config_repro();
    subs = cfg.subjects;

    rows = {};
    allD = [];
    for ii = 1:numel(subs)
        s = subs{ii};
        t0 = tic;

        Rt = load(fullfile(opt.root, sprintf('ret_%s.mat', s)), ...
                  'angle_adj','eccen','vexpl','nLH','nRH');
        Lb = load(fullfile(opt.root, sprintf('labels_%s.mat', s)), ...
                  'lh_V1_REmanual','rh_V1_REmanual');

        nLH = double(Rt.nLH);
        v1  = [double(Lb.lh_V1_REmanual(:)); double(Lb.rh_V1_REmanual(:)) + nLH];

        [aRep, veRep] = read_replicate(opt.deriv, s, opt.root);
        n = numel(Rt.angle_adj);
        if numel(aRep) ~= n
            error('diagnose_prf_angle_error:size', ...
                  '%s: replicate has %d vertices, analysis solution has %d.', ...
                  s, numel(aRep), n);
        end

        keep = false(n,1);  keep(v1) = true;
        keep = keep & double(Rt.eccen(:)) >= opt.eccRange(1) ...
                    & double(Rt.eccen(:)) <= opt.eccRange(2) ...
                    & double(Rt.vexpl(:)) > opt.veThresh ...
                    & isfinite(double(Rt.angle_adj(:))) & isfinite(aRep);
        if opt.requireBothVE
            keep = keep & veRep > opt.veThresh;
        end

        d = wrap180(double(Rt.angle_adj(keep)) - aRep(keep));
        allD = [allD; d]; %#ok<AGROW>
        st = circstats(d);
        rows(end+1,:) = {s, nnz(keep), st.bias, st.circSD, st.madSD, st.medAbs, ...
                         st.p90, toc(t0)}; %#ok<AGROW>
        if opt.verbose
            fprintf(['  %-14s n=%6d  bias %+6.2f  circSD %5.2f  robustSD %5.2f  ' ...
                     'median|d| %5.2f  (%4.1f s)\n'], s, nnz(keep), st.bias, ...
                     st.circSD, st.madSD, st.medAbs, toc(t0));
        end
    end

    E.perSubject = cell2table(rows, 'VariableNames', ...
        {'subject','n','biasDeg','circSDdeg','robustSDdeg','medAbsDeg','p90Deg','secs'});
    E.pooled   = circstats(allD);
    E.diffDeg  = allD;          % pooled between-solution differences, for plotting
    E.criteria = opt;

    % Per-subject sigma, then the across-subject summary. sigma for ONE solution is
    % the pairwise circSD over sqrt(2).
    sigPair = E.perSubject.circSDdeg;
    E.sigmaSingleDeg     = sigPair / sqrt(2);
    E.sigmaSingleMeanDeg = mean(E.sigmaSingleDeg);
    sr = deg2rad(E.sigmaSingleMeanDeg);
    E.lambda2 = exp(-2*sr^2);
    E.lambda4 = exp(-8*sr^2);

    % Robust variant, which is the one to quote if a few vertices flip hemifield.
    E.sigmaSingleRobustDeg = mean(E.perSubject.robustSDdeg) / sqrt(2);
    srr = deg2rad(E.sigmaSingleRobustDeg);
    E.lambda2Robust = exp(-2*srr^2);
    E.lambda4Robust = exp(-8*srr^2);

    if opt.verbose
        fprintf('\n%s\n', repmat('-',1,74));
        fprintf('Pooled over %d vertices: bias %+.2f deg, circSD %.2f, robustSD %.2f\n', ...
                numel(allD), E.pooled.bias, E.pooled.circSD, E.pooled.madSD);
        fprintf('sigma for ONE solution = circSD/sqrt(2):\n');
        fprintf('   circular  %.2f deg  -> lambda2 = %.3f, lambda4 = %.3f\n', ...
                E.sigmaSingleMeanDeg, E.lambda2, E.lambda4);
        fprintf('   robust    %.2f deg  -> lambda2 = %.3f, lambda4 = %.3f\n', ...
                E.sigmaSingleRobustDeg, E.lambda2Robust, E.lambda4Robust);
        fprintf(['\nFor comparison, reproducing the observed b1 ratio (0.395) from a\n' ...
                 'context-free world needs lambda2 = 0.395, i.e. sigma = %.0f deg.\n'], ...
                rad2deg(sqrt(-log(0.395)/2)));
    end
end

% ------------------------------------------------------------------------
function [a, ve] = read_replicate(deriv, subj, root)
% angle_adj and vexpl of the prfvista replicate, lh then rh.
%
% Prefers the LOCAL mirror written by COLLECT_PRF_REPLICATE, so repeat runs never
% touch the network. Falls back to reading the server directly if it is absent.
    local = fullfile(root, sprintf('ret_prfvista_%s.mat', subj));
    if isfile(local)
        S  = load(local, 'angle_adj', 'vexpl');
        a  = double(S.angle_adj(:));
        ve = double(S.vexpl(:));
        return
    end

    if ~isfolder(deriv)
        error('diagnose_prf_angle_error:mount', ...
              ['No local replicate at %s and %s is not mounted. Run ' ...
               '../server_extract/collect_prf_replicate.m with the volume mounted.'], ...
              local, deriv);
    end
    d = dir(fullfile(deriv, 'prfvista', subj, 'ses-*'));
    d = d([d.isdir]);
    if isempty(d)
        error('diagnose_prf_angle_error:noses', 'No ses-* under prfvista/%s.', subj);
    end
    f = fullfile(d(1).folder, d(1).name);
    a  = [grab(f,'lh.angle_adj.mgz'); grab(f,'rh.angle_adj.mgz')];
    ve = [grab(f,'lh.vexpl.mgz');     grab(f,'rh.vexpl.mgz')];
end

function v = grab(folder, name)
    m = MRIread(fullfile(folder, name));
    v = double(m.vol(:));
end

% ------------------------------------------------------------------------
function d = wrap180(d)
    d = mod(d + 180, 360) - 180;
end

% ------------------------------------------------------------------------
function st = circstats(d)
% Circular spread of a set of wrapped angular differences, in degrees.
    r  = deg2rad(d);
    Rv = mean(exp(1i*r), 'omitnan');
    st.bias   = rad2deg(angle(Rv));
    st.circSD = rad2deg(sqrt(max(-2*log(abs(Rv)), 0)));
    % Robust alternative: scaled MAD about the circular mean. Insensitive to the
    % handful of vertices whose two solutions land in opposite hemifields.
    dc        = wrap180(d - st.bias);
    st.madSD  = 1.4826 * median(abs(dc), 'omitnan');
    st.medAbs = median(abs(dc), 'omitnan');
    st.p90    = prctile(abs(dc), 90);
end

% ------------------------------------------------------------------------
function ensure_mriread()
    if exist('MRIread', 'file') == 2, return; end
    fsPath = ['/Users/jaw288/repos/Code/Projects/equivalent_input_noise_marc/' ...
              'external/freesurfer'];
    if isfolder(fsPath)
        addpath(fsPath);
    else
        error('diagnose_prf_angle_error:mriread', ...
              'MRIread not on the path and %s does not exist.', fsPath);
    end
end
