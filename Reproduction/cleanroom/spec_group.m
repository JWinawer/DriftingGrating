function G = spec_group(d, sigma, weighting)
% SPEC_GROUP  Combine per-observer values into a group estimate, one definition.
%
%   G = spec_group(d, sigma, 'equal')       % d nSubj x 1, sigma nSubj x 1
%   G = spec_group(d, sigma, 'precision')
%
% EQUAL is the ordinary across-observer t interval, and is the specification's primary
% (../SPECIFICATION.md section 6): with 8 runs, sigma-hat carries ~7 df and scatters even
% when every observer has identical true precision, so most of what you would weight by
% is estimation noise.
%
% PRECISION is the random-effects meta-analysis estimator, weights
%   w_i = 1 / (tau^2 + sigma_i^2),   tau^2 = max(0, var(y_i) - mean(sigma_i^2))
% with sigma_i MEASURED by resampling runs rather than inferred from a design that has
% no replication. tau^2 is common to every observer, so a large spread in RELIABILITY
% compresses into a small spread in WEIGHT -- see ../METHOD_DECISIONS.md section 4, which also
% explains why MATLAB's fitlme cannot express this and the closed form is used.
%
% NAN CONVENTION, shared with PRECISION_WEIGHTED_TABLE: if any observer is missing the
% whole summary is NaN, rather than an estimate over the observers that remain. An
% estimate silently computed over 7 of 8 is not comparable with one computed over 8.
%
% Returns G with .mean .lo .hi .p .n .tau .meanSigma .weightRatio .obsAgree

    d = d(:);  sigma = sigma(:);
    n = numel(d);
    G = struct('mean',NaN,'lo',NaN,'hi',NaN,'p',NaN,'n',nnz(isfinite(d)), ...
               'tau',NaN,'meanSigma',NaN,'weightRatio',NaN,'obsAgree',0);
    if G.n < n || n < 3 || any(~isfinite(sigma)), return; end

    s2   = sigma.^2;
    tau2 = max(0, var(d, 0) - mean(s2));
    tc   = tinv(0.975, n-1);

    switch lower(weighting)
        case 'equal'
            m  = mean(d);
            se = std(d) / sqrt(n);
            w  = ones(n,1) / n;
        case 'precision'
            w  = 1 ./ (tau2 + s2);  w = w / sum(w);
            m  = sum(w .* d);
            se = sqrt(sum(w.^2 .* (tau2 + s2)));
        otherwise
            error('spec_group:weighting', 'weighting must be ''equal'' or ''precision''.');
    end

    G.mean = m;  G.lo = m - tc*se;  G.hi = m + tc*se;
    G.p    = 2 * tcdf(-abs(m/se), n-1);
    G.tau  = sqrt(tau2);
    G.meanSigma   = sqrt(mean(s2));
    G.weightRatio = max(w) / min(w);
    G.obsAgree    = nnz(sign(d) == sign(m));
end
