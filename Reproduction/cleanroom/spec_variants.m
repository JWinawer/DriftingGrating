function V = spec_variants(tag)
% SPEC_VARIANTS  The three analysis variants that get drawn, and what each isolates.
%
%   V = spec_variants()          % all three, in reporting order
%   V = spec_variants('roipw')   % just one
%
% Route and weighting are ORTHOGONAL knobs, so the three variants are arranged to
% change one at a time. Comparing them pairwise is the point:
%
%   spec   vs  roi     isolates the ROUTE      (continuous thetaV vs eight wedges)
%   roi    vs  roipw   isolates the WEIGHTING  (equal vs inverse-variance)
%
% Changing both at once would confound them, which is why there is no
% harmonic-plus-precision variant here: it would leave no pair differing in one thing.
%
%   spec    harmonic model, continuous thetaV, equal weighting     <- PRIMARY
%           the settled specification, ../EXTRASTRIATE.md section 1
%   roi     eight polar-angle wedges, equal weighting
%           the manuscript route. Identical to spec on complete data -- they are one
%           estimator re-parameterised -- and diverging in proportion to empty cells,
%           because it averages over whichever wedges an observer happens to have.
%   roipw   eight polar-angle wedges, precision weighting
%           w_i = 1/(tau^2 + sigma_i^2) with sigma_i measured by resampling runs.
%
% Fields: .tag .route .weighting .label .short

    all = { ...
      'spec',  'harmonic', 'equal', ...
        'harmonic model, continuous \theta_V, equal weighting', 'continuous \theta_V, equal wt'; ...
      'roi',   'roi',      'equal', ...
        'ROI route, eight polar-angle wedges, equal weighting',  'ROI wedges, equal wt'; ...
      'roipw', 'roi',      'precision', ...
        'ROI route, eight polar-angle wedges, precision weighting', 'ROI wedges, precision wt'};

    V = struct('tag', all(:,1), 'route', all(:,2), 'weighting', all(:,3), ...
               'label', all(:,4), 'short', all(:,5));
    if nargin > 0 && ~isempty(tag)
        k = find(strcmpi({V.tag}, tag), 1);
        assert(~isempty(k), 'spec_variants:tag', ...
               'unknown variant ''%s''; expected one of %s', tag, strjoin({V.tag}, ', '));
        V = V(k);
    end
end
