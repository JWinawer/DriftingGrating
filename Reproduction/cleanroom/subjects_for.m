function subs = subjects_for(cfg, which)
% SUBJECTS_FOR  The observer list for one experiment. Use this, not cfg.subjects.
%
%   subs = subjects_for(cfg, 'dg')       % dg, honouring cfg.dgSubjectMode
%   subs = subjects_for(cfg, 'da')       % the 7 valid da observers
%   subs = subjects_for(cfg, 'matched')  % the only set valid for a dg-vs-da contrast
%   subs = subjects_for(cfg, expCfg)     % same, keyed off expCfg.name
%
% WHY THIS EXISTS. There used to be a single cfg.subjects shared by both
% experiments, so every array's subject dimension lined up by construction and
% cross-experiment code could pair observers by POSITION. That is no longer true:
% dg has 13 observers and da has 7 (sub-0395's da session used a pilot stimulus).
% Positional pairing across experiments now silently pairs the wrong people, which
% is why ASSERT_SAME_OBSERVERS exists and why nothing should hand-index the
% cfg.subjects_* fields.
%
% See CONFIG_REPRO and ../../AGENTS.md standing fact 7.

    if nargin < 2 || isempty(which), which = 'matched'; end
    if isstruct(which)
        if ~isfield(which, 'name')
            error('subjects_for:expCfg', 'expCfg struct has no .name field.');
        end
        which = which.name;
    end
    if ~(ischar(which) || isstring(which))
        error('subjects_for:which', 'WHICH must be ''dg'', ''da'', ''matched'', or an expCfg struct.');
    end

    switch lower(char(which))
        case 'dg'
            switch lower(cfg.dgSubjectMode)
                case 'matched', subs = cfg.subjects_matched;
                case 'all',     subs = cfg.subjects_dg;
                otherwise
                    error('subjects_for:mode', ...
                          'cfg.dgSubjectMode must be ''matched'' or ''all'', not ''%s''.', ...
                          cfg.dgSubjectMode);
            end
        case 'da'
            subs = cfg.subjects_da;
        case 'matched'
            subs = cfg.subjects_matched;
        otherwise
            error('subjects_for:which', ...
                  'Unknown experiment ''%s'' -- expected ''dg'', ''da'' or ''matched''.', char(which));
    end
end
