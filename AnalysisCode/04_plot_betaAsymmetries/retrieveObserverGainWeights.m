function gainWeights = retrieveObserverGainWeights(subjects, gainSummaryFile)
% retrieveObserverGainWeights - per-observer pRF gain, for gain-weighting BOLD asymmetries
%
%   gainWeights = retrieveObserverGainWeights(subjects, gainSummaryFile)
%
% subjects        - cell array of subject IDs, in the same order as the
%                    subject dimension of medianBOLDpa/medianBOLD.
% gainSummaryFile - path to gainSummary.mat (written by dg_computeGain.m,
%                    see AnalysisCode/01_calculate_observer_gain), holding a
%                    table with columns 'subject' and 'meanGain_avg' (each
%                    observer's mean pRF gain, averaged across the
%                    prfvista_mov and prfvista protocols).
%
% Returns a 1 x numel(subjects) row vector of meanGain_avg, in the same
% subject order, so callers can divide each observer's data by their own
% gain before averaging across observers.
%
% See also DG_COMPUTEGAIN

if ~exist(gainSummaryFile, 'file')
    error('gain summary file not found: %s', gainSummaryFile);
end
S = load(gainSummaryFile, 'summary');
summaryTable = S.summary;

gainWeights = nan(1, numel(subjects));
for si = 1:numel(subjects)
    row = strcmp(summaryTable.subject, subjects{si});
    if ~any(row)
        error('subject %s not found in gain summary table %s', ...
            subjects{si}, gainSummaryFile);
    end
    gainWeights(si) = summaryTable.meanGain_avg(row);
end

end
