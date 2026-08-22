function dg_inspectPrfResults(varargin)
% dg_inspectPrfResults - STEP 1: report what is inside the prfvista results
%
%   dg_inspectPrfResults()
%   dg_inspectPrfResults('root', '/path/to/derivatives/prfvista_mov')
%
% Run this FIRST, at the server, before computing anything. It only reads
% variable headers plus a couple of small fields, so it is fast and does not
% load the big arrays.
%
% It answers the one open question: do the prfvista results files carry the
% fitted betas and the stimulus? Those are what the gain calculation needs
% and they are NOT in the ret_*.mat summary maps that were copied to NY.
%
% Send the printed output back (a few KB). Easiest:
%   diary ~/dg_inspect.txt; dg_inspectPrfResults(); diary off
%
% See also DG_COMPUTEGAIN, RMMODELGAIN

p = inputParser;
p.addParameter('root', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
    'data_bids/derivatives/prfvista_mov'], @ischar);
p.parse(varargin{:});
root = p.Results.root;

fprintf('root: %s\n', root);
if ~exist(root, 'dir')
    fprintf('*** root does not exist. Pass the right path with ''root''.\n');
    return
end

subs = dir(fullfile(root, 'sub-*'));
subs = subs([subs.isdir]);
fprintf('found %d subject folders\n\n', numel(subs));

for ii = 1:numel(subs)
    sub  = subs(ii).name;
    hits = dir(fullfile(root, sub, 'ses-*', '*results.mat'));
    fprintf('===== %s : %d results file(s)\n', sub, numel(hits));

    for jj = 1:numel(hits)
        f = fullfile(hits(jj).folder, hits(jj).name);
        fprintf('  %s  (%.0f MB)\n', hits(jj).name, hits(jj).bytes/1e6);

        try
            w = whos('-file', f);
            for k = 1:numel(w)
                fprintf('     %-14s %-8s [%s]  %.1f MB\n', w(k).name, ...
                    w(k).class, num2str(w(k).size), w(k).bytes/1e6);
            end
        catch ME
            fprintf('     [could not read variable headers: %s]\n', ME.message);
        end

        % Look one level in, at whichever top-level variable holds the model
        try
            V = load(f);
            top = fieldnames(V);
            for k = 1:numel(top)
                describeNode(V.(top{k}), top{k}, 5);
            end
            clear V
        catch ME
            fprintf('     [could not load: %s]\n', ME.message);
        end
        fprintf('\n');
    end
end

fprintf(['\nWHAT MATTERS:\n' ...
    '  model.beta                       <- the amplitude, gain comes from this\n' ...
    '  params.analysis.allstimimages    <- stimulus, already HRF-convolved\n' ...
    '  params.analysis.X / .Y           <- sampling grid\n' ...
    'If all of those are present, dg_computeGain will work.\n']);

return


% ------------------------------------------------------------------------
function describeNode(node, name, indent)
pad = repmat(' ', 1, indent);
if iscell(node) && ~isempty(node)
    fprintf('%s%s{1}:\n', pad, name);
    describeNode(node{1}, [name '{1}'], indent + 2);
    return
end
if ~isstruct(node), return; end

f = fieldnames(node);
fprintf('%s%s fields: %s\n', pad, name, strjoin(f', ', '));

% Report the fields the gain calculation needs
wanted = {'beta','x0','y0','sigma','rss','rawrss','exponent','description','desc'};
for k = 1:numel(f)
    if ismember(f{k}, wanted)
        v = node.(f{k});
        if ischar(v)
            fprintf('%s  .%-12s ''%s''\n', pad, f{k}, v);
        elseif isnumeric(v)
            fprintf('%s  .%-12s %s [%s]\n', pad, f{k}, class(v), num2str(size(v)));
        elseif isstruct(v)
            fprintf('%s  .%-12s struct, fields: %s\n', pad, f{k}, ...
                strjoin(fieldnames(v)', ', '));
        end
    end
end

% Descend into params.analysis
if isfield(node, 'analysis') && isstruct(node.analysis)
    a = node.analysis;
    fprintf('%s  .analysis fields: %s\n', pad, strjoin(fieldnames(a)', ', '));
    for nm = {'allstimimages','allstimimages_unconvolved','X','Y','Hrf', ...
              'HrfMaxResponse','sampleRate','scan_number'}
        if isfield(a, nm{1})
            v = a.(nm{1});
            fprintf('%s    .%-26s %s [%s]\n', pad, nm{1}, class(v), ...
                num2str(size(v)));
        else
            fprintf('%s    .%-26s MISSING\n', pad, nm{1});
        end
    end
end

% Recurse into model/params one level
for nm = {'model','params'}
    if isfield(node, nm{1})
        describeNode(node.(nm{1}), [name '.' nm{1}], indent + 2);
    end
end
return
