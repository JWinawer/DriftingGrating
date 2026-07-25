% Run-count and per-run R2, per session. A misaligned/mispaired run shows up as
% one run with much lower R2 than its session's others (worst/median << 1).
COL = fullfile(getenv('HOME'),'dg_collect');
subs = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
projs = {'dg','da'};
fprintf('\n%-14s %-2s | nRuns | per-run median R2                       | worst/med\n','subject','ex');
fprintf('%s\n', repmat('-',1,92));
for si=1:numel(subs)
  for pj=1:numel(projs)
    m = matfile(fullfile(COL,sprintf('glm_%s_%s.mat',subs{si},projs{pj})));
    rr = squeeze(double(m.R2run));           % vertices x nRuns
    nruns = size(rr,2);
    runmed = median(rr,1,'omitnan');
    sessmed = median(runmed);
    worst = min(runmed)/sessmed;
    fprintf('%-14s %-2s | %4d  | %-38s | %.2f\n', ...
      subs{si}, projs{pj}, nruns, num2str(runmed,'%4.1f '), worst);
  end
end
fprintf('\n(worst/med near 1.0 = all runs consistent; << 1 would flag a mispaired/bad run)\n');
