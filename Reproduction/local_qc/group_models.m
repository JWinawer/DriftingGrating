% Per-session GLMsingle model-parameter sanity check across all 16 sessions.
COL = fullfile(getenv('HOME'),'dg_collect');
subs = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
projs = {'dg','da'};
fprintf('\n%-14s %-2s | pcnum  noise%%  medFRAC  %%shrunk(<.5)  modalHRF  HRF#  xvalPeak\n','subject','ex');
fprintf('%s\n', repmat('-',1,82));
PC=[]; MF=[]; SH=[]; lab={};
for si=1:numel(subs)
  for pj=1:numel(projs)
    m = matfile(fullfile(COL,sprintf('glm_%s_%s.mat',subs{si},projs{pj})));
    pcnum = double(m.pcnum);
    np = double(m.noisepool); np=np(:);
    fv = double(m.FRACvalue); fv=fv(:);
    hi = double(m.HRFindex); hi=hi(:);
    xv = double(m.xvaltrend); xv=xv(:);
    medfrac = median(fv,'omitnan');
    pShrunk = 100*mean(fv<0.5);
    modalHRF = mode(hi(~isnan(hi)));
    hrfN = numel(unique(hi(~isnan(hi))));
    fprintf('%-14s %-2s | %4d  %6.0f  %7.3f  %10.1f  %7d  %4d  %8.2f\n', ...
      subs{si}, projs{pj}, pcnum, 100*mean(np), medfrac, pShrunk, modalHRF, hrfN, max(xv));
    PC(end+1)=pcnum; MF(end+1)=medfrac; SH(end+1)=pShrunk; lab{end+1}=sprintf('%s/%s',subs{si},projs{pj});
  end
end
z=@(x)(x-mean(x))./std(x);
fprintf('\npcnum : min %d max %d median %d\n', min(PC),max(PC),median(PC));
fprintf('medFRAC range [%.3f %.3f]\n', min(MF), max(MF));
fprintf('\noutliers (|z|>2 on medFRAC or %%shrunk):\n');
zf=z(MF); zs=z(SH); any_out=false;
for i=1:numel(lab)
  if abs(zf(i))>2 || abs(zs(i))>2
    fprintf('  %-18s medFRAC=%.3f (z=%+.1f)  %%shrunk=%.1f (z=%+.1f)\n', lab{i}, MF(i), zf(i), SH(i), zs(i)); any_out=true;
  end
end
if ~any_out, fprintf('  none.\n'); end
