% Group GLM QC across all 8 subjects x {dg,da}, from local ~/dg_collect files.
COL = fullfile(getenv('HOME'),'dg_collect');
subs = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
projs = {'dg','da'};

fprintf('\n%-14s %-2s | wholeR2  V1all  V1patch  nPatch  noise%%  runR2[min-max]  V1def\n','subject','ex');
fprintf('%s\n', repmat('-',1,86));
for si=1:numel(subs)
  s = subs{si};
  Rt = load(fullfile(COL,sprintf('ret_%s.mat',s)),'nLH','nRH','eccen','vexpl');
  nLH = double(Rt.nLH); nRH = double(Rt.nRH);
  ecc = double(Rt.eccen(:)); vex = double(Rt.vexpl(:));
  L = load(fullfile(COL,sprintf('labels_%s.mat',s)));
  v1 = []; v1tag = 'NONE';
  if isfield(L,'lh_V1_REmanual') && isfield(L,'rh_V1_REmanual')
    v1 = [double(L.lh_V1_REmanual); nLH + double(L.rh_V1_REmanual)]; v1tag='REmanual';
  elseif isfield(L,'lh_V1') && isfield(L,'rh_V1')
    v1 = [double(L.lh_V1); nLH + double(L.rh_V1)]; v1tag='atlas';
  end
  for pj=1:numel(projs)
    p = projs{pj};
    m = matfile(fullfile(COL,sprintf('glm_%s_%s.mat',s,p)));
    r2 = double(m.R2); r2=r2(:);
    np = double(m.noisepool); np=np(:);
    rr = squeeze(double(m.R2run)); runmed = median(rr,1,'omitnan');
    wholeMed = median(r2,'omitnan');
    v1all=NaN; v1patch=NaN; npatch=0; lenwarn='';
    if numel(r2) ~= nLH+nRH, lenwarn = sprintf(' [!len %d/%d]',numel(r2),nLH+nRH); end
    if ~isempty(v1)
      v1 = v1(v1>=1 & v1<=numel(r2));
      v1all = median(r2(v1),'omitnan');
      inp = v1(ecc(v1)>=4 & ecc(v1)<=8 & vex(v1)>0.1);
      npatch = numel(inp); v1patch = median(r2(inp),'omitnan');
    end
    fprintf('%-14s %-2s | %6.2f  %5.2f  %6.2f  %5d  %5.0f  [%4.1f-%4.1f]   %s%s\n', ...
      s, p, wholeMed, v1all, v1patch, npatch, 100*mean(np), min(runmed), max(runmed), v1tag, lenwarn);
  end
end
fprintf('\n(V1patch = V1 label AND 4-8deg AND pRF-vexpl>0.1 -- the analysis patch)\n');
