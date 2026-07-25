% meanvol (EPI intensity -> dropout) and surface-alignment consistency.
COL = fullfile(getenv('HOME'),'dg_collect');
subs = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
projs = {'dg','da'};
fprintf('\n%-14s | surfLen  R2len(dg/da)  align | meanvolV1 dg/da | meanvolWhole dg/da\n','subject');
fprintf('%s\n', repmat('-',1,90));
for si=1:numel(subs)
  s = subs{si};
  Rt = load(fullfile(COL,sprintf('ret_%s.mat',s)),'nLH','nRH');
  nLH=double(Rt.nLH); nRH=double(Rt.nRH); surf=nLH+nRH;
  L = load(fullfile(COL,sprintf('labels_%s.mat',s)));
  v1=[]; if isfield(L,'lh_V1_REmanual'), v1=[double(L.lh_V1_REmanual); nLH+double(L.rh_V1_REmanual)]; end
  vals=struct(); r2len=[0 0];
  for pj=1:2
    m = matfile(fullfile(COL,sprintf('glm_%s_%s.mat',s,projs{pj})));
    mv = double(m.meanvol); mv=mv(:); r2len(pj)=numel(mv);
    vv = v1(v1>=1 & v1<=numel(mv));
    vals.(projs{pj}) = [median(mv(vv),'omitnan') median(mv,'omitnan')];
  end
  align = (r2len(1)==surf) && (r2len(2)==surf);
  fprintf('%-14s | %6d  %6d/%-6d %5s | %6.0f / %-6.0f | %6.0f / %-6.0f\n', ...
    s, surf, r2len(1), r2len(2), string(align), vals.dg(1), vals.da(1), vals.dg(2), vals.da(2));
end
fprintf('\n(align=true means GLM R2 length matches the retinotopy surface for both exps;\n meanvol is arbitrary EPI units -- watch for a session far below its peers = dropout)\n');
