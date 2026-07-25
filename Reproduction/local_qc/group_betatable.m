% Per-subject table: median R2 and mean beta (allsVblank = all-stationary - blank),
% for V1, V2 (manual ROIs) and whole surface, in dg and da.
COL = fullfile(getenv('HOME'),'dg_collect');
subs = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
projs = {'dg','da'};

% units sanity check on one session
S=load(fullfile(COL,'glm_sub-0037_dg.mat'),'top_contrasts'); c=S.top_contrasts;
chk = mean([c.s0Vb c.s90Vb c.s45Vb c.s135Vb],2);
fprintf('units check (sub-0037 dg): allsVblank vs mean(4 sVb): maxdiff=%.2e ; sample values %.3f %.3f %.3f\n\n',...
    max(abs(c.allsVblank-chk)), c.allsVblank(1000), c.allsVblank(5000), c.allsVblank(20000));

fprintf('%-14s | V1 R2 dg/da  V1 beta dg/da | V2 R2 dg/da  V2 beta dg/da | whole R2 dg/da  whole beta dg/da\n','subject');
fprintf('%s\n', repmat('-',1,104));
for si=1:numel(subs)
  s=subs{si};
  Rt=load(fullfile(COL,sprintf('ret_%s.mat',s)),'nLH','nRH'); nLH=double(Rt.nLH);
  L=load(fullfile(COL,sprintf('labels_%s.mat',s)));
  roi.V1=[double(L.lh_V1_REmanual); nLH+double(L.rh_V1_REmanual)];
  roi.V2=[double(L.lh_V2_REmanual); nLH+double(L.rh_V2_REmanual)];
  out=struct();
  for pj=1:2
    p=projs{pj};
    D=load(fullfile(COL,sprintf('glm_%s_%s.mat',s,p)),'R2','top_contrasts');
    r2=double(D.R2(:)); beta=double(D.top_contrasts.allsVblank(:));
    roi.whole=(1:numel(r2))';
    for rn = {'V1','V2','whole'}
      idx=roi.(rn{1}); idx=idx(idx>=1&idx<=numel(r2));
      out.(rn{1}).(p)=[median(r2(idx),'omitnan') mean(beta(idx),'omitnan')];
    end
  end
  fprintf('%-14s | %5.1f %5.1f  %6.3f %6.3f | %5.1f %5.1f  %6.3f %6.3f | %5.1f %5.1f  %6.3f %6.3f\n', s, ...
    out.V1.dg(1),out.V1.da(1),out.V1.dg(2),out.V1.da(2), ...
    out.V2.dg(1),out.V2.da(1),out.V2.dg(2),out.V2.da(2), ...
    out.whole.dg(1),out.whole.da(1),out.whole.dg(2),out.whole.da(2));
end
fprintf('\nR2 = median %%variance-explained; beta = mean of allsVblank (all-stationary minus blank).\n');
