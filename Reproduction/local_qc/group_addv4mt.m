% Add V4 (hV4) and MT (pMT) to the summary, WHOLE ROI (4-8 deg too sparse for these).
% V1/V2/V3 stay at 4-8 deg. Writes glm_summary.csv and prints a motion-focused V4/MT table.
COL=fullfile(getenv('HOME'),'dg_collect');
CSVp='/Users/jaw288/repos/Code/Projects/DriftingGrating/Support/allsubjectsTable.csv';
QC='/Users/jaw288/repos/Code/Projects/DriftingGrating/Reproduction/local_qc';
subs={'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
ns=numel(subs);
% per-area config: {outName, REmanual base, visual_area string, eccLo, eccHi}
cfg={ {'V1','V1','V1',4,8}, {'V2','V2','V2',4,8}, {'V3','V3','V3',4,8}, ...
      {'V4','hV4','hV4',0,inf}, {'MT','pMT','MT',0,inf}, {'wh','','',0,inf} };
na=numel(cfg);

%% R2 (median GLM) from glm files, per area
R2=nan(ns,na,2);
for si=1:ns
  s=subs{si}; Rt=load(fullfile(COL,sprintf('ret_%s.mat',s)),'nLH','eccen'); nLH=double(Rt.nLH); ecc=double(Rt.eccen(:));
  L=load(fullfile(COL,sprintf('labels_%s.mat',s)));
  for pj=1:2, p={'dg','da'}; p=p{pj};
    r2=double(getfield(load(fullfile(COL,sprintf('glm_%s_%s.mat',s,p)),'R2'),'R2')); r2=r2(:);
    for a=1:na
      c=cfg{a};
      if isempty(c{2}), R2(si,a,pj)=median(r2,'omitnan'); continue; end
      lf=['lh_' c{2} '_REmanual']; rf=['rh_' c{2} '_REmanual'];
      lab=[double(L.(lf)); nLH+double(L.(rf))]; lab=lab(lab>=1&lab<=numel(r2));
      idx=lab(ecc(lab)>=c{4}&ecc(lab)<=c{5}); R2(si,a,pj)=median(r2(idx),'omitnan');
    end
  end
end

%% raw % betas from CSV
cM={'cartexp_vertical_grating_rightwards_motion','cartexp_horizontal_grating_upwards_motion','cartexp_vertical_grating_leftwards_motion','cartexp_horizontal_grating_downwards_motion','cartexp_leftleaning_grating_upperrightwards_motion','cartexp_rightleaning_grating_upperleftwards_motion','cartexp_leftleaning_grating_lowerleftwards_motion','cartexp_rightleaning_grating_lowerrightwards_motion'};
cS={'cartexp_horizontal_stationary','cartexp_vertical_stationary','cartexp_rightleaning_grating_stationary','cartexp_leftleaning_grating_stationary'};
pM={'polexp_pinwheel_grating_clockwise_motion','polexp_annulus_grating_outwards_motion','polexp_pinwheel_grating_cclockwise_motion','polexp_annulus_grating_inwards_motion','polexp_ccspiral_grating_clockoutwards_motion','polexp_cspiral_grating_cclockoutwards_motion','polexp_ccspiral_grating_cclockinwards_motion','polexp_cspiral_grating_clockinwards_motion'};
pS={'polexp_annulus_grating_stationary','polexp_pinwheel_grating_stationary','polexp_cspiral_grating_stationary','polexp_ccspiral_grating_stationary'};
o=detectImportOptions(CSVp); o.VariableNamingRule='preserve';
o.SelectedVariableNames=[{'subject','visual_area','pRF_ecc'},cM,cS,{'cartexp_blank'},pM,pS,{'polexp_blank'}];
T=readtable(CSVp,o); gm=@(c) mean(T{:,c},2);
Stat=[gm(cS) gm(pS)]; Mot=[gm(cM) gm(pM)]; Blank=[T.cartexp_blank T.polexp_blank];
va=string(T.visual_area); ec=T.pRF_ecc; sj=string(T.subject);
S_=nan(ns,na,2); M_=nan(ns,na,2); B_=nan(ns,na,2); NV=nan(ns,na);
for si=1:ns, sm=sj==subs{si};
  for a=1:na
    c=cfg{a};
    if isempty(c{3}), m=sm; else, m=sm & va==c{3} & ec>=c{4} & ec<=c{5}; end
    NV(si,a)=nnz(m);
    for pj=1:2, S_(si,a,pj)=mean(Stat(m,pj),'omitnan'); M_(si,a,pj)=mean(Mot(m,pj),'omitnan'); B_(si,a,pj)=mean(Blank(m,pj),'omitnan'); end
  end
end

%% write full CSV
Out=table(); Out.subject=string(subs(:)); Out.subjectNum=(1:ns)';
exps={'dg','da'};
for a=1:na, nm=cfg{a}{1};
  for pj=1:2
    Out.(sprintf('%s_R_%s',nm,exps{pj}))=round(R2(:,a,pj),3);
    Out.(sprintf('%s_Stat_%s',nm,exps{pj}))=round(S_(:,a,pj),4);
    Out.(sprintf('%s_Mot_%s',nm,exps{pj}))=round(M_(:,a,pj),4);
    Out.(sprintf('%s_Blank_%s',nm,exps{pj}))=round(B_(:,a,pj),4);
    Out.(sprintf('%s_Beta_%s',nm,exps{pj}))=round(S_(:,a,pj)-B_(:,a,pj),4);
    Out.(sprintf('%s_BetaMot_%s',nm,exps{pj}))=round(M_(:,a,pj)-B_(:,a,pj),4);
  end
end
outCsv=fullfile(getenv('HOME'),'Downloads','glm_summary.csv');
writetable(Out,outCsv); copyfile(outCsv,fullfile(QC,'glm_summary.csv'));
fprintf('wrote %s (%d x %d)\n\n', outCsv, height(Out), width(Out));

%% motion-focused V4/MT table (whole ROI). Mot-Stat = motion selectivity (baseline-free)
for a=[5 4]  % MT then V4
  nm=cfg{a}{1};
  fprintf('=== %s (whole ROI) : Stat / Mot / Blank raw, Mot-Stat, R2 ===\n', nm);
  fprintf('%-14s | dg: Stat  Mot   Blk  |Mot-Stat| R2  | da: Stat  Mot   Blk  |Mot-Stat| R2  | nverts\n','subject');
  for si=1:ns
    fprintf('%-14s | %5.2f %5.2f %5.2f |  %5.2f | %4.1f | %5.2f %5.2f %5.2f |  %5.2f | %4.1f | %d\n', ...
      subs{si}, S_(si,a,1),M_(si,a,1),B_(si,a,1), M_(si,a,1)-S_(si,a,1), R2(si,a,1), ...
      S_(si,a,2),M_(si,a,2),B_(si,a,2), M_(si,a,2)-S_(si,a,2), R2(si,a,2), NV(si,a));
  end
  fprintf('\n');
end
