% Build glm_summary.csv: per subject, median GLM R2 + mean raw % beta for
% stationary-minus-blank (Beta) and motion-minus-blank (BetaMot), for V1/V2/V3
% (4-8 deg, no variance threshold) and whole surface. Columns named to match the
% user's plotting code: V<n>_R_<exp>, V<n>_Beta_<exp>, V<n>_BetaMot_<exp>, wh_*.
COL  = fullfile(getenv('HOME'),'dg_collect');
CSVp = '/Users/jaw288/repos/Code/Projects/DriftingGrating/Support/allsubjectsTable.csv';
subs  = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
areas = {'V1','V2','V3'};
ns=numel(subs);

%% R2 (median GLM) : REmanual label AND pRF-ecc in [4,8]; whole = all vertices
R2 = nan(ns,4,2);   % (subject, [V1 V2 V3 whole], [dg da])
for si=1:ns
  s=subs{si};
  Rt=load(fullfile(COL,sprintf('ret_%s.mat',s)),'nLH','eccen'); nLH=double(Rt.nLH); ecc=double(Rt.eccen(:));
  L=load(fullfile(COL,sprintf('labels_%s.mat',s)));
  for pj=1:2
    p={'dg','da'}; p=p{pj};
    r2=double(getfield(load(fullfile(COL,sprintf('glm_%s_%s.mat',s,p)),'R2'),'R2')); r2=r2(:);
    for a=1:3
      lab=[double(L.(['lh_' areas{a} '_REmanual'])); nLH+double(L.(['rh_' areas{a} '_REmanual']))];
      lab=lab(lab>=1&lab<=numel(r2)); idx=lab(ecc(lab)>=4 & ecc(lab)<=8);
      R2(si,a,pj)=median(r2(idx),'omitnan');
    end
    R2(si,4,pj)=median(r2,'omitnan');
  end
end

%% raw % betas from CSV
cartMot={'cartexp_vertical_grating_rightwards_motion','cartexp_horizontal_grating_upwards_motion','cartexp_vertical_grating_leftwards_motion','cartexp_horizontal_grating_downwards_motion','cartexp_leftleaning_grating_upperrightwards_motion','cartexp_rightleaning_grating_upperleftwards_motion','cartexp_leftleaning_grating_lowerleftwards_motion','cartexp_rightleaning_grating_lowerrightwards_motion'};
cartStat={'cartexp_horizontal_stationary','cartexp_vertical_stationary','cartexp_rightleaning_grating_stationary','cartexp_leftleaning_grating_stationary'};
polMot={'polexp_pinwheel_grating_clockwise_motion','polexp_annulus_grating_outwards_motion','polexp_pinwheel_grating_cclockwise_motion','polexp_annulus_grating_inwards_motion','polexp_ccspiral_grating_clockoutwards_motion','polexp_cspiral_grating_cclockoutwards_motion','polexp_ccspiral_grating_cclockinwards_motion','polexp_cspiral_grating_clockinwards_motion'};
polStat={'polexp_annulus_grating_stationary','polexp_pinwheel_grating_stationary','polexp_cspiral_grating_stationary','polexp_ccspiral_grating_stationary'};
opts=detectImportOptions(CSVp); opts.VariableNamingRule='preserve';
opts.SelectedVariableNames=[{'subject','visual_area','pRF_ecc'}, cartMot, cartStat, {'cartexp_blank'}, polMot, polStat, {'polexp_blank'}];
T=readtable(CSVp,opts);
getm=@(cols) mean(T{:,cols},2);
bStat=nan(height(T),2); bMot=nan(height(T),2);
bStat(:,1)=getm(cartStat)-T.cartexp_blank;  bMot(:,1)=getm(cartMot)-T.cartexp_blank;   % dg
bStat(:,2)=getm(polStat) -T.polexp_blank;   bMot(:,2)=getm(polMot) -T.polexp_blank;    % da
va=string(T.visual_area); ec=T.pRF_ecc; sj=string(T.subject);
BS=nan(ns,4,2); BM=nan(ns,4,2);
for si=1:ns
  sm=sj==subs{si};
  for a=1:3
    m=sm & va==areas{a} & ec>=4 & ec<=8;
    for pj=1:2, BS(si,a,pj)=mean(bStat(m,pj),'omitnan'); BM(si,a,pj)=mean(bMot(m,pj),'omitnan'); end
  end
  for pj=1:2, BS(si,4,pj)=mean(bStat(sm,pj),'omitnan'); BM(si,4,pj)=mean(bMot(sm,pj),'omitnan'); end
end

%% assemble output table with the requested column names
Out=table(); Out.subject=string(subs(:)); Out.subjectNum=(1:ns)';
regs={'V1','V2','V3','wh'}; exps={'dg','da'};
for r=1:4
  for pj=1:2
    Out.(sprintf('%s_R_%s',regs{r},exps{pj}))       = round(R2(:,r,pj),3);
    Out.(sprintf('%s_Beta_%s',regs{r},exps{pj}))    = round(BS(:,r,pj),4);
    Out.(sprintf('%s_BetaMot_%s',regs{r},exps{pj})) = round(BM(:,r,pj),4);
  end
end
outCsv=fullfile(getenv('HOME'),'Downloads','glm_summary.csv');
writetable(Out,outCsv);
copyfile(outCsv,'/Users/jaw288/repos/Code/Projects/DriftingGrating/Reproduction/local_qc/glm_summary.csv');
fprintf('wrote %s (%d rows x %d cols)\n', outCsv, height(Out), width(Out));
disp(Out(:,{'subject','V1_R_dg','V1_R_da','V1_Beta_dg','V1_Beta_da','V1_BetaMot_dg','V1_BetaMot_da'}));
