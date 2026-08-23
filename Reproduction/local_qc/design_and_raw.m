%% Part 1: example design matrix (dg, one run) + Part 2: no-blank-subtracted CSV
COL  = fullfile(getenv('HOME'),'dg_collect');
CSVp = '/Users/jaw288/repos/Code/Projects/DriftingGrating/Support/allsubjectsTable.csv';
QC   = '/Users/jaw288/repos/Code/Projects/DriftingGrating/Reproduction/local_qc';

%% -------- design matrix from a real dg design file --------
d = dir(fullfile(COL,'design','dg','_dg_01_Run1_S01_design_Run1.mat'));
e = getfield(load(fullfile(d(1).folder,d(1).name),'expDes'),'expDes');
tm = e.trialMat; nt = size(tm,1);
stimDur=e.stimDur_s; iti=e.itiDur_s; pad=e.runPadding_s; slot=stimDur+iti; total=e.total_s;
% map each trial to one of 13 condition rows: 1-8 motion, 9-12 stationary, 13 blank
oris = sort(unique(tm(tm(:,2)==1,3)));   dirs = sort(unique(tm(tm(:,2)==2,4)));
cond = zeros(nt,1);
for k=1:nt
  if tm(k,2)==0, cond(k)=13;
  elseif tm(k,2)==1, cond(k)=8+find(oris==tm(k,3),1);
  else, cond(k)=find(dirs==tm(k,4),1);
  end
end
dt=0.5; tvec=0:dt:total; DM=zeros(numel(tvec),13);
for k=1:nt
  on = pad + (k-1)*slot;                 % onset (s); grating shown for first stimDur
  dur = stimDur; if cond(k)==13, dur=iti* (slot/iti); end   % blank trial spans the slot
  m = tvec>=on & tvec<on+dur; DM(m,cond(k))=1;
end
f=figure('Visible','off','Position',[100 100 1100 500]);
imagesc(tvec, 1:13, DM'); colormap(flipud(gray)); hold on
yline(12.5,'r-','LineWidth',1.5);
set(gca,'YTick',1:13,'YTickLabel',[compose('mot%d',1:8), compose('stat%d',1:4), {'BLANK'}]);
xlabel('time in run (s)'); ylabel('GLM condition'); title('Example design matrix — dg, one run (52 trials, 8 runs/session)');
text(total/2, 15.2, sprintf(['pink noise is on screen the ENTIRE run: surround during gratings, ' ...
   'full-field during every ITI (%gs), every BLANK trial, and the %gs padding at run start/end. ' ...
   'No mean-luminance baseline exists.'], iti, pad), 'HorizontalAlignment','center','FontSize',9,'Color',[.5 0 0]);
ylim([0.5 15.8]);
% shade padding
patch([0 pad pad 0],[0.5 0.5 13.5 13.5],[1 .8 .8],'FaceAlpha',.4,'EdgeColor','none');
patch([total-pad total total total-pad],[0.5 0.5 13.5 13.5],[1 .8 .8],'FaceAlpha',.4,'EdgeColor','none');
pngp = fullfile(QC,'example_design_matrix.png');
print(f,'-dpng','-r150',pngp);
fprintf('nBlankTrials=%d of %d (%.0f%%); slot=%gs; run=%gs; wrote %s\n', nnz(cond==13), nt, 100*nnz(cond==13)/nt, slot, total, pngp);

%% -------- CSV incl. RAW (no blank subtracted) --------
subs={'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
areas={'V1','V2','V3'}; ns=numel(subs);
R2=nan(ns,4,2);
for si=1:ns
  s=subs{si}; Rt=load(fullfile(COL,sprintf('ret_%s.mat',s)),'nLH','eccen'); nLH=double(Rt.nLH); ecc=double(Rt.eccen(:));
  L=load(fullfile(COL,sprintf('labels_%s.mat',s)));
  for pj=1:2, p={'dg','da'}; p=p{pj};
    r2=double(getfield(load(fullfile(COL,sprintf('glm_%s_%s.mat',s,p)),'R2'),'R2')); r2=r2(:);
    for a=1:3, lab=[double(L.(['lh_' areas{a} '_REmanual'])); nLH+double(L.(['rh_' areas{a} '_REmanual']))];
      lab=lab(lab>=1&lab<=numel(r2)); idx=lab(ecc(lab)>=4&ecc(lab)<=8); R2(si,a,pj)=median(r2(idx),'omitnan'); end
    R2(si,4,pj)=median(r2,'omitnan');
  end
end
cM={'cartexp_vertical_grating_rightwards_motion','cartexp_horizontal_grating_upwards_motion','cartexp_vertical_grating_leftwards_motion','cartexp_horizontal_grating_downwards_motion','cartexp_leftleaning_grating_upperrightwards_motion','cartexp_rightleaning_grating_upperleftwards_motion','cartexp_leftleaning_grating_lowerleftwards_motion','cartexp_rightleaning_grating_lowerrightwards_motion'};
cS={'cartexp_horizontal_stationary','cartexp_vertical_stationary','cartexp_rightleaning_grating_stationary','cartexp_leftleaning_grating_stationary'};
pM={'polexp_pinwheel_grating_clockwise_motion','polexp_annulus_grating_outwards_motion','polexp_pinwheel_grating_cclockwise_motion','polexp_annulus_grating_inwards_motion','polexp_ccspiral_grating_clockoutwards_motion','polexp_cspiral_grating_cclockoutwards_motion','polexp_ccspiral_grating_cclockinwards_motion','polexp_cspiral_grating_clockinwards_motion'};
pS={'polexp_annulus_grating_stationary','polexp_pinwheel_grating_stationary','polexp_cspiral_grating_stationary','polexp_ccspiral_grating_stationary'};
o=detectImportOptions(CSVp); o.VariableNamingRule='preserve';
o.SelectedVariableNames=[{'subject','visual_area','pRF_ecc'},cM,cS,{'cartexp_blank'},pM,pS,{'polexp_blank'}];
T=readtable(CSVp,o); gm=@(c) mean(T{:,c},2);
% raw (no blank sub) and blank, per exp: [dg da]
Stat=[gm(cS) gm(pS)]; Mot=[gm(cM) gm(pM)]; Blank=[T.cartexp_blank T.polexp_blank];
va=string(T.visual_area); ec=T.pRF_ecc; sj=string(T.subject);
S_=nan(ns,4,2); M_=nan(ns,4,2); B_=nan(ns,4,2); BS=nan(ns,4,2); BM=nan(ns,4,2);
for si=1:ns, sm=sj==subs{si};
  for a=1:3, m=sm & va==areas{a} & ec>=4 & ec<=8;
    for pj=1:2, S_(si,a,pj)=mean(Stat(m,pj),'omitnan'); M_(si,a,pj)=mean(Mot(m,pj),'omitnan'); B_(si,a,pj)=mean(Blank(m,pj),'omitnan');
      BS(si,a,pj)=mean(Stat(m,pj)-Blank(m,pj),'omitnan'); BM(si,a,pj)=mean(Mot(m,pj)-Blank(m,pj),'omitnan'); end
  end
  for pj=1:2, S_(si,4,pj)=mean(Stat(sm,pj),'omitnan'); M_(si,4,pj)=mean(Mot(sm,pj),'omitnan'); B_(si,4,pj)=mean(Blank(sm,pj),'omitnan');
    BS(si,4,pj)=mean(Stat(sm,pj)-Blank(sm,pj),'omitnan'); BM(si,4,pj)=mean(Mot(sm,pj)-Blank(sm,pj),'omitnan'); end
end
Out=table(); Out.subject=string(subs(:)); Out.subjectNum=(1:ns)';
regs={'V1','V2','V3','wh'}; exps={'dg','da'};
for r=1:4, for pj=1:2
  Out.(sprintf('%s_R_%s',regs{r},exps{pj}))=round(R2(:,r,pj),3);
  Out.(sprintf('%s_Stat_%s',regs{r},exps{pj}))=round(S_(:,r,pj),4);      % raw stationary
  Out.(sprintf('%s_Mot_%s',regs{r},exps{pj}))=round(M_(:,r,pj),4);       % raw motion
  Out.(sprintf('%s_Blank_%s',regs{r},exps{pj}))=round(B_(:,r,pj),4);     % blank beta
  Out.(sprintf('%s_Beta_%s',regs{r},exps{pj}))=round(BS(:,r,pj),4);      % stationary - blank
  Out.(sprintf('%s_BetaMot_%s',regs{r},exps{pj}))=round(BM(:,r,pj),4);   % motion - blank
end, end
outCsv=fullfile(getenv('HOME'),'Downloads','glm_summary.csv');
writetable(Out,outCsv); copyfile(outCsv,fullfile(QC,'glm_summary.csv'));
fprintf('wrote %s (%d x %d)\n', outCsv, height(Out), width(Out));
disp(Out(:,{'subject','V1_Stat_dg','V1_Mot_dg','V1_Blank_dg','V1_Stat_da','V1_Mot_da','V1_Blank_da'}));
