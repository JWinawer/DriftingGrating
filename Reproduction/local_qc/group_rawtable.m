% Per-subject table: median GLM R2 (from extracted glm files) + mean RAW % beta
% (from allsubjectsTable.csv), for V1/V2/V3 restricted to 4-8 deg eccentricity
% (NO variance-explained threshold), and whole surface (no restriction).
COL  = fullfile(getenv('HOME'),'dg_collect');
CSVp = '/Users/jaw288/repos/Code/Projects/DriftingGrating/Support/allsubjectsTable.csv';
subs  = {'sub-0037','sub-0201','sub-0255','sub-wlsubj123','sub-wlsubj124','sub-0395','sub-0426','sub-0250'};
areas = {'V1','V2','V3'};
vk = @(s) matlab.lang.makeValidName(s);

%% R2 (GLM median) : REmanual label AND pRF-ecc in [4,8], no threshold
R2 = struct();
for si=1:numel(subs)
  s=subs{si};
  Rt=load(fullfile(COL,sprintf('ret_%s.mat',s)),'nLH','eccen','vexpl');
  nLH=double(Rt.nLH); ecc=double(Rt.eccen(:)); vex=double(Rt.vexpl(:));
  L=load(fullfile(COL,sprintf('labels_%s.mat',s)));
  for pj=1:2
    p={'dg','da'}; p=p{pj};
    r2=double(getfield(load(fullfile(COL,sprintf('glm_%s_%s.mat',s,p)),'R2'),'R2')); r2=r2(:);
    for a=1:numel(areas)
      lab=[double(L.(['lh_' areas{a} '_REmanual'])); nLH+double(L.(['rh_' areas{a} '_REmanual']))];
      lab=lab(lab>=1&lab<=numel(r2));
      idx=lab(ecc(lab)>=4 & ecc(lab)<=8);              % ecc only, NO vexpl
      R2.(vk(s)).(areas{a}).(p)=median(r2(idx),'omitnan');
      if strcmp(s,'sub-0037') && strcmp(areas{a},'V1')  % reconciliation vs 16.9/2.83
        idxv=lab(ecc(lab)>=4 & ecc(lab)<=8 & vex(lab)>0.1);
        R2.recon.(p)=[median(r2(lab),'omitnan') median(r2(idx),'omitnan') median(r2(idxv),'omitnan')];
      end
    end
    R2.(vk(s)).whole.(p)=median(r2,'omitnan');
  end
end

%% raw % beta from CSV : visual_area AND pRF_ecc in [4,8], no threshold
opts=detectImportOptions(CSVp); opts.VariableNamingRule='preserve';
opts.SelectedVariableNames={'subject','visual_area','pRF_ecc', ...
 'cartexp_horizontal_stationary','cartexp_vertical_stationary','cartexp_rightleaning_grating_stationary','cartexp_leftleaning_grating_stationary','cartexp_blank', ...
 'polexp_annulus_grating_stationary','polexp_pinwheel_grating_stationary','polexp_cspiral_grating_stationary','polexp_ccspiral_grating_stationary','polexp_blank'};
T=readtable(CSVp,opts);
cartS=[T.cartexp_horizontal_stationary T.cartexp_vertical_stationary T.cartexp_rightleaning_grating_stationary T.cartexp_leftleaning_grating_stationary];
polS =[T.polexp_annulus_grating_stationary T.polexp_pinwheel_grating_stationary T.polexp_cspiral_grating_stationary T.polexp_ccspiral_grating_stationary];
bDG=mean(cartS,2)-T.cartexp_blank; bDA=mean(polS,2)-T.polexp_blank;
va=string(T.visual_area); ec=T.pRF_ecc; sj=string(T.subject);
BETA=struct();
for si=1:numel(subs)
  s=subs{si}; sm=sj==s;
  for a=1:numel(areas)
    m=sm & va==areas{a} & ec>=4 & ec<=8;
    BETA.(vk(s)).(areas{a})=[mean(bDG(m),'omitnan') mean(bDA(m),'omitnan') nnz(m)];
  end
  BETA.(vk(s)).whole=[mean(bDG(sm),'omitnan') mean(bDA(sm),'omitnan') nnz(sm)];
end

%% print
fprintf('\n=== median GLM R2 (%%) and mean RAW %% beta, area in 4-8 deg (no var threshold) ===\n');
hdr={'V1','V2','V3','whole'};
fprintf('%-14s |', 'subject');
for h=hdr, fprintf(' %-27s |', [h{1} ': R2dg R2da  bDG    bDA']); end; fprintf('\n');
for si=1:numel(subs)
  s=subs{si}; fprintf('%-14s |', s);
  for a=1:numel(hdr)
    r=R2.(vk(s)).(hdr{a}); b=BETA.(vk(s)).(hdr{a});
    fprintf(' %4.1f %4.1f %6.3f %6.3f |', r.dg, r.da, b(1), b(2));
  end
  fprintf('\n');
end
fprintf('\nsub-0037 V1 R2 reconciliation (dg/da):\n');
fprintf('  whole V1 (no ecc)           : %.1f / %.1f\n', R2.recon.dg(1), R2.recon.da(1));
fprintf('  V1 & ecc4-8 (this table)    : %.1f / %.1f\n', R2.recon.dg(2), R2.recon.da(2));
fprintf('  V1 & ecc4-8 & pRF-R2>0.1    : %.1f / %.1f   <- the 16.9/2.83 patch\n', R2.recon.dg(3), R2.recon.da(3));
