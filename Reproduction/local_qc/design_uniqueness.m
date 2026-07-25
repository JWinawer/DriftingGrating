% Design-uniqueness check: fingerprint each per-run trialMat and find collisions.
% NOTE: with a FIXED rng seed (=0, as used here) all sessions are IDENTICAL by run
% position, so "duplicates" are EXPECTED, not a bug. See REPORT.md.
COL = fullfile(getenv('HOME'),'dg_collect','design');
for pj = {'dg','da'}
  d = dir(fullfile(COL,pj{1},'*_design_*.mat'));
  keys = containers.Map('KeyType','char','ValueType','any');
  for i=1:numel(d)
    S = load(fullfile(d(i).folder,d(i).name),'expDes');
    k = mat2str(S.expDes.trialMat);          % mat2str treats NaN==NaN (intended here)
    if isKey(keys,k), keys(k)=[keys(k) i]; else, keys(k)=i; end
  end
  fprintf('%s: %d files -> %d unique trialMats\n', pj{1}, numel(d), keys.Count);
end
