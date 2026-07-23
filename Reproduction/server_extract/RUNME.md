# For whoever has the data mounted — how to run this

Two commands. The first takes seconds and tells you whether the second will work.

## 1. Get the code

```bash
git pull
```

## 2. Check it will work (seconds, writes nothing)

```bash
matlab -batch "addpath(genpath('<repo>/Reproduction/server_extract')); collect_everything([],[],[],struct('preflightOnly',true))"
```

Replace `<repo>` with wherever this repository is on your machine.

**You can run this from any folder.** The script finds `AnalysisCode` on its own and
calls `setup_user` from there — you do not need to `cd` anywhere first.

It prints a block like:

```
=============== PREFLIGHT ===============
[ ok ] MRIread found (via setup_user('rania') in /.../AnalysisCode)
[ ok ] bidsDir reachable
[ ok ] expOutDir reachable
[ ok ] results.mat  16/16 found
[ ok ] retinotopy   8/8 found
[ ok ] label dirs   8/8 found
-----------------------------------------
Preflight OK.
=========================================
```

**If every line says `ok`, go to step 3.** Otherwise the `[FAIL]` lines say what to fix.
The two most likely:

- **`MRIread not on the MATLAB path`** — FreeSurfer is not where `setup_user.m` expects.
  Check the `freesurferDir` for your machine in `AnalysisCode/general/setup_user.m`.
  This is a hard stop, because without it no retinotopy maps are exported at all.
- **`bidsDir NOT reachable`** — the volume is not mounted, or it is at a different path.
  Pass yours: `collect_everything('/your/path/data_bids/')`.

## 3. Run it

```bash
matlab -batch "addpath(genpath('<repo>/Reproduction/server_extract')); collect_everything"
```

**This takes a while** — it reads 16 files of roughly 500 MB each, and over a network
mount that is the slow step. It writes about **1.2 GB** to `~/dg_collect/`.

Then send that folder back — rsync, ftp, whatever is easiest for something that size.

```bash
rsync -avz ~/dg_collect/ <destination>
```

## 4. Before you send it

Open `~/dg_collect/manifest.csv` and look for any row whose `status` is not `ok`. Rows
reading `missing`, `bad-structure`, or `read-failed` are worth mentioning when you send
the folder — they may be expected, but we cannot tell from here.

## Two questions we cannot answer from off-server

Worth a sentence in your reply if you know:

1. **Are `/Volumes/Vision/UsersShare/Rania/Project_dg/` and
   `/Volumes/server/Projects/Project_dg/` the same data, or two copies?** Different
   scripts in `AnalysisCode` point at each — `setup.json`, `meanWithinLabel.m` and
   `lme1_fit.m` use the first, while `roi2image*.m` and `analyzeROI_anotherMetric.m` use
   the second. If they are separate copies that have drifted apart, some figures were
   made from different data than others, and we would want to know which is authoritative.

2. **What is the actual server?** Everything in the code is a local mount point
   (`/Volumes/...`), so nothing records the machine or share. If direct access is ever
   worth setting up, we would need the hostname and share path.

## Options, if asked for

```matlab
% include the betas too -- adds ~8 GB, usually not needed
collect_everything([],[],[],struct('includeBetas',true))

% collect whatever is reachable even though preflight failed
collect_everything([],[],[],struct('force',true))

% somewhere other than ~/dg_collect
collect_everything([],[],'/path/to/output')
```

## What it does and does not do

It **only reads**. It writes nothing into the data directories and modifies nothing on
the server. Everything it writes goes to `~/dg_collect/`.

It deliberately **filters nothing** — whole surface, every retinotopy map, every label.
Three previous rounds of this were each cut short because something had been filtered out
before saving, and each one cost another trip. The output is v7.3, so individual
variables can be read later without loading whole files.
