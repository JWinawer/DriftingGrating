# Server extraction — for whoever has the data volume mounted

Self-contained MATLAB that reads data on the server and writes a small folder to send back.
**These scripts only read.** They write nothing into the data directories and modify nothing on the
server; everything they write goes to `~/dg_collect/` (or wherever you point them).

They deliberately **filter nothing** — whole surface, every retinotopy map, every label. Three
early rounds were each cut short because something had been filtered out before saving, and each
one cost another trip.

## The scripts

Each is a separate one-pass extraction; none depends on the others. Run the one that matches what
is needed.

| script | feeds | status |
|---|---|---|
| `collect_everything.m` | the GLM quality review — the main extraction | **done** 2026-07-24 (the 8 who did both experiments) and 2026-08-27 (the 5 `dg`-only observers) — **all 13** |
| `collect_prf_replicate.m` | the pRF polar-angle precision control (second, independent pRF solution) | **done** |
| `collect_runwise_betas.m` | per-run condition betas, V1, for the within-observer error estimate | **done** |
| `collect_runwise_betas_areas.m` | the same for eight maps — V1 V2 V3 V3a V3b hV4 MT MST | **done** |
| `collect_gain_areas.m` | per-observer × **map** pRF gain (the V1-only `gainSummary.csv` special-cased V1) | **done** |
| `collect_timeseries.m` | preprocessed fsnative BOLD, to calibrate the Fig-4A run-mismatch control | **done** for sub-0037/`dg` (8 runs, 2.35 GB); see Open items in `../../AGENTS.md` |
| `extract_for_transfer.m` | the original, narrower extraction | superseded by `collect_everything.m`; kept as a record |

Which analysis each one feeds is documented in the file header, and in
[`../SPECIFICATION.md`](../SPECIFICATION.md),
[`../supplement/SUPPLEMENT_harmonic_model.md`](../supplement/SUPPLEMENT_harmonic_model.md)
and [`../local_qc/DATA_QUALITY.md`](../local_qc/DATA_QUALITY.md).

## Running one

All of them take the same shape. Using `collect_everything` as the example:

**1. Preflight — seconds, writes nothing.**

```bash
matlab -batch "addpath(genpath('<repo>/Reproduction/server_extract')); collect_everything([],[],[],struct('preflightOnly',true))"
```

Run it from any folder; the script finds `AnalysisCode` on its own and calls `setup_user` from
there. It prints an `[ ok ]` / `[FAIL]` line per check. The two failures worth anticipating:

- **`MRIread not on the MATLAB path`** — FreeSurfer is not where `setup_user.m` expects. Check the
  `freesurferDir` for your machine in `AnalysisCode/general/setup_user.m`. Hard stop: without it no
  retinotopy maps are exported at all.
- **`bidsDir NOT reachable`** — the volume is not mounted, or is at a different path. Pass yours:
  `collect_everything('/your/path/data_bids/')`.

**2. Run it.**

```bash
matlab -batch "addpath(genpath('<repo>/Reproduction/server_extract')); collect_everything"
```

Expect it to be slow: it reads 16 `results.mat` files of roughly 500 MB each, and they are v7
format, so MATLAB cannot load parts of them selectively. Over a network mount that is the whole
cost. Measured mount throughput has ranged 0.8–3.4 MB/s (Abu Dhabi → New York), so treat any single
throughput measurement as unreliable for planning.

Run it under `caffeinate` — a sleep dropped the mount mid-run once, and **a dropped mount returns
empty rather than erroring**, so verify file sizes, not just presence.

**3. Check, then send.** Open `~/dg_collect/manifest.csv` and look for any row whose `status` is
not `ok`. `missing`, `bad-structure` and `read-failed` are worth mentioning when you send the
folder — they may be expected, but that cannot be told from off-server. Then `rsync -avz
~/dg_collect/ <destination>`; the folder is about 1.2 GB.

Common options:

```matlab
collect_everything([],[],[],struct('includeBetas',true))  % adds ~8 GB, usually not needed
collect_everything([],[],[],struct('force',true))         % collect what is reachable despite preflight
collect_everything([],[],'/path/to/output')               % somewhere other than ~/dg_collect
```

**Collecting a different observer set.** `dg` was run on 13 observers and `da` on 8; the script's
built-in list is the 8 who did both. Pass `subjects` and `projects` for any other set — this is how
the five `dg`-only observers were added on 2026-08-27:

```matlab
collect_everything([],[],[],struct( ...
    'subjects', {{'sub-0442','sub-wlsubj121','sub-wlsubj127','sub-0397','sub-0427'}}, ...
    'projects', {{'dg'}}, ...
    'collectDesign', false))
```

Three behaviours worth knowing:

- **`skipExisting` is on by default**, so an interrupted run resumes without re-reading what it
  already has. It treats a file at or below `emptyMB` (0.5 MB) as absent, because a dropped mount
  writes an empty file rather than erroring.
- **`collectDesign` is on by default but the design tree is subject-independent.** `collect_design`
  globs `expOutDir/<proj>/**` and never uses the subject in the path, so it is now collected **once
  per project** rather than once per observer. Set it `false` when the design tree is already
  collected.
- **`Total written` counts only rows marked `ok`.** The manifest also logs each `modelOutput.mat`
  as `present-not-copied` — on record, deliberately not copied — and those are reported on their own
  line. Before 2026-08-27 they were added into the total, which reported ~10 GB written for a run
  that wrote 452 MB, and made five ordinary rows look like five failures.

**`manifest.csv` is overwritten by each run.** It describes only the observers in that call, so copy
it aside first if it holds a record you want to keep. The current file is the merge of the two runs
and covers all 13; the individual runs are kept as `manifest_8subj_2026-07-24.csv` and
`manifest_5subj_2026-08-27.csv` in `~/dg_collect/`.

**`collect_gain_areas` rewrites `gain_areas_summary.csv` from whatever subjects it is given**, so
run it on the **whole** set, not just the new observers, or the summary silently loses the rest. Its
per-subject `.mat` files are cached and reloaded, so already-collected observers cost no server reads.

## Deliberately avoided: polar angle

The V1 patch is built from the label plus **eccentricity and `vexpl` only** — never polar angle.
Polar angle in this project has a Benson-vs-conventional convention difference that has already
caused one false-alarm bug report
([`../STIMULUS_CONVENTIONS.md`](../STIMULUS_CONVENTIONS.md)). Nothing here needs it, so
nothing here touches it.

## Two questions still unanswered from off-server

Worth a sentence in your reply if you know. Also tracked in
[`../../AGENTS.md`](../../AGENTS.md) §5:

1. **Are `/Volumes/Vision/UsersShare/Rania/Project_dg/` and `/Volumes/server/Projects/Project_dg/`
   the same data, or two copies?** Different scripts in `AnalysisCode` point at each — `setup.json`,
   `meanWithinLabel.m` and `lme1_fit.m` use the first; `roi2image*.m` and
   `analyzeROI_anotherMetric.m` use the second. If they are separate copies that have drifted,
   some figures were made from different data than others.
2. **What is the actual server?** Everything in the code is a local mount point (`/Volumes/...`),
   so nothing records the machine or share.
