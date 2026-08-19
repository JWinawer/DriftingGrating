# Server extraction — for whoever has access to the data volume

Self-contained MATLAB that reads data on the server and writes a small folder to send back.
**It only reads. It writes nothing into the data directories and modifies nothing.**

## Which script

Four scripts have accumulated here, each a separate one-pass extraction. Run the one that
matches what is needed; none depends on the others.

| script | feeds | status |
|---|---|---|
| `collect_everything.m` | the GLM quality review (`../local_qc/REPORT.md`) — the main extraction | done 2026-07-24, all 8 observers |
| `collect_prf_replicate.m` | the pRF polar-angle precision control (`../HARMONIC_MODEL.md`, supplement §S6) | done |
| `collect_runwise_betas.m` | per-run condition betas for the within-observer error estimate (supplement §S5.5) | done |
| `extract_for_transfer.m` | the original narrower extraction, superseded by `collect_everything.m` | kept as a record |

`RUNME.md` is the step-by-step for `collect_everything.m` specifically.

## Why this exists

Two open questions about Figures 5–8 both need data that is not in `allsubjectsTable.csv`, and
this collects both in one pass so there is no second round trip.

1. **Are the GLM fits sound for every observer?** No GLMsingle quality metric enters the analysis
   pipeline at any stage — the only quality filter is on the *pRF* fit (`pRF_r2 > 0.1`). Two
   observers look anomalous in the polar experiment: **sub-0201**, whose blank beta is larger than
   all 12 stimulus betas in *both* experiments, and **sub-0037**, which shows no differentiation
   between any conditions in `da` despite responding strongly in `dg`. Those patterns are what a
   bad run or a motion artefact produces, and `R2run` would show it.
2. **Is there a usable per-observer BOLD gain?** We want to normalize observers into commensurate
   units, but the divisor has to be independent of the conditions being analysed. The retinotopy
   scan is a separate session, so `prfvista_mov` is the natural source — we just don't know from
   off-server which maps that pipeline saves. The script inventories them.

## How to run

```matlab
cd <wherever you put this file>
extract_for_transfer
```

If the data is not at the default path, pass it explicitly:

```matlab
extract_for_transfer('/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/')
```

It needs `setup_user`, `read_label`, `get_surfsize` and `MRIread` on the path — the same helpers
`meanWithinLabel.m` uses. Nothing else.

**Expect it to take a while.** It reads 16 `results.mat` files of roughly 500 MB each. They are v7
format, so MATLAB cannot load parts of them selectively — each is read whole and then discarded.
Over a network mount this is the slow step.

## What to send back

The folder it prints at the end (`glm_qc_for_transfer/` by default). Zip it.

- With the V1 restriction working: a few MB per file at most, so tens of MB total. (The first
  run of this script also restricted to 4–8°, giving ~50 kB files. It now keeps all of V1 and
  saves each vertex's eccentricity, so the eccentricity band can be chosen after the fact —
  bigger files, but no further round trips.)
- If the V1 restriction fails (it warns), files are ~10 MB each and the folder is ~160 MB. Still
  far smaller than the 8 GB of `results.mat`, but too big to email — use a file transfer.
- **If all else fails, `summary.csv` alone answers most of question 1.** It is a few kB.

## What to check before sending

The script prints a line per subject × experiment. Worth a glance:

- Any row saying **`MISSING results.mat`** — is the path right for that subject?
- Any row saying **`UNEXPECTED structure`** — that subject was run with a different
  `hRF_setting` and has no GLMsingle TYPED output. Worth knowing.
- The **`V1 patch UNAVAILABLE`** warning — everything still gets saved, but the files are much
  bigger and we lose the restriction to the analysed patch. The error message says why.

## A note on one thing that is deliberately avoided

The patch is built from the V1 label plus **eccentricity and `vexpl` only** — never polar angle.
Polar angle in this project has a Benson-vs-conventional convention difference that has already
caused one false-alarm bug report (see `../AUDIT.md`). Nothing here needs polar angle, so it
does not touch it.
