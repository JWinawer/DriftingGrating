# Reproduction — independent recomputation of Figures 5–8

This folder recomputes the manuscript's Figures 5–8 (V1 orientation asymmetries, Cartesian vs.
polar gratings) **starting from a single tidy table**, `../Support/allsubjectsTable.csv`
(git-ignored), and is where every follow-up analysis has since been built.

**Read [`../AGENTS.md`](../AGENTS.md) first** — it holds the standing facts, what is settled, what
is open, and an index of every document here. This file covers only *how the code is organised and
run*.

## The four routes

- **`cleanroom/`** — Path A, the primary one. Standalone MATLAB: CSV → filter → per-vertex
  orientation−blank contrasts → polar-angle-wedge aggregates → the 4 asymmetries → Figs 5–8. It
  does **not** call the original `AnalysisCode` analysis functions, so it is an independent check.
  The per-vertex harmonic model, the precision-weighted tables, and every diagnostic live here too.
- **`bridge/`** — Path B. Regenerates the group `.mat` arrays the existing `AnalysisCode` stage-04
  scripts expect, then runs that original code. Built after Path A validated; cross-path arrays are
  identical (0.0 difference).
- **`server_extract/`** — self-contained, read-only extractions to run on the machine with
  `/Volumes/Vision` mounted. Deliberately filter nothing. See its `README.md`.
- **`local_qc/`** — the data-quality review built on those extractions.

Outputs land in `figures/` and `cleanroom/_cache/` (both git-ignored); tables that are meant to be
kept are committed under `supplement/`.

## Data facts (from `allsubjectsTable.csv`)

- ~2.22M rows = 8 subjects × 2 hemispheres × all visual areas. One row per surface vertex.
- Per-vertex: `pRF_angle`, `pRF_angle_bin` (nearest of 0/45/…/315), `pRF_ecc`, `pRF_r2` (a
  fraction), `pRF_sigma`, `visual_area`, `included` (a loose prefilter — **not** the figure set).
- Per experiment, 13 condition betas: `cartexp_*` (Cartesian / `dg`) and `polexp_*` (polar / `da`),
  each = 8 motion + 4 stationary + 1 blank. Plus `dg_beta_mean/std`, `da_beta_mean/std` = the
  mean/(N−1)-std across that experiment's 13 conditions. (Those `*_beta_std` columns supply the
  z-scoring divisor, so the z-scored sensitivity variants need no GLM re-run — but z-scoring is not
  used; see [`WHY_NOT_ZSCORE.md`](WHY_NOT_ZSCORE.md).)

**Subjects, fixed order** (those who did both experiments):
`sub-0037, sub-0201, sub-0255, sub-wlsubj123, sub-wlsubj124, sub-0395, sub-0426, sub-0250`

## Inclusion filter

`visual_area == "V1"`, `4 ≤ pRF_ecc ≤ 8`, `pRF_r2 > 0.1`. Applied explicitly — do **not** rely on
the `included` column, which is a superset.

This is the filter for the **figure/asymmetry** analysis, where 4–8° is a stimulus-matching
constraint (the two grating sets are spatial-frequency-matched only near 6°). It is *not* the right
filter for **fit-quality** questions, and the server extractions deliberately apply neither
restriction. Don't re-apply it when asking whether a session produced a usable response.

The extrastriate supplement runs the same filter over eight maps, at 4–8° and 2–10°, with a
measured coverage criterion deciding which are reportable — [`EXTRASTRIATE.md`](EXTRASTRIATE.md) §6.

## How to run

**Path A (cleanroom)** — from `cleanroom/` in MATLAB:

| command | what it does |
|---|---|
| `run_all_repro` | Figs 5–8 (both variants) + the validation table |
| `run_fig5_6`, `run_fig7`, `run_fig8` | individual figures |
| `validate_against_manuscript` | the eight asymmetries against the manuscript values |
| `test_harmonic_model` | assertions guarding the harmonic model — **run before interpreting it** |
| `run_harmonic_model(false)` | the per-vertex harmonic model, raw variant (`true` = z-scored sensitivity check) |
| `precision_weighted_table` | the per-observer precision-weighted tables |
| `asymmetry_tables` | the asymmetry and context-effect tables on the current route |
| **`run_spec_outputs`** | **Figures 5/6 and every table under the settled specification** — see [`SPEC_FIGURES.md`](SPEC_FIGURES.md). Needs `~/dg_collect/`; ~5 min. |
| `spec_profiles`, `spec_tables`, `spec_areas_summary` | the pieces that driver runs, callable per map |

The first call builds the V1 cache (`_cache/v1.mat`, ~10 s); later calls reuse it.

**Path B (bridge)** — from `bridge/` in MATLAB (`AnalysisCode/` is added to the path
automatically):

| command | what it does |
|---|---|
| `run_pathB_values` | all 8 asymmetries through the real repo functions; prints existing-code vs. clean-room vs. manuscript |
| `run_pathB_figures` | regenerates Figs 5/6 through the real `plot1`/`plot2` |
| `resolve_da_HV` | runs the real `compute_derivativeDirections.m`; prints the per-θ `da` H−V comparison |

## One thing not to repeat

This folder once reported a polar-angle-ordering bug in
`AnalysisCode/.../compute_derivativeDirections.m`. **That report is retracted.** The original code
is correct and the manuscript reproduces on all eight asymmetries; the apparent discrepancy came
from two bugs *inside this reproduction* — swapped c-/cc-spirals in `cleanroom/config_repro.m` and
a Benson-vs-conventional polar-angle frame mismatch in `bridge/`. Both are fixed.

[`AUDIT.md`](AUDIT.md) is the standing account and owns every angle convention in the chain.
