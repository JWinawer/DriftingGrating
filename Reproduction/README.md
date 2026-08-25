# Reproduction — independent recomputation, and where the analysis now lives

This folder recomputes the draft's Figures 5–8 (V1 orientation asymmetries, Cartesian vs. polar
gratings) **starting from a single tidy table**, `../Support/allsubjectsTable.csv` (git-ignored), and
is where every follow-up analysis has since been built.

**Read [`../AGENTS.md`](../AGENTS.md) first** — it holds the standing facts, what is settled, what is
open, and an index of every document here. This file covers only *how the code is organised and run*.

| document | what it owns |
|---|---|
| [`SPECIFICATION.md`](SPECIFICATION.md) | the settled analysis specification, and the evidence for each choice |
| [`RESULTS.md`](RESULTS.md) | every current number |
| [`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) | five closed methodological choices |
| [`MISSING_DATA.md`](MISSING_DATA.md) | what empty cells do, and why some maps are not reported by polar angle |
| [`STIMULUS_CONVENTIONS.md`](STIMULUS_CONVENTIONS.md) | stimuli, condition indices, angle frames |
| [`local_qc/DATA_QUALITY.md`](local_qc/DATA_QUALITY.md) | GLM data quality |
| [`local_qc/RELIABILITY.md`](local_qc/RELIABILITY.md) | reliability of the analysed measurements |

## The four code paths

- **`cleanroom/`** — the primary one. Standalone MATLAB: CSV → filter → per-vertex
  orientation−blank contrasts → polar-angle aggregates → the four asymmetries → the figures. It does
  **not** call the original `AnalysisCode` analysis functions, so it is an independent check. The
  settled specification, the per-vertex harmonic model, the precision-weighted tables and every
  diagnostic live here.
- **`bridge/`** — regenerates the group `.mat` arrays the original `AnalysisCode` stage-04 scripts
  expect, then runs that original code. Built after the clean-room validated; cross-path arrays are
  identical (0.0 difference).
- **`server_extract/`** — self-contained, read-only extractions to run on the machine with
  `/Volumes/Vision` mounted. Deliberately filter nothing. See its [`README.md`](server_extract/README.md).
- **`local_qc/`** — the data-quality review built on those extractions.

Outputs land in `figures/` and `cleanroom/_cache/` (both git-ignored); tables and figures that are
meant to be kept are committed under `supplement/`.

## The data table

`../Support/allsubjectsTable.csv`, ~2.22M rows = 8 subjects × 2 hemispheres × all visual areas, one
row per surface vertex.

- Per vertex: `pRF_angle`, `pRF_angle_bin` (nearest of 0/45/…/315), `pRF_ecc`, `pRF_r2` (a fraction),
  `pRF_sigma`, `visual_area`, `included` (a loose prefilter — **not** the figure set).
- Per experiment, 13 condition betas: `cartexp_*` (Cartesian / `dg`) and `polexp_*` (polar / `da`),
  each 8 motion + 4 stationary + 1 blank. Plus `dg_beta_mean/std`, `da_beta_mean/std` = the
  mean and (N−1)-std across that experiment's 13 conditions. (The `*_beta_std` columns supply the
  z-scoring divisor, so the z-scored sensitivity variants need no GLM re-run — but z-scoring is not
  used; see [`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) §1.)

**Subjects, fixed order** (those who did both experiments):
`sub-0037, sub-0201, sub-0255, sub-wlsubj123, sub-wlsubj124, sub-0395, sub-0426, sub-0250`

## The inclusion filter

`visual_area == "V1"`, `4 ≤ pRF_ecc ≤ 8`, `pRF_r2 > 0.1`. Applied explicitly — do **not** rely on the
`included` column, which is a superset.

This is the filter for the **figure/asymmetry** analysis, where 4–8° is a stimulus-matching
constraint (the two grating sets are spatial-frequency-matched only near 6°). It is *not* the right
filter for **fit-quality** questions, and the server extractions deliberately apply neither
restriction. Don't re-apply it when asking whether a session produced a usable response.

The extrastriate supplement runs the same filter over eight maps, at 4–8° and 2–10°, with a measured
coverage criterion deciding which are reportable — [`SPECIFICATION.md`](SPECIFICATION.md) §7.

## How to run

**The settled specification** — everything currently reported. From `cleanroom/` in MATLAB:

| command | what it does |
|---|---|
| **`run_spec_outputs`** | **Figures 5/6, the profile, the hierarchy figure, and every table, in all three variants.** Needs `~/dg_collect/`; ~3 min |
| `spec_profiles`, `spec_tables`, `spec_areas_summary` | the pieces that driver runs, callable per map |

→ [`SPECIFICATION.md`](SPECIFICATION.md) §8 for the variants, the outputs and the code map.

**The earlier reproduction and the diagnostics** — still runnable, and the source of several results
the specification builds on:

| command | what it does |
|---|---|
| `run_all_repro` | Figs 5–8 (both variants) + the validation table |
| `run_fig5_6`, `run_fig7`, `run_fig8` | individual figures |
| `validate_against_manuscript` | the eight asymmetries against the draft's values |
| `test_harmonic_model` | assertions guarding the harmonic model — **run before interpreting it** |
| `run_harmonic_model(false)` | the per-vertex harmonic model (`true` = z-scored sensitivity check) |
| `precision_weighted_table` | the per-observer precision-weighted tables |
| `asymmetry_tables` | the asymmetry and context-effect tables on the pre-specification ROI route |
| `splithalf_reliability` | the split-half reliability of the analysed profile |

**The missing-data diagnostics** — run on demand, none of them part of `run_spec_outputs`.
→ [`MISSING_DATA.md`](MISSING_DATA.md) §8 for what each one answers.

| command | what it does |
|---|---|
| `cell_occupancy('area','MT')` | vertices per (observer × polar-angle ROI) for any map, from labels alone |
| `diagnose_cell_loss('donor','MT')` | delete a map's empty cells from V1 and refit, both routes |
| `diagnose_loss_structure('donor','MT')` | is the resulting shift systematic, and how much data each map actually has |
| `diagnose_pooled_fit('donor','MT')` | one pooled group fit versus averaging per-observer fits |

The first call builds the V1 cache (`_cache/v1.mat`, ~10 s); later calls reuse it.

**The bridge** — from `bridge/` in MATLAB (`AnalysisCode/` is added to the path automatically):

| command | what it does |
|---|---|
| `run_pathB_values` | all 8 asymmetries through the real repo functions; prints existing-code vs. clean-room vs. draft |
| `run_pathB_figures` | regenerates Figs 5/6 through the real `plot1`/`plot2` |
| `resolve_da_HV` | runs the real `compute_derivativeDirections.m`; prints the per-θ `da` H−V comparison |

## One thing not to repeat

This folder once reported a polar-angle-ordering bug in
`AnalysisCode/.../compute_derivativeDirections.m`. **That report is retracted.** The original code is
correct; the apparent discrepancy came from two bugs *inside this reproduction*, both fixed.
[`STIMULUS_CONVENTIONS.md`](STIMULUS_CONVENTIONS.md) is the standing account and owns every angle
convention in the chain.
