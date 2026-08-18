# Reproduction of Figures 5–8 from `allsubjectsTable.csv`

This folder reproduces the manuscript's Figures 5–8 (V1 orientation asymmetries, Cartesian vs.
polar gratings) **starting from a single tidy table**, `../Support/allsubjectsTable.csv`
(git-ignored). See `../AGENTS.md` for what those figures are, and the plan the work follows.

> **Newest first: [`HARMONIC_MODEL.md`](HARMONIC_MODEL.md) (2026-08-18)** and its reader-facing
> version [`supplement/SUPPLEMENT_harmonic_model.md`](supplement/SUPPLEMENT_harmonic_model.md).
> A per-vertex harmonic model that replaces the eight wedges with each vertex's own pRF polar
> angle, to test whether the Cartesian-vs-polar difference is within-ROI geometry rather than a
> context effect. It is not: geometry accounts for 6–8%. Also establishes that the Fig-7 LME and
> the wedge-average route give **identical** estimates (the design is balanced and the codes
> orthogonal), so the choice between them is only a choice of error bar.
>
> **Before that: [`local_qc/REPORT.md`](local_qc/REPORT.md) (2026-07-24)** is the current
> authority on data quality and on the z-scoring decision, and supersedes parts of the docs in
> this folder. Headline: the "blank" condition is **full-field pink noise, not a baseline**
> (so a low GLM R² means weak differentiation among stimuli, not an absent response), there are
> **no bad observers and no processing errors**, and the recommendation is to **use the raw,
> non-z-scored analyses**. The orientation asymmetries themselves are unaffected.

Two independent routes:

- **`cleanroom/`** — Path A. Standalone MATLAB that goes CSV → filter → per-vertex
  orientation−blank contrasts → polar-angle-wedge medians → the 4 asymmetries → Figs 5–8. It does
  **not** call the original `AnalysisCode` analysis functions, so it is an independent check.
- **`bridge/`** — Path B. Regenerates the group `.mat` arrays the existing `AnalysisCode`
  stage-04 scripts expect, then runs that original code. (Built after Path A validates.)

Two further folders, added later as the work moved from reproduction to data quality:

- **`server_extract/`** — a self-contained, read-only extraction (`collect_everything.m`) to run
  on the machine with `/Volumes/Vision` mounted. Deliberately filters nothing. See `RUNME.md`.
- **`local_qc/`** — the GLM data-quality review built on that extraction: `REPORT.md`,
  `glm_summary.csv`, per-question scripts, and a draft manuscript caveat paragraph.

Both produce **two variants** of each figure: **z-scored** (each vertex's betas standardized
across its 13 conditions before analysis — the manuscript's default) and **non-z-scored**. Since
the per-vertex mean cancels in an (orientation − blank) subtraction, the z-scored contrast is
just `(orientation − blank) / beta_std` and the raw contrast is `(orientation − blank)`; the
`*_beta_std` columns in the CSV supply the divisor (no GLM re-run needed).

## Data facts (from `allsubjectsTable.csv`)
- ~2.22M rows = 8 subjects × 2 hemispheres × all visual areas. One row per surface vertex.
- Per-vertex: `pRF_angle`, `pRF_angle_bin` (nearest of 0/45/…/315), `pRF_ecc`, `pRF_r2` (a
  fraction), `pRF_sigma`, `visual_area`, `included` (a loose prefilter — **not** the figure set).
- Per experiment, 13 condition betas: `cartexp_*` (Cartesian / "dg") and `polexp_*` (polar /
  "da"), each = 8 motion + 4 stationary + 1 blank. Plus `dg_beta_mean/std`, `da_beta_mean/std`
  = the mean/(N−1)-std across that experiment's 13 conditions.

## Analysis inclusion filter (applied explicitly)
`visual_area == "V1"`, `4 ≤ pRF_ecc ≤ 8`, `pRF_r2 > 0.1`. Do **not** rely on the `included`
column (it is a superset of this set).

This is the filter for the **figure/asymmetry** analysis, where the 4–8° band is a
stimulus-matching constraint (Cartesian and polar gratings are spatial-frequency-matched only
near 6°). It is *not* the right filter for **fit-quality** questions, and the server extraction
deliberately applies neither restriction — see `GLM_QUALITY.md` §6a. Don't re-apply these when
asking whether a session produced a usable response.

## Subjects (fixed order, those who did both experiments)
`sub-0037, sub-0201, sub-0255, sub-wlsubj123, sub-wlsubj124, sub-0395, sub-0426, sub-0250`

## How to run

**Path A (cleanroom)** — from `cleanroom/` in MATLAB: `run_all_repro` regenerates Figs 5–8
(both variants) and prints the validation table. Or run pieces: `run_fig5_6`, `run_fig7`,
`run_fig8`, `validate_against_manuscript`. First call builds the V1 cache (`_cache/v1.mat`,
~10 s); later calls reuse it. Outputs → `figures/cleanroom/`.

**Path B (bridge)** — from `bridge/` in MATLAB (needs `AnalysisCode/` on the path, added
automatically):
- `resolve_da_HV` — runs the real `compute_derivativeDirections.m`; prints the per-θ da H-V
  comparison that pins the artifact.
- `run_pathB_values` — computes all 8 asymmetries through the real repo functions; prints
  existing-code vs clean-room vs manuscript.
- `run_pathB_figures` — regenerates Figs 5/6 through the real `plot1/plot2` (→ `figures/bridge/`).

## Key result — ⚠️ RETRACTED 2026-07-22

This folder previously reported a polar-angle-ordering bug in
`AnalysisCode/.../compute_derivativeDirections.m`. **That finding is wrong.** The original
code is correct and the manuscript reproduces exactly on all eight asymmetries. The apparent
discrepancy came from two bugs in this reproduction — swapped spirals in
`cleanroom/config_repro.m` and a Benson/conventional frame mismatch in `bridge/`.

See [`../AUDIT.md`](../AUDIT.md). Fix those two before trusting any *derived* asymmetry
computed here; the direct asymmetries were unaffected.

## Status
- [x] Path A: data path validated — direct asymmetries + `dg` derived reproduce the manuscript
- [x] Path A: Figs 5, 6, 7, 8 generated (z-scored + raw)
- [x] Path B: CSV bridged into the existing pipeline; cross-path arrays identical (0.0 diff)
- [x] Path B: Figs 5/6 regenerated through the real `plot1/plot2`; da H-V artifact resolved
- [x] Finding documented in [FINDINGS.md](FINDINGS.md) — subsequently **retracted**, see AUDIT.md
- [x] z-scored vs non-z-scored comparison — mechanism in [ZSCORE_FIG7.md](ZSCORE_FIG7.md),
      **decided** in [local_qc/REPORT.md](local_qc/REPORT.md) §4: use the raw variants
- [x] GLM fit quality, all 8 observers × both experiments — [GLM_QUALITY.md](GLM_QUALITY.md),
      then the full unfiltered review in [local_qc/REPORT.md](local_qc/REPORT.md): no errors,
      no bad data, all observers retained
- [x] Geometry-vs-context test — per-vertex harmonic model in
      [HARMONIC_MODEL.md](HARMONIC_MODEL.md), written up in
      [supplement/](supplement/SUPPLEMENT_harmonic_model.md): within-ROI geometry explains only
      6–8%, pRF angle error is an order of magnitude too small, the Cartesian-frame context
      effect stands, the polar-frame one is uninformative (not absent)
- [ ] **Open:** add a GLM-`R2` column to `allsubjectsTable.csv` so the vertex filter can screen
      on it alongside `pRF_r2` ([NEXT_STEPS.md](NEXT_STEPS.md) step 5)
- [ ] **Open:** confirm the fixed stimulus `rng` seed was intended
      ([local_qc/REPORT.md](local_qc/REPORT.md) §3)
