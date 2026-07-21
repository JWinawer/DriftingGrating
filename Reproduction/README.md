# Reproduction of Figures 5–8 from `allsubjectsTable.csv`

This folder reproduces the manuscript's Figures 5–8 (V1 orientation asymmetries, Cartesian vs.
polar gratings) **starting from a single tidy table**, `../Support/allsubjectsTable.csv`
(git-ignored). See `../AGENTS.md` for what those figures are, and the plan the work follows.

Two independent routes:

- **`cleanroom/`** — Path A. Standalone MATLAB that goes CSV → filter → per-vertex
  orientation−blank contrasts → polar-angle-wedge medians → the 4 asymmetries → Figs 5–8. It does
  **not** call the original `AnalysisCode` analysis functions, so it is an independent check.
- **`bridge/`** — Path B. Regenerates the group `.mat` arrays the existing `AnalysisCode`
  stage-04 scripts expect, then runs that original code. (Built after Path A validates.)

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

## Subjects (fixed order, those who did both experiments)
`sub-0037, sub-0201, sub-0255, sub-wlsubj123, sub-wlsubj124, sub-0395, sub-0426, sub-0250`

## Path A (cleanroom) — how to run
From `cleanroom/` in MATLAB: `run_all_repro` regenerates Figs 5–8 (both variants) and prints the
validation table. Or run pieces: `run_fig5_6`, `run_fig7`, `run_fig8`,
`validate_against_manuscript`. First call builds the V1 cache (`_cache/v1.mat`, ~10 s); later
calls reuse it. Outputs land in `figures/cleanroom/`.

## Status
- [x] Path A: data path validated — 7/8 asymmetries reproduce the manuscript exactly
- [x] Path A: Figs 5, 6, 7, 8 generated (z-scored + raw variants)
- [x] Finding: `da` horizontal-vertical diverges (clean ~0 vs manuscript −0.45) — see
      [FINDINGS.md](FINDINGS.md)
- [ ] Path B: bridge CSV → existing code; resolve `da` H-V per-θ definitively
