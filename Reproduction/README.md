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

## Key result

The reproduction confirms the analysis **except** for the two first-harmonic *derived* (reference-
frame-mismatched) asymmetries, which carry a polar-angle-ordering bug (a cardinal-meridian
horizontal/vertical and radial/tangential swap) in `AnalysisCode/.../compute_derivativeDirections.m`.
Corrected values: `da` horizontal-vertical ≈ 0 (not −0.45) and `dg` radial-tangential ≈ 0.23.
Full write-up and evidence in [FINDINGS.md](FINDINGS.md).

## Status
- [x] Path A: data path validated — direct asymmetries + `dg` derived reproduce the manuscript
- [x] Path A: Figs 5, 6, 7, 8 generated (z-scored + raw)
- [x] Path B: CSV bridged into the existing pipeline; cross-path arrays identical (0.0 diff)
- [x] Path B: Figs 5/6 regenerated through the real `plot1/plot2`; da H-V artifact resolved
- [x] Finding documented in [FINDINGS.md](FINDINGS.md)
- [ ] **Next:** z-scored vs non-z-scored comparison — see [NEXT_STEPS.md](NEXT_STEPS.md)
