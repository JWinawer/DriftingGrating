# AGENTS.md — Target analyses (Figures 5–8)

This file lets a future agent understand the analyses behind **Figures 5–8** of the
manuscript *"Local orientation asymmetries in V1 depend on global stimulus properties"*
(Ezzo, Carrasco, Rokers, Winawer) **without reading the manuscript**. The draft itself is
in `Manuscript draft/` (git-ignored, not tracked).

These four figures are the analyses of active interest.

## Issues reviewed, and where each was resolved

Six questions have been worked through. **All are resolved**; the documents below are the standing
accounts, and the working documents they replace are in `Reproduction/_archive/`.

| # | issue | question | resolution | document |
|---|---|---|---|---|
| 1 | **Z-scoring** | Should the per-vertex betas be z-scored before analysis? | **No.** The blank is pink noise, so `beta_std` is not a gain. | [`Reproduction/WHY_NOT_ZSCORE.md`](Reproduction/WHY_NOT_ZSCORE.md) |
| 2 | **Data quality** | Are there outlier observers or datapoints to exclude? | **No.** No processing errors, no bad runs; all 8 observers retained. Measurement reliability of the analysed quantities quantified from runs (§2.7). | [`Reproduction/local_qc/REPORT.md`](Reproduction/local_qc/REPORT.md) |
| 3 | **Angle conventions** | Are the absolute and location-dependent angle codings correct? | **Yes.** Code is consistent and correct; an earlier bug report is retracted. | [`Reproduction/AUDIT.md`](Reproduction/AUDIT.md) |
| 4 | **Harmonic model** | Does ROI binning manufacture a context effect that is really local stimulus geometry relative to each vertex's pRF? | **No.** Geometry accounts for 6–8%; pRF angle error is an order of magnitude too small. | [`Reproduction/HARMONIC_MODEL.md`](Reproduction/HARMONIC_MODEL.md), [supplement](Reproduction/supplement/SUPPLEMENT_harmonic_model.md) |
| 5 | **Context effects** | What does a within-subject assessment support? | **Cartesian-frame asymmetries yes; polar-frame uninformative** — absence of evidence, not evidence of absence. | [`Reproduction/HARMONIC_MODEL.md`](Reproduction/HARMONIC_MODEL.md) Result 3 |
| 6 | **LME** | Does the Fig-7 mixed model add anything? | **No.** Identical estimates to the subtraction route; its SE is anti-conservative. Recommend omitting from the manuscript. | [`Reproduction/LME.md`](Reproduction/LME.md) |
| 7 | **Extrastriate** | Does the V1 route work in maps with empty polar-angle ROIs? | **No, and it is fixed.** Empty cells are not random in polar angle, so the ROI average changes estimand; de-trending by the fitted profile repairs it and is exact on complete data. | [`Reproduction/EXTRASTRIATE.md`](Reproduction/EXTRASTRIATE.md) |

Two facts from issue 2 change how everything else should be read, and are worth stating up front:

1. **The "blank" condition is full-field pink noise, not a mean-luminance baseline.** Pink noise is
   on the screen for the entire run, so every reported response is grating − pink-noise and **there
   is no true baseline anywhere**. A low GLM R² therefore means weak *differentiation among* the
   stimuli, **not** an absent response. The orientation asymmetries are unaffected — the common
   blank term cancels in each.
2. **There are no bad observers and no processing errors.** The two observers previously flagged as
   anomalous (sub-0037, sub-0201) show normal MT motion-selectivity and V4 grating preference in the
   same sessions. The exclusion once argued for is **withdrawn**.

---

## Study in one paragraph

fMRI (n=8 observers, V1 defined by pRF retinotopy) measuring whether V1 BOLD **orientation
asymmetries are context-dependent on the global stimulus**. Two experiments with *locally
matched* stimuli (1 cpd at 6° eccentricity):

- **Cartesian gratings** — horizontal, vertical, two obliques. Project code name: **`dg`**
  (drifting grating).
- **Polar gratings** — pinwheel, annulus, CW spiral, CCW spiral. Project code name: **`da`**
  (drifting annulus).

Only **stationary** gratings are analyzed here (the GLM also estimates 8 motion directions,
reported in a future study). V1 vertices are binned into **8 polar-angle wedges** (centered
0,45,…,315°, each ±22.5°), restricted to **4–8° eccentricity** and **pRF R² > 0.1**.

**Four orientation asymmetries**, in two reference frames:

| Asymmetry | Reference frame | Harmonic | LME predictor coding |
|---|---|---|---|
| horizontal vs. vertical | Cartesian (absolute) | 1st | −1 / 0 / +1 |
| cardinal vs. oblique | Cartesian (absolute) | 2nd | −1 / +1 |
| radial vs. tangential | polar (gaze-relative) | 1st | −1 / 0 / +1 |
| polar-cardinal vs. polar-oblique | polar (gaze-relative) | 2nd | −1 / +1 |

**Core claim** (what the figures are meant to show): an asymmetry is **strong when its
reference frame matches the global stimulus, weak/absent when it does not**. Cartesian
asymmetries (BOLD *lower* for horizontal and for cardinal — the "inverted" direction) appear
for Cartesian gratings but vanish for polar; polar asymmetries (radial > tangential,
polar-cardinal > polar-oblique) appear for polar gratings but not Cartesian.

---

## The figures

In every figure, panels **A–D** are the four asymmetries in this order:
**A** = horizontal vs. vertical, **B** = cardinal vs. oblique, **C** = radial vs. tangential,
**D** = polar-cardinal vs. polar-oblique.

### Figure 5 — Independent asymmetries, **Cartesian** gratings (`dg`)
Each asymmetry quantified **independently** (not jointly modeled). Per polar-angle wedge, take
the median z-scored % signal change (orientation minus blank) across vertices, then average
across observers. **Top:** polar plot (radius = V1 % change; angle = visual-field polar angle),
one colored trace per orientation class. **Bottom:** per-observer pairwise plot with a
bootstrap 95% CI on the difference (1000 resamples of observers). Result: A and B significant
(H<V, cardinal<oblique); C and D absent for Cartesian gratings.

### Figure 6 — Independent asymmetries, **Polar** gratings (`da`)
Identical analysis and layout to Fig 5, on the polar-grating experiment. Result: Cartesian
asymmetries (A, B) weakened/absent; polar asymmetries present — C (radial > tangential)
significant, D (polar-cardinal > polar-oblique) a positive trend.

### Figure 7 — Joint **linear mixed-effects** model (bar plots)
All four asymmetries entered into **one** LME per experiment, to estimate each asymmetry's
weight while controlling for observer variability. 256 data points (8 observers × 4
orientations × 8 polar angles). Model:

```
y_im = β0 + β1·(horizVvert) + β2·(cardVobl) + β3·(radVtang) + β4·(polarcardVpolarobl)
       + b_0m (observer random intercept) + ε_im
```

Fit by maximum likelihood; fixed-effect coefficients = each asymmetry's weight (in
standardized units, relative to the mean). 68% CIs from a bootstrap (1000 iterations,
resampling 8 observers with replacement); global intercept removed to compare across
experiments. Run separately for `dg` and `da` → **panel A = Cartesian, panel B = polar**.
**This is the figure whose significance pattern differs between the z-scored and non-z-scored
versions.**

> **The random intercept does nothing to these estimates.** The design is balanced (4 ori × 8
> wedges × 8 observers, no empty cells) and the four codes are exactly orthogonal to each other
> and to the intercept, so each asymmetry contrast sums to zero within an observer and the
> observer term — shrunk or not, random or fixed — cancels out of it. The LME fixed effects equal
> the plain mean of the per-observer wedge-averaged contrasts to < 2e-16, on all four asymmetries
> in both experiments. So Fig 7 "controls for observer variability" only in the sense of changing
> the *standard error*, not the estimate — and it changes it in the wrong direction (DF = 502, the
> wedge-level count, not 7). See `Reproduction/HARMONIC_MODEL.md` Result 3.

### Figure 8 — Per-location context effects (data + model)
Summarizes the LME's predicted response for **each orientation at every polar-angle location**,
irrespective of the asymmetry groupings. One small polar plot per visual-field location: local
orientation on the angular axis, response magnitude as distance from origin; **black = observer
data, red = LME model estimate** (the 4 orientations are duplicated to close the plot).
Panel A = Cartesian (`dg`), panel B = polar (`da`). Highlights that the largest context effect
is on the **horizontal meridian** (suppressed horizontal + enhanced radial oppose each other).

---

## Code map

Analysis is MATLAB (R2022b). Paths in the scripts point at a mounted data volume
(`/Volumes/Vision/UsersShare/Rania/Project_dg/...`) and `~/Documents/GitHub/DriftingGrating`;
the fMRI data itself is **not** in the repo. `projectName` selects the experiment:
**`'dg'` = Cartesian, `'da'` = polar.**

**Data pipeline (bottom → top):**

1. **`AnalysisCode/01_process_singlesubjectGLM/main_singlesub.m`** — per subject, GLMsingle
   produces `betamaps` (vertices × 13 conditions: 8 motion + 4 orientations + 1 blank). Saves
   `betas_nonzscored.mat`; if the `normalize` flag is set, z-scores **per vertex across the 13
   conditions** (`zscore(betamaps_condOnly, 0, 2)`, ~line 158) and saves `betas_zscored.mat`.
   Blank-subtracted orientation contrasts (e.g. `orientation_minus_baseline`) are computed and
   stored in `results.mat` → `results.contrasts`. **This is where z-scoring is applied and
   propagates downstream.**
2. **`AnalysisCode/03_process_groupBetas/meanWithinLabel.m`** — loads each subject's
   `results.mat`, reads `results.contrasts.<name>` per vertex, bins vertices into V1 × 8
   polar-angle wedges (4–8°, R² > 0.1), takes the median within each wedge, and saves group
   matrices **`medianBOLDpa`** (contrasts × polarAngles × ROIs × subjects) and **`medianBOLD`**
   (contrasts × ROIs × subjects). Whether these are z-scored is inherited from stage 1.
3. **`AnalysisCode/04_plot_betaAsymmetries/`** — consumes `medianBOLDpa` / `medianBOLD` to make
   the figures.

**Figure → script:**

| Figure | Script | Notes |
|---|---|---|
| **5** (Cartesian) | `04_plot_betaAsymmetries/plot_NeuralAsymmetries.m` with `projectName='dg'` | calls `plot1_experimentalCond` (polar plots) + `plot2_experimentalCond` (pairwise / mean-across-PA) |
| **6** (Polar) | same driver with `projectName='da'` | same two helpers |
| **7** (LME weights) | `04_plot_betaAsymmetries/lme1_fit.m` (also `lme1_fit_plotNewFormat.m`) | z-score vs. not currently toggled by hand here — see `ylabel('zscored BOLD psc')` (~line 557) and the commented `ylim([-0.4 0.4]) % for non zscore` (~line 569) |
| **8** (per-location) | `04_plot_betaAsymmetries/lme2_ploteachDirLoc.m` | data (black) + model (red) per location |

---

## Normalization — we do not z-score (decided 2026-07-24)

All analyses use **raw beta weights in percent signal change**. The standing account is
[`Reproduction/WHY_NOT_ZSCORE.md`](Reproduction/WHY_NOT_ZSCORE.md); the working documents it
replaces are in `Reproduction/_archive/`.

One line of reasoning: the "blank" is full-field pink noise, not a mean-luminance baseline, so
`beta_std` is not a gain — it conflates BOLD gain with tuning strength, and dividing the
(blank-independent) asymmetries by it puts the blank back in. No valid substitute divisor exists in
these 13 conditions.

Two things to carry forward:

- **Do not frame Fig 7 around which asymmetry is largest.** Subject bootstraps run 0.10–0.31 at
  n = 6 and 0.28–0.63 at n = 8 — unanimous in direction, never decisive. The cross-experiment
  comparison is what the data support.
- **The radial/tangential context effect is the one conclusion that turns on this choice** — it
  reaches significance z-scored and does not raw (`Reproduction/HARMONIC_MODEL.md` Result 3). State
  the dependence rather than burying it.

The earlier argument for excluding sub-0037 and sub-0201 as having "no measurable gain" is
**withdrawn**: `local_qc/REPORT.md` §2.5–2.6 clears both with positive evidence, and §1 removes the
premise. All 8 observers are retained.

## GLM fit quality — checked twice, sound
All 8 observers × both experiments were extracted from the server (2026-07-23, then again
unfiltered on 2026-07-24) and audited: no coding or processing error, uniform model parameters, no
bad run, no dropout, correct surface co-registration. See `Reproduction/local_qc/REPORT.md` §2 and
`Reproduction/_archive/GLM_QUALITY.md`. What remains genuinely open is that **no GLMsingle metric enters the
pipeline** — the only quality filter is still on the *pRF* fit (`pRF_r2 > 0.1`), and neither
`allsubjectsTable.csv` variant carries a GLM column. Adding one is step 5 of
`Reproduction/_archive/NEXT_STEPS.md`.

---

## Reproduction & the retracted "bug"
`Reproduction/` reproduces Figs 5–8 from `Support/allsubjectsTable.csv` two ways (a clean-room
MATLAB pipeline and a bridge into this repo's own code). It once reported a
polar-angle-ordering bug in `compute_derivativeDirections.m`. **That report is retracted** —
see `Reproduction/AUDIT.md`, a full experiment-code→figure audit of the stimulus conventions.

The original `AnalysisCode` pipeline is **correct**, and an independent recomputation from the
CSV reproduces the manuscript on all eight asymmetries, including `da` horizontal−vertical
= −0.446 (manuscript −0.45). The discrepancy was caused by two bugs inside `Reproduction/`:
swapped c-/cc-spirals in `cleanroom/config_repro.m` (flipping the four oblique wedges) and a
Benson-vs-conventional polar-angle frame mismatch in `bridge/` (flipping the four cardinals).

**Key convention** (established in `Reproduction/AUDIT.md`): the shared condition index 26–29 is each
stimulus's *local orientation at the upper vertical meridian* (0°/90°/45°/135°), and the wedge
dimension of `medianBOLDpa` is in **Benson** order — conventional `[90 45 0 315 270 225 180 135]`.
Do not "align" that array with `[0 45 90 …]`; that would introduce the bug into working code.

## Related docs

Listed newest first — later documents supersede earlier ones where they overlap.

- **`Reproduction/HARMONIC_MODEL.md` (2026-08-17).** A per-vertex harmonic model that replaces
  the eight polar-angle wedges with each vertex's own pRF polar angle, to separate within-ROI
  local-orientation geometry from genuine context effects. Vertices are weighted for equal
  polar-angle coverage, which is what the published ROI analysis does implicitly and what makes
  the four predictors orthogonal. Findings: geometry accounts for only
  **5.6% / 7.6%** of the Cartesian-vs-polar gap in the horizontal/vertical and cardinal/oblique
  asymmetries, so **the context claim survives**; but in raw (% signal change) units the
  **radial/tangential asymmetry shows no DETECTABLE difference between experiments** (dg 0.119,
  da 0.162). Do NOT upgrade that to "the polar-frame asymmetries are context-invariant": the
  interval admits an effect larger than the card-obl context effect, dropping sub-0395 makes it
  significant, and the Cartesian-vs-polar difference of differences is n.s. (p=0.26). Absence of
  evidence only. Also note the radial/tangential comparison reverses under z-scoring, which is why
  it is stated for the raw analysis only. All context tests are WITHIN SUBJECT (per-observer
  difference first); an LME with experiment x asymmetry interactions is anti-conservative here
  (DF=502 not 7, p smaller by 5-25x) and must not be quoted — see
  `Reproduction/cleanroom/diagnose_context_asymmetry.m`. pRF polar-angle error was measured (not
  assumed) from the two independent pRF fits at **σ = 3.9°**, against the ~39° that would be needed
  for measurement error to explain the result away. Code in `Reproduction/cleanroom/harmonic_*`
  and `run_harmonic_model.m`; `test_harmonic_model.m` asserts the convention against `lme_codes`.
  The summary-statistic route was then defended by MEASURING within-observer error rather than
  assuming it away (split-half over the 35 balanced run splits, plus a bootstrap over runs;
  `diagnose_within_observer_error.m`, fed by `server_extract/collect_runwise_betas.m`): it is
  23-39% of the across-observer variance, and disattenuating changes no conclusion. The binding
  limitation is between-observer variability at n = 8, not measurement noise — so a mixed model
  would not have rescued anything.
- **`Reproduction/supplement/SUPPLEMENT_harmonic_model.md` (2026-08-17).** The same work written
  for a paper supplement — motivation, model design, results, interpretation, with four committed
  figures. Written for readers, not agents; `HARMONIC_MODEL.md` remains the fuller internal account
  (z-scored sensitivity, weighting diagnostics, implementation notes).
- **`Reproduction/local_qc/REPORT.md` (2026-07-24) — start here.** GLM data-quality review of
  all 8 observers × both experiments from the unfiltered server extraction. The pink-noise
  reference, the clearing of both flagged observers, the fixed-`rng` design finding, and the
  recommendation to drop z-scoring. Draft manuscript caveat text alongside it in
  `manuscript_caveat_paragraph.md`; scripts and `glm_summary.csv` in the same folder.
- `Reproduction/_archive/GLM_QUALITY.md` (2026-07-23) — **archived.** The first GLM fit-quality
  audit, on the 4–8° pRF-filtered extraction, superseded by `local_qc/REPORT.md` (the unfiltered
  re-extraction its own §6a called for, with a richer §2.1 table). Its measurement-reliability
  section now lives at `local_qc/REPORT.md` §2.7.
- `Reproduction/_archive/NEXT_STEPS.md` — the task setups those two audits came from. The z-scoring and
  fit-quality tasks are now closed; step 5 (a GLM-`R2` column in `allsubjectsTable.csv`) is the
  live remainder.
- **`Reproduction/LME.md`** — issue 6. Why the Fig-7 mixed model returns the identical estimates to
  the subtraction route, what a trial-level or precision-weighted version would and would not add,
  and the recommendation to omit it from the manuscript.
- **`Reproduction/WHY_NOT_ZSCORE.md`** — the standing account of why analyses are raw, not z-scored.
  Read this rather than the archived working documents; it carries the decision, its one live
  dependency (the radial/tangential context result), and pointers to the diagnostics that still run.
- `Reproduction/_archive/ZSCORE_FIG7.md` (2026-07-22) — why z-scoring reverses the radTan/H−V rank order
  in Fig 7B: observer reweighting, and the algebra of the reversal (§1–§5), which remains the
  best account of the mechanism. **Its §7–§8 conclusions — the missing GLM check and the
  two-observer exclusion — are superseded.** Scripts: `cleanroom/diagnose_zscore_fig7.m`,
  `diagnose_response_signs.m`, `compare_subject_weighting.m`.
- `Reproduction/_archive/FINDINGS.md` — **fully retracted**; kept only as a record of the reproduction's
  own two bugs. Do not act on it; `Reproduction/AUDIT.md` is the correct account.
- `Reproduction/server_extract/` — the read-only server extraction (`collect_everything.m`) that
  produced the data behind `local_qc/REPORT.md`, plus `RUNME.md` for whoever has the volume
  mounted.
- `README.md` — experiment (stimulus-presentation) overview.
- `AnalysisCode/README.rtf` — detailed, near function-by-function description of the full
  analysis pipeline (includes stages upstream and downstream of Figs 5–8).
