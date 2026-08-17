# AGENTS.md — Target analyses (Figures 5–8)

This file lets a future agent understand the analyses behind **Figures 5–8** of the
manuscript *"Local orientation asymmetries in V1 depend on global stimulus properties"*
(Ezzo, Carrasco, Rokers, Winawer) **without reading the manuscript**. The draft itself is
in `Manuscript draft/` (git-ignored, not tracked).

These four figures are the analyses of active interest. The open question driving the work
is **whether each vertex's GLM beta weights should be z-scored before analysis** — two
versions of every one of Figs 5–8 currently exist (z-scored vs. not), and we need to decide
which to adopt. See [The z-scoring question](#the-z-scoring-question) below.

> ## Read this first
> **[`Reproduction/local_qc/REPORT.md`](Reproduction/local_qc/REPORT.md) (2026-07-24) is the
> current authority** on data quality and on the z-scoring decision. It supersedes parts of
> this file, of `Reproduction/ZSCORE_FIG7.md`, `GLM_QUALITY.md`, and `NEXT_STEPS.md`, all of
> which predate it. Two findings there change how everything else should be read:
>
> 1. **The "blank" condition is full-field pink noise, not a mean-luminance baseline.** Pink
>    noise is on the screen for the entire run, so every reported response is
>    grating − pink-noise and **there is no true baseline anywhere**. A low GLM R² therefore
>    means weak *differentiation among* the stimuli, **not** an absent response. The
>    orientation asymmetries are unaffected — the common blank term cancels in each.
> 2. **There are no bad observers and no processing errors.** The two observers previously
>    flagged as anomalous (sub-0037, sub-0201) show normal MT motion-selectivity and V4
>    grating preference in the same sessions. The exclusion argued for in `ZSCORE_FIG7.md` §8
>    is **withdrawn**.
>
> Its recommendation: **use the non-z-scored analyses and remove z-scoring**, including the
> "beta weights were standardized" Methods language and the σ-unit statistics.

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

## The z-scoring question

Two versions of Figs 5–8 exist, differing only in whether each vertex's GLM betas were
**z-scored across the 13 conditions before analysis** (applied at
`01_process_singlesubjectGLM/main_singlesub.m` via the `normalize` flag).

- The manuscript **Methods currently commit to z-scoring** ("beta weights for each vertex were
  standardized"; the LME uses "the z-scored, blank-subtracted GLM coefficients"), and in-text
  statistics are reported in σ units.
- **Figs 5, 6, and 8** look qualitatively the same between the two versions (asymmetry shapes
  and significance unchanged; only axis units/scale differ).
- **Fig 7 differs substantively**: which asymmetries reach significance in the joint LME
  changes between the z-scored and non-z-scored fits. This is the crux of the decision.

When working on these analyses, keep the two variants consistent across all four figures and
across the Methods prose, and treat Fig 7 as the case where the choice actually matters.

> **Settled 2026-07-24 — use the raw (non-z-scored) analyses.**
> [`Reproduction/local_qc/REPORT.md`](Reproduction/local_qc/REPORT.md) §4 resolves this
> against z-scoring, on grounds that were not available when the account below was written:
> the blank is **pink noise, not a baseline**, so `beta_std` is the spread among
> contrast-pattern responses rather than a gain. It conflates BOLD gain with orientation-tuning
> strength and motion sensitivity, and dividing the (blank-independent) raw asymmetries by it
> *reintroduces* a dependence on the blank. In the raw analysis the largest polar asymmetry is
> horizontal−vertical — the same ordering the analysis below arrived at, but without needing to
> exclude any observer. The account below remains a correct description of *how* z-scoring
> changes Fig 7; its §7–§8 conclusions do not survive.

**What the Fig 7 difference turned out to be** (2026-07-22, full account in
`Reproduction/ZSCORE_FIG7.md`): z-scoring changes which polar-grating asymmetry is largest
(radial−tangential z-scored, horizontal−vertical raw). The cause is **observer reweighting**, not
vertex-level gain control — a single scalar per subject reproduces the whole effect. At the
subject level `beta_std` is overall response amplitude (r = +0.94), *not* data quality, so
z-scoring weights observers by 1/responsiveness; the two it up-weights most have anomalous polar
data. Neither ordering survives a subject bootstrap (P = 0.69 z-scored, 0.28 raw, 0.43 under an
uncontaminated divisor). Two consequences for anyone picking this up:

- **~~Normalising observers is justified, but it has a precondition.~~ Withdrawn.** This bullet
  argued for excluding sub-0037 and sub-0201 as having no measurable gain in the polar
  experiment. `local_qc/REPORT.md` §2.5–2.6 rules that out with positive evidence — both
  observers' MT is motion-selective and their V4 prefers gratings to pink noise in the same
  sessions — and §1 removes the premise, since a "blank-referenced gain" measured against pink
  noise is not a gain at all. **All 8 observers are retained.**
- **Do not frame Fig 7 around which asymmetry is largest** (bootstraps 0.10–0.31 at n=6;
  unanimous in direction, never decisive). The cross-experiment comparison is what the data
  support. *(This bullet stands.)*
- **The GLM fits have been quality-checked — twice, and they are sound.** All 8 observers ×
  both experiments were extracted from the server (2026-07-23, then again unfiltered on
  2026-07-24) and audited: no coding or processing error, uniform model parameters, no bad run,
  no dropout, correct surface co-registration. See `Reproduction/local_qc/REPORT.md` §2 and
  `Reproduction/GLM_QUALITY.md`. What remains genuinely open is that **no GLMsingle metric
  enters the pipeline** — the only quality filter is still on the *pRF* fit (`pRF_r2 > 0.1`),
  and neither `allsubjectsTable.csv` variant carries a GLM column. Adding one is step 5 of
  `Reproduction/NEXT_STEPS.md`.

---

## Reproduction & the retracted "bug"
`Reproduction/` reproduces Figs 5–8 from `Support/allsubjectsTable.csv` two ways (a clean-room
MATLAB pipeline and a bridge into this repo's own code). It once reported a
polar-angle-ordering bug in `compute_derivativeDirections.m`. **That report is retracted** —
see `AUDIT.md`, a full experiment-code→figure audit of the stimulus conventions.

The original `AnalysisCode` pipeline is **correct**, and an independent recomputation from the
CSV reproduces the manuscript on all eight asymmetries, including `da` horizontal−vertical
= −0.446 (manuscript −0.45). The discrepancy was caused by two bugs inside `Reproduction/`:
swapped c-/cc-spirals in `cleanroom/config_repro.m` (flipping the four oblique wedges) and a
Benson-vs-conventional polar-angle frame mismatch in `bridge/` (flipping the four cardinals).

**Key convention** (established in `AUDIT.md`): the shared condition index 26–29 is each
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
  **radial/tangential asymmetry does not differ between the two experiments** (dg 0.119, da 0.162), so it is the
  Cartesian-frame asymmetries that are context-dependent, not the polar-frame ones — a sharper,
  more asymmetric statement than the current framing. That last point reverses under z-scoring,
  which is why it is stated for the raw analysis only. pRF polar-angle error was measured (not
  assumed) from the two independent pRF fits at **σ = 3.9°**, against the ~39° that would be needed
  for measurement error to explain the result away. Code in `Reproduction/cleanroom/harmonic_*`
  and `run_harmonic_model.m`; `test_harmonic_model.m` asserts the convention against `lme_codes`.
- **`Reproduction/supplement/SUPPLEMENT_harmonic_model.md` (2026-08-17).** The same work written
  for a paper supplement — motivation, model design, results, interpretation, with four committed
  figures. Written for readers, not agents; `HARMONIC_MODEL.md` remains the fuller internal account
  (z-scored sensitivity, weighting diagnostics, implementation notes).
- **`Reproduction/local_qc/REPORT.md` (2026-07-24) — start here.** GLM data-quality review of
  all 8 observers × both experiments from the unfiltered server extraction. The pink-noise
  reference, the clearing of both flagged observers, the fixed-`rng` design finding, and the
  recommendation to drop z-scoring. Draft manuscript caveat text alongside it in
  `manuscript_caveat_paragraph.md`; scripts and `glm_summary.csv` in the same folder.
- `Reproduction/GLM_QUALITY.md` (2026-07-23) — the first GLM fit-quality audit, on the 4–8°
  pRF-filtered extraction. Its R² tables stand; its §5 exclusion support and its §6/§6a open
  questions are answered by `local_qc/REPORT.md`.
- `Reproduction/NEXT_STEPS.md` — the task setups those two audits came from. The z-scoring and
  fit-quality tasks are now closed; step 5 (a GLM-`R2` column in `allsubjectsTable.csv`) is the
  live remainder.
- `Reproduction/ZSCORE_FIG7.md` (2026-07-22) — why z-scoring reverses the radTan/H−V rank order
  in Fig 7B: observer reweighting, and the algebra of the reversal (§1–§5), which remains the
  best account of the mechanism. **Its §7–§8 conclusions — the missing GLM check and the
  two-observer exclusion — are superseded.** Scripts: `cleanroom/diagnose_zscore_fig7.m`,
  `diagnose_response_signs.m`, `compare_subject_weighting.m`.
- `Reproduction/FINDINGS.md` — **fully retracted**; kept only as a record of the reproduction's
  own two bugs. Do not act on it; `AUDIT.md` is the correct account.
- `Reproduction/server_extract/` — the read-only server extraction (`collect_everything.m`) that
  produced the data behind `local_qc/REPORT.md`, plus `RUNME.md` for whoever has the volume
  mounted.
- `README.md` — experiment (stimulus-presentation) overview.
- `AnalysisCode/README.rtf` — detailed, near function-by-function description of the full
  analysis pipeline (includes stages upstream and downstream of Figs 5–8).
