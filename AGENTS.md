# AGENTS.md — Target analyses (Figures 5–8)

This file lets a future agent understand the analyses behind **Figures 5–8** of the
manuscript *"Local orientation asymmetries in V1 depend on global stimulus properties"*
(Ezzo, Carrasco, Rokers, Winawer) **without reading the manuscript**. The draft itself is
in `Manuscript draft/` (git-ignored, not tracked).

These four figures are the analyses of active interest. The open question driving the work
is **whether each vertex's GLM beta weights should be z-scored before analysis** — two
versions of every one of Figs 5–8 currently exist (z-scored vs. not), and we need to decide
which to adopt. See [The z-scoring question](#the-z-scoring-question) below.

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

---

## Reproduction & a known bug
`Reproduction/` reproduces Figs 5–8 from `Support/allsubjectsTable.csv` two ways (a clean-room
MATLAB pipeline and a bridge into this repo's own code). It reproduces every direct asymmetry
exactly, and surfaced a **polar-angle-ordering bug in `compute_derivativeDirections.m`**: it
swaps the horizontal/vertical (and radial/tangential) labels at the four cardinal meridians,
corrupting the two *first-harmonic derived* asymmetries — `da` horizontal-vertical (reported
−0.45, correct ≈ 0) and `dg` radial-tangential (current code ≈ 0, correct 0.23). Direct
asymmetries and 2nd-harmonic derived ones are unaffected. See `Reproduction/FINDINGS.md`.

## Related docs
- `Reproduction/FINDINGS.md` — reproduction results and the derived-direction bug (read this
  before trusting any *derived* asymmetry: `da` H-V/card-obl, `dg` rad-tang/polar-card).
- `README.md` — experiment (stimulus-presentation) overview.
- `AnalysisCode/README.rtf` — detailed, near function-by-function description of the full
  analysis pipeline (includes stages upstream and downstream of Figs 5–8).
