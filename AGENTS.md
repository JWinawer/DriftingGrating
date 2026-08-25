# AGENTS.md — DriftingGrating

Entry point for this repository. It tells you what the analyses are, what is settled, what is
still open, and which document owns each topic — **without reading the manuscript**.

**Nothing here is published.** This is a **manuscript in preparation** — no journal version, no
preprint. Where the documents say "the manuscript route", "the manuscript values" or "the
manuscript asymmetries", they mean *the analysis and numbers in the current draft*, which is still
moving. Do not describe any of it as published, established, or final.

The draft is a Google Doc, owned by Rania Ezzo, reachable through the shortcut in `Manuscript/`
(git-ignored). Authors: Ezzo, Carrasco, Rokers, Winawer. Two titles are in circulation and it is
worth knowing which you are looking at — the document's own heading reads *"Local orientation
asymmetries in V1 depend on global stimulus properties"*, while the Drive filename reads *"Local
orientation asymmetries in fMRI BOLD magnitude depend on the global stimulus orientation"*.

The analyses of active interest are **Figures 5–8** plus the **extrastriate supplement**.

**The draft lags these analyses, and that is expected.** This repository is where the analysis work
happens and it has moved well ahead of the text; the two are expected to converge, not to match at
any given moment. So a disagreement between a number here and a number in the draft is normally the
draft not having caught up — **not** an error in either place, and not something to "fix" by editing
the analysis to match. Rania Ezzo, the first author, manages the manuscript; the convergence happens
there, on their schedule. What is useful from this side is keeping the list of divergences accurate
and specific (§7), so that catching up is a matter of transcription rather than re-derivation.

Two repository conventions:

- **`_archive/` folders hold retired material.** Check the superseding document named in the
  folder's `README.md` before acting on anything there.
- **Corrections replace the superseded text.** A withdrawn or re-measured claim carries a short
  dated note saying what changed; the old prose and numbers are deleted rather than kept below it
  (git history is the record). Warnings about what *not* to do are current guidance, not history —
  keep those. See [`CLAUDE.md`](CLAUDE.md).

---

## 1. Standing facts

Six facts change how everything else should be read. Each is settled; each has an owning document.

1. **The "blank" is full-field pink noise, not a mean-luminance baseline.** It is on screen for the
   whole run, so every reported response is grating − pink-noise and **there is no true baseline
   anywhere**. A low GLM R² means weak *differentiation among* stimuli, not an absent response.
   The orientation asymmetries are unaffected — the common blank term cancels in each.
   → [`local_qc/REPORT.md`](Reproduction/local_qc/REPORT.md) §1
2. **No bad observers, no processing errors; all 8 observers are retained.** The two once flagged
   as anomalous (sub-0037, sub-0201) show normal MT motion-selectivity and V4 grating preference in
   the same sessions. The exclusion once argued for is **withdrawn**.
   → [`local_qc/REPORT.md`](Reproduction/local_qc/REPORT.md) §2
3. **Analyses are raw % signal change — not z-scored** (decided 2026-07-24). `beta_std` is not a
   gain when there is no baseline, and dividing the blank-independent asymmetries by it puts the
   blank back in. → [`WHY_NOT_ZSCORE.md`](Reproduction/WHY_NOT_ZSCORE.md)
4. **The across-vertex aggregate is the MEAN, with per-observer pRF-gain rescaling** (decided
   2026-08-19). That route reproduces all eight manuscript asymmetries to ±0.003; the median route
   misses `dg` horiz−vert by 0.07. → [`local_qc/GLM_SUMMARY_SECTION.md`](Reproduction/local_qc/GLM_SUMMARY_SECTION.md)
5. **Polar-angle wedges are stored in Benson order** — conventional `[90 45 0 315 270 225 180 135]`.
   Shared condition index 26–29 is each stimulus's *local orientation at the upper vertical
   meridian* (0°/90°/45°/135°). Do **not** "align" the wedge dimension of `medianBOLDpa` with
   `[0 45 90 …]`; that introduces a bug into working code.
   → [`AUDIT.md`](Reproduction/AUDIT.md)
6. **All context tests are within subject** — form the per-observer difference first, then test
   across the 8 observers. The Fig-7 LME is not an alternative: it returns the identical estimate
   and only changes the standard error, in the anti-conservative direction (DF = 502, not 7).
   → [`LME.md`](Reproduction/LME.md)

---

## 2. Study in one paragraph

fMRI (n = 8 observers, V1 defined by pRF retinotopy) measuring whether V1 BOLD **orientation
asymmetries are context-dependent on the global stimulus**. Two experiments with *locally matched*
stimuli (1 cpd at 6° eccentricity):

- **Cartesian gratings** — horizontal, vertical, two obliques. Project code name **`dg`**.
- **Polar gratings** — pinwheel, annulus, CW spiral, CCW spiral. Project code name **`da`**.

Only **stationary** gratings are analyzed here (the GLM also estimates 8 motion directions,
reported in a future study). V1 vertices are binned into **8 polar-angle wedges** (centered
0,45,…,315°, each ±22.5°), restricted to **4–8° eccentricity** and **pRF R² > 0.1**. The 4–8° band
is where spatial frequency is most closely matched between the two experiments — it is *not* the
stimulated extent, which was roughly 0.5–12°.

**Four orientation asymmetries**, in two reference frames:

| Asymmetry | Reference frame | Harmonic | LME predictor coding |
|---|---|---|---|
| horizontal vs. vertical | Cartesian (absolute) | 1st | −1 / 0 / +1 |
| cardinal vs. oblique | Cartesian (absolute) | 2nd | −1 / +1 |
| radial vs. tangential | polar (gaze-relative) | 1st | −1 / 0 / +1 |
| polar-cardinal vs. polar-oblique | polar (gaze-relative) | 2nd | −1 / +1 |

**Core claim:** an asymmetry is **strong when its reference frame matches the global stimulus, weak
or absent when it does not**. Cartesian asymmetries (BOLD *lower* for horizontal and for cardinal —
the "inverted" direction) appear for Cartesian gratings and weaken for polar; polar asymmetries
(radial > tangential, polar-cardinal > polar-oblique) appear for polar gratings.

The claim survives every test run against it, with one qualification: the **Cartesian-frame**
context effects are established, and the **polar-frame** ones are *uninformative rather than
absent* — absence of evidence, not evidence of absence
([`HARMONIC_MODEL.md`](Reproduction/HARMONIC_MODEL.md) Result 3).

---

## 3. The figures

In every figure, panels **A–D** are the four asymmetries in the order above:
**A** = horizontal vs. vertical, **B** = cardinal vs. oblique, **C** = radial vs. tangential,
**D** = polar-cardinal vs. polar-oblique.

**Figure 5 — independent asymmetries, Cartesian gratings (`dg`).** Each asymmetry quantified
independently. Per polar-angle wedge, take the mean % signal change (orientation minus blank)
across vertices, gain-rescale per observer, then average across observers. *Top:* polar plot
(radius = V1 % change, angle = visual-field polar angle), one colored trace per orientation class.
*Bottom:* per-observer pairwise plot with a bootstrap 95% CI on the difference (1000 resamples of
observers). Result: A and B significant (H < V, cardinal < oblique); C and D weak for Cartesian
gratings.

**Figure 6 — independent asymmetries, polar gratings (`da`).** Identical analysis and layout.
Result: Cartesian asymmetries weakened; radial > tangential present; polar-cardinal vs.
polar-oblique a positive trend.

**Figure 7 — joint LME (bar plots). Being removed** (JW, 2026-08-19). It is exactly redundant with
Figures 5/6 once the same across-vertex aggregate is used: `fit_lme_fig7` returns
`[−0.547, −0.221, 0.104, 0.039]` for `dg`, identical to the independent asymmetry means. The
apparent difference between routes was an aggregate mismatch (`lme1_fit.m` read `meanBOLDpa`,
`plot_NeuralAsymmetries.m` passed `medianBOLDpa`), not anything the model added.
→ [`LME.md`](Reproduction/LME.md)

**Figure 8 — per-location context effects.** One small polar plot per visual-field location: local
orientation on the angular axis, response magnitude as distance from origin; **black = observer
data, red = LME model estimate**. Panel A = `dg`, panel B = `da`. Highlights that the largest
context effect is on the **horizontal meridian** (suppressed horizontal + enhanced radial oppose
each other). *Fig 8's red curve is the same LME that Fig 7 drops — see Open items.*

**Extrastriate supplement.** The same four asymmetries across eight visual maps, under one settled
specification. → [`EXTRASTRIATE.md`](Reproduction/EXTRASTRIATE.md)

---

## 4. The settled analysis specification

Arrived at 2026-08-24 and applied unchanged to every map, V1 included: per observer × map gain
applied at the observer boundary; four-term harmonic model fitted per vertex with **continuous**
θv; equal-coverage weighting at 45°; fit per observer then averaged; equal weighting primary with
precision weighting reported alongside. No per-map customisation anywhere.

Full statement, the evidence for each of the four decisions, and the coverage criterion that
decides which maps can be reported at polar-angle resolution:
[`EXTRASTRIATE.md`](Reproduction/EXTRASTRIATE.md) §1–§6. Adopting it moved V1 by at most 0.040 and
changed no significance.

---

## 5. Code map

Analysis is MATLAB (R2022b). Paths in the scripts point at a mounted data volume
(`/Volumes/Vision/UsersShare/Rania/Project_dg/...`) and `~/Documents/GitHub/DriftingGrating`; the
fMRI data itself is **not** in the repo. `projectName` selects the experiment: **`'dg'` =
Cartesian, `'da'` = polar.**

**Data pipeline (bottom → top):**

1. **`AnalysisCode/01_process_singlesubjectGLM/main_singlesub.m`** — per subject, GLMsingle
   produces `betamaps` (vertices × 13 conditions: 8 motion + 4 orientations + 1 blank). Saves
   `betas_nonzscored.mat`. The `normalize` flag (~line 158) would z-score per vertex across the 13
   conditions and save `betas_zscored.mat`; **it is off, and z-scoring is not used** (standing fact
   3). Blank-subtracted orientation contrasts land in `results.mat` → `results.contrasts`.
2. **`AnalysisCode/01_calculate_observer_gain/`** — per-vertex pRF gain from the retinotopy model,
   the divisor behind the per-observer rescaling. See its `README.md`.
3. **`AnalysisCode/03_process_groupBetas/meanWithinLabel.m`** — bins vertices into V1 × 8
   polar-angle wedges (4–8°, R² > 0.1) and saves group matrices **`meanBOLDpa`** / **`medianBOLDpa`**
   (contrasts × polarAngles × ROIs × subjects) and **`meanBOLD`** / **`medianBOLD`**. Downstream
   callers use the **mean** (standing fact 4).
4. **`AnalysisCode/04_plot_betaAsymmetries/`** — consumes those arrays to make the figures.

**Figure → script:**

| Figure | Script |
|---|---|
| **5** (Cartesian) | `04_plot_betaAsymmetries/plot_NeuralAsymmetries.m`, `projectName='dg'` → `plot1_experimentalCond` (polar plots) + `plot2_experimentalCond` (pairwise) |
| **6** (Polar) | same driver, `projectName='da'` |
| **7** (LME weights) | `04_plot_betaAsymmetries/lme1_fit.m` — figure being removed |
| **8** (per-location) | `04_plot_betaAsymmetries/lme2_ploteachDirLoc.m` |

**Known cosmetic wart, left alone deliberately:** `plot1_/plot2_experimentalCond.m` and the
`compute_derivative*` functions still name their first parameter `medianBOLDpa` while receiving
mean data. Renaming was skipped because upstream also edits those files.

---

## 6. Questions worked through

All resolved. The document listed is the standing account; working documents they replace are in
[`Reproduction/_archive/`](Reproduction/_archive/README.md).

| # | issue | question | resolution | document |
|---|---|---|---|---|
| 1 | **Z-scoring** | Should per-vertex betas be z-scored before analysis? | **No.** The blank is pink noise, so `beta_std` is not a gain. | [`WHY_NOT_ZSCORE.md`](Reproduction/WHY_NOT_ZSCORE.md) |
| 2 | **Data quality** | Are there outlier observers or datapoints to exclude? | **No.** No processing errors, no bad runs; all 8 observers retained. | [`local_qc/REPORT.md`](Reproduction/local_qc/REPORT.md) |
| 3 | **Angle conventions** | Are the absolute and location-dependent angle codings correct? | **Yes.** Code is consistent and correct; an earlier bug report is retracted. | [`AUDIT.md`](Reproduction/AUDIT.md) |
| 4 | **Harmonic model** | Does ROI binning manufacture a context effect that is really local stimulus geometry relative to each vertex's pRF? | **No.** Geometry accounts for 6–8%; pRF angle error (σ = 3.9°) is an order of magnitude too small. | [`HARMONIC_MODEL.md`](Reproduction/HARMONIC_MODEL.md), [supplement](Reproduction/supplement/SUPPLEMENT_harmonic_model.md) |
| 5 | **Context effects** | What does a within-subject assessment support? | **Cartesian-frame asymmetries yes; polar-frame uninformative** — absence of evidence, not evidence of absence. | [`HARMONIC_MODEL.md`](Reproduction/HARMONIC_MODEL.md) Result 3 |
| 6 | **LME** | Does the Fig-7 mixed model add anything? | **No.** Identical estimates to the subtraction route; anti-conservative SE. Figure being removed. | [`LME.md`](Reproduction/LME.md) |
| 7 | **Which aggregate / which gain** | Mean or median across vertices; gain rescaling or not? | **Mean, gain-rescaled** — reproduces the manuscript values to ±0.003. | [`local_qc/GLM_SUMMARY_SECTION.md`](Reproduction/local_qc/GLM_SUMMARY_SECTION.md) |
| 8 | **Extrastriate** | Does the V1 route work in maps with empty polar-angle ROIs? | **No, and it is replaced.** One specification now applies unchanged to every map (§4). | [`EXTRASTRIATE.md`](Reproduction/EXTRASTRIATE.md) |

---

## 7. Open items

Everything still live, in one place. The topic documents point here rather than keeping their own
lists.

**Where the draft and the analysis have diverged** — Rania owns these, in the manuscript
and in the upstream `AnalysisCode` files

- **The abstract's polar-frame claim has moved.** The draft reads "the radial asymmetry is 50%
  larger for polar gratings and the polar cardinal asymmetry is the same." The direction still
  holds — `da` rad−tang exceeds `dg` rad−tang on both routes (0.150 vs 0.104 on the manuscript
  route; 0.155 vs 0.119 under the settled spec) — but on the current analysis **the difference is
  not detectable**: *p* = 0.56 paired across observers, *p* = .65 under the spec. So the ratio
  should not be quoted as a quantity. The supporting statement to replace it with is in
  [`HARMONIC_MODEL.md`](Reproduction/HARMONIC_MODEL.md) Result 3: polar-frame context effects are
  **uninformative rather than absent**, and the interval admits an effect larger than the
  cardinal/oblique context effect that *is* significant. → also
  [`EXTRASTRIATE.md`](Reproduction/EXTRASTRIATE.md) §7
- **Bootstrap vs. *t* intervals disagree on two polar asymmetries** (`da` card−obl, `da` rad−tang).
  The supplement currently calls polar rad−tang "clearly non-zero"; *t* on 7 df gives *p* = 0.081.
  Context effects are unaffected. → [`LME.md`](Reproduction/LME.md) §7
- **Geometric vs. arithmetic mean gain.** The Methods say geometric; `lme1_fit.m`,
  `plot1_experimentalCond.m` and `plot2_experimentalCond.m` use `mean(gainWeights)` and should be
  `exp(mean(log(gainWeights)))`. Worth ~1% on every effect size; changes no *t* or *p*. The
  clean-room already defaults to geometric.

**Figure decisions on this side**

- **Figure 8's model overlay.** Fig 8 plots the LME estimate in red against the data in black. With
  Fig 7 removed, whether Fig 8 keeps that overlay, or is redrawn without it, has not been decided.
- **Figure 4C.** Either drop it, or present it as an illustration for one observer with no
  quantitative claim attached — the R²-by-eccentricity profile is flat, for reasons that are
  properties of the measure and the design. → [`local_qc/GLM_SUMMARY_SECTION.md`](Reproduction/local_qc/GLM_SUMMARY_SECTION.md)

**Analysis still to do**

- **Test the hierarchy trend.** The claim that the context effect declines monotonically V1 → V2 →
  V3 rests on six individually significant cells falling in the same order. The trend itself — a
  within-observer V1 − V3 difference — has not been computed. → [`EXTRASTRIATE.md`](Reproduction/EXTRASTRIATE.md) §8
- **Finish the run-mismatch calibration.** The sign result is settled and measured on 56 run pairs;
  what a further pull adds is absolute ROI-level magnitude. Check **sub-0201/da** (ρ = 0.36, only
  66% of pairs negative) before stating the claim over all run pairs.
  → [`local_qc/GLM_SUMMARY_SECTION.md`](Reproduction/local_qc/GLM_SUMMARY_SECTION.md) Result 5

**Loose ends, low priority**

- **A GLM-`R2` column in `allsubjectsTable.csv`**, so the vertex filter can screen on it alongside
  `pRF_r2`. Currently no GLMsingle metric enters the pipeline at any stage. Split-half reliability
  now answers the underlying question directly, which lowers the urgency.
- **Confirm the fixed stimulus `rng` seed was intended** — every observer's run *K* is
  byte-identical to every other's. → [`local_qc/REPORT.md`](Reproduction/local_qc/REPORT.md) §3
- **sub-0201 per-condition confirmation** — needs the raw single-trial betas (`modelmd`), ~1 GB,
  not extracted.
- **Cross-run design order** was not byte-matched against each session's data-run order. Run counts
  and per-run R² already rule out a mispaired run.

---

## 8. Document index

**Reproduction and analysis** — [`Reproduction/README.md`](Reproduction/README.md) explains the two
routes, the data table, and how to run them.

| document | what it owns |
|---|---|
| [`Reproduction/EXTRASTRIATE.md`](Reproduction/EXTRASTRIATE.md) | **The settled specification** and the results across eight visual maps. The most recent word on method; where V1 and the supplement now agree. |
| [`Reproduction/HARMONIC_MODEL.md`](Reproduction/HARMONIC_MODEL.md) | The per-vertex harmonic model — geometry vs. context in V1, the pRF-angle-error control, cross-experiment prediction. The internal account. |
| [`Reproduction/supplement/SUPPLEMENT_harmonic_model.md`](Reproduction/supplement/SUPPLEMENT_harmonic_model.md) | The same work written for readers, with four committed figures. |
| [`Reproduction/AUDIT.md`](Reproduction/AUDIT.md) | Stimulus-design audit, experiment code → GLM → CSV → figures. Owns every angle convention, and the retraction of the reported polar-angle bug. |
| [`Reproduction/WHY_NOT_ZSCORE.md`](Reproduction/WHY_NOT_ZSCORE.md) | Why analyses are raw, and the one conclusion that turns on it. |
| [`Reproduction/LME.md`](Reproduction/LME.md) | Why the Fig-7 mixed model adds nothing, precision weighting done properly, and the bootstrap-vs-*t* disagreement. |
| [`Reproduction/local_qc/REPORT.md`](Reproduction/local_qc/REPORT.md) | GLM data quality, all 8 observers × both experiments. The pink-noise finding, the clearing of both flagged observers, measurement reliability. |
| [`Reproduction/local_qc/GLM_SUMMARY_SECTION.md`](Reproduction/local_qc/GLM_SUMMARY_SECTION.md) | Split-half reliability, which aggregation route the manuscript uses, the Fig-4 controls, and draft manuscript text. |
| [`Reproduction/local_qc/manuscript_caveat_paragraph.md`](Reproduction/local_qc/manuscript_caveat_paragraph.md) | Draft Methods/Discussion paragraph for the pink-noise caveat. |
| [`Reproduction/server_extract/README.md`](Reproduction/server_extract/README.md) | The read-only server extractions, one per script, for whoever has `/Volumes/Vision` mounted. |
| [`Reproduction/_archive/README.md`](Reproduction/_archive/README.md) | Retired working documents and what superseded each. **Nothing there is current.** |

**Elsewhere in the repo**

- [`README.md`](README.md) — experiment (stimulus-presentation) overview and the `Models/` notebook.
- `AnalysisCode/README.rtf` — near function-by-function description of the full analysis pipeline,
  including stages upstream and downstream of Figs 5–8.
- [`AnalysisCode/01_calculate_observer_gain/README.md`](AnalysisCode/01_calculate_observer_gain/README.md)
  — what pRF gain means here and how it is computed.
