# AGENTS.md — DriftingGrating

Entry point for this repository: what the study is, what the analyses found, what is still open, and
which document owns each topic — **without reading the manuscript**.

**Nothing here is published.** This is a **manuscript in preparation** — no journal version, no
preprint. The draft is a Google Doc owned by Rania Ezzo, reachable through the shortcut in
`Manuscript/` (git-ignored). Authors: Ezzo, Carrasco, Rokers, Winawer. Where a document says "the
draft", it means that Google Doc, which is still moving. Do not describe any of this as published,
established, or final. (Two titles are in circulation, and it is worth knowing which one you are
looking at: the document's own heading reads *"Local orientation asymmetries in V1 depend on global
stimulus properties"*, while the Drive filename reads *"Local orientation asymmetries in fMRI BOLD
magnitude depend on the global stimulus orientation"*.)

**The draft lags these analyses, and that is expected.** This repository is where the analysis work
happens and it has moved ahead of the text; the two are expected to converge, not to match at any
given moment. A disagreement between a number here and a number in the draft is normally the draft
not having caught up — **not** an error in either place, and not something to "fix" by editing the
analysis to match. Rania manages the manuscript; convergence happens there, on their schedule. What
is useful from this side is keeping §5 accurate and specific, so catching up is transcription rather
than re-derivation.

---

## 1. The study in one page

fMRI, **n = 8 observers**, V1 defined by pRF retinotopy, measuring whether V1 BOLD **orientation
asymmetries depend on the global structure of the stimulus**. Two experiments with *locally matched*
stimuli (1 cpd at 6° eccentricity):

- **Cartesian gratings** — horizontal, vertical, two obliques. Project code name **`dg`**.
- **Polar gratings** — pinwheel, annulus, CW spiral, CCW spiral. Project code name **`da`**.

Only **stationary** gratings are analysed here; the GLM also estimates 8 motion directions, reported
in a future study. Vertices are restricted to **4–8° eccentricity** and **pRF R² > 0.1**. The 4–8°
band is where spatial frequency is most closely matched between the two experiments — it is *not*
the stimulated extent, which was roughly 0.5–12°.

**Polar angle enters continuously, not in bins.** The settled specification fits each vertex at its
own pRF polar angle (§3). The eight 45° wedges (centred 0, 45, …, 315°, each ±22.5°) are still
everywhere in the project — they are the earlier route, they are the `roi` variant kept for
comparison, they are the unit the coverage criterion counts, and they are what the original
`AnalysisCode` computes — so "wedge" and "ROI" in these documents mean that binning, and it is
never the primary analysis.

**Four orientation asymmetries**, in two reference frames:

| asymmetry | reference frame | harmonic |
|---|---|---|
| horizontal vs. vertical | Cartesian (absolute) | 1st |
| cardinal vs. oblique | Cartesian (absolute) | 2nd |
| radial vs. tangential | polar (gaze-relative) | 1st |
| polar-cardinal vs. polar-oblique | polar (gaze-relative) | 2nd |

**Core claim.** An asymmetry is **strong when its reference frame matches the global stimulus, weak
or absent when it does not**. Cartesian asymmetries — BOLD *lower* for horizontal and for cardinal,
the "inverted" direction — appear for Cartesian gratings and weaken for polar; polar asymmetries
(radial > tangential, polar-cardinal > polar-oblique) appear for polar gratings.

**Where the claim stands.** The **Cartesian-frame** context effects are established and survive every
test run against them. The **polar-frame** ones are *uninformative rather than absent* — absence of
evidence, not evidence of absence. → [`Reproduction/RESULTS.md`](Reproduction/RESULTS.md)

---

## 2. Six standing facts

Each changes how everything else should be read. Each is settled.

1. **The "blank" is full-field pink noise, not a mean-luminance baseline.** It is on screen for the
   whole run, so every reported response is grating − pink-noise and **there is no true baseline
   anywhere**. A low GLM R² means weak *differentiation among* stimuli, not an absent response. The
   orientation asymmetries are unaffected — the common blank term cancels in each.
   → [`local_qc/DATA_QUALITY.md`](Reproduction/local_qc/DATA_QUALITY.md) §1
2. **No bad observers, no processing errors; all 8 observers are retained.** The two once flagged as
   anomalous (sub-0037, sub-0201) show normal MT motion-selectivity and V4 grating preference in the
   same sessions. The exclusion once argued for is **withdrawn**.
   → [`local_qc/DATA_QUALITY.md`](Reproduction/local_qc/DATA_QUALITY.md) §2
3. **Analyses are raw % signal change — not z-scored.** `beta_std` is not a gain when there is no
   baseline, and dividing the blank-independent asymmetries by it puts the blank back in.
   → [`METHOD_DECISIONS.md`](Reproduction/METHOD_DECISIONS.md) §1
4. **The across-vertex aggregate is the MEAN, with per-observer pRF-gain rescaling.** That route
   reproduces all eight of the draft's asymmetries to ±0.003; the median route misses `dg` horiz−vert
   by 0.07. → [`METHOD_DECISIONS.md`](Reproduction/METHOD_DECISIONS.md) §2
5. **Polar-angle wedges are stored in Benson order** — conventional `[90 45 0 315 270 225 180 135]`.
   Shared condition index 26–29 is each stimulus's *local orientation at the upper vertical
   meridian*. Do **not** "align" the wedge dimension of `meanBOLDpa`/`medianBOLDpa` with
   `[0 45 90 …]`; that introduces a bug into working code.
   → [`STIMULUS_CONVENTIONS.md`](Reproduction/STIMULUS_CONVENTIONS.md) §3
6. **All context tests are within subject** — form the per-observer difference first, then test
   across the 8 observers. The Figure-7 LME is not an alternative: it returns the identical estimate
   and only changes the standard error, in the anti-conservative direction (DF = 502, not 7).
   → [`METHOD_DECISIONS.md`](Reproduction/METHOD_DECISIONS.md) §3

---

## 3. The analysis, in short

**One specification, applied unchanged to V1 and to all eight supplement maps** (settled 2026-08-24):
per observer × map gain applied at the observer boundary; a four-term harmonic model fitted per
vertex with **continuous** pRF polar angle; equal-coverage weighting at 45°; fit per observer then
averaged; equal weighting primary with precision weighting reported alongside; *t* intervals primary;
polar-angle-resolved results reported only where measured coverage passes a fixed criterion. No
per-map customisation anywhere.

Everything — Figures 5 and 6, the polar-angle profile, the extrastriate supplement figure, and eight
CSV tables in three variants — is regenerated by one command:

```bash
matlab -batch "cd('Reproduction/cleanroom'); run_spec_outputs"
```

→ [`SPECIFICATION.md`](Reproduction/SPECIFICATION.md) for the method and the evidence behind each
choice; [`RESULTS.md`](Reproduction/RESULTS.md) for every number.

**The figures.** Panels A–D are always the four asymmetries in the order above.
**Figure 5** = Cartesian gratings, **Figure 6** = polar gratings — polar plots on top (fitted model
plus observed wedge means), each observer's difference below.
**Figure 7** (joint LME bar plots) **is being removed**: it is exactly redundant with Figures 5/6
once the same aggregate is used.
**Figure 8** (per-location context effects) is untouched; its red model overlay is the Fig-7 LME, so
its fate is open (§5).
**Extrastriate supplement** = the same four asymmetries across V1, V2, V3 and V3a.

---

## 4. Code map

Analysis is MATLAB (R2022b). Paths in the original scripts point at a mounted data volume
(`/Volumes/Vision/UsersShare/Rania/Project_dg/...`); the fMRI data itself is **not** in the repo.
`projectName` selects the experiment: **`'dg'` = Cartesian, `'da'` = polar.**

**Original pipeline (`AnalysisCode/`), bottom → top:**

1. **`01_process_singlesubjectGLM/main_singlesub.m`** — per subject, GLMsingle produces `betamaps`
   (vertices × 13 conditions: 8 motion + 4 orientations + 1 blank), saved as `betas_nonzscored.mat`.
   The `normalize` flag would z-score per vertex across the 13 conditions; **it is off, and z-scoring
   is not used** (standing fact 3).
2. **`01_calculate_observer_gain/`** — per-vertex pRF gain from the retinotopy model, the divisor
   behind the per-observer rescaling. See its `README.md`.
3. **`03_process_groupBetas/meanWithinLabel.m`** — bins vertices into V1 × 8 polar-angle wedges and
   saves **`meanBOLDpa`** / `medianBOLDpa` (contrasts × polarAngles × ROIs × subjects). Downstream
   callers use the **mean** (standing fact 4).
4. **`04_plot_betaAsymmetries/`** — consumes those arrays to make the figures:
   `plot_NeuralAsymmetries.m` (Figs 5 and 6), `lme1_fit.m` (Fig 7, being removed),
   `lme2_ploteachDirLoc.m` (Fig 8).

**Reproduction (`Reproduction/`)** — an independent recomputation from a single tidy table, and where
every follow-up analysis now lives. `cleanroom/` is the primary route and holds the settled
specification; `server_extract/` holds read-only extractions for whoever has the data volume mounted;
`local_qc/` holds the data-quality review. (A fourth path, `bridge/`, ran the original
`AnalysisCode` functions for cross-checking; it agreed exactly, nothing depended on it afterwards,
and it was deleted 2026-08-25.) → [`Reproduction/README.md`](Reproduction/README.md)

---

## 5. Open items

**Everything still live, in one place.** The topic documents do not keep parallel lists; each item
below points at the document that has the detail. If something is open and is not here, it is a gap
in this list, not a second list somewhere else.

### Where the draft and the analysis have diverged — Rania owns these

- **The abstract's polar-frame claim has moved.** The draft reads "the radial asymmetry is 50% larger
  for polar gratings and the polar cardinal asymmetry is the same." The direction holds (`da`
  rad−tang 0.155 vs `dg` 0.119) but **the difference is not detectable** (*p* = .65), so the ratio
  should not be quoted as a quantity. The statement to replace it with is about *evidence*: the
  Cartesian-frame asymmetries show robust context dependence; for the polar-frame ones these data are
  uninformative. → [`RESULTS.md`](Reproduction/RESULTS.md) §4
- **Report *t* intervals, or say which method is used.** One V1 cell still disagrees between methods
  — `da` card−obl, bootstrap [−0.083, −0.002] excluding zero against *t* [−0.092, 0.011], *p* = .105.
  *t* is primary. Context effects are unaffected.
  → [`METHOD_DECISIONS.md`](Reproduction/METHOD_DECISIONS.md) §5
- **Geometric vs. arithmetic mean gain.** The Methods say geometric; `lme1_fit.m`,
  `plot1_experimentalCond.m` and `plot2_experimentalCond.m` use `mean(gainWeights)` and should be
  `exp(mean(log(gainWeights)))`. Worth ~1% on every effect size; changes no *t* or *p*. The clean-room
  already defaults to geometric.
  → [`METHOD_DECISIONS.md`](Reproduction/METHOD_DECISIONS.md) §2
- **"Both Cartesian effects decline monotonically up the hierarchy" is established for one of the
  two.** The within-observer V1 − V3 difference of the context effect is significant for horiz−vert
  (−0.215, *p* = .028) and not for card−obl (−0.066, *p* = .153). The cells fall in order in each
  map, but the difference *between* maps is resolved for only one asymmetry.
  → [`RESULTS.md`](Reproduction/RESULTS.md) §5

### Figure decisions on this side

- **Figure 8's model overlay.** Fig 8 plots the LME estimate in red against the data in black. With
  Figure 7 removed, whether Fig 8 keeps that overlay, or is redrawn without it, has not been decided.
- **Figure 4C.** Either drop it, or present it as an illustration for one observer with no
  quantitative claim attached — the R²-by-eccentricity profile is flat, for reasons that are
  properties of the measure and the design.
  → [`local_qc/RELIABILITY.md`](Reproduction/local_qc/RELIABILITY.md) §4
- **A Figure 4D is proposed and not drawn.** A split-half scatter — the 32 orientation-differential
  profile values from one half of the runs against the other, one colour per observer — would put
  the reliability claim in the same "show the data" register as 4A–C. Nothing depends on it.
  → [`local_qc/RELIABILITY.md`](Reproduction/local_qc/RELIABILITY.md) §5

### Analysis still to do

- **Calibrate the §7 coverage thresholds against the thing they proxy.** The criterion is written as
  thresholds on empty cells, median vertices and weight ratio; the quantity behind it — the width of
  the band on the estimate under a map's actual coverage — is now measured directly, but the
  thresholds have not been re-derived from it. Nothing reported depends on this: V3 and MT are
  separated by a factor of seven in band width, nowhere near the boundary.
  → [`MISSING_DATA.md`](Reproduction/MISSING_DATA.md) §7

- **Finish the run-mismatch calibration.** The sign result is settled and measured on 56 run pairs;
  what a further pull adds is absolute ROI-level magnitude. Check **sub-0201/da** (ρ = 0.36, only 66%
  of pairs negative) before stating the claim over all run pairs.
  → [`local_qc/RELIABILITY.md`](Reproduction/local_qc/RELIABILITY.md) §4

- **Does `equalcell` weighting help in the maps that *are* reported?** In the missing-data
  simulation, declining to compensate an observer for their missing ROIs beat compensating them on
  the Cartesian terms. That is a lead from one simulation, untested on V1, V2, V3 and V3a as
  analysed, and it is not a recommendation until it is.
  → [`MISSING_DATA.md`](Reproduction/MISSING_DATA.md) §5

### Loose ends, low priority

- **A GLM-`R2` column in `allsubjectsTable.csv`**, so the vertex filter could screen on it alongside
  `pRF_r2`. No GLMsingle metric enters the pipeline at any stage. Split-half reliability now answers
  the underlying question directly, which lowers the urgency.
- **Confirm the fixed stimulus `rng` seed was intended** — every observer's run *K* is byte-identical
  to every other's. → [`local_qc/DATA_QUALITY.md`](Reproduction/local_qc/DATA_QUALITY.md) §3
- **sub-0201 per-condition confirmation** — needs the raw single-trial betas (`modelmd`), ~1 GB, not
  extracted.
- **Cross-run design order** was not byte-matched against each session's data-run order. Run counts
  and per-run R² already rule out a mispaired run.
- **The normalization modelling in `Models/` is a separate, unfinished strand** with its own
  open questions and TO-DO list, kept in the first markdown cell of
  `Models/model_DriftingGratings.ipynb`. Nothing in `Reproduction/` depends on it, so it is not
  itemised here. → [`README.md`](README.md)
- **The Benson-order wedge centres are re-derived by hand in three downstream files.** Storing them
  alongside the data, or converting once at the source, would remove the trap that caused the
  retracted bug report. The current code's *output* is correct, so this is optional.
  → [`STIMULUS_CONVENTIONS.md`](Reproduction/STIMULUS_CONVENTIONS.md) §7
- **Two questions only answerable at the server.** Are
  `/Volumes/Vision/UsersShare/Rania/Project_dg/` and `/Volumes/server/Projects/Project_dg/` the same
  data or two copies — different `AnalysisCode` scripts point at each, and if they have drifted then
  some figures were made from different data than others. And what machine actually hosts them;
  nothing in the code records anything but a local mount point.
  → [`server_extract/README.md`](Reproduction/server_extract/README.md)

---

## 6. Document index

**Start here, then read what you need.** Every document below states the current position; none of
them stacks superseded versions.

| document | what it owns |
|---|---|
| [`Reproduction/README.md`](Reproduction/README.md) | How the reproduction code is organised, the data table, and how to run each route. |
| [`Reproduction/SPECIFICATION.md`](Reproduction/SPECIFICATION.md) | **The settled analysis specification** — the model, the four decisions and the evidence for each, the coverage criterion, how to regenerate everything, and the traps in the code. |
| [`Reproduction/RESULTS.md`](Reproduction/RESULTS.md) | **Every current number** — V1 asymmetries and context effects, results across the visual hierarchy, the tested hierarchy trend, cautions, and what the figures show. |
| [`Reproduction/METHOD_DECISIONS.md`](Reproduction/METHOD_DECISIONS.md) | Five closed choices and why: no z-scoring, mean not median, no mixed model, precision weighting as a check not a default, *t* intervals not bootstrap. |
| [`Reproduction/MISSING_DATA.md`](Reproduction/MISSING_DATA.md) | **What empty cells do**, simulated in V1 — the bias from holes and why it is systematic, why a pooled group fit does not help, and why sparsity rather than holes is what keeps MT out of the supplement figure. |
| [`Reproduction/STIMULUS_CONVENTIONS.md`](Reproduction/STIMULUS_CONVENTIONS.md) | What each stimulus was, what each condition index means, and which polar-angle frame each part of the pipeline uses. Owns every angle convention. |
| [`Reproduction/local_qc/DATA_QUALITY.md`](Reproduction/local_qc/DATA_QUALITY.md) | GLM data quality, all 8 observers × both experiments — the pink-noise finding, the clearing of both flagged observers, the fixed-seed designs. |
| [`Reproduction/local_qc/RELIABILITY.md`](Reproduction/local_qc/RELIABILITY.md) | Split-half reliability of the analysed measurements, each effect against its own measurement error, the Figure 4 controls, and draft manuscript text. |
| [`Reproduction/supplement/FIGURE_VARIANTS.md`](Reproduction/supplement/FIGURE_VARIANTS.md) | **All twelve figures `run_spec_outputs` produces**, laid out figure by figure so the three variants can be compared by eye, each with the numbers behind it and what to look for. |
| [`Reproduction/supplement/SUPPLEMENT_harmonic_model.md`](Reproduction/supplement/SUPPLEMENT_harmonic_model.md) | The per-vertex harmonic model written for readers, with four committed figures — geometry vs. context, the cross-experiment prediction, and the pRF-angle-error control. **On the model's own route, not the settled specification** — see the note at its head. |
| [`Reproduction/server_extract/README.md`](Reproduction/server_extract/README.md) | The read-only server extractions, one per script, for whoever has `/Volumes/Vision` mounted. |
| [`Reproduction/_archive/README.md`](Reproduction/_archive/README.md) | Two retired working documents, kept only because current code cites them. **Nothing there is current.** |

**Elsewhere in the repo**

- [`README.md`](README.md) — experiment (stimulus-presentation) overview and the `Models/` notebook.
- `AnalysisCode/README.rtf` — near function-by-function description of the full analysis pipeline,
  including stages upstream and downstream of Figures 5–8.
- [`AnalysisCode/01_calculate_observer_gain/README.md`](AnalysisCode/01_calculate_observer_gain/README.md)
  — what pRF gain means here and how it is computed.
- [`Simulations/README.md`](Simulations/README.md) — two optic-flow illustration scripts. Nothing in
  the analyses depends on them.

---

## 7. Two repository conventions

- **`_archive/` folders hold retired material.** Check the superseding document named in the folder's
  `README.md` before acting on anything there.
- **Corrections replace the superseded text.** A withdrawn or re-measured claim carries a short dated
  note saying what changed; the old prose and numbers are deleted rather than kept below it (git
  history is the record). Warnings about what *not* to do are current guidance, not history — keep
  those. See [`CLAUDE.md`](CLAUDE.md).
