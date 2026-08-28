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

**Acquisition and design.** 3T Siemens Prisma, multiband EPI, TR = 1 s, 2 mm isotropic; fMRIPrep and
FreeSurfer, so every analysis is vertex-wise on the cortical surface rather than voxel-wise. Each run
is 52 trials — 13 conditions × 4 repeats, 3 s stimulus and 2 s inter-trial interval — at 50% Michelson
contrast in an aperture of 12.2° radius. pRF parameters come from a separate retinotopy protocol fitted
with `vistasoft`. **Thirteen** observers were run on `dg` and 8 on `da`; which of them enters which
analysis is standing fact 7, and it is no longer the plain "8 in both" that the hardcoded `[1:8]`
index in the older scripts assumes. (Recorded 2026-08-27, merged from upstream; these protocol details
are not stated anywhere else in this repository.)

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

**Where the claim stands** (re-derived on the corrected 7 observers, 2026-08-27). **Three of the four
asymmetries now show detectable context dependence**: horizontal−vertical (−0.302, *p* = .021),
cardinal−oblique (−0.149, *p* = .032) and **radial−tangential (−0.107, *p* = .001, 7/7 observers)**.
Only **polar-cardinal versus polar-oblique** remains uninformative (+0.006, *p* = .81) — absence of
evidence, not evidence of absence. The radial−tangential result is new: on the previous 8-observer set
it was *p* = .65 and was reported as uninformative. What changed is the observer set, not the method —
sub-0395's `da` session used the wrong stimulus (fact 7), and that observer was the only one whose
`da` radial−tangential ran the other way.

**One conclusion went backwards in the same re-run.** The V1 − V3 attenuation of the
horizontal−vertical context effect is now *p* = .053, against *p* = .028 before; the estimate barely
moved (−0.215 → −0.211) but a degree of freedom was lost. It is **suggestive, not established** — do
not quote it as demonstrated. → [`Reproduction/RESULTS.md`](Reproduction/RESULTS.md)

---

## 2. Seven standing facts

Each changes how everything else should be read. Each is settled.

1. **The "blank" is full-field pink noise, not a mean-luminance baseline.** It is on screen for the
   whole run, so every reported response is grating − pink-noise and **there is no true baseline
   anywhere**. A low GLM R² means weak *differentiation among* stimuli, not an absent response. The
   orientation asymmetries are unaffected — the common blank term cancels in each.
   → [`local_qc/DATA_QUALITY.md`](Reproduction/local_qc/DATA_QUALITY.md) §1
2. **No observer is excluded for data quality.** The two once flagged as anomalous (sub-0037,
   sub-0201) show normal MT motion-selectivity and V4 grating preference in the same sessions. The
   exclusion once argued for is **withdrawn**. This is about *data quality only*; sub-0395 is
   separately excluded from `da` because it was run on the wrong stimulus, which is a fact about the
   session and not a judgement about its data (fact 7). (Reworded 2026-08-27: this used to read "all 8
   observers are retained", which is no longer true of `da`.)
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
7. **The observer set is no longer "the same 8 in both experiments".** Settled 2026-08-27.
   **sub-0395's `da` session used a pilot stimulus** whose annuli did **not** scale spatial period
   with eccentricity, so those data are not comparable with the other observers' and the observer is
   **excluded from `da`**. Separately, **13 observers were run on `dg`**, not 8. The agreed scheme is:
   **`dg` = all 13, `da` = the valid 7, and any contrast between the two experiments = those same 7.**
   Upstream implements this (see §4); **`Reproduction/` on this side does not yet** — every number in
   these documents is still the old 8-and-8, including sub-0395 in `da` (§5).

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

**New upstream material, merged 2026-08-27.** Rania's side has added a large body of code that the
analyses above do not depend on, but that anyone reading `AnalysisCode/` will now meet:

- **`02_ttave/`** — run time-series QC. Observed against GLMsingle-predicted BOLD for each run
  (`run_runTimeseries.m`), the cross-subject group average of the same (`run_groupAverageRunTimeseries.m`,
  which refits a canonical HRF to the group-mean trace because per-vertex HRF fits cannot be averaged),
  a within-run against across-run comparison, and full-model GLM R² on the surface
  (`computeGLMsingleR2.m`).
- **`01_calculate_observer_gain/computeObserverGainWeightsByROI.m`** — pRF gain per observer × cortical
  area rather than V1 only, on the same 4–8° and R² ≥ 0.1 vertex criteria.
- **`01_calculate_observer_precisionWeights/`** — per-observer measurement reliability, computed by the
  same across-run resampling as this side's precision weighting, which its header cites.
- **`04_plot_betaAsymmetries/fitAsymmetryRegression.m`** and the plots that read its cache — a weighted
  least-squares successor to `lme1_fit.m`, one joint fit per experiment × cortical area with a fixed
  sum-to-zero observer intercept. It is parallel to `lme1_fit.m`, not a replacement: neither reads the
  other. → [`REGRESSION.md`](AnalysisCode/04_plot_betaAsymmetries/REGRESSION.md)

Note that `REGRESSION.md` uses **"ROI" for cortical area and "location" for polar-angle wedge** — the
code's convention, which is the opposite of the manuscript's. These documents use "wedge" for the
polar-angle bin (§1).

**`da`/sub-wlsubj124 presented the same trial sequence twice — do not "fix" it.** Runs 1 and 2 of that
one observer share a byte-identical `expDes.trialMat`, and this is baked into `matrices_onset{1..2}`
in their `rawInfo.mat`. No other observer × experiment pair is affected. Upstream checked this against
the data (each run's BOLD correlates best with its own run's model) and concluded the sequence really
was presented twice and was recorded correctly — a quirk of the session, **not** a mislabelling bug, so
GLMsingle should not be re-run and `modelOutput.mat` should not be edited. It matters only where runs
are averaged as if independent: `findDuplicateDesignRuns.m` detects the case and
`run_groupAverageRunTimeseries.m` drops the later run of a duplicate pair. Nothing in the asymmetry
analyses averages runs this way, so none of the results in `Reproduction/` are affected.

**Figure-export loops need `drawnow` before `print()`.** Two figures made back-to-back in one loop can
otherwise be written with identical content, because `print()` grabs a stale render. Any new
figure-generating loop should call `drawnow;` immediately before each `print(...)`.

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

- **The abstract's polar-frame claim: half of it is now supported, half is not.** The draft reads
  "the radial asymmetry is 50% larger for polar gratings and the polar cardinal asymmetry is the
  same." Re-derived on the corrected 7 observers 2026-08-27: the **radial** half is **supported and
  understated** — `da` 0.198 against `dg` 0.091 is a factor of 2.2, with the context effect
  significant at *p* = .001 and all 7 observers agreeing. The **polar-cardinal** half should still not
  be quoted: "the same" claims equality, and +0.006 [−0.053, 0.065] *p* = .81 is an absence of
  evidence, not a demonstration that the two are equal. (This item previously said the radial claim
  was not quotable either; that was the 8-observer result, *p* = .65.)
  → [`RESULTS.md`](Reproduction/RESULTS.md) §4
- **Report *t* intervals, or say which method is used.** One V1 cell still disagrees between methods
  — `da` card−obl, bootstrap [−0.083, −0.002] excluding zero against *t* [−0.092, 0.011], *p* = .105.
  *t* is primary. Context effects are unaffected.
  → [`METHOD_DECISIONS.md`](Reproduction/METHOD_DECISIONS.md) §5
- **Geometric vs. arithmetic mean gain — closed.** `lme1_fit.m`, `plot1_experimentalCond.m` and
  `plot2_experimentalCond.m` now use `exp(mean(log(gainWeights)))`, matching the Methods and the
  clean-room. Worth ~1% on every effect size and no change to any *t* or *p*, as expected. (Was an open
  request to Rania; fixed upstream and merged here 2026-08-27.)
  → [`METHOD_DECISIONS.md`](Reproduction/METHOD_DECISIONS.md) §2
- **"Both Cartesian effects decline monotonically up the hierarchy" is now established for
  *neither*.** Re-run on 7 observers 2026-08-27: V1 − V3 horiz−vert is −0.211, *p* = **.053**, and
  card−obl −0.074, *p* = .150. The horiz−vert estimate barely moved (−0.215 → −0.211); it lost a
  degree of freedom. It also straddles the routes — *p* = .053 under `spec`, .044 under `roi` — which
  is a sign of an underpowered comparison rather than a route effect, and neither side of .05 should
  be quoted. The **V2 → V3 step** is solid for both (*p* = .016 and .009); it is the V1-to-V3
  difference that is unresolved. **The 13-observer `dg` set does not fix this** (corrected
  2026-08-27; an earlier version of this item said it would). The quantity the trend test tracks is
  the *context effect*, `dg − da`, so it is a paired difference and needs the matched 7 in both maps
  no matter how many observers `dg` has. Getting more power here needs more observers who did **both**
  experiments. A 13-observer `dg`-only attenuation across maps is answerable, but it is a different
  question from the one this table tests.
  → [`RESULTS.md`](Reproduction/RESULTS.md) §5
- **The observer scheme is implemented and the specification outputs have been regenerated on the
  corrected 7; `dg`-on-13 is still blocked.** 2026-08-27. All twelve figures and every CSV under
  `Reproduction/supplement/` are now the 7-observer result, and `RESULTS.md` has been rewritten around
  them — see the top of that document for what moved. What exists in the code:
  - `cfg.subjects` is split into `cfg.subjects_dg` (13), `cfg.subjects_da` (7) and
    `cfg.subjects_matched` (**derived**, never typed twice), with `cfg.dgSubjectMode`
    (`'matched'` default, `'all'` for the 13) mirroring the repo's `dg_subjectMode`.
    `cfg.subjects` survives as a legacy alias resolving to the **matched 7**, so the ~104
    existing readers get the conservative set rather than a silent change of meaning.
  - [`subjects_for.m`](Reproduction/cleanroom/subjects_for.m) resolves the list for an
    experiment; [`assert_same_observers.m`](Reproduction/cleanroom/assert_same_observers.m)
    refuses to pair two results that are not the same people in the same order.
  - `bin_and_aggregate` and `harmonic_vertex_data` resolve per experiment and **rebind
    `cfg.subjects` locally**, which keeps positional consumers such as
    `observer_gain_weights` aligned without editing every call site. `harmonic_vertex_data`
    now carries `D.subjects` out with the data, and `fit_harmonic_vertex` indexes by it.
  - The guard is wired into `harmonic_crossexp` (the site that pairs `si` across both
    experiments) and `spec_profiles` (one list, experiment on the inner loop).
  - Asking for the 13 against the current table now **errors** with what to do about it,
    rather than aggregating five observers to silent NaN rows.

  **`dg` on 13 now runs — the table was rebuilt 2026-08-27.**
  [`build_table_rows.m`](Reproduction/server_extract/build_table_rows.m) reproduces
  `createTables.m` for one observer from the already-collected local files, so no mount is
  needed, and the five `dg`-only observers were appended to
  `Support/allsubjectsTable.csv` (now 13 observers, 3.79M rows, 2.0 GB) with the `polexp_*`
  columns NaN. `cfg.dgSubjectMode = 'all'` then gives a 13-observer `dg` fit on both routes.

  **It is validated, not assumed.** [`validate_table_rows.m`](Reproduction/server_extract/validate_table_rows.m)
  rebuilds an observer already in the table and compares all 38 columns. Three observers —
  sub-0037, sub-0250, sub-wlsubj123 — reproduce the shipped rows exactly (worst difference
  5e-13, CSV text precision). After appending, the original 8 observers' V1 rows are
  **bit-identical** to the previous cache, so nothing already published moved.

  Two traps that a count-based check would have passed, both now covered by the validator:
  `collect_labels` stores label indices **already** 1-based while `createTables.m` adds the 1
  itself, so adding it again shifts every vertex by one — every label keeps its exact size and
  only boundary vertices move to the neighbouring area. And comparing NaN columns with `max()`
  silently reports a perfect match, because `max` ignores NaN.

  **`V2d` is deliberately not applied, and should never be added.** `left_V2d` is not an
  analysis ROI — it is a **diagnostic**, requested to check the polar-angle convention. V1, V2
  and V3 are bilateral and span dorsal and ventral, so together they cover the whole visual
  field and cannot tell one angle convention from another; **left** V2d is dorsal and
  one-hemisphere, so it should represent only the **lower right** visual field — 270–360° if
  angles are conventional, 90–180° if they are still Benson. It comes out 96–99% in 270–360°
  with a circular mean near 315°, which is the confirmation recorded in
  [`STIMULUS_CONVENTIONS.md`](Reproduction/STIMULUS_CONVENTIONS.md) §3 point 3. It answered its
  question. It is also probably **atlas-derived rather than hand-drawn**, unlike every other
  label used here, and `createTables.m` lists it last — so applying it would overwrite
  hand-drawn V2 with atlas boundaries (1,190 V2 vertices and 66 V1 vertices for sub-0037 alone).
  The shipped table's zero `left_V2d` rows are the **correct** state, not an omission being
  copied for consistency.

  **First result on 13.** All four `dg` asymmetries replicate with tighter intervals: horiz−vert
  −0.465 *p*<.001, card−obl −0.186 *p*=.003, rad−tang 0.095 *p*=.001, polc−polo 0.035 *p*=.007
  (wedge route, equal weighting; against −0.520/−0.193/0.080/0.025 on the matched 7). Six added
  observers, every asymmetry confirmed, every *p* smaller.

  **Done on 13 (2026-08-27).** `spec_profiles` now resolves the observer list **per experiment**
  and carries `S.dg.subjects` / `S.da.subjects` with the data;
  `diagnose_within_observer_error` takes `subjects` and `experiments` so σᵢ is measured on each
  experiment's own list. `run_spec_outputs` fits `dg` **twice** per variant — 13 observers for
  Figure 5, the matched 7 for everything paired — and the axis limits span both sets so Figures 5
  and 6 stay on one scale. `spec_tables` gained `context` (emit per-experiment asymmetries from a
  13/7 profile) and `fileSuffix` (the `*_dg13.csv` tables behind Figure 5); it still **refuses** a
  paired contrast whenever the two observer lists differ.

  | figure | observers |
  |---|--:|
  | Figure 5 (`dg`) | **13** |
  | Figure 6 (`da`) | 7 |
  | profile figure | 7 |
  | Figure S5 (hierarchy) | 7 |
  | every paired table | 7 |

  **All four `dg` asymmetries hold on 13 with smaller *p* than on 7**: horiz−vert −0.467, card−obl
  −0.189, rad−tang +0.106, polc−polo +0.066, all *p* ≤ .003, against −0.514/−0.193/+0.091/+0.056 on
  the matched 7. Six added observers, every asymmetry confirmed — the clearest replication evidence
  in the project. → [`FIGURE_VARIANTS.md`](Reproduction/supplement/FIGURE_VARIANTS.md)

  **`SPECIFICATION.md`, `FIGURE_VARIANTS.md` and the `local_qc/` documents have NOT been updated** —
  their numbers are still the 8-observer set and they contradict `RESULTS.md` until someone goes
  through them. `RESULTS.md` and this file are current.

  Two properties of the re-run worth carrying: dropping sub-0395 changes each *per-observer* value by
  a single **uniform factor of 1.0215** (the group gain is a property of the set, so that part is unit
  restoration and moves no *t*, *p* or CI); and **τ̂² pins to exactly 0 for both polar-frame
  asymmetries at n = 7**, with the precision weight ratio reaching 15:1. The τ̂² = 0 rows are a
  small-sample artifact, not a finding, and are why precision weighting gives *p* = .024 rather
  than .001 for the V1 radial−tangential context effect. Equal weighting is primary and is unaffected.

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
- [`AnalysisCode/04_plot_betaAsymmetries/REGRESSION.md`](AnalysisCode/04_plot_betaAsymmetries/REGRESSION.md)
  — upstream's weighted-least-squares successor to the Fig-7 LME: the model, the weighting, and its
  subject selection. Nothing in `Reproduction/` depends on it. Read the terminology note first (§4).
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
