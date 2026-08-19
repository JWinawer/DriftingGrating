# GLM data-quality review — DriftingGrating (Ezzo et al.)

*Quality-control pass on the GLMsingle output for the manuscript "Local orientation asymmetries
in V1 depend on global stimulus properties." 2026-07-24.*

**Data reviewed.** Locally extracted GLMsingle output for all 8 observers × 2 experiments
(Cartesian `dg`, polar `da`) — whole-surface `R2`, `R2run`, `FRACvalue`, `noisepool`, `HRFindex`,
`pcnum`, `meanvol`, cross-validation curves, the ~31 blank-referenced contrast maps, the
retinotopy maps, every FreeSurfer label, and the per-run design/timing files. Copied read-only from
`/Volumes/Vision` to `~/dg_collect/` (1.2 GB); the raw single-trial betas (`modelmd`) were not
copied. Scripts, tables, and figures referenced below are in `Reproduction/local_qc/`.

---

## Summary

1. **No coding or processing error** in the GLM analyses. Model parameters are consistent and
   sensible across all 16 sessions; runs are internally consistent; surfaces are correctly
   co-registered; there is no signal dropout (§2).
2. **No genuine non-responders and no bad data** — with *positive* evidence, not merely absence of
   error. The two observers that had looked anomalous (sub-0037, sub-0201) show normal, expected
   responses in higher areas: MT is motion-selective and V4 prefers gratings to pink noise in both,
   even where V1 does not differentiate (§2.5–2.6).
3. **The reference ("blank") condition is a full-field pink-noise stimulus, not a mean-luminance
   baseline (§1).** This changes how response magnitude and variance-explained should be
   interpreted, but leaves the orientation-asymmetry results — the paper's claims — intact.
4. **Recommendation: use the non-z-scored analyses; remove the z-scored ones** (§4), including the
   "beta weights were standardized" language and σ-unit statistics in the Methods. The pink-noise
   reference removes the justification for z-scoring.

A separate, non-analysis observation: the stimulus used a fixed random seed, so every observer and
both experiments share identical run designs (§3).

---

## 1. The reference ("blank") condition is full-field pink noise

This is the central finding; it reframes the interpretation of everything else.

**What "blank" is.** In both experiments the "blank" condition is not a mean-luminance screen but a
**full-field 1/f (pink) noise** image — the same texture that fills the surround around the grating
aperture. This is explicit in the stimulus code:

- `ExperimentCode/Config/constConfig.m` — `pinknoise = oneoverf(1.1, <full screen>)` builds
  `const.pinknoiseTex`.
- `ExperimentCode/Stim/my_stim.m` — each grating trial draws the pink noise, then the grating
  within a circular aperture, then inner/outer masks, so the pink noise shows through everywhere
  outside the aperture.
- `ExperimentCode/Stim/my_blank.m` — draws only the full-field pink noise.
- `ExperimentCode/Stim/my_padding.m` and the shared inter-trial code in `my_stim.m` — every 2 s
  inter-trial interval and the 10 s run-flanking padding are also full-field pink noise.

So **pink noise is on the screen for the entire run; there is no mean-luminance baseline anywhere**
(`example_design_matrix.png`: the 4 blank trials per run, ~8% of trials, are short and
interspersed). The manuscript Methods already describe this — *"blank (pink noise background)"* and
*"blank stimuli contained only the pink noise background that was present throughout the run."*

**Consequences.**

1. Every reported response is **grating − pink-noise**; within the aperture it is essentially a
   grating-vs-pink-noise contrast, never grating-vs-nothing.
2. GLMsingle R² is variance **across conditions**. If pink noise, stationary gratings, and moving
   gratings evoke similar BOLD, there is little variance for the model to fit, so **R² and beta
   weights are near zero even when the region responds strongly to all of them.** Low estimates
   therefore indicate *weak differentiation among the stimuli*, not an absence of response.
3. **No true baseline exists, so absolute response magnitude cannot be recovered.** Comparisons of
   overall amplitude (across areas or observers) are differential responses relative to pink noise.

**Clean demonstration (sub-0037, V1, raw beta relative to the GLM baseline).** Holding observer,
area, and eccentricity band (4–8°) fixed, the only difference is the experiment:

| sub-0037, V1 | stationary | moving | blank (pink noise) | spread | median R² |
|---|--:|--:|--:|--:|--:|
| Cartesian (dg) | −0.59 | 0.36 | −0.74 | 1.10 | 16.6 |
| polar (da) | −0.19 | −0.15 | −0.20 | 0.05 | 2.9 |

In dg the three conditions differ (V1 distinguishes the Cartesian gratings, especially moving, from
pink noise) → variance to fit → R² 16.6. In da they are essentially identical (moving ≈ stationary ≈
pink noise) → no variance → R² floors at 2.9. This is variance-explained reflecting condition
*differences*, not response magnitude. (At the high end, sub-0395 has large condition differences —
blank −1.0, moving +1.4 — and correspondingly high R² ≈ 26.)

**What is unaffected.** The orientation **asymmetries** (horizontal vs vertical, cardinal vs
oblique, radial vs tangential, polar-cardinal vs polar-oblique) are contrasts *among* the grating
orientations, and the common blank term cancels in each. The pink-noise reference therefore affects
**overall magnitude and variance-explained**, but **not the asymmetry results.**

---

## 2. GLM data quality: no errors, no bad data

### 2.1 Fit quality and response by area

Median GLM R² and mean raw % response (stationary − blank), for V1/V2/V3 restricted to the 4–8°
stimulated eccentricity band (no variance-explained threshold on the pRF or the GLM) and the whole
surface unrestricted. Source: extracted GLM files (R²) and `allsubjectsTable.csv` (raw % betas);
full machine-readable values incl. motion and un-subtracted responses are in `glm_summary.csv`.

| observer | V1 R²dg | V1 R²da | V1 βdg | V1 βda | V2 R²dg | V2 R²da | V2 βdg | V2 βda | V3 R²dg | V3 R²da | V3 βdg | V3 βda | wh R²dg | wh R²da | wh βdg | wh βda |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| sub-0037 | 16.6 | 2.9 | 0.151 | 0.011 | 3.8 | 3.3 | 0.029 | −0.012 | 3.5 | 3.3 | 0.271 | 0.081 | 3.4 | 3.1 | 0.037 | 0.006 |
| sub-0201 | 3.5 | 4.5 | −0.224 | −0.349 | 3.1 | 3.9 | −0.178 | −0.062 | 3.0 | 3.7 | −0.004 | 0.290 | 3.0 | 3.9 | 0.003 | 0.007 |
| sub-0255 | 3.6 | 4.2 | 0.041 | 0.280 | 3.5 | 3.5 | −0.146 | 0.195 | 3.3 | 3.4 | 0.117 | 0.290 | 3.3 | 3.3 | 0.025 | 0.045 |
| sub-wlsubj123 | 4.8 | 4.4 | 0.543 | 0.424 | 3.4 | 3.7 | 0.401 | 0.216 | 3.4 | 3.7 | 0.523 | 0.354 | 3.4 | 3.5 | 0.063 | 0.016 |
| sub-wlsubj124 | 4.5 | 3.4 | 1.137 | 0.563 | 2.5 | 2.4 | 0.511 | 0.258 | 2.5 | 2.3 | 0.537 | 0.345 | 2.7 | 2.9 | 0.058 | 0.015 |
| sub-0395 | 26.4 | 19.8 | 1.558 | 1.304 | 8.0 | 6.2 | 0.839 | 0.990 | 7.1 | 7.1 | 0.887 | 1.132 | 3.0 | 3.9 | 0.083 | 0.095 |
| sub-0426 | 4.6 | 9.9 | −0.096 | 0.575 | 4.5 | 4.6 | −0.085 | 0.103 | 3.9 | 4.2 | 0.159 | 0.342 | 3.2 | 4.3 | 0.011 | 0.042 |
| sub-0250 | 6.3 | 3.6 | 0.353 | 0.378 | 8.1 | 3.7 | 0.421 | 0.265 | 4.9 | 3.6 | 0.504 | 0.455 | 4.0 | 3.2 | 0.087 | 0.050 |

*R² = median % variance explained; β = mean % signal change, stationary − blank. V1 responses are
concentrated in the stimulated 4–8° band: over the whole V1 label sub-0037's dg R² is only 5.4,
but within 4–8° it is 16.6 — the eccentricity restriction, not any variance threshold, is what
isolates the responsive patch.*

### 2.2 Model parameters — uniform (`group_models.m`)

Number of GLMdenoise noise regressors 2–9 (median 5, all in the normal range); fractional-ridge
shrinkage effectively identical across sessions (median fraction 0.05, as expected when most of the
surface is noise pool); HRF selection uniform. No session was fit with anomalous settings.

### 2.3 Runs — consistent (`group_runs.m`)

Run counts vary (most 8; sub-0255/dg has 9, sub-0395/da has 6), but each GLM used its own count (a
mismatch would have failed at fit time), and every session's per-run R² is internally consistent
(worst-run / session-median 0.82–0.98). No mispaired or dead run.

### 2.4 EPI intensity and surface alignment — clean (`group_extra.m`)

V1 mean EPI intensity is uniform across sessions (no dropout), and the GLM surface length matches
the retinotopy surface for every observer (retinotopy, labels, and GLM are correctly co-registered).

### 2.5 The two previously-flagged observers, explained

- **sub-0037** is the one observer with a large Cartesian-vs-polar difference in V1: in dg its V1
  differentiates the Cartesian gratings from pink noise (R² 16.6), in da it does not (R² 2.9). No
  error accounts for it — no bad run, no dropout, no misalignment, and the design is identical
  across observers (§3) so a design mix-up would have no effect. Given the pink-noise reference
  (§1), this is **weak differentiation of the polar gratings from pink noise, not an absence of
  response**. Higher areas confirm the session is sound: in the *same* da session, this observer's
  MT is motion-selective and V4 prefers gratings to pink noise (§2.6). So it is a genuine property
  of V1, not an error or bad data.
- **sub-0201** is the only observer whose mean V1 stationary−blank is negative in both experiments
  (−0.22 dg, −0.35 da): the windowed grating drove V1 no more than full-field pink noise. Under §1
  this is expected for an observer in whom pink noise and grating evoke comparable V1 responses; it
  is not bad data, and its GLM fit quality (R², model parameters) is normal — and its MT is
  strongly motion-selective (da moving−stationary +1.02, R² 19.2; §2.6), so the session clearly
  contains good, well-differentiated responses. Confirming the earlier per-condition "blank-beta"
  observation fully would require the raw single-trial betas, which were not extracted (Open items).

### 2.6 Higher areas (V4, MT) — positive evidence of good data

Because V1's small receptive fields respond similarly to gratings and to pink noise, a weak V1
grating−blank contrast is uninformative about data quality (§1). Higher areas provide a stronger
test. V4 and MT are reported over the **whole ROI** (their pRFs extend well beyond 4–8°, and the
4–8° band leaves too few vertices — as few as 12 for MT; whole-ROI counts are ~900–2000 for V4,
~370–730 for MT). All values are raw % signal change (`group_addv4mt.m`, `glm_summary.csv`).

**MT is motion-selective (moving − stationary > 0) in every observer and both experiments** — direct
evidence the motion system responded:

| observer | dg moving−stationary | dg R² | da moving−stationary | da R² |
|---|--:|--:|--:|--:|
| sub-0037 | +0.28 | 4.3 | +0.25 | 3.5 |
| sub-0201 | +0.32 | 3.1 | +1.02 | 19.2 |
| sub-0255 | +0.57 | 5.8 | +0.93 | 14.5 |
| sub-wlsubj123 | +0.73 | 13.9 | +0.74 | 12.6 |
| sub-wlsubj124 | +0.27 | 2.6 | +0.64 | 2.7 |
| sub-0395 | +0.78 | 10.6 | +1.02 | 17.1 |
| sub-0426 | +1.01 | 15.6 | +0.95 | 12.3 |
| sub-0250 | +0.75 | 9.9 | +1.10 | 20.8 |

**V4 prefers gratings to pink noise (stationary − blank > 0) in all 16 sessions** (range +0.06 to
+0.82; V4 shows little motion preference, as expected for a form area). For sub-0037, V4
stationary−blank is +0.38 (dg) / +0.31 (da) even though its V1 barely differentiates in da (+0.01).

The picture is physiologically coherent — MT prefers motion, V4 prefers coherent gratings over
noise, V1 responds similarly to both — and it establishes that for **both** previously-flagged
observers the visual system responded normally. sub-0037's da/V1 result is therefore a genuine V1
property (no differentiation of polar gratings from pink noise), not weak or bad data.

### 2.7 Measurement reliability of the analysed quantities, from runs

R² and response magnitude say how well the GLM fits. They do not say how reliable the **quantities
actually analysed** are. That is measured directly by resampling the measurement:
`../cleanroom/diagnose_within_observer_error.m` uses per-run condition betas
(`../server_extract/collect_runwise_betas.m`) and takes both a split-half over all 35 balanced 4-vs-4
run splits and a bootstrap over runs. The design is balanced — 8 runs, 4 trials per condition per
run — so both estimators are balanced by construction.

Within-observer SE of each asymmetry (% BOLD, bootstrap over runs):

| asymmetry | dg | da |
|---|---|---|
| horiz−vert | 0.121 | 0.042 |
| card−obl | 0.076 | 0.024 |
| rad−tang | 0.042 | 0.121 |
| polc−polo | 0.024 | 0.067 |

Note the structure: the SE is ~3× larger for the asymmetry **matched** to each experiment (dg
horiz−vert, da rad−tang) than for the derived one. A matched contrast uses the same two stimulus
conditions in every wedge, so its noise does not average across wedges; a derived contrast rotates
which stimuli it draws on and averages more effectively. Worth knowing when designing a similar
comparison.

Against the across-observer spread in the context difference:

| asymmetry | mean diff | SD across observers | within-observer SE | within/total |
|---|---|---|---|---|
| horiz−vert | −0.268 | 0.265 | 0.128 | 23% |
| card−obl | −0.175 | 0.126 | 0.079 | 39% |
| rad−tang | −0.034 | 0.245 | 0.128 | 28% |
| polc−polo | −0.004 | 0.137 | 0.071 | 27% |

**Two things follow, and they are the reason this belongs in a data-quality section.**

First, the measurements are good: measurement error is a minority of the observed spread for every
asymmetry, so most of what separates observers is genuine individual variation, not noise.
Disattenuating — removing the measurement variance and re-testing, the ceiling on what any mixed
model could recover — changes no conclusion (horiz−vert p 0.024 → 0.014, card−obl 0.006 → 0.0015,
rad−tang 0.70 → 0.66). **The binding limitation is between-observer variability at n = 8, not
measurement noise.**

Second, **precision weighting is not needed.** Weighting observers by 1/(τ² + σᵢ²) — where σᵢ² is
observer *i*'s measurement variance and τ² the between-observer variance of the true effects —
rather than equally shifts no group estimate materially and changes no *p*-value's interpretation.
The weights spread by only ~3× despite a 14× spread in reliability, because τ² is common to every
observer and dominates. Details and the full table are in [`../LME.md`](../LME.md) §5.

> **An earlier estimate of this quantity was withdrawn.** It resampled *vertices*, which is invalid:
> it holds the GLM betas fixed and only reshuffles which vertices enter the wedge median, so it
> characterises which patch of V1 was sampled rather than the reliability of the measurement — and
> vertex responses are not independent. It understated the error by 4–5×. Reliability has to come
> from resampling the measurement, i.e. across runs.

---

## 3. The stimulus used a fixed random seed (all designs identical)

Fingerprinting every per-run design (`design_uniqueness.m`): the 110 Cartesian design files reduce
to 12 distinct designs and the 68 polar files to 8, clustered exactly by run position. Every
observer's run *K* is byte-identical to every other observer's run *K* (verified with `isequaln`;
`rng` seed 0 with identical generator state), and Cartesian run *K* equals polar run *K*. Runs
differ within a session (run 1 ≠ run 2), but the whole sequence is deterministic from the fixed
seed, so it is the same for everyone.

This is a property of the stimulus code (the fixed `rng(0)`), not an analysis bug. Its one
consequence is that the stimulus order was not randomized per observer, so any stimulus-order
confound is shared across observers rather than averaging out. Worth confirming with the lead author
whether the fixed seed was intended. (A useful corollary: because all designs are identical, any
accidental use of one observer's or one experiment's design for another would have no effect —
which is part of why the earlier "wrong design matrix" concern is ruled out.)

---

## 4. Recommendation: use the non-z-scored analyses; remove z-scoring

The pink-noise reference (§1) settles the outstanding z-scoring question (`../_archive/NEXT_STEPS.md`,
`../_archive/ZSCORE_FIG7.md`) **against** z-scoring:

- z-scoring divides each vertex by `beta_std` (the SD across its 13 conditions). This is a valid
  gain normalizer only if the responses are anchored to a true baseline. With no baseline,
  `beta_std` is the spread among contrast-pattern responses — it conflates overall BOLD gain with
  orientation-tuning strength and motion sensitivity. Two observers with equal gain but different
  tuning have different `beta_std`; z-scoring forces them to match, erasing a genuine difference as
  if it were nuisance gain.
- The raw orientation asymmetries are blank-independent, but dividing them by `beta_std`
  reintroduces a dependence on the pink-noise blank and the motion conditions. z-scoring can only
  add distortion here.
- No valid alternative normalizer is available (the retinotopy model stores no gain; no divisor
  from the 13 conditions qualifies). The appropriate default is therefore not to normalize.

Practical consequences of removing z-scoring: the Methods drop the "beta weights were standardized"
statement and the σ-unit in-text statistics; Figures 5–8 use the raw (% signal change) variants,
which already exist; and in Figure 7 the largest polar asymmetry is horizontal−vertical (the raw
result). The framing rests on the cross-experiment context effect. The absolute-magnitude caveat
(§1) applies to both variants equally, so it is not a mark against the raw analysis.

A draft manuscript paragraph capturing the pink-noise-baseline caveat is in
`manuscript_caveat_paragraph.md`.

---

## Open items

- **sub-0201, full per-condition confirmation.** The negative mean V1 response is consistent with
  the earlier "blank beta exceeds the stimulus betas" observation, but a per-condition check (blank
  vs each of the 13) needs the raw single-trial betas (`modelmd`), which were not extracted. A
  targeted re-pull of that one field for sub-0201 (~1 GB) would settle it.
- **Fixed random seed** — confirm with the lead author whether it was intended (§3).
- **Cross-run design order.** Run counts and per-run R² are consistent (rules out a *mispaired*
  run), but each session's data-run order was not byte-matched against its design-file order; low
  priority given §2.3.

## Files (in `Reproduction/local_qc/`)

- `REPORT.md` — this document.
- `glm_summary.csv` — per-observer table: R², and raw and blank-subtracted responses (stationary,
  moving, blank) for V1/V2/V3 (4–8°), V4 and MT (whole ROI), and whole surface. Also at
  `~/Downloads/`.
- `example_design_matrix.png` — one run's design (blank frequency, pink-noise timing).
- `manuscript_caveat_paragraph.md` — draft caveat text.
- Scripts, re-runnable against `~/dg_collect/`: `group_qc.m` (R² by area), `group_models.m` (model
  parameters), `group_runs.m` (run consistency), `group_extra.m` (EPI/alignment),
  `design_uniqueness.m` (design fingerprinting), `group_rawtable.m` / `group_summary_csv.m` /
  `design_and_raw.m` / `group_addv4mt.m` (response tables incl. V4/MT + design figure).

*Extraction note: the copy runs read-only over a high-latency mount; run it under `caffeinate` (a
sleep dropped the mount mid-run once, and a dropped mount returns empty rather than erroring — so
verify file sizes, not just presence).*
