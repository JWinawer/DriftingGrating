# GLM data quality

*Quality-control pass on the GLMsingle output for all 8 observers × 2 experiments, 2026-07-24.
Reviewed: whole-surface `R2`, `R2run`, `FRACvalue`, `noisepool`, `HRFindex`, `pcnum`, `meanvol`,
cross-validation curves, the ~31 blank-referenced contrast maps, the retinotopy maps, every
FreeSurfer label, and the per-run design/timing files — copied read-only from `/Volumes/Vision` to
`~/dg_collect/` (1.2 GB). The raw single-trial betas (`modelmd`) were not copied.*

**Three conclusions.**

1. **No coding or processing error.** Model parameters are consistent and sensible across all 16
   sessions; runs are internally consistent; surfaces are correctly co-registered; there is no
   signal dropout (§2).
2. **No bad observers and no bad data**, with *positive* evidence rather than absence of error. The
   two observers that had looked anomalous (sub-0037, sub-0201) show normal responses in higher
   areas: MT is motion-selective and V4 prefers gratings to pink noise in both, in the same sessions
   (§2.4). All 8 observers are retained; the exclusion once argued for is **withdrawn**.
3. **The reference ("blank") condition is full-field pink noise, not a mean-luminance baseline**
   (§1). This changes how response magnitude and variance-explained should be read, but leaves the
   orientation asymmetries — the paper's claims — intact.

How reliable the *analysed quantities* are is a separate question, answered by resampling runs:
[`RELIABILITY.md`](RELIABILITY.md).

---

## 1. The blank is full-field pink noise

This is the central finding; it reframes everything else.

**What "blank" is.** In both experiments the blank condition is a **full-field 1/f (pink) noise**
image — the same texture that fills the surround around the grating aperture. This is explicit in the
stimulus code: `constConfig.m` builds `const.pinknoiseTex`; `my_stim.m` draws the pink noise, then
the grating within a circular aperture, then inner/outer masks; `my_blank.m` draws only the pink
noise; and `my_padding.m` plus the shared inter-trial code make every 2 s inter-trial interval and
the 10 s run-flanking padding full-field pink noise too.

So **pink noise is on the screen for the entire run; there is no mean-luminance baseline anywhere.**
The draft's Methods already describe this — *"blank stimuli contained only the pink noise background
that was present throughout the run."*

**Three consequences.**

1. Every reported response is **grating − pink-noise**; within the aperture it is essentially a
   grating-vs-pink-noise contrast, never grating-vs-nothing.
2. GLMsingle R² is variance **across conditions**. If pink noise, stationary gratings and moving
   gratings evoke similar BOLD, there is little variance to fit, so **R² and beta weights are near
   zero even where the region responds strongly to all of them.** Low estimates indicate *weak
   differentiation among the stimuli*, not absence of response.
3. **No true baseline exists, so absolute response magnitude cannot be recovered.** Comparisons of
   overall amplitude across areas or observers are differential responses relative to pink noise.

**Clean demonstration (sub-0037, V1, 4–8°).** Holding observer, area and eccentricity fixed, the only
difference is the experiment:

| sub-0037, V1 | stationary | moving | blank (pink noise) | spread | median R² |
|---|--:|--:|--:|--:|--:|
| Cartesian (dg) | −0.59 | 0.36 | −0.74 | 1.10 | 16.6 |
| polar (da) | −0.19 | −0.15 | −0.20 | 0.05 | 2.9 |

In `dg` the three conditions differ → variance to fit → R² 16.6. In `da` they are essentially
identical → no variance → R² floors at 2.9. That is variance-explained reflecting condition
*differences*, not response magnitude.

**What is unaffected.** The orientation **asymmetries** are contrasts *among* the grating
orientations, and the common blank term cancels in each. The pink-noise reference affects overall
magnitude and variance-explained, but **not the asymmetry results**.

### Draft paragraph for the Methods or Discussion

> In both experiments, the reference ("blank") condition was not a mean-luminance screen but a
> full-field 1/f (pink) noise image — the same background that surrounded the grating aperture — and
> this background was present continuously throughout each run, during the inter-trial intervals and
> run-flanking periods as well as the blank trials. The experiment therefore contains no true
> no-stimulus baseline, and all responses are estimated relative to this broadband-noise condition.
> Because the GLM partitions variance across conditions, a region that responds similarly to the
> gratings and to the pink-noise background yields small beta weights and low variance explained even
> if it responds robustly to all of these patterns; low response estimates thus indicate weak
> *differentiation* among the stimuli rather than an absence of response. Comparisons of overall
> response amplitude — across visual areas or observers — should accordingly be read as differential
> responses relative to the pink-noise background rather than as absolute BOLD magnitudes.
> Critically, the orientation asymmetries that are the focus of this study are contrasts *among* the
> grating orientations, within which the common reference term cancels; these measures are therefore
> unaffected by the nature of the reference condition.

---

## 2. No errors, no bad data

### 2.1 Fit quality and response by area

Median GLM R² and mean raw % response (stationary − blank), for V1/V2/V3 restricted to 4–8° (no
variance-explained threshold on the pRF or the GLM) and the whole surface unrestricted. Full
machine-readable values, including motion and un-subtracted responses, in `glm_summary.csv`.

| observer | V1 R²dg | V1 R²da | V1 βdg | V1 βda | V2 R²dg | V2 R²da | V3 R²dg | V3 R²da | wh R²dg | wh R²da |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| sub-0037 | 16.6 | 2.9 | 0.151 | 0.011 | 3.8 | 3.3 | 3.5 | 3.3 | 3.4 | 3.1 |
| sub-0201 | 3.5 | 4.5 | −0.224 | −0.349 | 3.1 | 3.9 | 3.0 | 3.7 | 3.0 | 3.9 |
| sub-0255 | 3.6 | 4.2 | 0.041 | 0.280 | 3.5 | 3.5 | 3.3 | 3.4 | 3.3 | 3.3 |
| sub-wlsubj123 | 4.8 | 4.4 | 0.543 | 0.424 | 3.4 | 3.7 | 3.4 | 3.7 | 3.4 | 3.5 |
| sub-wlsubj124 | 4.5 | 3.4 | 1.137 | 0.563 | 2.5 | 2.4 | 2.5 | 2.3 | 2.7 | 2.9 |
| sub-0395 | 26.4 | 19.8 | 1.558 | 1.304 | 8.0 | 6.2 | 7.1 | 7.1 | 3.0 | 3.9 |
| sub-0426 | 4.6 | 9.9 | −0.096 | 0.575 | 4.5 | 4.6 | 3.9 | 4.2 | 3.2 | 4.3 |
| sub-0250 | 6.3 | 3.6 | 0.353 | 0.378 | 8.1 | 3.7 | 4.9 | 3.6 | 4.0 | 3.2 |

> ⚠️ **Draw no spatial-specificity claim from this table.** 4–8° is the band where **spatial
> frequency is most closely matched between the two experiments** — it is *not* the stimulated
> extent, which was roughly **0.5–12°**. So a 4–8°-vs-whole-label contrast does not show that 4–8° is
> "the responsive patch"; the whole label also covers eccentricities the stimulus never reached.
> Measured directly, the group-median R²-by-eccentricity profile across V1 is roughly flat, for
> reasons that are properties of the measure and the design —
> [`RELIABILITY.md`](RELIABILITY.md) §4.

### 2.2 Model parameters, runs, alignment — all uniform

- **Model parameters** (`group_models.m`): 2–9 GLMdenoise noise regressors (median 5, normal range);
  fractional-ridge shrinkage effectively identical across sessions (median fraction 0.05, as expected
  when most of the surface is noise pool); HRF selection uniform. No session fit with anomalous
  settings.
- **Runs** (`group_runs.m`): counts vary (most 8; sub-0255/dg has 9, sub-0395/da has 6), but each GLM
  used its own count, and every session's per-run R² is internally consistent (worst-run /
  session-median 0.82–0.98). No mispaired or dead run.
- **EPI and surfaces** (`group_extra.m`): V1 mean EPI intensity uniform across sessions (no dropout);
  GLM surface length matches the retinotopy surface for every observer.

### 2.3 The two previously-flagged observers, explained

- **sub-0037** is the one observer with a large Cartesian-vs-polar difference in V1: `dg` R² 16.6,
  `da` R² 2.9. No error accounts for it — no bad run, no dropout, no misalignment, and the design is
  identical across observers (§3). Given the pink-noise reference, this is **weak differentiation of
  the polar gratings from pink noise, not absence of response**, and §2.4 confirms the session is
  sound. A genuine property of V1.
- **sub-0201** is the only observer whose mean V1 stationary−blank is negative in both experiments
  (−0.22 dg, −0.35 da): the windowed grating drove V1 no more than full-field pink noise. Expected
  under §1 for an observer in whom the two evoke comparable V1 responses. Its GLM fit quality is
  normal and its MT is strongly motion-selective (da moving−stationary +1.02, R² 19.2).

### 2.4 Higher areas (V4, MT) — positive evidence of good data

Because V1's small receptive fields respond similarly to gratings and to pink noise, a weak V1
grating−blank contrast is uninformative about data quality. Higher areas are the stronger test. V4
and MT are reported over the **whole ROI** (their pRFs extend well beyond 4–8°, and that band leaves
as few as 12 vertices for MT). Raw % signal change, from `group_addv4mt.m` / `glm_summary.csv`.

**MT is motion-selective (moving − stationary > 0) in every observer and both experiments:**

| observer | dg mov−stat | dg R² | da mov−stat | da R² |
|---|--:|--:|--:|--:|
| sub-0037 | +0.28 | 4.3 | +0.25 | 3.5 |
| sub-0201 | +0.32 | 3.1 | +1.02 | 19.2 |
| sub-0255 | +0.57 | 5.8 | +0.93 | 14.5 |
| sub-wlsubj123 | +0.73 | 13.9 | +0.74 | 12.6 |
| sub-wlsubj124 | +0.27 | 2.6 | +0.64 | 2.7 |
| sub-0395 | +0.78 | 10.6 | +1.02 | 17.1 |
| sub-0426 | +1.01 | 15.6 | +0.95 | 12.3 |
| sub-0250 | +0.75 | 9.9 | +1.10 | 20.8 |

**V4 prefers gratings to pink noise (stationary − blank > 0) in all 16 sessions** (+0.06 to +0.82; V4
shows little motion preference, as expected for a form area). For sub-0037, V4 stationary−blank is
+0.38 (dg) / +0.31 (da) even though its V1 barely differentiates in `da` (+0.01).

The picture is physiologically coherent — MT prefers motion, V4 prefers coherent gratings over noise,
V1 responds similarly to both — and establishes that **both** previously-flagged observers had a
normally responding visual system.

---

## 3. The stimulus used a fixed random seed

Fingerprinting every per-run design (`design_uniqueness.m`): the 110 Cartesian design files reduce to
12 distinct designs and the 68 polar files to 8, clustered exactly by run position. **Every
observer's run *K* is byte-identical to every other observer's run *K***, and Cartesian run *K*
equals polar run *K* (`rng` seed 0 with identical generator state). Runs differ within a session, but
the whole sequence is deterministic.

This is a property of the stimulus code, not an analysis bug. Its one consequence is that stimulus
order was not randomised per observer, so any stimulus-order confound is shared across observers
rather than averaging out. Worth confirming with the lead author that the fixed seed was intended. A
useful corollary: because all designs are identical, any accidental use of one observer's or one
experiment's design for another would have had no effect — which is part of why the "wrong design
matrix" concern is ruled out.

---

## Files (in `Reproduction/local_qc/`)

- `DATA_QUALITY.md` — this document.
- `RELIABILITY.md` — split-half reliability of the analysed quantities, measurement error against the
  effect, the Figure 4 controls, and draft manuscript text.
- `glm_summary.csv` — per-observer R², raw and blank-subtracted responses (stationary, moving, blank)
  for V1/V2/V3 (4–8°), V4 and MT (whole ROI), and whole surface.
- `splithalf_reliability.csv` — per-observer split-half reliabilities (the `.mat` is git-ignored).
- `example_design_matrix.png` — one run's design, showing blank frequency and pink-noise timing.
- Scripts, re-runnable against `~/dg_collect/`: `group_qc.m` (R² by area), `group_models.m` (model
  parameters), `group_runs.m` (run consistency), `group_extra.m` (EPI/alignment),
  `design_uniqueness.m` (design fingerprinting), `group_rawtable.m` / `group_summary_csv.m` /
  `design_and_raw.m` / `group_addv4mt.m` (response tables incl. V4/MT + design figure).

*Extraction note: the copy runs read-only over a high-latency mount; run it under `caffeinate` — a
sleep dropped the mount mid-run once, and a dropped mount returns empty rather than erroring, so
verify file sizes, not just presence.*
