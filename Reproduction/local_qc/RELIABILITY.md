# Reliability of the analysed measurements

**How reproducible are the quantities Figures 5–8 are built from, and how large is each effect
relative to the noise in its own measurement?** Both are answered by resampling the measurement —
splitting the runs — rather than by variance explained. Computed by
[`../cleanroom/splithalf_reliability.m`](../cleanroom/splithalf_reliability.m) and
[`../cleanroom/diagnose_within_observer_error.m`](../cleanroom/diagnose_within_observer_error.m);
per-observer values in `splithalf_reliability.csv`.

**Headline.** The orientation-differential ROI profile has a split-half reliability of **0.87**
(range 0.71–0.98) for Cartesian gratings and **0.62** (0.50–0.88) for polar gratings, exceeding 0.50
in all 16 sessions. Every asymmetry the paper reports is 2.1–5.3× its own within-observer standard
error, and disattenuating changes no conclusion. **The binding limitation is between-observer
variability at n = 8, not measurement noise.**

---

## 1. Why reliability rather than R²

R² is the wrong currency for this dataset. The reference condition is full-field pink noise present
throughout every run ([`DATA_QUALITY.md`](DATA_QUALITY.md) §1), so GLM R² measures differentiation
*among* conditions, not responsiveness. Median V1 R² in the analysed 4–8° band is 4.7% (dg) and 4.4%
(da) against 3.3% for the whole surface — a reader who takes that at face value concludes the data
are bad.

The draft's "GLM summary" section currently reports three different quantities all called R² —
variance explained in a group-mean ROI time series (0.42 / 0.28), fit to a trial-triggered average
(0.37–0.21), and vertexwise GLM R² (Fig 4C) — each from single examples, and none of them is the
quantity the claims rest on. What Figures 5–8 require is that the analysed profile be
**reproducible**.

## 2. Method

The design is balanced — 8 runs, 4 trials per condition per run — so a 4-vs-4 run split is balanced by
construction. For each observer and experiment: recompute the analysed profile (4 orientations × 8
polar-angle wedges) from each half, correlate the two halves, Fisher-z average over **all 35**
balanced splits, and Spearman–Brown correct to the full run count. For the two observers who depart
from 8 runs (sub-0255 dg has 9, sub-0395 da has 6) each run is dropped in turn and the even procedure
averaged over the drops, correcting to 8 rather than 9 — conservative.

Reliability is reported at two spatial levels × two contrast levels:

| | includes overall response level | orientation-differential only |
|---|---|---|
| **per vertex** (nVert 1126–1786) | `vertexFull` | `vertexDiff` |
| **ROI profile** (32 values — what Figs 5–8 analyse) | `roiFull` | `roiDiff` |

For `roiFull`, the 32 numbers are the response to each of the 4 orientations in each of the 8 ROIs,
for one observer in one experiment — the same 32 numbers Figures 5/6 plot. Each half of the runs
yields its own 4 × 8 matrix; *r* is the Pearson correlation of the flattened matrices. For `roiDiff`
each ROI's mean over the 4 orientations is first removed, so the 32 values carry 24 independent
degrees of freedom. That strips response-level differences which are reliable but carry no
orientation information, leaving exactly the subspace the four asymmetries live in. **`roiDiff` is
the strict test and the one to report.**

**Validation.** This route (run-averaged GLMsingle single-trial betas) reproduces the ROI pipeline,
and its within-observer SEs agree with the independent bootstrap-over-runs estimator to ±0.005. The
vertex counts match `gainSummary.csv`'s independently-computed `nVertices` exactly for all 8
observers, confirming the same V1 / 4–8° / pRF R² > 0.1 selection.

## 3. Results

### 3.1 Split-half reliability

Spearman–Brown corrected, per observer. **Invariant to the gain rescaling** — a per-observer scalar
cancels in a within-observer correlation — so these do not depend on that choice at all.

| observer | nRun dg/da | vertexFull dg | da | vertexDiff dg | da | roiFull dg | da | **roiDiff dg** | **da** |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| sub-0037 | 8 / 8 | 0.734 | 0.121 | 0.797 | 0.301 | 0.899 | 0.328 | **0.916** | **0.500** |
| sub-0201 | 8 / 8 | 0.485 | 0.693 | 0.404 | 0.202 | 0.675 | 0.578 | **0.713** | **0.546** |
| sub-0255 | 9 / 8 | 0.631 | 0.714 | 0.748 | 0.707 | 0.823 | 0.894 | **0.854** | **0.877** |
| sub-wlsubj123 | 8 / 8 | 0.863 | 0.753 | 0.811 | 0.395 | 0.916 | 0.833 | **0.944** | **0.691** |
| sub-wlsubj124 | 8 / 8 | 0.924 | 0.602 | 0.710 | 0.144 | 0.931 | 0.649 | **0.880** | **0.556** |
| sub-0395 | 8 / 6 | 0.944 | 0.903 | 0.894 | 0.668 | 0.976 | 0.893 | **0.976** | **0.800** |
| sub-0426 | 8 / 8 | 0.707 | 0.734 | 0.713 | 0.482 | 0.752 | 0.756 | **0.845** | **0.507** |
| sub-0250 | 8 / 8 | 0.615 | 0.718 | 0.690 | 0.651 | 0.663 | 0.775 | **0.797** | **0.823** |
| **median** | | 0.721 | 0.716 | 0.731 | 0.438 | 0.861 | 0.766 | **0.867** | **0.623** |

**The ROI level is the right level, and the data say so.** Per-vertex orientation-differential
reliability is much lower (median 0.73 dg, **0.44** da), but each analysed value is a mean over
~1100–1800 vertices, so per-vertex reliability is not what any claim depends on; reporting it would
understate the reliability of the actual measurements.

### 3.2 Measurement error against the effect

Across-vertex mean, gain-rescaled. `SNR` = |group mean| / within-observer SE; `p disatt` removes the
measurement variance and re-tests.

| | asymmetry | group mean | SD across obs | within-obs SE | SNR | *p* observed | *p* disatt |
|---|---|--:|--:|--:|--:|--:|--:|
| **dg** | horiz−vert | −0.547 | 0.158 | 0.127 | 4.3 | <0.0001 | <0.0001 |
| | card−obl | −0.221 | 0.158 | 0.073 | 3.0 | 0.0056 | 0.0030 |
| | rad−tang | 0.104 | 0.080 | 0.037 | 2.8 | 0.0081 | 0.0044 |
| | polc−polo | 0.039 | 0.044 | 0.018 | 2.1 | 0.0396 | 0.0274 |
| **da** | horiz−vert | −0.218 | 0.176 | 0.041 | 5.3 | 0.0098 | 0.0086 |
| | card−obl | −0.029 | 0.045 | 0.020 | 1.4 | 0.1129 | 0.0811 |
| | rad−tang | 0.150 | 0.170 | 0.136 | 1.1 | 0.0416 | 0.0043 |
| | polc−polo | 0.040 | 0.105 | 0.073 | 0.5 | 0.3185 | 0.1795 |

*(These are the pre-specification ROI-route values, which is the route the split-half work was built
on. The specification's values are within 0.03 of them and change nothing here —
[`../RESULTS.md`](../RESULTS.md) §2.)*

**What survives as a claim.** The effects are large relative to their own measurement error (SNR ≥ 2.8
for all four Cartesian-experiment asymmetries and for `da` horiz−vert); disattenuation changes no
conclusion — every asymmetry significant at *p* < 0.05 stays significant and every null stays null;
and the only genuinely noise-limited cells are `da` rad−tang and `da` polc−polo, both of which a
better measurement would *strengthen*.

> ⚠️ **Do not headline the "fraction of between-observer variance due to measurement error".** It is
> large either when measurement error is large **or when observers genuinely agree with one another**
> — opposite situations. `dg` horiz−vert reads 65% because the effect is large and near-identical in
> every observer, so there is little true variation for measurement error to be a small fraction *of*
> (per observer, |effect| / own SE runs 1.9–7.9, median 5.0; the effect is negative in all 8 observers
> and all 8 ROIs). `da` rad−tang reads 64% because it is genuinely noise-limited (SNR 1.1).
> **Report SNR and the consistency counts instead** — SNR separates the two cases, the ratio
> conflates them.

> ⚠️ **Reliability must come from resampling the *measurement*, i.e. runs.** An earlier estimate
> resampled *vertices*, which holds the GLM betas fixed and only reshuffles which vertices enter the
> wedge median — it characterises which patch of V1 was sampled, not the reliability of the
> measurement, and vertex responses are not independent. It understated the error by 4–5×.

### 3.3 A structural fact worth knowing

The within-observer SE is roughly **3× larger for the asymmetry matched to each experiment**
(`dg` horiz−vert 0.121, `da` rad−tang 0.121) than for the derived one (0.042 in each case):

| asymmetry | dg | da |
|---|--:|--:|
| horiz−vert | 0.121 | 0.042 |
| card−obl | 0.076 | 0.024 |
| rad−tang | 0.042 | 0.121 |
| polc−polo | 0.024 | 0.067 |

A matched contrast uses the same two stimulus conditions in every wedge, so its noise does not
average across wedges; a derived contrast rotates which stimuli it draws on and averages more
effectively. Worth knowing when designing a similar comparison.

---

## 4. Figure 4 — the controls

### 4A, the run-mismatch control

Figure 4A's mismatched-run comparison was a single anecdote (two runs, one ROI). For a **fixed** (not
refit) prediction — which is why R² can go negative — writing ρ = corr(ŷ_c, ŷ_m) and *k* = ‖ŷ_m‖/‖ŷ_c‖:

> R²_m = *Q* · (2ρ*k* − *k*²),  where *Q* = ‖ŷ_c‖²/‖y‖²

**The sign result is measured and settled.** R²_m < 0 whenever *k* > 2ρ. The fixed-seed designs are
near-orthogonal — median ρ = 0.006 (dg), 0.020 (da), *k* ≈ 1.04 — so the ratio is negative for every
run pair in all 16 sessions, and direct measurement on sub-0037/dg confirms **measured R²_m negative
in 56 of 56 run pairs**, at both per-vertex and ROI-mean level. That is the claim the manuscript
needs. The check also *explains* the figure: the difference between −0.25 and −0.53 in the two panels
is almost entirely *k*, the ratio of predicted amplitudes, not a difference in data quality.

> ⚠️ **Use *Q*, not R²_c, in that formula, and use it only for a group statement.** The original form
> used R²_c and assumed ⟨e, ŷ_c⟩ = 0; measurement showed that term carries the whole bias (−0.079)
> while ⟨e, ŷ_m⟩ — the term originally flagged as the risk — is negligible (−0.010). Substituting *Q*
> cuts the mean gap from −0.087 to −0.014 and the accounting closes. But the correlation between
> measured and predicted R²_m across the 56 pairs is only 0.15, and the *Q* rescaling *increases*
> per-pair scatter even as it fixes the mean. **Never use it to predict a specific run pair.**

**What averaging over vertices buys, measured:** per-vertex R²_c 0.051 against ROI-mean 0.181 — only
**3.5×** from 1232 vertices, nowhere near what independent noise would give, because the mean
cross-vertex residual correlation is **0.273**. This is the quantity that blocks converting the
scale-free ratio into an absolute ROI-level R², and it could not be derived from predictions.

> **Note on comparability.** Measured per-vertex R²_c (0.051) is well below GLMsingle's own `R2run`
> (0.174), for two reasons, neither a bug: GLMsingle fits 52 **single-trial** betas per run where this
> uses 13 condition betas fit across runs, and the extraction projects out only polynomials while
> GLMsingle also projects noise-pool PCs. The values are not interchangeable.

**Still open:** check **sub-0201/da** (ρ = 0.36, only 66% of pairs negative — some of its run-pair
designs are genuinely correlated) before stating the claim over all run pairs. Code:
[`../cleanroom/run_mismatch_local.m`](../cleanroom/run_mismatch_local.m),
[`check_fig4a.m`](../cleanroom/check_fig4a.m),
[`calibrate_run_mismatch.m`](../cleanroom/calibrate_run_mismatch.m).

### 4C, the R²-by-eccentricity map — no quantitative claim

**Recommendation: either drop Figure 4C, or present it as an illustration for one observer with no
quantitative claim attached.**

4–8° is the window where spatial frequency is matched between the experiments, **not the stimulated
extent**, which was roughly 0.5–12°. Measured directly (V1 median GLM R² by eccentricity, `dg`, all 8
observers), the group profile is essentially flat — 3.6 to 6.5 across 0–25°, with no drop beyond 12°.

That flatness is what *should* be expected, for three reasons:

1. **R² is sign-blind.** A vertex *suppressed* by the stimulus has as much cross-condition variance as
   a driven one and earns the same R². Variance explained cannot localise a stimulus even in
   principle.
2. **pRF extent, not just centre.** A pRF whose centre lies beyond the stimulus can still overlap it.
   Measured: at 14–18°, median σ = 2.17 and 59% of V1 vertices have (ecc − 2σ) < 12°; at 18–25°,
   only 18% do. Extent accounts for the near periphery, not the far.
3. **Surround suppression, with no way to see it.** Vertices well outside the stimulus would normally
   show a *negative* response, but the rapid event-related design plus the absence of any true blank
   makes that transition invisible. Measured: median stationary−blank stays positive at every
   eccentricity and *rises* (+0.17 at 4–8°, +0.42 at 18–25°) — no crossing at all.

> ⚠️ **Vertices with pRF centres beyond 12° are not an unstimulated control.** The pRF *mapping*
> stimulus only extended to ~12°, so such a vertex must nonetheless have responded to a stimulus
> filling the same aperture as the gratings — it is an extrapolated fit from a vertex that *did*
> respond. No eccentricity contrast built from them is a valid test. (They also rest on few vertices
> — median n ≈ 100–120 beyond 14° against ≈1350 at 4–8° — and their median pRF σ is non-monotonic,
> pointing to poorly-constrained fits at the label edge.)

**A split-half reliability map would not do better.** Reliability is *also* sign-blind, so a reliably
suppressed vertex scores as high as a reliably driven one. Where reliability wins is
*interpretability of magnitude*: in the responsive patch it reads 0.86–0.92, numbers a reader can
evaluate, where R² reads 16.6% and 26.4%, which look poor to anyone unaware of the pink-noise
reference. But that is an argument for the reliability **statistic** (§3.1), not for rendering it on a
surface — and building the map would need `modelmd` re-pulled (~6.8 GB), since
`collect_runwise_betas.m` saved V1 only. **Do not spend a pull on it.**

*(A moving-vs-stationary map would be genuinely diagnostic — MT+ positive in 8/8 observers, hV4 near
zero, and the pink-noise reference cancels exactly since both conditions share it. It is set aside
only because this paper analyses the stationary conditions. Worth remembering for a paper where it is
not.)*

---

## 5. Draft text for the "GLM summary" section

*Keep the existing Fig 4A–C paragraphs as illustration; these follow them.*

> **Before transcribing:** the SNRs and *p* values in the third paragraph are §3.2's, i.e. the
> **pre-specification ROI route** the split-half work was built on, not the settled specification.
> They are within 0.03 of the specification's values and change nothing qualitatively
> ([`../RESULTS.md`](../RESULTS.md) §2), but re-read them off the current numbers before they go
> into the manuscript. The reliabilities in §3.1 are unaffected — a per-observer scalar cancels in a
> within-observer correlation, so they do not depend on the route or the gain at all.

> Because the reference condition was full-field pink noise rather than a mean-luminance screen, and
> this background was present continuously throughout each run, the variance explained by the GLM
> reflects differentiation among the stimulus conditions rather than overall responsiveness: a region
> that responds strongly but similarly to the gratings and to the pink-noise background yields low
> variance explained. We therefore assessed the quality of the first-stage estimates by their
> reliability — the reproducibility of the response estimates that enter the asymmetry analyses —
> rather than by variance explained.
>
> Each observer's runs were split in half. Because the design was balanced (four trials per condition
> in every run), a 4-vs-4 split of the eight runs is balanced by construction. We recomputed the
> analysed profile — the response to each of the four orientations in each of the eight polar-angle
> ROIs — separately from each half, correlated the two halves, averaged over all 35 balanced splits,
> and applied the Spearman–Brown correction to the full run count. Because differences in overall
> responsiveness among ROIs are reliable but carry no orientation information, we report the
> reliability of the profile after removing each ROI's mean response across the four orientations,
> which leaves exactly the component that carries the four asymmetries. This orientation-differential
> profile had a split-half reliability of r = 0.87 (range across observers 0.71–0.98) for the
> Cartesian experiment and r = 0.62 (range 0.50–0.88) for the polar experiment, exceeding 0.5 in all
> 16 sessions (Table SX). We report reliability at the level of the ROI profile rather than of
> individual vertices because each analysed value is a summary over more than a thousand vertices;
> the reliability of a single vertex is correspondingly lower (median r = 0.73 and 0.44 for the two
> experiments) but is not the quantity on which any of the following analyses depends.
>
> The same run splits give the measurement error of each asymmetry directly, so each effect can be
> expressed relative to the noise in its own measurement. In the Cartesian experiment the
> horizontal–vertical asymmetry was 4.3 times its within-observer standard error, and the
> cardinal–oblique, radial–tangential and polar cardinal–polar oblique asymmetries were 3.0, 2.8 and
> 2.1 times theirs; in the polar experiment the horizontal–vertical asymmetry was 5.3 times its
> standard error and the radial–tangential asymmetry 1.1 times. Individual observers therefore
> resolve the larger asymmetries on their own, which is why the horizontal–vertical difference in the
> Cartesian experiment was negative in all eight observers and at all eight locations. Removing the
> measurement variance entirely and re-testing changes no conclusion: every asymmetry significant at
> p < 0.05 remains so and every null remains null. The exception in the other direction is the
> radial–tangential asymmetry in the polar experiment, whose measurement error is comparable to the
> effect itself (p = 0.042 observed, 0.004 disattenuated); that effect is limited by measurement
> precision rather than by the number of observers, and we treat it as marginal throughout.
>
> As an independent check that the GLM recovered known response properties, we examined two regions
> about which this study makes no claims. Area MT responded more strongly to moving than to
> stationary gratings in all 16 sessions (+0.25 to +1.10% signal change), and area V4 responded more
> strongly to gratings than to the pink-noise background in all 16 sessions (+0.06 to +0.82%)
> (Table SY). The pattern is physiologically coherent — MT prefers motion, V4 prefers coherent
> gratings over noise — and confirms that the sessions contain well-differentiated responses,
> including in the two observers whose V1 responses were weakest.

**Fig 4D suggestion.** A split-half scatter — the 32 orientation-differential profile values from one
half against the other, one colour per observer — puts the reliability claim in the same "show the
data" register as 4A–C. Tracked in [`../../AGENTS.md`](../../AGENTS.md) §5.

**Supplementary tables.** Table SX = `splithalf_reliability.csv`. Table SY = the MT / V4 columns of
`glm_summary.csv` ([`DATA_QUALITY.md`](DATA_QUALITY.md) §2.4).
