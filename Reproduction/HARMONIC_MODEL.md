# Per-vertex harmonic model — geometry vs context (2026-08-17)

> **Revised 2026-08-17b.** Vertices are now weighted for **equal polar-angle coverage**, which is
> what the published ROI analysis does implicitly and what makes the four predictors orthogonal.
> This is the primary specification; the previous natural-vertex-density numbers are kept
> alongside as a clearly labelled secondary column. See *Vertex weighting* below. Conclusions are
> unchanged; the individual numbers moved by up to ~0.03.

A model-based test of whether the Cartesian-vs-polar difference in the V1 orientation
asymmetries is explained by **within-ROI local-orientation geometry** rather than by
longer-range **context effects**.

Every analysis behind Figs 5–8 is ROI-based: vertices are binned into eight ±22.5° polar-angle
wedges, and each stimulus gets one label per wedge — a horizontal grating counts as "horizontal"
for the Cartesian asymmetries and as "radial" for the polar ones. But it is only exactly radial
for a vertex whose pRF sits on the horizontal meridian; across a 45°-wide wedge the local offset
runs over ±22.5°. The two stimulus sets are only *approximately* matched within those wedges.
This model removes the binning entirely and regresses on each vertex's own pRF polar angle.

> **Conclusion.** Within-ROI local-orientation geometry accounts for **6–8%** of the
> Cartesian-vs-polar gap in the horizontal/vertical and cardinal/oblique asymmetries (5.6% and
> 7.6% under the primary equal-coverage weighting; 7.6% and 8.4% at natural vertex density). The
> remaining ~93% survives, and survives equalising overall response gain and tightening the pRF
> quality floor. The manuscript's context claim stands on this test.
>
> **Corrected 2026-08-17c.** An earlier version of this document claimed that the radial/tangential
> asymmetry is *not* context-dependent. That is not supported. It shows no DETECTABLE difference
> between experiments, but the interval still admits an effect larger than the cardinal/oblique
> context effect, one observer drives the null, and the Cartesian-frame effect is not reliably
> larger than the polar-frame one. Report absence of evidence, not evidence of absence. See
> `cleanroom/diagnose_context_asymmetry.m`.

---

## The model

Per vertex *v*, take the four stationary-orientation betas, subtract the vertex's mean across
those four, concatenate across vertices, and predict with four global weights:

```
y_vk = b1*cos(2*theta) + b2*cos(4*theta)
     + b3*cos(2*(theta - thetaV)) + b4*cos(4*(theta - thetaV))
```

| term | asymmetry | +1 means |
|---|---|---|
| `b1` `cos(2θ)` | horizontal vs vertical | horizontal |
| `b2` `cos(4θ)` | cardinal vs oblique | cardinal |
| `b3` `cos(2(θ−θV))` | radial vs tangential | radial |
| `b4` `cos(4(θ−θV))` | polar-cardinal vs polar-oblique | polar-cardinal |

Demeaning removes the blank along with everything else common to the four conditions, so the
pink-noise-baseline problem (`local_qc/REPORT.md` §1) does not touch this model — only
orientation *differences* enter. All four predictors sum to zero across the four conditions in
both experiments, so no intercept is needed and demeaning does not distort the design.

Coefficients are reported as **2·b**, the pro-minus-con difference, matching the manuscript's
deltas, in **percent signal change**. Fitted per subject by weighted least squares (see *Vertex
weighting*), then averaged, with a 1000-resample bootstrap over the 8 observers — a pooled
vertex-level fit would weight observers by V1 surface area (their vertex counts range 3.4k–8.1k).

### Angle convention

`theta` is the orientation of the **bars**, in conventional visual-field degrees (0° = rightward
horizontal meridian, CCW positive): **0° is a horizontal grating, 90° a vertical one**. This
matches `AUDIT.md` §2 and `config_repro.oriAngle`, so the coefficients read directly against the
Fig 7 LME weights.

This is *not* the "direction of luminance variation" convention, which is rotated 90° and would
flip the sign of **both** first-harmonic terms: `b1 > 0` would mean vertical > horizontal, and
`b3 > 0` would mean tangential > radial, since a radial grating has its bars *along* the radius
(bar orientation = `thetaV`). The second harmonics are unaffected either way.

### Vertex weighting — equal polar-angle coverage (primary)

The published ROI analysis aggregates each of eight polar-angle wedges and then weights the eight
wedges **equally**. Under that weighting the four asymmetry predictors of `lme_codes` are exactly
orthogonal: their Gram matrix is `diag(16, 32, 16, 32)`, every off-diagonal zero.

A naive per-vertex fit does not inherit that. It weights by actual vertex density, and V1
over-represents the horizontal meridian. The consequence is not cosmetic. For this design the
uncentered correlation between the `b1` column, `cos(2θ)`, and the `b3` column, `cos(2(θ−θV))`, is
*analytically* the weighted mean of `cos(2θV)`:

```
r(b1,b3) = sum(w .* cosd(2*thetaV)) / sum(w)
```

At natural density that is **+0.35** (per-observer 0.20–0.46). `b1` and `b3` are precisely the pair
whose separation this model exists to adjudicate — how much apparent horizontal-vs-vertical is
radial-vs-tangential leaking through geometry — so a 0.35 design correlation is a confound
imported from cortical magnification, not from the question.

`harmonic_weights.m` fixes it: bin `thetaV` into 24 bins of 15° over [0,360), give each vertex
`w = 1/count(its bin)`, rescale to `mean(w) == 1`. Every occupied bin then carries the same total
weight. This drives the correlation to **+0.016**, restoring the published design's orthogonality
**while keeping θV continuous** — so the only remaining difference from the published analysis is
the thing actually under study, true θV versus wedge centre.

| | natural density | equal coverage |
|---|---|---|
| `r(b1,b3)` | +0.349 | +0.016 |
| `cond(X′WX)` | 3.31 | 2.07 |
| max VIF | 1.16 | 1.00 |
| effective n (Kish) | 100% | 52% of vertices |

The other four cross-block correlations (`b1`–`b2`, `b1`–`b4`, `b2`–`b3`, `b3`–`b4`) pair a 2nd
with a 4th harmonic and are zero *at every vertex individually*, hence zero under any weighting.
`b2`–`b4` remains the genuinely collinear pair (§ *Why the model is interpretable*), and equal
coverage improves it too.

Both the identity above and the four exact zeros are **asserted** to 1e-12 in
`run_harmonic_model` §2, for both experiments and both weightings. Checking the design diagnostic
against the closed form is a strong test that the same weights reached the design *and* the
summary: a `sqrt(w)`-vs-`w` slip, or weights applied to the fit but not the diagnostics, breaks it
immediately. (Weighted least squares must minimise `sum(w.*(y-X*b).^2)`, i.e.
`b = (X.*sqrt(w)) \ (y.*sqrt(w))`; each vertex contributes four stacked rows all carrying the same
weight.)

The price is precision: equalising coverage discards the redundancy of the over-sampled horizontal
meridian and leaves an effective n of about half the vertices. In practice the bootstrap CIs are
comparable in width, because they are taken over the 8 observers, not over vertices.

Two consequences worth knowing:

- **Fit A now reproduces the ROI pipeline exactly.** With θV quantised to the wedge centres, the
  weights are binned from those same centres, so every wedge carries equal weight — literally the
  published weighting. `FitA·2` now equals `LMEmean` to the printed precision, where the
  natural-density fit differed by up to 0.03.
- **For the Fit A vs Fit B comparison the weights are held fixed** (binned from the true pRF
  angle for both), so A − B is a pure θV manipulation. Without that pin, roughly a third of the
  `dg` horizontal−vertical A→B shift would be the weighting moving from 8 wedges to 24 bins rather
  than geometry. `fit_harmonic_vertex`'s `opts.weightSource` controls this.

`'natural'` remains available and is reported alongside throughout §2.

### Why the model is interpretable, not a black box

Four orientations at 45° spacing give each vertex's demeaned response vector exactly **three**
degrees of freedom: the 4th-harmonic pair collapses to one dimension because `sin(4θ)` and
`cos(4θ)` are proportional at the sampled points — degenerately so for `dg`, where `sin(4θ)` is
identically zero, and with a per-vertex ratio for `da`, where it is not (see the structural facts
below).

That count is **per vertex**, and it is not a cap on how many parameters can be estimated: `b1`–`b4`
are shared across vertices whose `thetaV` differs, so identification comes from *across*-vertex
variation in `thetaV`, not from the three within-vertex dimensions. The four-column design is full
rank. What the per-vertex structure does determine is *how* each coefficient is identified, and how
precisely.

Writing the three per-vertex amplitudes on `cos2`, `sin2`, `cos4` as *A, B, C*, the model is
**exactly equivalent** to three scalar regressions across vertices:

```
Cartesian (absolute frame)          Polar (radial-relative frame)
A = b1 + b3*cos(2*thetaV)           A = b3 + b1*cos(2*thetaV)
B =      b3*sin(2*thetaV)           B =    - b1*sin(2*thetaV)
C = b2 + b4*cos(4*thetaV)           C = b4 + b2*cos(4*thetaV)
```

The two experiments are exact mirror images, with the absolute and polar terms swapping roles.
Read off the first line: **`b1` and `b3` are separated only by the `thetaV` modulation of A, plus
the B channel** — which is precisely the question of how much apparent horizontal-vs-vertical is
radial-vs-tangential leaking through geometry. This makes `b1`/`b3` the **strong** pair: because
`b3` alone drives the B channel, the two separate even within a single vertex, provided it is off
the meridians where `sin(2*thetaV) = 0`. `b2`/`b4` have no such channel and need across-vertex
spread in `cos(4*thetaV)` — see the structural facts below.
`figures/cleanroom/harmonic_decomposition_raw.png`
plots these six panels with the fitted curves; the functional form holds well. It plots **2A, 2B
and 2C**, because each of those is exactly a pro-minus-con response difference between two
conditions — `dg`: horizontal−vertical, 45°−135°, cardinal−oblique; `da`: radial−tangential,
ccspiral−cspiral, polar-cardinal−polar-oblique — so the panels show the raw measured quantity in
percent signal change, and each curve's offset and amplitude are the same `2b` values reported
everywhere else. Its grey points are
the mean **across observers** of each observer's 15°-bin mean, with SEM across observers (n = 8).
The unit of inference is the observer, not the vertex: a vertex-level SEM over the ~400–500
vertices per bin would be roughly an order of magnitude too small, since vertices within an
observer are nothing like independent. Bins with fewer than 3 contributing observers are dropped.

Two structural facts, both asserted in `cleanroom/test_harmonic_model.m`:

- The two 4th-harmonic columns are proportional at each vertex with ratio `cos(4θV)`
  (`cos(4(θ−θV)) = cos(4θV)·cos(4θ)` for `dg`; the reverse for `da`). `b2` and `b4` are therefore
  identified only by the meridian-vs-oblique contrast. In practice this is benign here —
  cond(X′WX) = 2.1 and max VIF = 1.00 under equal coverage (3.3 and 1.16 at natural density),
  because polar-angle coverage is good.
- For `da`, `cos(2(θ−θV))` is **vertex-independent**. So the context-free null predicts `da`'s
  radial−tangential asymmetry to equal `dg`'s exactly. That is the sharpest form of the null, and
  the data do not reject it (below).

---

## Validation

The per-vertex machinery reproduces the published ROI pipeline before any of its departures mean
anything (§1 of the run output). Raw variant, % signal change:

| | | ref | pubMed | rtMed | rtMean | LMEmean | FitA·2 |
|---|---|---|---|---|---|---|---|
| **dg** | horiz−vert | −0.480 | −0.480 | −0.489 | −0.554 | −0.554 | −0.554 |
| | card−obl | −0.204 | −0.204 | −0.207 | −0.237 | −0.237 | −0.237 |
| | rad−tang | 0.116 | 0.116 | 0.118 | 0.109 | 0.109 | 0.109 |
| | polc−polo | 0.029 | 0.029 | 0.044 | 0.042 | 0.042 | 0.042 |
| **da** | horiz−vert | −0.211 | −0.211 | −0.209 | −0.215 | −0.215 | −0.215 |
| | card−obl | −0.030 | −0.030 | −0.033 | −0.032 | −0.032 | −0.032 |
| | rad−tang | 0.150 | 0.150 | 0.161 | 0.151 | 0.151 | 0.151 |
| | polc−polo | 0.034 | 0.034 | 0.044 | 0.040 | 0.040 | 0.040 |

`rtMean` equals `LMEmean` exactly (the four asymmetries are zero-sum contrasts, so demeaning
cancels under a linear aggregator); `rtMed` differs slightly from `pubMed` because the median is
not linear. **`FitA·2` now equals `LMEmean` to the printed precision** — under equal coverage the
per-vertex fit at the wedge centres weights each wedge equally, which is exactly what the ROI
route does. The natural-density fit differed by up to 0.03 (it weighted each observer's wedges by
vertex count); that residual gap is now closed, and the per-vertex machinery is an exact
re-expression of the published analysis before any of its departures.

### A note on units

The four constants per experiment transcribed in `cleanroom/validate_against_manuscript.m`
(`[-1.155 -0.40 0.23 0.06]`, `[-0.45 -0.06 0.60 0.17]`) are reproduced to three decimals by the
**z-scored** pipeline and not by the raw one, so as transcribed they are in σ units. Z-scoring
has been dropped (`local_qc/REPORT.md` §4; confirmed against the revised manuscript 2026-08-17),
so those are superseded. The percent-signal-change equivalents from the same published route are
the `ref` column above — **these are the values a non-z-scored manuscript should carry.**

---

## Result 1 — how much is within-ROI geometry? (6–8%)

Fit A quantises `thetaV` to the wedge centre; Fit B uses the vertex's true pRF polar angle.
Identical vertices, identical weights (both binned from the true pRF angle), **only `thetaV`
differs** — so A − B isolates the within-wedge local-orientation artifact, with no confound from
mean-vs-median, inclusion, or the weighting itself.

2·b, raw, with 95% bootstrap CI. **Primary — equal polar-angle coverage:**

| | Fit A (wedge centre) | Fit B (true pRF angle) | B − A |
|---|---|---|---|
| **dg** horiz−vert | −0.577 [−0.687 −0.458] | −0.577 [−0.686 −0.459] | 0.000 |
| **dg** card−obl | −0.240 [−0.352 −0.124] | −0.240 [−0.351 −0.124] | 0.000 |
| **dg** rad−tang | 0.114 [0.057 0.180] | 0.119 [0.061 0.194] | +0.005 |
| **dg** polc−polo | 0.043 [0.016 0.074] | 0.077 [0.040 0.121] | +0.034 |
| **da** horiz−vert | −0.208 [−0.305 −0.103] | −0.228 [−0.330 −0.117] | −0.020 |
| **da** card−obl | −0.032 [−0.066 −0.004] | −0.048 [−0.102 −0.004] | −0.015 |
| **da** rad−tang | 0.163 [0.025 0.266] | 0.162 [0.025 0.264] | −0.001 |
| **da** polc−polo | 0.033 [−0.034 0.097] | 0.034 [−0.033 0.097] | 0.000 |

**Secondary — natural vertex density** (the previous specification, one vertex one vote):

| | Fit A (wedge centre) | Fit B (true pRF angle) | B − A |
|---|---|---|---|
| **dg** horiz−vert | −0.540 [−0.661 −0.404] | −0.546 [−0.669 −0.409] | −0.006 |
| **dg** card−obl | −0.230 [−0.351 −0.116] | −0.234 [−0.357 −0.117] | −0.004 |
| **dg** rad−tang | 0.137 [0.074 0.207] | 0.145 [0.079 0.224] | +0.008 |
| **dg** polc−polo | 0.043 [0.018 0.071] | 0.078 [0.040 0.123] | +0.035 |
| **da** horiz−vert | −0.224 [−0.308 −0.129] | −0.254 [−0.344 −0.149] | −0.030 |
| **da** card−obl | −0.043 [−0.087 −0.012] | −0.062 [−0.129 −0.014] | −0.019 |
| **da** rad−tang | 0.148 [0.024 0.238] | 0.164 [0.047 0.250] | +0.016 |
| **da** polc−polo | 0.034 [−0.027 0.094] | 0.035 [−0.027 0.095] | +0.001 |

The quantity the manuscript's claim rests on is the **gap between experiments**:

| | Fit A gap | Fit B gap | explained by geometry | *(natural density)* |
|---|---|---|---|---|
| horiz−vert | −0.369 | −0.349 | **5.6%** | *7.6%* |
| card−obl | −0.208 | −0.192 | **7.6%** | *8.4%* |

**Correcting for exact per-vertex local orientation shrinks the Cartesian-vs-polar gap by 6–8%.**
It does not come close to abolishing it, under either weighting. Note also that the correction
moves the *polar* asymmetries slightly **away** from zero (da horiz−vert goes from −0.208 to
−0.228), i.e. the wedge binning was mildly *understating* the polar experiment's Cartesian-frame
asymmetries, not manufacturing the difference.

Equal coverage makes the mechanism cleaner to state. For `dg` the `cos(2θ)` and `cos(4θ)` columns
are vertex-independent, and once the design is orthogonal their coefficients no longer depend on
the `thetaV` specification at all — `dg` horizontal−vertical and cardinal−oblique are identical to
three decimals under Fit A and Fit B. **The entire geometry correction therefore lands on the
polar experiment's Cartesian-frame coefficients**, plus `dg`'s 4th-harmonic polar term. At natural
density the same correction was smeared across all eight coefficients by the 0.35 `b1`/`b3`
collinearity.

The one coefficient that moves materially is `dg` polar-cardinal−polar-oblique, 0.043 → 0.077
(+79%) — a 4th-harmonic term in `thetaV`, which is exactly the term most blurred by 45°-wide
wedges. It is small in absolute terms.

## Result 2 — cross-experiment prediction

Fit on one experiment, predict the other's per-vertex responses. Variance explained in the
held-out experiment, raw, 95% CI:

| | R² | ceiling | % of ceiling |
|---|---|---|---|
| dg → da | −0.313 [−0.709 0.032] | 0.262 | — |
| dg → da, free gain | 0.140 [0.067 0.217] (gain 0.46) | 0.262 | 53% |
| da → dg | 0.167 [0.064 0.272] | 0.443 | 38% |
| da → dg, free gain | 0.202 [0.112 0.306] (gain 1.32) | 0.443 | 46% |

`dg → da` without a gain is *worse than predicting zero*: the Cartesian coefficients applied to
polar geometry overshoot, because polar responses are about half the amplitude. With a single
free gain the transfer recovers roughly half the ceiling in both directions — the shared
structure is real but incomplete.

Nested models on the concatenated data:

| level | R² |
|---|---|
| (i) shared b1..b4 | 0.300 [0.250 0.345] |
| (ii) shared b1..b4 + one gain on da | 0.352 [0.289 0.414] (gain 0.51) |
| (iii) separate b1..b4 | 0.386 [0.315 0.462] |

Of the shared→separate gap, **61% is overall gain and 39% is genuine reference-frame difference**.
The (ii)→(iii) step is the context effect proper.

## Result 3 — which asymmetries are context-dependent

Coefficient differences `b_dg − b_da` (* = 95% CI excludes 0):

| | raw | gain-equalised |
|---|---|---|
| horiz−vert | −0.174 [−0.262 −0.100] * | −0.213 [−0.303 −0.126] * |
| card−obl | −0.096 [−0.147 −0.042] * | −0.107 [−0.160 −0.049] * |
| rad−tang | −0.022 [−0.081 0.072] | +0.014 [−0.035 0.076] |
| polc−polo | +0.022 [−0.017 0.070] | +0.027 [−0.000 0.059] |

**The two Cartesian-frame asymmetries are significantly larger in the Cartesian experiment**,
after accounting for exact per-vertex local orientation and after equalising overall gain — and
gain equalisation makes the difference *larger*, not smaller. This is the manuscript's claim, and
it survives.

**The two polar-frame asymmetries show no DETECTABLE difference between experiments.**
Radial−tangential is 0.119 [0.061 0.194] for `dg` and 0.162 [0.025 0.264] for `da`, difference
−0.043 [−0.161 0.143]. **Do not upgrade this to "radial/tangential is context-invariant."**
`diagnose_context_asymmetry.m` shows three reasons the null is not a result:

- the interval admits an effect as large as |0.158|, which exceeds the card−obl context effect
  (−0.175) that *is* significant, and is 59% of the horiz−vert one (−0.268);
- 6 of 8 observers show `da` > `dg` (sign test p = 0.29), and dropping sub-0395 alone — the only
  observer with a negative `da` radial−tangential value, −0.244 — makes the difference significant
  (−0.113 [−0.182 −0.036]). No other leave-one-out does;
- the within-subject difference of differences is n.s.: horiz−vert vs rad−tang 0.106 [−0.033 0.275],
  p = 0.26; card−obl vs rad−tang −0.004, p = 0.92. So the Cartesian-frame context effect is *not*
  reliably larger than the polar-frame one.

Report absence of evidence, not evidence of absence. Separating a genuinely one-sided context effect
from a two-sided one of unequal size needs more observers.

**Within-observer measurement error was measured, by resampling runs** (`diagnose_within_observer_error.m`,
on GLMsingle single-trial betas pulled by `server_extract/collect_runwise_betas.m`). It is 0.07–0.13
for the context differences, i.e. **23–39% of the across-observer variance**. Disattenuating changes
no conclusion — horiz−vert p 0.024 → 0.014, card−obl 0.006 → 0.0015, rad−tang 0.70 → 0.66 — so the
binding limitation is between-observer variability at n = 8, not measurement noise. (An earlier
vertex-resampling estimate of 0.02–0.03 was invalid and understated the error by 4–5×.)

**All context tests are within subject**, forming the per-observer difference first. An LME with
experiment × asymmetry interactions is **anti-conservative here and must not be quoted**: only the
intercept is random, so `fitlme` tests the interaction on DF = 502 (wedge-level observations) rather
than 7, with no Satterthwaite/Kenward-Roger correction, giving p smaller by 5–25× (horiz−vert
p = 0.0009 vs the paired 0.024). Adding random slopes, including on the interaction terms, leaves
DF at 502. Because the 4 × 8 design is balanced and orthogonal, the LME fixed effect is *identical*
to the mean of the per-observer contrasts (0.268459 either way), so the mixed model adds nothing to
the estimate and only misstates its uncertainty.

That equivalence is not confined to one contrast. Running both routes over the same `M` from
`bin_and_aggregate`, **all four asymmetries in both experiments agree to < 2e-16** (`dg` −0.479829,
−0.204375, 0.115887, 0.029248; `da` −0.211370, −0.029590, 0.150079, 0.033570 — identical either
way). Two conditions make it exact, and both hold here:

- the four codes from `lme_codes` are **exactly orthogonal**, to each other and to the intercept
  (cross-product matrix diagonal, 16/32/16/32; every column sums to 0, both experiments);
- the design is **completely balanced** — no empty wedge/subject cells in either experiment, minimum
  22 vertices per wedge.

Three consequences, all worth knowing before anyone proposes "fixing" the LME:

- **The Gaussian prior cannot touch the four asymmetries.** With `(1|subject)` only the subject
  intercept is shrunk, and each asymmetry contrast sums to zero within a subject, so the intercept
  cancels out of it exactly — shrunk or not. Random *slopes* would shrink each observer's BLUP
  toward the group mean, but the fixed effect would still be the unweighted mean of the per-observer
  contrasts.
- **Treating subject as a FIXED effect changes nothing.** Fixed-subject OLS returns the same four
  estimates to ~1e-16; it merely spends 7 more DF on intercepts (DFE 244 vs 251) while still fitting
  one asymmetry slope for everyone and testing it against wedge-level residuals. Random-vs-fixed is
  not the lever.
- **The lever is whether the asymmetry SLOPE varies by observer, and how DF is computed.** Fitting
  random slopes and then querying `fixedEffects(lme,'DFMethod','satterthwaite')` gives DF ≈ n−1 and
  very nearly reproduces the paired test. `fitlme`'s default *residual* DF reports the wedge-level
  number whatever the random structure — which is what the DF = 502 above refers to. The blocker is
  the DF method as much as the random-effects specification.

So the choice between the ROI route and the LME is **only** a choice of standard error; the
estimates are the same object. Note also that the equivalence depends on balance: deleting a single
wedge for a single observer makes the two routes diverge by ~1.2e-2. Not a concern at the published
inclusion criteria, but it would become one if a stricter filter emptied a wedge.

**This particular conclusion depends on the z-scoring decision.** In the z-scored variant the
radial−tangential difference between experiments is −0.204 [−0.341 −0.031], i.e. significant, and
the "polar asymmetries are stronger for polar gratings" reading returns. Since z-scoring has been
dropped, the raw result above is the one that stands — but the dependence should be known, and it
is the same mechanism `_archive/ZSCORE_FIG7.md` §1–§5 describes.

## Result 4 — is the four-term model adequate?

The complete harmonic basis at harmonics 2 and 4 adds four `sin` columns, which should vanish
under left–right visual-field symmetry. They do, with one exception (raw, Fit B):

- `dg`: `sin2ori` 0.016, `sin2rad` 0.001, `sin4rad` −0.000 — all CIs include zero.
- `da`: `sin2ori` −0.005, `sin2rad` −0.030 — CIs include zero; `sin4ori` −0.0116
  [−0.0216 −0.0029] excludes zero but is an order of magnitude below the core terms.

Core coefficients barely move when the `sin` columns are added (`dg` b1 −0.288453 → −0.288454;
`da` b1 −0.114082 → −0.114207).
The four-term model is adequate. The small `da sin4ori` is a chirality term — a CW-vs-CCW spiral
asymmetry the model does not capture — worth a footnote, not a fifth parameter.

---

## Caveats

### pRF polar-angle error — measured, and far too small to matter

This was the most serious threat to the result, because the two experiments depend on `thetaV`
through **different** terms:

| | `dg` (Cartesian) | `da` (polar) |
|---|---|---|
| `b1`, `b2` | regressor is `thetaV`-free → **unattenuated** | depends on `thetaV` → **attenuated** |
| `b3`, `b4` | depends on `thetaV` → **attenuated** | regressor is `thetaV`-free → **unattenuated** |

So angle error deflates `b1(da)` but not `b1(dg)`, inflating the very gap the context claim rests
on — it biases *toward* the conclusion. Forward-simulating a genuinely context-free world (one
shared coefficient set, plus angle error σ) reproduces the observed `b1` ratio of 0.395 at
**σ ≈ 39°**. So the question is entirely quantitative: how big is σ?

**It was measured, not assumed.** Every observer was fitted twice and independently, from two pRF
runs with different stimuli, held in sibling derivative trees: `prfvista_mov` — which is the
solution the analysis uses, verified identical to the local `ret_<subject>.mat` to
`max|difference| = 0` — and `prfvista`. `collect_prf_replicate.m` mirrors the second solution into
`~/dg_collect`; `diagnose_prf_angle_error.m` compares them over the published vertex set
(V1_REmanual, 4–8°, `vexpl > 0.1` in **both** solutions; 10,956 vertices):

```
pooled bias  +0.19°     circular SD 5.70°     robust SD 3.75°
sigma for ONE solution = circSD/sqrt(2):  3.90°  (robust 2.80°)
  -> lambda2 = exp(-2 sigma^2) = 0.991      lambda4 = exp(-8 sigma^2) = 0.964
```

**σ ≈ 4°, against the ~39° required.** Second-harmonic attenuation is under 1%. Three confirmations:

- Adding the measured error again as fresh jitter to the real data moves every coefficient by
  ≤ 0.009 (`dg`) and ≤ 0.003 (`da`).
- Attenuation-correcting every coefficient leaves the gaps essentially untouched: horiz−vert
  −0.3487 → −0.3466, card−obl −0.1922 → −0.1904.
- Given the measured σ, the context-free null predicts specific `da/dg` coefficient ratios. Three
  of four are badly violated, one of them in the *opposite direction*:

  | | null predicts | observed | |
  |---|---|---|---|
  | horiz−vert | 0.991 | 0.395 | ✗ |
  | card−obl | 0.964 | 0.199 | ✗ |
  | rad−tang | 1.009 | 1.361 | ✓ |
  | polc−polo | 1.038 | 0.435 | ✗ (sign of the discrepancy reversed) |

**One qualification.** The two fits share a scanning session (mov and static are both ses-02, or
ses-05/ses-03 for sub-0201/wlsubj124), so they differ in stimulus, not session. The disagreement
therefore captures noise- and stimulus-driven error but *misses* anything shared within a session —
surface registration, distortion correction, head position. Those would produce a spatially smooth
offset in `thetaV` rather than independent noise, biasing rather than attenuating. Reaching σ = 39°
would require ~38° of purely shared systematic error, which is not compatible with the clean
retinotopic organisation of these maps or with the 2–3° median between-solution agreement.

**Residual limits of this argument.** `b4` is the weakest leg of the ratio table — `da`'s `b4` CI
includes zero, so the reversed sign there is suggestive rather than decisive. And a smooth
registration bias remains formally unbounded by this test.

Raising the pRF R² floor bounds it a second way. Every coefficient grows with the floor, and the
`thetaV`-dependent ones grow most (`dg` rad−tang 0.119 → 0.178 from R² > 0.1 to > 0.5). But the
cross-experiment conclusion is unchanged at R² > 0.3: horiz−vert −0.184 [−0.276 −0.110] *,
card−obl −0.096 [−0.148 −0.039] *, rad−tang −0.018 [−0.074 0.076] n.s. (`polc−polo` becomes
marginally significant gain-equalised at that floor, 0.040 [0.019 0.065] — a small 4th-harmonic
term, not part of the main claim.)

**Other limits.** No per-vertex gain term — the model assumes one global set of weights; adding a
data-estimated gain is effectively the z-scoring that `local_qc/REPORT.md` §4 rejected. Only 3 of
4 harmonic degrees of freedom are observable per vertex, so the 4th-harmonic pair rests on
across-vertex `cos(4θV)` variation alone. The gain-equalised column in the **z-scored** variant is
numerically unstable (per-subject gains vary wildly; CIs run to ±5) and should not be read.

---

## How to run

```bash
matlab -batch "cd('Reproduction/cleanroom'); test_harmonic_model; run_harmonic_model(false)"
```

`test_harmonic_model` must pass first — its §1 asserts that the four predictors reduce exactly to
the repo's `lme_codes` at the eight wedge centres, which is the guard against the sign/frame slips
that produced both retracted bugs in `_archive/FINDINGS.md`. `run_harmonic_model` §2 then asserts the
weighting identities (see *Vertex weighting*) before printing any coefficients.
`run_harmonic_model(true)` runs the z-scored sensitivity variant.

Outputs land in `Reproduction/figures/cleanroom/` (git-ignored):
`harmonic_coefficients_<variant>.csv` (now carrying a `weighting` column, `equalcoverage` and
`natural` rows for both fits), `harmonic_decomposition_<variant>.png/pdf`,
`harmonic_coefficients_<variant>.png/pdf`, `harmonic_crossprediction_<variant>.png/pdf`; full
results struct in `cleanroom/_cache/harmonic_<variant>.mat`.

## Files

| File | Role |
|---|---|
| `cleanroom/harmonic_predictors.m` | the design matrix; owns the angle convention |
| `cleanroom/harmonic_vertex_data.m` | inclusion filter + per-vertex demeaned responses |
| `cleanroom/harmonic_weights.m` | equal-polar-angle-coverage vertex weights; owns the weighting |
| `cleanroom/fit_harmonic_vertex.m` | per-subject weighted LS, bootstrap, collinearity diagnostics |
| `cleanroom/predict_harmonic.m` | apply coefficients to an experiment's geometry |
| `cleanroom/harmonic_decompose.m` | the A/B/C reformulation |
| `cleanroom/harmonic_crossexp.m` | cross-prediction, gain, nested levels |
| `cleanroom/harmonic_roi_roundtrip.m` | push per-vertex values through the published ROI pipeline |
| `cleanroom/test_harmonic_model.m` | assertions — run first |
| `cleanroom/diagnose_prf_angle_error.m` | measures pRF polar-angle σ from the two independent fits |
| `cleanroom/diagnose_context_asymmetry.m` | within-subject context tests, leave-one-out, equivalence bound, and why the LME is not used |
| `cleanroom/diagnose_within_observer_error.m` | within-observer SE from split-half and bootstrap over runs |
| `server_extract/collect_runwise_betas.m` | per-run condition betas for V1 from GLMsingle single-trial fits (needs the mount, once) |
| `server_extract/collect_prf_replicate.m` | mirrors the second pRF solution to `~/dg_collect` (needs the mount, once) |
| `cleanroom/run_harmonic_model.m` | driver, prints the full report |
| `cleanroom/plot_harmonic.m` | the three figures |

Reuses `config_repro`, `load_allconditions`, `compute_vertex_contrasts`, `bin_and_aggregate`,
`compute_asymmetries`, `lme_codes`, `fit_lme_fig7`, `bootstrap_ci` unchanged.
