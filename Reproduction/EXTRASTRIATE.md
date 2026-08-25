# Orientation asymmetries across the visual hierarchy

**Status: specification settled 2026-08-24.** One analysis, applied unchanged to all eight
supplement maps — V1, V2, V3, V3a, V3b, hV4, MT, MST —
arrived at by working through four decisions in turn. Section 1 is the specification; sections 2–5
are the evidence for each choice; sections 6–8 are the results and what qualifies them. Section 9
records two traps in the code.

Numbers: [`supplement/precision_weighted_areas.csv`](supplement/precision_weighted_areas.csv) (eight
maps x two eccentricity bands, 192 rows, each with its coverage and a `reportable` flag), [`supplement/gain_areas_summary.csv`](supplement/gain_areas_summary.csv).
Code: [`cleanroom/diagnose_within_observer_error.m`](cleanroom/diagnose_within_observer_error.m),
[`cleanroom/precision_weighted_table.m`](cleanroom/precision_weighted_table.m),
[`cleanroom/precision_weighted_cells.m`](cleanroom/precision_weighted_cells.m),
[`server_extract/collect_gain_areas.m`](server_extract/collect_gain_areas.m),
[`server_extract/collect_runwise_betas_areas.m`](server_extract/collect_runwise_betas_areas.m).

---

## 1. The specification

| | choice |
|---|---|
| **Gain** | on, per observer × **map**, applied at the observer boundary |
| **Model** | four-term harmonic, fitted per vertex, **continuous θv** |
| **Weighting** | equal coverage at **45°** — the polar-angle ROIs themselves |
| **Fitting** | per observer, then averaged across observers |
| **Precision** | σ from bootstrapping runs; **equal weighting primary**, precision-weighted reported alongside |
| **Diagnostics** | VIF and design correlations per map; polar-angle resolution reported only where coverage passes the §6 criterion |

No per-map customisation anywhere. Which maps get reported at polar-angle resolution is a
reporting decision taken on a measured criterion (§6), not a change of method.

## 2. Why the V1 procedure needed replacing

The V1 route took each observer's asymmetry as the mean over their eight polar-angle ROIs. That is
only meaningful if every observer contributes the same eight. In V1 they do — the sparsest cell
holds 22 vertices. Outside V1 they do not, and two facts make the gaps dangerous rather than
inconvenient:

- **The asymmetry varies strongly with polar angle.** For dg polar-cardinal−polar-oblique the
  fitted profile alternates ±0.13–0.26 with cardinal/oblique ROI, against a mean effect of 0.04.
- **The loss is not random with respect to polar angle.** Map boundaries lie on the meridians, and
  vertices near a boundary have the poorest pRF fits, so cells empty along exactly the axis the
  asymmetry varies along.

Deleting meridian ROIs from V1 in simulation biases the ROI route by **+0.071** on dg polc−polo, an
effect whose true value is 0.040 — cell loss nearly triples it.

Note what is *not* broken. The inclusion mask is orientation-independent, so all four orientations
of a cell vanish together and the within-ROI subtraction is never partially missing. Only the
average over ROIs breaks.

## 3. Decision 1 — gain, per observer × map

The pRF-gain rescaling is a per-observer scalar applied before observers are combined. It touches
no within-observer quantity, so it is separable from every other decision — but it does **not**
commute with precision weighting (scaling `y_i` and `σ_i` by `c_i` changes τ̂² and hence the
weights), so it is applied at the observer boundary, never to a finished group estimate.

It was briefly turned off, because the only gain available was computed in V1 alone and applying a
V1 number to V2/V3/hV4 special-cases V1 as the source. [`collect_gain_areas.m`](server_extract/collect_gain_areas.m)
removed that objection by computing gain over all eight maps — the V1 restriction in the original
`dg_computeGain` was for speed, not design, so widening it needed no transfer beyond the same
`results.mat`.

Gain falls up the hierarchy while observer heterogeneity stays flat:

| | V1 | V2 | V3 | hV4 (2–10°) |
|---|---|---|---|---|
| group geometric mean, % BOLD | 4.369 | 3.843 | 3.188 | 3.654 |
| across-observer range | 1.67× | 1.75× | 1.80× | 1.77× |

**The rescaling is amplitude-neutral per map** — the geometric mean of the scale factors is exactly
1 (verified to 12 decimals in every map), so it equalises observers *within* a map without altering
that map's response level. The V1 > V2 > V3 decline, and the asymmetry attenuation reported below,
are therefore not artefacts of the normalisation.

Per-map scale factors correlate with the V1-derived one at r = 0.88 (V2), 0.73 (V3), 0.66 (hV4), so
a V1 scalar is a fair proxy for V2 and a poor one for hV4 — but the whole thing is worth at most
0.0073 in the group estimate (0.0127 for gain versus no gain), the smallest of the four decisions.

Validation: V1 4–8° reproduces the existing `gainSummary.csv` for all eight observers to **5.3e-15**
via the independent `dg_computeGain` path. The mov/stat protocols are combined by **geometric**
mean, superseding that file's arithmetic average.

## 4. Decision 2 — equal-coverage weighting at 45°

Vertex density across the eight ROIs spans 5.4× in V1, 3.4× V2, 2.4× V3, 95× hV4. Natural weighting
therefore estimates an average over **cortex**; equal coverage estimates an average over the
**visual field**. When the question is polar-angle dependence, the cortical average confounds the
asymmetry with cortical magnification. **This is an estimand choice, not numerical hygiene** — the
proof being V2, where `r(b1,b3) = 0.001` under natural weighting (no collinearity to fix) and the
estimate still moves by 0.0275.

That it also orthogonalises the design is a second, weaker, map-specific benefit: mean `r(b1,b3)`
under natural weighting is 0.349 in V1, 0.182 in V3, 0.375 in hV4 — but 0.001 in V2. What matters
is not the density range but its alignment with cos(2θv): V1's excess sits at 0° and 180° where
cos(2θv) = +1 and does not cancel; V2's sits at obliques and does.

**45°, not 15°.** The bin width controls leverage, since `w = 1/count` gives a bin holding one
vertex the same total weight as a bin holding two hundred:

| map | width | min count | max weight ratio | eff N | r(b1,b3) |
|---|---|---|---|---|---|
| V1 | 45° | 33 | **11×** | 1005 | 0.075 |
| V1 | 15° | 2 | 101× | 742 | 0.016 |
| V3 | 45° | 6 | 27× | 433 | 0.086 |
| V3 | 15° | 1 | 89× | 249 | 0.035 |

A residual `r(b1,b3)` of 0.075 is a variance inflation of 1.006 — nothing. A 100× leverage on a
single vertex is a real hazard. 45° also keeps **one** binning in the pipeline and makes the
estimand exactly "equal weight per polar-angle ROI", which is what the ROI route computes.

Finer binning does not help: with continuous θv in V1 the residual runs natural 0.349, 8 bins 0.075,
24 bins 0.016, 48 bins 0.009, **96 bins 0.023, 360 bins 0.080** — inverse-count weighting on fine
bins is a noisy density estimate, so it turns back up. A von Mises kernel gives 0.026, no better.

**Empty weighting bins are a non-event.** They are not analysis units; only vertices that exist are
indexed, so an empty bin never divides and never creates a missing value (verified: 0 non-finite
weights across all maps and observers). The bin width has nothing to do with the missing-data
problem.

## 5. Decisions 3 and 4 — precision weighting, and identifiability

**Equal weighting is primary.** With 8 runs, σ̂ carries ~7 df and scatters even when every observer
has identical true precision: the null max/min over 8 observers has median 2.20× and 95th percentile
3.55×, against an observed median of 3.37×. Only 14 of 32 rows exceed the null, so most of what you
would weight by is estimation noise. Simulation at this dataset's operating point (τ median 0.049,
true spread ~2–2.5× after discounting the chance component) puts precision weighting at **4% worse**
than equal weighting; it only pays when τ is small *and* the spread is large.

Empirically it changes little — max |equal − precision| is 0.029 (V1), 0.016 (V2), 0.017 (V3),
0.024 (hV4) — and where it does most work is the τ̂² = 0 rows, which are the least trustworthy. Both
columns are reported; this is only about which is quoted.

A real and predictable structure sits underneath: **matched-frame contrasts are ~3× noisier**
(mean σ 0.081 versus 0.027) because the contrast projects onto correlated trial-wise noise, and they
carry more genuine between-observer spread (9/16 rows exceed the null, versus 5/16). But that
variation is *between asymmetries*, and precision weighting only ever uses the between-observer part.

**Identifiability is continuous, and nothing is near degenerate at real coverage** — max |r| = 0.467,
max VIF 1.28 across all maps, weightings and θv choices. But the second-harmonic pair deserves care.
Per cell, polc−polo is exactly ± card−obl with the sign alternating cardinal/oblique (maximum
difference over all cells: 0.000, both experiments), so the two are one measurement separated only
by how cells combine across ROIs — which is what ROI loss damages. **In a map with class-structured
ROI loss they are not two findings**, under either route.

That is also why **θv is continuous**. Binned θv puts cos(4θv) at only two values, making
`b4 = b2 × cos(4θv)` a two-point design — exactly orthogonal at full coverage, exactly degenerate
(VIF = ∞) once one ROI class is lost. Continuous θv never degenerates: VIF 4.12 with all four
cardinal ROIs removed.

## 6. Coverage, and which maps can be reported

The supplement covers eight maps. All eight were extracted and all have per-map gain; whether each
can support a *polar-angle-resolved* analysis is a separate, measurable question. Criterion: at most
2 empty (observer × ROI) cells of 64, a median of at least 20 vertices per cell, and a maximum
precision-weight ratio below 25.

| map | band | empty cells | median vertices/cell | max weight ratio | reportable |
|---|---|---|---|---|---|
| V1 | 4–8° | 0 | 150 | 17.1 | **yes** |
| V1 | 2–10° | 0 | 326 | 18.1 | **yes** |
| V2 | 4–8° | 0 | 136 | 15.8 | **yes** |
| V2 | 2–10° | 0 | 325 | 8.6 | **yes** |
| V3 | 4–8° | 1 | 80 | 10.3 | **yes** |
| V3 | 2–10° | 0 | 202 | 10.1 | **yes** |
| V3a | 4–8° | 4 | 18 | 15.4 | no |
| V3a | 2–10° | 1 | 44 | 9.5 | **yes** |
| hV4 | 4–8° | 24 | 4 | 22.2 | no |
| hV4 | 2–10° | 10 | 34 | 8.8 | no |
| V3b | 4–8° | 18 | 10 | 170.2 | no |
| V3b | 2–10° | 10 | 25 | 176.9 | no |
| MT | 4–8° | 30 | 1 | ∞ | no |
| MT | 2–10° | 19 | 9 | 48.2 | no |
| MST | 4–8° | 28 | 2 | 146.8 | no |
| MST | 2–10° | 21 | 9 | 123.6 | no |

So **V1, V2 and V3 at 4–8°, plus V3a at 2–10°**. hV4, V3b, MT and MST do not qualify at either
band: MT has a median of *one* vertex per cell at 4–8°, and MT at 4–8° returns all-NaN because one
observer has no surviving vertices at all. Their weight ratios above 100× are the same failure seen
from the other side — with one or two vertices in a cell, σ̂ is not a measurement.

This is not a claim that those maps have no asymmetries. It is a claim that **this design cannot
resolve them by polar angle**, which is what `local_qc/group_addv4mt.m` already concluded when it
went whole-ROI for hV4 and MT. Whole-ROI analyses of those maps remain available and are unaffected.

Note that hV4 at 2–10° now falls below the line on the empty-cell count, where earlier passes
reported it with caveats. The criterion is stricter than the prose caveat was, and consistently
applied.

All sixteen runs are in `supplement/precision_weighted_areas.csv` with their coverage and a
`reportable` flag, so the excluded numbers are inspectable rather than absent.

## 7. Results

Equal-weighted, spec settings. ⚠ marks τ̂² pinned at zero.

**Context effects (dg − da), 4–8°.** Both Cartesian effects decline monotonically up the hierarchy
and remain significant throughout. V3a is shown at 2–10°, the only band where it qualifies, so it is
not on the same footing as the other three.

| asymmetry | V1 | V2 | V3 | V3a (2–10°) |
|---|---|---|---|---|
| horiz−vert | −0.316 [−0.520, −0.113] *p*=.008 | −0.260 [−0.369, −0.152] *p*=.001 | −0.102 [−0.186, −0.018] *p*=.024 | −0.084 *p*=.074 |
| card−obl | −0.181 [−0.311, −0.050] *p*=.014 | −0.161 [−0.256, −0.066] *p*=.005 | −0.115 [−0.220, −0.010] *p*=.036 | −0.152 [−0.259, −0.046] *p*=.012 |
| rad−tang | −0.036 *p*=.647 | +0.045 *p*=.304 | −0.060 *p*=.055 ⚠ | −0.154 [−0.207, −0.101] *p*<.001 ⚠ |
| polc−polo | +0.039 *p*=.354 | +0.046 *p*=.259 | +0.035 *p*=.380 | +0.037 *p*=.046 ⚠ |

**V3a, 2–10°, per experiment.** Cartesian asymmetries continue to attenuate (dg horiz−vert −0.099
*p*=.027, card−obl −0.150 *p*=.010); radial−tangential in the polar experiment stays large
(0.189 [0.116, 0.262], *p*<.001), consistent with V1–V3. Its two significant polar-frame *context*
effects both sit in τ̂² = 0 rows and should not be read as findings.

**Per-experiment asymmetries, 4–8°.**

| exp | asymmetry | V1 | V2 | V3 |
|---|---|---|---|---|
| dg | horiz−vert | −0.548 *p*<.001 | −0.343 *p*<.001 | −0.143 *p*=.004 |
| dg | card−obl | −0.221 *p*=.006 | −0.167 *p*=.007 | −0.113 *p*=.022 |
| dg | rad−tang | 0.119 *p*=.006 | 0.193 *p*<.001 | 0.074 *p*=.002 |
| dg | polc−polo | 0.072 *p*=.005 | 0.072 *p*=.019 | 0.027 *p*=.347 |
| da | horiz−vert | −0.232 *p*=.008 | −0.082 *p*=.006 | −0.041 *p*=.059 |
| da | card−obl | −0.041 *p*=.105 | −0.006 *p*=.438 ⚠ | 0.002 *p*=.901 |
| da | rad−tang | 0.155 *p*=.028 | 0.148 *p*<.001 ⚠ | 0.134 *p*<.001 ⚠ |
| da | polc−polo | 0.033 *p*=.352 | 0.026 *p*=.236 ⚠ | −0.007 *p*=.661 ⚠ |

The Cartesian asymmetries attenuate up the hierarchy in both experiments; radial−tangential in the
polar experiment stays large and significant in all three maps. **The context dependence is a
V1-weighted phenomenon that weakens but does not disappear through V2 and V3.**

### What adopting the spec cost in V1

Against the legacy numbers in `supplement/precision_weighted.csv`: max difference **0.0396**, and
**zero of twelve rows change significance**. The largest movers are dg polc−polo (0.039 → 0.072)
and the dg−da polc−polo context effect (−0.001 → +0.039) — both polc−polo, as expected, since b4
has the least redundancy under binning and so moves most when binning is dropped.

### The two routes agree

At matched θv the harmonic and ROI-space routes are the same estimator: **7.5e-16** (V1) and
**9.4e-16** (V2) where every cell is populated, diverging only in proportion to empty cells —
1.7e-3 in V3 (2 empty), 3.3e-2 in hV4 (20). Compared at their own defaults they differ by
0.04–0.05, and that gap is entirely the within-wedge local-orientation term, not a disagreement
about missing data.

### hV4, 2–10°

Reported for completeness and **not comparable to the maps above** — different eccentricity band, 20
empty cells, and a route disagreement of 0.053.

| exp | asymmetry | estimate |
|---|---|---|
| dg | horiz−vert | −0.158 [−0.285, −0.030] *p*=.022 |
| dg | card−obl | −0.094 [−0.201, 0.014] *p*=.079 |
| dg | rad−tang | 0.041 [−0.019, 0.102] *p*=.150 |
| dg | polc−polo | −0.002 [−0.048, 0.044] *p*=.911 |
| da | horiz−vert | 0.024 [−0.031, 0.079] *p*=.336 |
| da | card−obl | −0.004 [−0.081, 0.074] *p*=.916 |
| da | rad−tang | 0.233 [0.124, 0.341] *p*=.001 |
| da | polc−polo | 0.048 [−0.053, 0.150] *p*=.299 |
| dg−da | horiz−vert | −0.181 [−0.315, −0.048] *p*=.015 |
| dg−da | card−obl | −0.090 [−0.245, 0.065] *p*=.211 |
| dg−da | rad−tang | −0.192 [−0.340, −0.044] *p*=.018 |
| dg−da | polc−polo | −0.050 [−0.175, 0.074] *p*=.371 |

## 8. Cautions

**τ̂² is pinned at zero in six rows** — three in V2, three in V3, none in V1 or hV4. There the
across-observer spread is at or below the measured within-observer noise, so the interval is driven
entirely by measured σ rather than by observer variability. These are not ordinary random-effects
intervals.

**card−obl and polc−polo are not two findings in V3 or hV4.** See §5: they are one measurement per
cell, and class-structured ROI loss is what separates them.

**The hierarchy trend has not been tested.** The claim that the context effect declines
monotonically V1 → V2 → V3 rests on six individually significant cells falling in the same order in
both asymmetries. That is more than any single test, but the trend itself — a within-observer
V1 − V3 difference — has not been computed. This is the first thing to do next.

**Multiplicity.** The tables hold 48 tests. Nothing here is corrected for that; the Cartesian
context effects survive comfortably, the marginal cells do not.

## 9. Traps in this code

Two of these cost real time and would recur silently. Both are the reason an independent check
path is worth keeping.

- **Bootstrap draws must be generated one row at a time.** Pre-generating them as a matrix fills
  column-major and consumes the random stream in a different order than the original per-iteration
  call, **silently changing every bootstrap SE** — no error, no warning, different numbers. Draws do
  have to be generated up front, since the fitting path reseeds the global stream; generate them a
  row at a time.
- **Do not mask a full-surface vector with a mask defined over `vertIndex`.** An indexing bug in the
  gain summary did exactly that and read the wrong elements. It was caught only because an
  independent path (`dg_computeGain`) existed to check against.

And one reading trap: **read significance off the equal-weighted column, which is primary.**
Adopting continuous θv was once flagged as moving dg polc−polo across *p* = .05; that was read off
the precision-weighted column, and under equal weighting the status holds.

Three claims from earlier passes were wrong and are withdrawn; the correct positions are stated in
§§1–7 above. For the record: the harmonic route's inputs were said to be V1-only (they are not — the
cached per-vertex table covers all four maps and the fitting code is area-agnostic; the only real
gap was precision weighting); the harmonic route was said to be most exposed to the second-harmonic
degeneracy (backwards — the degeneracy is in the data, and the harmonic route makes it exact and
visible); and hV4's positive `da` horiz−vert was flagged as needing scrutiny (it got it — under the
spec it is 0.024 [−0.031, 0.079], *p* = .34, an ROI-binning artefact).
