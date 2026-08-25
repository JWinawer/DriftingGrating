# Results

Every number here is produced by one command (`run_spec_outputs`) under the settled specification —
[`SPECIFICATION.md`](SPECIFICATION.md). Estimates are 2·b in **percent signal change**, equal-weighted
across the 8 observers, with ***t* intervals on 7 df**. `dg` = Cartesian gratings, `da` = polar
gratings.

---

## 1. The claim, and where it stands

An orientation asymmetry is **strong when its reference frame matches the global stimulus, weak or
absent when it does not**.

- **Established.** The two *Cartesian-frame* asymmetries — horizontal vs. vertical, and cardinal vs.
  oblique — are substantially larger with Cartesian gratings than with polar gratings. This survives
  computing each vertex's exact local orientation, equalising response gain, tightening the pRF
  quality floor, changing the aggregation route, and changing the observer weighting.
- **Uninformative, not absent.** The two *polar-frame* asymmetries — radial vs. tangential, and
  polar-cardinal vs. polar-oblique — show **no detectable difference** between experiments. This is
  absence of evidence, not evidence of absence: the interval admits a polar-frame context effect
  larger than the cardinal/oblique one that *is* significant. §4 gives the three reasons.

---

## 2. V1, 4–8° — the four asymmetries in each experiment

`spec_asymmetries_spec_v1_4-8.csv`. These are Figures 5 and 6.

| exp | asymmetry | estimate | *t* 95% CI | *p* | observers agreeing |
|---|---|--:|---|--:|--:|
| dg | horiz−vert | **−0.548** | [−0.685, −0.411] | <.001 | 8/8 |
| dg | card−obl | **−0.221** | [−0.356, −0.086] | .006 | 7/8 |
| dg | rad−tang | **0.119** | [0.046, 0.192] | .006 | 8/8 |
| dg | polc−polo | **0.072** | [0.030, 0.113] | .005 | 8/8 |
| da | horiz−vert | **−0.232** | [−0.382, −0.081] | .008 | 7/8 |
| da | card−obl | −0.041 | [−0.092, 0.011] | .105 | 6/8 |
| da | rad−tang | **0.155** | [0.023, 0.287] | .028 | 7/8 |
| da | polc−polo | 0.033 | [−0.045, 0.110] | .352 | 4/8 |

Note the direction of the Cartesian-frame effects: BOLD is **lower** for horizontal than vertical,
and **lower** for cardinal than oblique — the "inverted" direction relative to behavioural
oblique effects.

## 3. V1, 4–8° — the four context effects

`spec_context_spec_v1_4-8.csv`. Each observer's `dg − da` difference is formed first, then tested
across the 8 observers — the correct analysis for a balanced within-subject design, and the one
reported everywhere ([`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) §3).

| asymmetry | estimate | *t* 95% CI | *p* | observers agreeing |
|---|--:|---|--:|--:|
| horiz−vert | **−0.316** | [−0.520, −0.113] | **.008** | 8/8 |
| card−obl | **−0.181** | [−0.311, −0.050] | **.014** | 7/8 |
| rad−tang | −0.036 | [−0.215, 0.143] | .647 | 7/8 |
| polc−polo | +0.039 | [−0.054, 0.132] | .354 | 5/8 |

## 4. Why the polar-frame null is not a result

Three checks, all on the harmonic model (`cleanroom/diagnose_context_asymmetry.m`; full account in
[`supplement/SUPPLEMENT_harmonic_model.md`](supplement/SUPPLEMENT_harmonic_model.md) §S5.3):

- **The interval is wide.** A radial/tangential context effect as large as |0.158| remains
  compatible with the data — larger than the cardinal/oblique context effect that *is* reported as
  significant, and 59% of the horizontal/vertical one.
- **It rests on one observer.** Six of eight show `da` rad−tang exceeding `dg` (sign test *p* = 0.29).
  Dropping sub-0395 — the only observer with a negative `da` radial−tangential value, −0.244 — makes
  the difference significant (−0.113 [−0.182, −0.036]). No other leave-one-out does.
- **The frames are not reliably different from each other.** A within-subject difference of
  differences comparing the two frames' context effects is not significant: horiz−vert vs rad−tang
  0.106 [−0.033, 0.275], *p* = 0.26; card−obl vs rad−tang −0.004, *p* = 0.92.

**Consequence for the draft.** The abstract's "the radial asymmetry is 50% larger for polar gratings"
should not be quoted as a quantity. The direction holds (`da` 0.155 vs `dg` 0.119), but the
difference is not detectable (*p* = .65). Replace it with the statement about *evidence*: the
Cartesian-frame asymmetries show robust context dependence, while for the polar-frame ones these
data are uninformative. Distinguishing a genuinely one-sided context effect from a two-sided one of
unequal size needs more observers.

---

## 5. Across the visual hierarchy

### Which maps can be reported

Coverage measured, not assumed, on the [`SPECIFICATION.md`](SPECIFICATION.md) §7 criterion.
`spec_areas_coverage_spec.csv`.

| map | band | empty cells /64 | median vertices/cell | max weight ratio | reportable |
|---|---|--:|--:|--:|---|
| V1 | 4–8° | 0 | 150 | 17.1 | **yes** |
| V1 | 2–10° | 0 | 326 | 18.1 | **yes** |
| V2 | 4–8° | 0 | 136 | 15.8 | **yes** |
| V2 | 2–10° | 0 | 325 | 8.6 | **yes** |
| V3 | 4–8° | 1 | 80 | 10.3 | **yes** |
| V3 | 2–10° | 0 | 202 | 10.1 | **yes** |
| V3a | 4–8° | 4 | 18 | 15.4 | no |
| V3a | 2–10° | 1 | 44 | 9.5 | **yes** |
| V3b | 4–8° | 18 | 10 | 170.2 | no |
| V3b | 2–10° | 10 | 25 | 176.9 | no |
| hV4 | 4–8° | 24 | 4 | 22.2 | no |
| hV4 | 2–10° | 10 | 34 | 8.8 | no |
| MT | 4–8° | 30 | 1 | ∞ | no |
| MT | 2–10° | 19 | 9 | 48.2 | no |
| MST | 4–8° | 28 | 2 | 146.8 | no |
| MST | 2–10° | 21 | 9 | 123.6 | no |

So **V1, V2 and V3 at 4–8°, plus V3a at 2–10°**. MT at 4–8° has a median of *one* vertex per cell and
one observer with no surviving vertices at all, so its estimates are NaN by construction. Excluded
numbers are still in `spec_areas_*.csv` behind a `reportable` flag, so they are inspectable rather
than absent.

### Per-experiment asymmetries, 4–8°

⚠ marks a row where τ̂² is pinned at zero (see §7).

| exp | asymmetry | V1 | V2 | V3 |
|---|---|--:|--:|--:|
| dg | horiz−vert | −0.548 *p*<.001 | −0.343 *p*<.001 | −0.143 *p*=.004 |
| dg | card−obl | −0.221 *p*=.006 | −0.167 *p*=.007 | −0.113 *p*=.022 |
| dg | rad−tang | 0.119 *p*=.006 | 0.193 *p*<.001 | 0.074 *p*=.002 |
| dg | polc−polo | 0.072 *p*=.005 | 0.072 *p*=.019 | 0.027 *p*=.347 |
| da | horiz−vert | −0.232 *p*=.008 | −0.082 *p*=.006 | −0.041 *p*=.059 |
| da | card−obl | −0.041 *p*=.105 | −0.006 *p*=.438 ⚠ | 0.002 *p*=.901 |
| da | rad−tang | 0.155 *p*=.028 | 0.148 *p*<.001 ⚠ | 0.134 *p*<.001 ⚠ |
| da | polc−polo | 0.033 *p*=.352 | 0.026 *p*=.236 ⚠ | −0.007 *p*=.661 ⚠ |

**The Cartesian asymmetries attenuate up the hierarchy in both experiments; radial−tangential in the
polar experiment stays large and significant in all three maps.**

### Context effects (dg − da), 4–8°

| asymmetry | V1 | V2 | V3 | V3a (2–10°) |
|---|---|---|---|---|
| horiz−vert | −0.316 [−0.520, −0.113] *p*=.008 | −0.260 [−0.369, −0.152] *p*=.001 | −0.102 [−0.186, −0.018] *p*=.024 | −0.084 *p*=.074 |
| card−obl | −0.181 [−0.311, −0.050] *p*=.014 | −0.161 [−0.256, −0.066] *p*=.005 | −0.115 [−0.220, −0.010] *p*=.036 | −0.152 [−0.259, −0.046] *p*=.012 |
| rad−tang | −0.036 *p*=.647 | +0.045 *p*=.304 | −0.060 *p*=.055 ⚠ | −0.154 [−0.207, −0.101] *p*<.001 ⚠ |
| polc−polo | +0.039 *p*=.354 | +0.046 *p*=.259 | +0.035 *p*=.380 | +0.037 *p*=.046 ⚠ |

V3a is shown at 2–10°, the only band where it qualifies, so it is not on the same footing as the
other three. Its two significant polar-frame *context* effects both sit in τ̂² = 0 rows and should not
be read as findings. Its per-experiment values continue the pattern: `dg` horiz−vert −0.099
*p*=.027, `dg` card−obl −0.150 *p*=.010, `da` rad−tang 0.189 [0.116, 0.262] *p*<.001.

**The context dependence is a V1-weighted phenomenon that weakens but does not disappear through V2
and V3.**

### The hierarchy trend, tested

The monotonic decline was originally read off six individually significant cells falling in order.
It was **computed** on 2026-08-24, *within observer* — for each observer, the difference of the
context effect between two maps, then a *t* test across the 8 — and it holds for one asymmetry, not
both. `spec_areas_trend_spec.csv`.

| comparison | asymmetry | mean | *t* 95% CI | *p* | obs |
|---|---|--:|---|--:|--:|
| V1 − V2 | horiz−vert | −0.056 | [−0.211, 0.099] | .422 | 5/8 |
| V1 − V2 | card−obl | −0.019 | [−0.092, 0.053] | .547 | 5/8 |
| V2 − V3 | horiz−vert | **−0.159** | [−0.256, −0.061] | **.006** | 8/8 |
| V2 − V3 | card−obl | −0.046 | [−0.092, −0.001] | .047 | 6/8 |
| V2 − V3 | rad−tang | 0.105 | [0.008, 0.202] | .037 | 7/8 |
| **V1 − V3** | **horiz−vert** | **−0.215** | [−0.398, −0.031] | **.028** | 7/8 |
| **V1 − V3** | **card−obl** | −0.066 | [−0.162, 0.031] | .153 | 5/8 |

(The four polc−polo rows and the two remaining rad−tang rows are all *p* > .3.)

**What this establishes.** The context effect on **horizontal−vertical** genuinely attenuates from V1
to V3 (−0.215, *p* = .028, 7/8 observers), with the V2 → V3 step carrying most of it (*p* = .006,
8/8).

**What it does not.** For **cardinal−oblique** the V1 − V3 trend is not significant (−0.066,
*p* = .153); only the V2 − V3 step reaches *p* = .047, one of twelve uncorrected tests. So "both
Cartesian effects decline monotonically up the hierarchy" is established for **horizontal−vertical
only**. The individual cells are each significant in each map and do fall in order; the *difference
between maps* is resolved for one asymmetry.

The trend is **always equal-weighted**, whatever variant is requested, so the `roi` and `roipw` trend
tables are byte-identical (asserted, not assumed). Precision weighting it would need the
within-observer covariance of one observer's estimate in two maps, and those come from the same runs,
so their errors are correlated by an unmeasured amount; assuming independence would understate the SE.

---

## 6. Geometry does not explain the context effect

The alternative to a context effect is that the ROI binning mislabels local orientation: a horizontal
grating counts as "radial" for the whole 45°-wide wedge, though it is exactly radial only on the
meridian. The harmonic model tests this by discarding the binning and using each vertex's own pRF
polar angle.

**Within-ROI geometry accounts for 6–8% of the Cartesian-versus-polar gap** (5.6% for
horizontal−vertical, 7.6% for cardinal−oblique). The remaining ~93% survives, and survives equalising
overall response gain and tightening the pRF quality floor to R² > 0.3.

Two supporting results:

- **pRF polar-angle error is an order of magnitude too small to matter.** Angle error would deflate
  the polar experiment's Cartesian-frame coefficients but not the Cartesian experiment's, biasing
  *toward* the conclusion; reproducing the observed ratio in a context-free world would need
  σ ≈ 39°. Measured from two independent pRF solutions per observer: **σ = 3.9°**, giving
  second-harmonic attenuation under 1%.
- **Overall gain explains part of the difference but not the reference-frame structure.** Polar
  gratings elicit roughly half the response amplitude. Of the improvement from letting the two
  experiments have separate coefficients, 61% is that single scale factor and 39% is a genuine
  reference-frame difference. Equalising gain makes the Cartesian-frame context effects *larger*, not
  smaller.

Full account, with figures and the model-adequacy tests:
[`supplement/SUPPLEMENT_harmonic_model.md`](supplement/SUPPLEMENT_harmonic_model.md).

---

## 7. Cautions

- **τ̂² is pinned at zero in six rows** — three in V2, three in V3, none in V1. There the
  across-observer spread is at or below the measured within-observer noise, so the interval is driven
  entirely by measured σ rather than by observer variability. These are not ordinary random-effects
  intervals.
- **card−obl and polc−polo are not two findings where ROI coverage is class-structured.** Per cell
  they are one measurement with an alternating sign ([`SPECIFICATION.md`](SPECIFICATION.md) §6);
  only how cells combine across ROIs separates them, and that is exactly what ROI loss damages.
- **Multiplicity.** The tables hold 48 asymmetry/context tests plus 12 trend tests, uncorrected. The
  Cartesian context effects survive comfortably; the marginal cells do not.
- **One interval-method disagreement remains.** For `da` card−obl the percentile bootstrap
  [−0.083, −0.002] excludes zero and the *t* interval [−0.092, 0.011] does not (*p* = .105). **Report
  the *t* reading**: at n = 8 the percentile bootstrap has poor coverage and does not account for
  uncertainty in the spread. Both columns are in `spec_asymmetries_spec_v1_4-8.csv` with a
  `ci_methods_disagree` flag, so the choice stays visible rather than being made silently.
  (`da` rad−tang used to disagree too; under the specification both intervals now exclude zero.)
- **This is not a claim that hV4, V3b, MT and MST have no asymmetries** — only that this design
  cannot resolve them by polar angle. Whole-ROI analyses of those maps are unaffected.

---

## 8. Robustness: neither the route nor the weighting changes a conclusion

Same vertices, same gain, same observers; only the named thing differs.

**The route** (`spec` → `roi`, i.e. continuous θ_V → wedge-centre θ_V, both equal-weighted). Largest
movement in V1 is 0.033 (`dg` polc−polo, 0.072 → 0.039); `dg` horiz−vert and card−obl are identical
to three decimals because their regressors are θ_V-free, so the route cannot touch them. Context
effects: horiz−vert −0.316 → −0.330, card−obl −0.181 → −0.192, rad−tang −0.036 → −0.046, polc−polo
+0.039 → −0.001. **No cell changes significance.**

**The weighting** (`roi` → `roipw`, equal → precision). Largest change 0.025, on `da` rad−tang
(0.150 → 0.175). Every context effect unchanged to within 0.015. **No cell changes significance.**
That is what [`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) §4 predicts: τ² is common to every observer
and here exceeds the mean σᵢ², so a 14× spread in *reliability* compresses into a ~3× spread in
*weight*.

**The hierarchy trend is unaffected by either**: V1 − V3 horiz−vert −0.215 *p*=.028 (`spec`) vs
−0.221 *p*=.030 (`roi`); card−obl −0.066 *p*=.153 vs −0.084 *p*=.060.

**The two routes are the same estimator at matched θ_V**, and they agree to **7.5e-16** (V1) and
**9.4e-16** (V2) where every cell is populated, diverging only in proportion to empty cells —
1.7e-3 in V3 (2 empty), 3.3e-2 in hV4 (20). At their own defaults they differ by 0.04–0.05, and that
gap is entirely the within-wedge local-orientation term, not a disagreement about missing data. The
specification is preferable because the ROI route degrades where cells go empty, which is the
extrastriate case.

---

## 9. What the figures show

**Figures 5 (`dg`) and 6 (`da`).** Panels A–D are the four asymmetries in the standard order:
**A** horizontal vs. vertical, **B** cardinal vs. oblique, **C** radial vs. tangential,
**D** polar-cardinal vs. polar-oblique.

- **Row 1, polar plots.** Lines are the fitted harmonic model; markers are the observed wedge means
  over the same vertices. Where they separate — visibly in the radial/tangential and polar-cardinal
  panels — the gap *is* the within-wedge local-orientation term that binning θ_V conflates with
  context. It is shown rather than argued.
- **The radius is demeaned**, so it is a deviation from each vertex's mean over the four
  orientations, not a raw percent signal change. That is deliberate: demeaning removes the
  pink-noise blank, so only orientation differences enter.
- **Row 2 plots each observer's difference**, with the group mean and its interval — not the old
  pro-vs-con paired plot. Under this parameterisation the two classes are demeaned and partition the
  four orientations, so `con = −pro` exactly for the second-harmonic asymmetries and a paired plot
  would show the same information mirrored. The difference is the estimand, so the difference is
  what is shown.

**Axis scales are shared across both figures**, not per panel: every polar panel at ±0.62, every
difference panel at [−0.79, 0.38], all profile panels at [−1.10, 0.58]. `dg`-versus-`da` is the
paper's claim, so drawing the two on different scales would undercut the one comparison a reader most
needs to make by eye. A consequence worth stating: on the shared scale the two second-harmonic panels
are small. **That is the finding, not a drawing problem** — `card−obl` and `polc−polo` *are* about a
quarter the size of the first-harmonic contrasts, and per-panel autoscaling concealed it. One
definition, in `cleanroom/spec_axis_limits.m`.

**Why the polar plots stay at 45°.** The pro/con classes only exist at multiples of 45° for half the
panels. For `dg`, *horizontal* and *cardinal* are the same classes at every polar angle, so A/B could
bin arbitrarily finely — but *radial* means offset 0, and at θ_V = 15° the four offsets are
165°/75°/30°/120°, so no stimulus is radial and C/D have nothing to average. For `da` it is the
mirror image. Binning finely would give 24 markers in four panels and 8 in the other four.

**`Figure_5_6_spec_profile`** is the finely-binned version of exactly the panels that permit it,
which is why it has three columns rather than four: four orientations at 45° spacing give each
vertex's demeaned response exactly three degrees of freedom, so a profile against continuous θ_V can
show three curves and no more. The fourth coefficient is identified *across* vertices, from the θ_V
modulation of the first and third — which is what the tilt of these curves is. In the `roi` variant
the curve becomes a step function, showing what wedge-binning assumes about polar-angle structure.

**Figure 7 (joint LME) is being removed** — it is exactly redundant with Figures 5/6
([`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) §3).

**Figure 8 (per-location context effects)** is untouched by this work. One small polar plot per
visual-field location, black = observer data, red = LME model estimate; it highlights that the
largest context effect is on the horizontal meridian, where suppressed horizontal and enhanced radial
oppose each other. Its red curve is the Fig-7 LME, so whether it keeps the overlay is open
([`../AGENTS.md`](../AGENTS.md) §5).

---

## 10. Not done

- **The 2–10° band is computed but not used for the main figures.** Figures 5/6 are V1 4–8°, the band
  where spatial frequency is matched between the experiments. The 2–10° numbers are in the CSVs, and
  V3a is reportable only there.
- **Whole-ROI analyses of hV4 and MT** were not rerun; the coverage criterion is about polar-angle
  resolution only.
