# Results

Every number here is under the settled specification — [`SPECIFICATION.md`](SPECIFICATION.md).
Estimates are 2·b in **percent signal change**, equal-weighted across the **7** observers, with
***t* intervals on 6 df**. `dg` = Cartesian gratings, `da` = polar gratings.

**The observer set changed on 2026-08-27, and every number below moved with it.** sub-0395's `da`
session used a pilot stimulus whose annuli did not scale spatial period with eccentricity, so that
observer is excluded from `da` and from every `dg`-vs-`da` contrast — see `../AGENTS.md` standing
fact 7. The previous version of this document reported 8 observers. Two conclusions changed, in
opposite directions: the radial−tangential context effect became detectable (§4), and the
V1−V3 hierarchy difference stopped being (§5).

**Where each number comes from.** §2, §3, §5, §7 and §8 are regenerated wholesale by one command,
`run_spec_outputs`, and are transcribed from its CSVs. §4 is derived from one of those CSVs
(`spec_perobserver_spec_v1_4-8.csv`) by tests the driver does not itself run. §6 comes from the
harmonic-model diagnostics, which are run on demand and are **not** part of `run_spec_outputs` —
those numbers are on the model's own route, not the specification's, and the difference is stated
where it matters.

---

## 1. The claim, and where it stands

An orientation asymmetry is **strong when its reference frame matches the global stimulus, weak or
absent when it does not**.

- **Established.** The two *Cartesian-frame* asymmetries — horizontal vs. vertical, and cardinal vs.
  oblique — are substantially larger with Cartesian gratings than with polar gratings. This survives
  computing each vertex's exact local orientation, equalising response gain, tightening the pRF
  quality floor, changing the aggregation route, and changing the observer weighting.
- **Established for radial vs. tangential too** (changed 2026-08-27). It is larger with polar
  gratings than with Cartesian ones — context effect **−0.107** [−0.154, −0.061], *p* = .001, with
  **all 7** observers in the same direction. On the previous 8-observer set this was −0.036,
  *p* = .65, and was reported as uninformative; §4 says what changed and why the change is credible.
- **Still nothing for polar-cardinal vs. polar-oblique.** +0.006 [−0.053, 0.065], *p* = .81, 4 of 7
  observers. This one remains absence of evidence rather than evidence of absence, and the caution
  in §4 continues to apply to it.

---

## 2. V1, 4–8° — the four asymmetries in each experiment

`spec_asymmetries_spec_v1_4-8.csv`. These are Figures 5 and 6.

| exp | asymmetry | estimate | *t* 95% CI | *p* | observers agreeing |
|---|---|--:|---|--:|--:|
| dg | horiz−vert | **−0.514** | [−0.662, −0.366] | <.001 | 7/7 |
| dg | card−obl | **−0.193** | [−0.337, −0.050] | .016 | 6/7 |
| dg | rad−tang | **0.091** | [0.044, 0.138] | .003 | 7/7 |
| dg | polc−polo | **0.056** | [0.027, 0.085] | .003 | 7/7 |
| da | horiz−vert | **−0.212** | [−0.384, −0.041] | .023 | 6/7 |
| da | card−obl | −0.044 | [−0.103, 0.015] | .119 | 5/7 |
| da | rad−tang | **0.198** | [0.119, 0.277] | .001 | 7/7 |
| da | polc−polo | 0.050 | [−0.023, 0.124] | .146 | 4/7 |

Note the direction of the Cartesian-frame effects: BOLD is **lower** for horizontal than vertical,
and **lower** for cardinal than oblique — the "inverted" direction relative to behavioural
oblique effects.

Note also `da` rad−tang (**0.198**) against `dg` rad−tang (0.091): the radial asymmetry is more than
twice as large for polar gratings, and every observer shows it in both experiments.

## 3. V1, 4–8° — the four context effects

`spec_context_spec_v1_4-8.csv`. Each observer's `dg − da` difference is formed first, then tested
across the 7 observers — the correct analysis for a balanced within-subject design, and the one
reported everywhere ([`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) §3).

| asymmetry | estimate | *t* 95% CI | *p* | observers agreeing |
|---|--:|---|--:|--:|
| horiz−vert | **−0.302** | [−0.538, −0.065] | **.021** | 7/7 |
| card−obl | **−0.149** | [−0.281, −0.018] | **.032** | 6/7 |
| rad−tang | **−0.107** | [−0.154, −0.061] | **.001** | 7/7 |
| polc−polo | +0.006 | [−0.053, 0.065] | .815 | 4/7 |

## 4. The radial−tangential context effect, and the one that is still null

**Changed 2026-08-27.** This section used to be titled "Why the polar-frame null is not a result" and
argued that *both* polar-frame asymmetries were uninformative. That is no longer the position for
radial−tangential. The change is the observer set, not the method: sub-0395's `da` session used a
pilot stimulus (standing fact 7), and with that observer excluded the effect is
**−0.107 [−0.154, −0.061], *p* = .001, 7 of 7 observers**.

**The previous version of this document had already computed this number.** Its leave-one-out check
read: "Dropping sub-0395 — the only observer with a negative `da` radial−tangential value, −0.180 —
makes the difference significant (−0.110 [−0.157, −0.062], *p* = .001). No other leave-one-out comes
close: the next smallest *p* is .68." That is the result now obtained as the primary analysis; the
remaining difference (−0.110 against −0.107) is the group-gain constant, which is computed over
whichever observers are in the set.

**Why this is not a result found by dropping an inconvenient observer.** Three things have to line up,
and they do. The exclusion was decided on **stimulus grounds and independently of any result** — the
pilot annuli did not scale spatial period with eccentricity. The asymmetry that changed is exactly
the one that error would corrupt: radial versus tangential depends directly on the annulus geometry
that was wrong. And sub-0395 was not adding symmetric noise — its `da` rad−tang was **−0.180**, the
only negative value among the eight, pulling against every other observer.

**Still null: polar-cardinal versus polar-oblique.** +0.006 [−0.053, 0.065], *p* = .815, 4 of 7
observers. The cautions below were written about both polar-frame asymmetries and **continue to apply
to this one**:

- **The interval admits a real effect.** A polc−polo context effect of |0.065| stays compatible with
  the data, which is not negligible against the card−obl effect of −0.149 that *is* reported.
- **Absence of evidence is not evidence of absence.** Distinguishing a genuinely one-sided context
  effect from a two-sided one of unequal size needs more observers, and n = 7 is fewer than before.
- **`tau` estimates as exactly 0 for both polar-frame asymmetries at n = 7**, and the precision
  weight ratio reaches 15:1. That affects the precision-weighted *variant* only — equal weighting is
  primary and is what the tables above report — but a zero between-observer variance component is a
  small-sample artifact, not a finding, and should not be quoted as one.

**Consequence for the draft.** The abstract's "the radial asymmetry is 50% larger for polar gratings"
is now **supported, and understated**: `da` 0.198 against `dg` 0.091 is a factor of 2.2, with the
context effect itself significant at *p* = .001. What should *not* be carried over is the companion
claim that the polar-cardinal asymmetry is unchanged between experiments — that one is still
uninformative rather than shown to be equal.

*(Earlier history: the three checks in this section were re-derived under the specification on
2026-08-25, having been carried over from the pre-specification harmonic route. The argument was
first made on the harmonic model's own route in
[`supplement/SUPPLEMENT_harmonic_model.md`](supplement/SUPPLEMENT_harmonic_model.md) §S5.3, whose
numbers differ — see the note at the head of that document.)*

---

## 5. Across the visual hierarchy

### Which maps can be reported

Coverage measured, not assumed, on the [`SPECIFICATION.md`](SPECIFICATION.md) §7 criterion.
`spec_areas_coverage_spec.csv`.

| map | band | empty cells /64 | median vertices/cell | max **precision**-weight ratio | reportable |
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

**Why the excluded maps are excluded, measured 2026-08-25.** Reproducing a map's whole coverage
profile in V1 — its empty cells and its per-cell vertex counts together — gives a 90% band on `dg`
rad−tang of 0.018 under V3's coverage and **0.131 under MT's**, against an effect of 0.119. The
binding problem is how few vertices these maps have, not the shape of the gaps: V1 with MT's 30
cells deleted still holds 7045 vertices and the specification absorbs that loss to within 0.066,
while MT holds 545. → [`MISSING_DATA.md`](MISSING_DATA.md)

### Per-experiment asymmetries, 4–8°

⚠ marks a row where τ̂² is pinned at zero (see §7).

| exp | asymmetry | V1 | V2 | V3 |
|---|---|--:|--:|--:|
| dg | horiz−vert | −0.514 *p*<.001 | −0.330 *p*<.001 | −0.128 *p*=.010 |
| dg | card−obl | −0.193 *p*=.016 | −0.137 *p*=.012 | −0.084 *p*=.036 |
| dg | rad−tang | 0.091 *p*=.003 | 0.184 *p*=.002 | 0.068 *p*=.006 |
| dg | polc−polo | 0.056 *p*=.003 | 0.055 *p*=.027 | 0.007 *p*=.738 |
| da | horiz−vert | −0.212 *p*=.023 | −0.085 *p*=.011 | −0.037 *p*=.118 |
| da | card−obl | −0.044 *p*=.119 | −0.002 *p*=.765 ⚠ | −0.009 *p*=.445 |
| da | rad−tang | 0.198 *p*=.001 ⚠ | 0.146 *p*<.001 ⚠ | 0.132 *p*<.001 ⚠ |
| da | polc−polo | 0.050 *p*=.146 | 0.033 *p*=.173 ⚠ | −0.001 *p*=.959 ⚠ |

**The Cartesian asymmetries attenuate up the hierarchy in both experiments; radial−tangential in the
polar experiment stays large and significant in all three maps.**

### Context effects (dg − da), 4–8°

| asymmetry | V1 | V2 | V3 | V3a (2–10°) |
|---|---|---|---|---|
| horiz−vert | −0.302 [−0.538, −0.065] *p*=.021 | −0.245 [−0.369, −0.122] *p*=.003 | −0.090 [−0.185, 0.004] *p*=.058 | −0.077 *p*=.126 |
| card−obl | −0.149 [−0.281, −0.018] *p*=.032 | −0.135 [−0.224, −0.046] *p*=.010 | −0.075 [−0.139, −0.012] *p*=.028 | −0.132 [−0.248, −0.016] *p*=.031 |
| rad−tang | **−0.107 [−0.154, −0.061] *p*=.001 ⚠** | +0.038 [−0.074, 0.151] *p*=.434 | −0.064 [−0.134, 0.006] *p*=.067 ⚠ | −0.156 [−0.208, −0.105] *p*<.001 ⚠ |
| polc−polo | +0.006 [−0.053, 0.065] *p*=.815 ⚠ | +0.021 [−0.058, 0.101] *p*=.534 | +0.008 [−0.063, 0.080] *p*=.791 | +0.041 [0.004, 0.077] *p*=.035 ⚠ |

V3a is shown at 2–10°, the only band where it qualifies, so it is not on the same footing as the
other three. Its two significant polar-frame *context* effects both sit in τ̂² = 0 rows and should not
be read as findings.

**A caution that now applies to the V1 radial−tangential result too.** That cell is significant
(*p* = .001) but its row is also τ̂² = 0 — the between-observer variance component is pinned at zero at
n = 7, which narrows the interval. Equal weighting, which is primary, does not use τ̂², so the estimate
and *t* interval in the table stand; but the τ̂² = 0 flag is why the precision-weighted variant gives
*p* = .024 rather than .001 for the same cell, and that gap is a small-sample property rather than a
disagreement about the effect. Read the V1 rad−tang context effect as established in **direction and
sign** (7/7 observers) and as **less precisely bounded** than its interval alone suggests.

**Note the V2 reversal.** rad−tang context at V2 is +0.038 (*p* = .43) — the opposite sign to V1 and
V3, and not significant in either direction. The polar-frame context effect is not a smooth gradient
up the hierarchy the way the Cartesian ones are.

**The context dependence is a V1-weighted phenomenon that weakens but does not disappear through V2
and V3.**

### The hierarchy trend, tested

The monotonic decline was originally read off six individually significant cells falling in order.
It was **computed** on 2026-08-24, *within observer* — for each observer, the difference of the
context effect between two maps, then a *t* test across observers. `spec_areas_trend_spec.csv`.

**Re-run on the corrected 7 observers, 2026-08-27. The V1 − V3 result no longer reaches
significance**, and this is the one conclusion that the observer correction weakened rather than
strengthened.

| comparison | asymmetry | mean | *t* 95% CI | *p* | obs |
|---|---|--:|---|--:|--:|
| V1 − V2 | rad−tang | −0.146 | [−0.264, −0.028] | .023 | 6/7 |
| V2 − V3 | horiz−vert | **−0.155** | [−0.269, −0.041] | **.016** | 7/7 |
| V2 − V3 | card−obl | −0.059 | [−0.098, −0.021] | .009 | 6/7 |
| V2 − V3 | rad−tang | 0.102 | [−0.011, 0.216] | .070 | 6/7 |
| **V1 − V3** | **horiz−vert** | −0.211 | [−0.426, 0.003] | **.053** | 6/7 |
| **V1 − V3** | card−obl | −0.074 | [−0.184, 0.036] | .150 | 5/7 |
| V1 − V3 | rad−tang | −0.044 | [−0.109, 0.022] | .155 | 6/7 |

(polc−polo rows are all *p* > .9.)

**What this establishes.** The **V2 → V3** step is solid for horizontal−vertical (−0.155, *p* = .016,
7/7 observers) and for cardinal−oblique (−0.059, *p* = .009).

**What it does not.** The **V1 − V3** difference for horizontal−vertical is now *p* = **.053** and its
interval includes zero by a hair — on the previous 8-observer set it was −0.215, *p* = .028. The
estimate barely moved (−0.215 → −0.211); what changed is the degrees of freedom, from 7 to 6. So the
claim "the horizontal−vertical context effect attenuates from V1 to V3" is now **suggestive rather
than established**, and should not be stated as demonstrated. Nothing about V1 − V3 for
cardinal−oblique changed: it was not significant before and is not now.

**Do not read the *p* = .053 as "nearly significant" and report it as a decline.** The honest summary
is that a real V1→V3 attenuation is likely — the individual cells do fall in order in each map, and
the V2 → V3 step is resolved — but the V1-to-V3 *difference* is not resolved by these data at n = 7.
Recovering it is one of the concrete gains that the 13-observer `dg` set could deliver, since the
hierarchy comparison is within-experiment and does not need `da`.

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
  the *t* reading**: at n = 7 the percentile bootstrap has poor coverage and does not account for
  uncertainty in the spread. Both columns are in `spec_asymmetries_spec_v1_4-8.csv` with a
  `ci_methods_disagree` flag, so the choice stays visible rather than being made silently.
  (`da` rad−tang used to disagree too; under the specification both intervals now exclude zero.)
- **This is not a claim that hV4, V3b, MT and MST have no asymmetries** — only that this design
  cannot resolve them by polar angle. Whole-ROI analyses of those maps are unaffected.

---

## 8. Robustness: neither the route nor the weighting changes a conclusion

Same vertices, same gain, same observers; only the named thing differs.

**The route** (`spec` → `roi`, i.e. continuous θ_V → wedge-centre θ_V, both equal-weighted). Largest
movement in V1 is 0.031 (`dg` polc−polo); `dg` horiz−vert and card−obl are identical to three decimals
because their regressors are θ_V-free, so the route cannot touch them. Context effects:
horiz−vert −0.302 → −0.321, card−obl −0.149 → −0.165, rad−tang −0.107 → −0.114, polc−polo
+0.006 → −0.035. **No cell changes significance.**

**The weighting** (`roi` → `roipw`, equal → precision). Largest change 0.017, on `da` rad−tang
(0.195 → 0.212). **No cell changes significance.**

**The radial−tangential context effect survives both**, which is the check that matters most for the
conclusion that changed on 2026-08-27: −0.107 *p*=.001 (`spec`), −0.114 *p*=.007 (`roi`),
−0.129 *p*=.021 (`roipw`). Significant on all three, same sign, same order of magnitude. The *p*
weakens under precision weighting for the reason given in §5 — τ̂² = 0 at n = 7 — but the effect does
not depend on the choice.

**The hierarchy trend is now the one place a variant crosses *p* = .05** (changed 2026-08-27, with
the observer correction). V1 − V3 horiz−vert is −0.211 *p*=**.053** under `spec` and −0.231
*p*=**.044** under `roi`; card−obl is −0.074 *p*=.150 vs −0.090 *p*=.072. The estimates agree to
within 0.02 and the conclusion should not be taken from whichever side of .05 a given route lands on
— **the honest reading is that this comparison is underpowered at n = 7**, not that the route
matters. Previously, at n = 8, both routes agreed the V1 − V3 horiz−vert trend was significant.

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

**`Figure_S5_spec_hierarchy`, the extrastriate supplement figure.** Same four asymmetries, same four
panel order, across the maps §5 says can be resolved by polar angle. Top row is each asymmetry in
each experiment; bottom row is the context effect (`dg − da`) over the same maps. Filled markers are
V1, V2 and V3 at 4–8°; the open marker is V3a, which qualifies only at 2–10° and so is not on the
same footing. **Scales are shared within each row and deliberately not between them** — within a row
the panels are the same quantity, across rows they are an asymmetry and a difference of asymmetries,
and one scale there would assert a comparison nobody is making. Maps that fail the criterion are
absent by design, not oversight; their numbers stay in `spec_areas_*.csv` behind the `reportable`
flag, which is what `plot_spec_hierarchy.m` filters on.

**Figure 7 (joint LME) is being removed** — it is exactly redundant with Figures 5/6
([`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) §3).

**Figure 8 (per-location context effects)** is untouched by this work. One small polar plot per
visual-field location, black = observer data, red = LME model estimate; it highlights that the
largest context effect is on the horizontal meridian, where suppressed horizontal and enhanced radial
oppose each other. Its red curve is the Fig-7 LME, so whether it keeps the overlay is open
([`../AGENTS.md`](../AGENTS.md) §5).

---

## 10. Scope — what these results deliberately do not cover

Not a TODO list; open items are in [`../AGENTS.md`](../AGENTS.md) §5.

- **The 2–10° band is computed but not used for the main figures.** Figures 5/6 are V1 4–8°, the band
  where spatial frequency is matched between the experiments. The 2–10° numbers are in the CSVs, and
  V3a is reportable only there.
- **Whole-ROI analyses of hV4 and MT** were not rerun; the coverage criterion is about polar-angle
  resolution only.
