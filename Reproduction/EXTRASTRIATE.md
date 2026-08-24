# Orientation asymmetries beyond V1

**Status: revised 2026-08-24.** The recommended route changed after the first pass. It is now the
**harmonic route**, applied identically to every map including V1, with the ROI-space de-trending
route kept as a cross-check. The two agree exactly wherever coverage is adequate and disagree only
in hV4, which is itself the finding about hV4. Section 8 records what the first pass got wrong.

Code: [`cleanroom/diagnose_within_observer_error.m`](cleanroom/diagnose_within_observer_error.m)
(`'route'`), [`cleanroom/precision_weighted_table.m`](cleanroom/precision_weighted_table.m),
[`cleanroom/precision_weighted_cells.m`](cleanroom/precision_weighted_cells.m),
[`cleanroom/fit_cell_meta.m`](cleanroom/fit_cell_meta.m),
[`server_extract/collect_runwise_betas_areas.m`](server_extract/collect_runwise_betas_areas.m).
Numbers: [`supplement/precision_weighted_areas.csv`](supplement/precision_weighted_areas.csv), both
routes, 96 rows.

---

## 1. Why the V1 procedure does not transfer

The V1 route takes each observer's asymmetry as the mean over their eight polar-angle ROIs, then
combines observers by precision weighting ([`LME.md`](LME.md) §5). That mean is only meaningful if
every observer contributes the same eight ROIs. In V1 they do — the sparsest cell still holds 22
vertices. Outside V1 they do not.

Two facts make the resulting gaps dangerous rather than merely inconvenient:

- **The asymmetry varies strongly with polar angle.** For dg polar-cardinal−polar-oblique the fitted
  profile alternates ±0.13–0.26 with cardinal/oblique ROI, against a mean effect of 0.039.
- **The loss is not random with respect to polar angle.** Map boundaries lie on the meridians, and
  vertices near a boundary have the poorest pRF fits, so cells empty out along exactly the axis the
  asymmetry varies along.

Averaging "whichever ROIs survived" therefore estimates a different quantity for each observer.
Deleting meridian ROIs from V1 in simulation biases the V1 route by **+0.071** on dg polc−polo — an
effect whose true value is 0.040. Cell loss nearly triples it.

Note what is *not* broken. The inclusion mask is orientation-independent
([`bin_and_aggregate.m:16`](cleanroom/bin_and_aggregate.m)), so all four orientations of a cell
vanish together and the within-ROI subtraction is never partially missing. Only the average over
ROIs breaks.

## 2. The procedure: one model for every map

The fix is to stop binning. Polar angle enters as a continuous covariate in a four-term harmonic
model fitted per vertex, so there are no cells to go empty; a coverage gap costs precision instead
of changing what is being estimated.

**Step 1 — extract run-wise betas for eight areas.**
[`collect_runwise_betas_areas.m`](server_extract/collect_runwise_betas_areas.m) writes the mean
GLMsingle single-trial beta per (vertex, condition, run) for V1, V2, V3, V3a, V3b, hV4, MT and MST.
All eight in one pass: the cost is streaming `modelmd` (425 MB per subject-experiment, 6.8 GB over
16 files) and it is identical for one area or eight, because the whole variable must cross the wire
before anything can be subset.

**Step 2 — fit the harmonic model, inside a bootstrap over runs.** For each observer, experiment and
map, resample the 8 runs 500×; on each draw, average the resampled runs, demean each vertex across
the four orientations, and fit

```
y_vk = b1*cos(2*theta) + b2*cos(4*theta) + b3*cos(2*(theta-thetaV)) + b4*cos(4*(theta-thetaV))
```

by weighted least squares, with `harmonic_weights`' equal-coverage weighting. Runs, not vertices:
an earlier vertex-resampling estimate was withdrawn because it holds the GLM betas fixed and
characterises which patch of cortex was sampled, not the reliability of the measurement.

**Step 3 — read off the four asymmetries.** Evaluate the fitted coefficients at the eight canonical
ROI centres and push that through the unmodified `compute_asymmetries`, so the output is in exactly
the units of the published numbers. The spread over the 500 draws is that observer's σ.

**Step 4 — combine observers exactly as before.** Feed the per-observer values and their σ into the
unchanged precision-weighted estimator from `LME.md` §5.

### θv is binned by default, and that is what makes V1 safe

Quantising the *regressor* to the eight ROI centres makes the fit algebraically the published ROI
analysis. It does **not** reintroduce the empty-cell problem: binning a predictor creates no cell
that can go empty — the fit simply uses whatever vertices exist. Consequently

> **V1 does not move at all.** Point estimates and bootstrap SEs agree with the ROI route to
> **1.5e-14** and **1.1e-14**, and the V1 table still reproduces
> [`supplement/precision_weighted.csv`](supplement/precision_weighted.csv) to **5.0e-16**.

So the same model can be used throughout, main text included, at zero cost to the existing numbers.
`'thetaV','continuous'` — each vertex's own pRF angle — is the scientifically preferable model and is
available, but it would move V1 by the within-wedge local-orientation term quantified in
[`HARMONIC_MODEL.md`](HARMONIC_MODEL.md), so it is offered for comparison rather than defaulted to.

## 3. The ROI-space cross-check

[`precision_weighted_cells.m`](cleanroom/precision_weighted_cells.m) implements an independent
route that stays in ROI space: fit the polar-angle profile `alpha_p` from the incomplete design with
[`fit_cell_meta.m`](cleanroom/fit_cell_meta.m) — a multilevel random-effects meta-analysis with the
**known** per-cell sampling covariance from the run bootstrap — then take each observer's number as
the mean over their *present* ROIs of `d_ip - alpha_p`. Because the `alpha_p` sum to zero, complete
data returns the existing numbers exactly.

Two things learned building it are worth keeping:

- **The cell-level GLS estimate `mu` should not be reported.** On complete V1 data it moves by up to
  0.05 depending on how the within-observer covariance is modelled, because GLS re-weights ROIs
  *within* an observer using a covariance estimated from 8 runs, and sampling noise is correlated
  with effect size (r up to 0.52 after removing observer and ROI means). Five of eight estimates move
  toward zero and three away, so it is not a correctable attenuation. This is `LME.md` §5's own
  caveat about estimated-weight GLS at n = 8, in practice.
- **The within-observer sampling covariance is structured, not noise.** Mean between-ROI correlation
  is **+0.6 when the asymmetry's reference frame matches the experiment** and **−0.1 when it does
  not**. A matched frame uses the same two stimulus conditions in every ROI, so a run-level
  fluctuation moves all eight together; an unmatched frame rotates which stimulus counts as "radial"
  with polar angle, so the same fluctuation pushes some ROIs up and others down.

## 4. Validation

| check | result |
|---|---|
| Harmonic and ROI routes agree on V1 and V2 | point estimates **1.5e-14**, bootstrap SEs **1.1e-14** |
| V1 table still reproduces the stored numbers | **5.0e-16** |
| Multi-area extraction preserves the V1-only extraction | V1 subset **bit-identical**, max diff 0 |
| Cell-level covariance consistent with the observer-level SE | aggregates back to **2.8e-17** |
| Recovery under simulated meridian cell loss (200 draws) | mean RMSE: ROI route 0.025, de-trend 0.016, **harmonic 0.014** |

Under simulated loss the harmonic route is the most accurate of the three, and the de-trend route's
benefit is less a lower average error than a narrower range of it — the ROI route's RMSE spans
0.009–0.073 across the eight cells (8×), de-trending 0.009–0.021 (2.3×).

## 5. Coverage

Empty (observer × ROI) cells, of 128 per map, at 4–8° with pRF R² > 0.1:

| map | empty | note |
|---|---|---|
| V1 | 0 | sparsest cell 22 vertices |
| V2 | 0 | |
| V3 | 2 | one cell per experiment |
| hV4 | 48 | median 4 vertices per cell; the 270° ROI is empty for **all** observers |

hV4 at 4–8° is not analysable at polar-angle resolution — a *structural* gap, not a sampling gap, and
no model recovers an ROI no observer has. Widening to 2–10° leaves 20 empty cells and no
all-observer gap, which is the band reported below.

## 6. Results

Harmonic route, precision-weighted. ⚠ marks τ̂² pinned at zero.

**Context effects (dg − da), 4–8°.** Both Cartesian effects decline monotonically up the hierarchy
and remain significant throughout.

| asymmetry | V1 | V2 | V3 |
|---|---|---|---|
| horiz−vert | −0.335 [−0.530, −0.140] *p*=.005 | −0.260 [−0.357, −0.163] *p*<.001 | −0.102 [−0.191, −0.013] *p*=.030 |
| card−obl | −0.183 [−0.303, −0.062] *p*=.009 | −0.163 [−0.254, −0.072] *p*=.004 | −0.107 [−0.200, −0.015] *p*=.028 |
| rad−tang | −0.061 *p*=.393 | +0.042 *p*=.170 | −0.067 [−0.119, −0.015] *p*=.019 ⚠ |
| polc−polo | −0.001 *p*=.974 | +0.010 *p*=.679 | +0.023 *p*=.456 |

**Per-experiment asymmetries, 4–8°.**

| exp | asymmetry | V1 | V2 | V3 |
|---|---|---|---|---|
| dg | horiz−vert | −0.549 *p*<.001 | −0.356 *p*<.001 | −0.136 *p*=.004 |
| dg | card−obl | −0.213 *p*=.006 | −0.167 *p*=.006 | −0.110 *p*=.024 |
| dg | rad−tang | 0.112 *p*=.005 | 0.173 *p*<.001 | 0.068 *p*=.001 |
| dg | polc−polo | 0.040 *p*=.035 | 0.037 *p*=.003 | 0.018 *p*=.271 |
| da | horiz−vert | −0.214 *p*=.011 | −0.092 *p*=.002 | −0.030 *p*=.026 |
| da | card−obl | −0.032 *p*=.079 | −0.005 *p*=.523 | −0.006 *p*=.439 |
| da | rad−tang | 0.176 *p*=.007 | 0.158 *p*<.001 ⚠ | 0.134 *p*<.001 ⚠ |
| da | polc−polo | 0.040 *p*=.300 | 0.019 *p*=.355 ⚠ | −0.007 *p*=.722 |

The Cartesian asymmetries attenuate up the hierarchy in both experiments; radial−tangential in the
polar experiment stays large and significant in all three maps. The headline reading is that the
context dependence is a V1-weighted phenomenon that weakens but does not disappear through V2 and V3.

### The two routes agree where coverage is adequate

| map | max \|de-trend − harmonic\| | rows changing significance |
|---|---|---|
| V1 | 0.0000 | 0 |
| V2 | 0.0000 | 0 |
| V3 | 0.0021 | 0 |
| hV4 | 0.0375 | 3 |

All disagreement is in hV4 — including the result the first pass flagged as suspicious:

| | de-trend | harmonic |
|---|---|---|
| hV4 da horiz−vert | +0.082 [0.053, 0.111] **sig** | +0.045 [−0.008, 0.098] **ns** |

That sign flip against V1/V2/V3 was largely an artifact of ROI binning, and the unbinned route
dissolves it. **The routes agree wherever coverage is adequate and disagree only where it is not**,
which is evidence about hV4's coverage rather than about hV4.

### hV4, 2–10°, harmonic route

Reported for completeness and **not comparable to the maps above** — different eccentricity band,
20 empty cells, and a route disagreement of 0.0375.

| exp | asymmetry | estimate |
|---|---|---|
| dg | horiz−vert | −0.146 [−0.225, −0.067] *p*=.003 |
| dg | card−obl | −0.071 [−0.175, 0.032] *p*=.148 |
| dg | rad−tang | 0.039 [0.005, 0.072] *p*=.029 |
| dg | polc−polo | 0.009 [−0.016, 0.034] *p*=.431 |
| da | horiz−vert | 0.045 [−0.008, 0.098] *p*=.084 |
| da | card−obl | 0.006 [−0.037, 0.049] *p*=.741 |
| da | rad−tang | 0.223 [0.087, 0.359] *p*=.006 |
| da | polc−polo | 0.059 [−0.051, 0.169] *p*=.242 |
| dg−da | horiz−vert | −0.171 [−0.268, −0.075] *p*=.004 |
| dg−da | card−obl | −0.074 [−0.188, 0.039] *p*=.167 |
| dg−da | rad−tang | −0.179 [−0.340, −0.018] *p*=.034 |
| dg−da | polc−polo | −0.057 [−0.172, 0.057] *p*=.274 |

## 7. Cautions

**τ̂² is pinned at zero in four rows** — two in V2, two in V3, none in V1 or hV4. There the
across-observer spread is at or below the measured within-observer noise, so the weights spread the
full σ range (up to 16×) and the interval is driven entirely by measured σ rather than by observer
variability. These are not ordinary random-effects intervals and should not be read as such.

**The one new significant cell sits in that regime.** The V3 rad−tang context effect, *p* = .019, has
τ̂² = 0. It is also one of 48 tests and survives no correction for multiplicity. Treat it as noise
unless something independent supports it.

**card−obl and polc−polo are not two findings in a map with class-structured ROI loss.** Per cell,
polc−polo is exactly ± card−obl, sign alternating with cardinal/oblique ROI (maximum difference over
all cells: 0.000, both experiments) — at cardinal ROIs radial and tangential *are* horizontal and
vertical. The two are one measurement, separated only by how cells combine across ROIs, which is
what ROI loss damages. Under 150 meridian-deletion draws on V1 (dg):

| route | error correlation | max difference | rmse card−obl | rmse polc−polo |
|---|---|---|---|---|
| ROI / de-trend | +0.465 | 5.1e-02 | 0.0100 | 0.0229 |
| harmonic | +1.0000 | 1.9e-15 | 0.0096 | 0.0096 |

The harmonic errors coincide because its 4th-harmonic part is `d_p = A + B*s_p` with `s_p = ±1`, so
with two parameters and two ROI classes the fit is saturated and a deletion touching only cardinals
leaves both errors equal. It is nonetheless **more accurate on both**, markedly so on polc−polo. It
does not create the dependence; it makes it exact and visible where the ROI route leaves it blurred
behind noisier estimates that merely look more independent. The constraint applies to V3 and
especially hV4 under either route — not to V1 or V2.

## 8. What the first pass got wrong

Recorded rather than rewritten, per the repository's correction convention.

- **"The harmonic route's inputs are still V1-only."** False. `_cache/areas.mat` supplies the
  per-vertex table for all four maps and the fitting code is area-agnostic; fitting V2 needed no
  code change. The only real gap was precision weighting, which is what §2 step 2 now provides.
- **"It is the route most exposed to the identifiability problem."** Backwards — see §7. The
  degeneracy is in the data; the harmonic route makes it visible and is more accurate regardless.
- **hV4's positive da horiz−vert was presented as needing scrutiny.** It got it: the effect is an
  ROI-binning artifact and is not significant under the harmonic route.
- **A bootstrap regression was introduced and fixed.** Pre-generating draws as
  `randi(nRun,[nBoot nRun])` fills column-major and consumes the random stream in a different order
  than the original per-iteration `randi(nRun,[1 nRun])`, silently changing every bootstrap SE in
  the repository (the V1 table drifted 2.8e-3 from the stored CSV). Draws must still be generated up
  front — the fitting path reseeds the global stream — but one row at a time.

## 9. Open

**hV4 should not be reported at polar-angle resolution.** At 4–8° it has 48 empty cells of 128 and a
structural gap at 270°; at 2–10° the two routes disagree by 0.0375 with three significance changes.
Whole-ROI, as [`local_qc/group_addv4mt.m`](local_qc/group_addv4mt.m) already does, is the defensible
treatment.

**The hierarchy trend has not been tested.** The claim that the context effect declines
monotonically V1 → V2 → V3 rests on six individually significant cells that fall in the same order
in both asymmetries. That is more than any single test, but the trend itself — a within-observer
V1 − V3 difference — has not been computed.

**`'thetaV','continuous'` has not been run across maps.** It is the better model, and the difference
from the binned default is exactly the within-wedge local-orientation term of `HARMONIC_MODEL.md`.
Running it would move the V1 numbers, so it is a decision rather than a default.
