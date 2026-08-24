# Orientation asymmetries beyond V1

**Status: first pass, 2026-08-24.** The procedure is validated; the numbers are a first run and
carry three cautions listed at the end. Code: [`cleanroom/precision_weighted_cells.m`](cleanroom/precision_weighted_cells.m),
[`cleanroom/fit_cell_meta.m`](cleanroom/fit_cell_meta.m), [`cleanroom/validate_cell_meta.m`](cleanroom/validate_cell_meta.m),
[`server_extract/collect_runwise_betas_areas.m`](server_extract/collect_runwise_betas_areas.m).
Numbers: [`supplement/precision_weighted_areas.csv`](supplement/precision_weighted_areas.csv).

---

## 1. Why the V1 procedure does not transfer

The V1 route takes each observer's asymmetry as the mean over their eight polar-angle ROIs, then
combines observers by precision weighting ([`LME.md`](LME.md) §5). That mean is only meaningful if
every observer contributes the same eight ROIs. In V1 they do — the sparsest cell still holds 22
vertices. Outside V1 they do not.

Two facts make the resulting gaps dangerous rather than merely inconvenient:

- **The asymmetry varies strongly with polar angle.** Fitted ROI profiles are large: for dg
  polar-cardinal−polar-oblique the profile alternates ±0.13–0.26 with cardinal/oblique ROI, against
  a mean effect of 0.039.
- **The loss is not random with respect to polar angle.** Visual map boundaries lie on the
  meridians, and vertices near a boundary have the poorest pRF fits, so cells empty out along
  exactly the axis the asymmetry varies along.

Averaging "whichever ROIs survived" therefore estimates a different quantity for each observer, and
the group mean is a composition-weighted mixture. Deleting meridian ROIs from V1 in simulation
biases the V1 route by **+0.071** on dg polc−polo — an effect whose true value is 0.040. Cell loss
nearly triples it.

Note what is *not* broken. The inclusion mask is orientation-independent
([`bin_and_aggregate.m:16`](cleanroom/bin_and_aggregate.m)), so all four orientations of a cell
vanish together and the within-ROI subtraction is never partially missing. Only the average over
ROIs breaks.

## 2. Procedure

Write `d_ip` for the asymmetry of observer *i* in polar-angle ROI *p* — one number per cell, 8 × 8
per map and experiment.

**Step 1 — extract run-wise betas for eight areas.** `collect_runwise_betas_areas.m` writes the mean
GLMsingle single-trial beta per (vertex, condition, run) for V1, V2, V3, V3a, V3b, hV4, MT and MST.
All eight areas are taken in one pass because the cost is streaming `modelmd` (425 MB per
subject-experiment, 6.8 GB over 16 files) and that cost is identical for one area or eight — the
whole variable must cross the wire before anything can be subset.

**Step 2 — measure each cell's sampling error.** Resample the observer's runs 500× and recompute,
giving `d_ip`'s variance *and* the full 8 × 8 covariance across ROIs within each observer. Runs, not
vertices: an earlier vertex-resampling estimate was withdrawn because it holds the GLM betas fixed
and characterises which patch of cortex was sampled, not the reliability of the measurement.

The covariance is emphatically not diagonal, and the structure is interpretable rather than noise —
mean between-ROI correlation is **+0.6 when the asymmetry's reference frame matches the experiment**
and **−0.1 when it does not**. A matched frame uses the same two stimulus conditions in every ROI,
so a run-level fluctuation moves all eight ROIs together; an unmatched frame rotates which stimulus
counts as "radial" with polar angle, so the same fluctuation pushes some ROIs up and others down.

**Step 3 — estimate the polar-angle profile.** Fit, per map, experiment and asymmetry,

```
d_ip = mu + alpha_p + u_i + v_ip + e_ip
```

by REML, with `alpha_p` fixed and sum-to-zero (the profile is systematic, not a draw from a
population), `u_i` an observer offset, `v_ip` observer-by-ROI heterogeneity, and `e_ip` sampling
error with the **known** covariance from step 2. Missing cells are simply absent rows. This is a
multilevel random-effects meta-analysis — `metafor::rma.mv` in R; MATLAB's `fitlme` cannot express
known per-cell sampling variances, which is why the fit is written out directly.

**Step 4 — de-trend, then average.** Take each observer's number as

```
y_i = mean over present p of ( d_ip - alpha_p )
```

After removing the polar-angle profile, which ROIs an observer happens to have is no longer
confounded with where the asymmetry is large.

**Step 5 — combine observers exactly as before.** Feed `y_i` and its bootstrap variance into the
unchanged precision-weighted estimator from `LME.md` §5.

### Why the model supplies only the profile, not the estimate

The fit in step 3 also produces `mu` directly, and it is tempting to report it. It should not be
used. On **complete** V1 data — where the answer is known — `mu` moves by up to 0.05 depending on
how the within-observer covariance is modelled, because GLS re-weights ROIs *within* an observer
using a covariance estimated from 8 runs, and sampling noise is correlated with effect size
(r up to 0.52 after removing observer and ROI means). Five of eight estimates move toward zero and
three away, so it is not a correctable attenuation. This is `LME.md` §5's own caveat — that
estimated-weight GLS at n = 8 can add variance rather than remove it — appearing in practice.

De-trending uses the model only for the profile, which is what the missing cells actually require,
and leaves the weighting at the observer level where it is already validated.

### Why this is safe to adopt

The `alpha_p` sum to zero over the ROI set, so **with complete data `y_i` is the plain ROI mean and
the procedure returns the existing numbers exactly.** Adopting it costs nothing in V1.

## 3. Validation

| check | result |
|---|---|
| Complete data reduces to the existing estimator | V1 and V2 reproduce `precision_weighted.csv` to **5.6e-17** |
| Multi-area extraction preserves the V1-only extraction | V1 subset **bit-identical** — same indices, `isequaln`, max diff 0 |
| Cell-level covariance is consistent with the observer-level SE | aggregates back to **2.8e-17** |
| Recovery under simulated meridian cell loss (200 draws) | mean RMSE **0.025 → 0.016** |

On the simulation, the benefit is not a lower average error so much as a narrower range of it: the
V1 route's RMSE spans 0.009–0.073 across the eight cells (8×), de-trending spans 0.009–0.021 (2.3×).
It removes the catastrophic case rather than improving the typical one. De-trending wins 4 of 8
cells outright and is slightly worse in the others, which is the expected cost of estimating seven
extra profile parameters.

## 4. Coverage

Empty (observer × ROI) cells, of 128 per map:

| map | 4–8° | note |
|---|---|---|
| V1 | 0 | sparsest cell 22 vertices |
| V2 | 0 | |
| V3 | 2 | one cell per experiment |
| hV4 | 48 | median 4 vertices per cell; the 270° ROI is empty for **all** observers |

hV4 at 4–8° is not analysable at polar-angle resolution — that is a structural gap, not a sampling
gap, and no model recovers an ROI no observer has. Widening to 2–10° leaves 20 empty cells and no
all-observer gap, which is what is reported below; `local_qc/group_addv4mt.m` had already gone
whole-ROI for hV4 on the same grounds.

## 5. Results

**Context effects (dg − da), 4–8°, precision-weighted.** Both Cartesian effects decline
monotonically up the hierarchy and remain significant throughout.

| asymmetry | V1 | V2 | V3 |
|---|---|---|---|
| horiz−vert | −0.335 [−0.530, −0.140] *p* = .005 | −0.260 [−0.357, −0.163] *p* < .001 | −0.101 [−0.189, −0.013] *p* = .030 |
| card−obl | −0.183 [−0.303, −0.062] *p* = .009 | −0.163 [−0.254, −0.072] *p* = .004 | −0.107 [−0.198, −0.015] *p* = .029 |
| rad−tang | −0.061 *p* = .393 | +0.042 *p* = .170 | −0.067 [−0.119, −0.015] *p* = .019 ⚠ |
| polc−polo | −0.001 *p* = .974 | +0.010 *p* = .679 | +0.021 *p* = .493 |

**Per-experiment asymmetries, 4–8°, precision-weighted.** ⚠ marks τ̂² pinned at zero.

| exp | asymmetry | V1 | V2 | V3 |
|---|---|---|---|---|
| dg | horiz−vert | −0.549 *p* < .001 | −0.356 *p* < .001 | −0.135 *p* = .004 |
| dg | card−obl | −0.213 *p* = .006 | −0.167 *p* = .006 | −0.110 *p* = .024 |
| dg | rad−tang | 0.112 *p* = .005 | 0.173 *p* < .001 | 0.068 *p* = .001 |
| dg | polc−polo | 0.040 *p* = .035 | 0.037 *p* = .003 | 0.017 *p* = .318 |
| da | horiz−vert | −0.214 *p* = .011 | −0.092 *p* = .002 | −0.030 *p* = .024 |
| da | card−obl | −0.032 *p* = .079 | −0.005 *p* = .523 | −0.007 *p* = .410 |
| da | rad−tang | 0.176 *p* = .007 | 0.158 *p* < .001 ⚠ | 0.134 *p* < .001 ⚠ |
| da | polc−polo | 0.040 *p* = .300 | 0.019 *p* = .355 ⚠ | −0.006 *p* = .730 |

The Cartesian asymmetries attenuate up the hierarchy in both experiments; radial−tangential in the
polar experiment stays large and significant in all three maps. The headline reading is that the
context dependence is a V1-weighted phenomenon that weakens but does not disappear through V2 and V3.

## 6. Three cautions

**τ̂² is pinned at zero in four rows** — two in V2, two in V3, none in V1. There the across-observer
spread is at or below the measured within-observer noise, so the weights spread the full σ range
(up to 16×) and the interval is driven entirely by measured σ rather than by observer variability.
Those are not ordinary random-effects intervals and should not be read as such.

**The one new significant cell sits in that regime.** V3 rad−tang context effect, *p* = .019, with
τ̂² = 0. It is also one of 48 tests in this table and survives no correction for multiplicity. Treat
it as noise unless something independent supports it.

**hV4 is reported but is not comparable.** Different eccentricity band (2–10° against 4–8°), 20
empty cells, and da horiz−vert comes out **+0.082 [0.053, 0.111]** — positive, against negative in
all three of V1/V2/V3. A clean sign flip with a tight interval, in the map with the worst coverage,
is more likely a coverage artifact than a result. It needs scrutiny before it is used.

## 7. Open

**Second-harmonic identifiability has not been checked per map.** Under meridian ROI loss, card−obl
and polc−polo become perfectly confounded — in the harmonic parameterisation their errors correlate
at exactly +1.0000, and the design correlation r(b2, b4) goes to −1 when a whole ROI class is lost.
V3 has empty cells and hV4 has 20, so those two columns may be partly the same measurement in those
maps while the table presents them as separate numbers. This should be checked before either is
reported.

**The harmonic route is the better estimator and is not yet wired to the extrastriate path.** On the
same simulated cell loss it reaches mean RMSE 0.014 against de-trending's 0.016 and the V1 route's
0.025, and it reproduces the ROI route *exactly* on complete data. But `fit_harmonic_vertex` and its
inputs are still V1-only, and it is the route most exposed to the identifiability problem above.
