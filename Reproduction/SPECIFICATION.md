# The analysis specification

**Settled 2026-08-24.** One analysis, applied unchanged to V1 and to all eight supplement maps
(V1, V2, V3, V3a, V3b, hV4, MT, MST), with no per-map customisation anywhere.

This document is *how the analysis works and why*. The numbers it produces are in
[`RESULTS.md`](RESULTS.md).

---

## 1. The specification in one table

| | choice |
|---|---|
| **Vertices** | one visual map, 4–8° eccentricity (2–10° reported alongside), pRF R² > 0.1 |
| **Units** | raw % signal change; each vertex demeaned across its four orientations. No z-scoring |
| **Model** | four-term harmonic, fitted per vertex, with **continuous** pRF polar angle θ_V |
| **Vertex weighting** | equal coverage over the eight 45° polar-angle ROIs |
| **Gain** | on, per observer × **map**, applied at the observer boundary |
| **Combining observers** | fit per observer, then average. **Equal weighting is primary**; precision weighting reported alongside |
| **Intervals** | *t* on n−1 = 7 df, primary. Percentile bootstrap reported alongside |
| **Reporting** | polar-angle-resolved results only where measured coverage passes the §7 criterion |

Adopting this specification moved the V1 numbers by at most 0.040 against the earlier ROI route and
changed no significance.

---

## 2. The model

Per vertex *v*, take the four stationary-orientation betas, subtract that vertex's mean across the
four, and predict with four weights shared across vertices:

```
y_vk = b1·cos(2θ) + b2·cos(4θ) + b3·cos(2(θ − θ_V)) + b4·cos(4(θ − θ_V))
```

| term | asymmetry | +1 means |
|---|---|---|
| `b1` | horizontal vs. vertical | horizontal |
| `b2` | cardinal vs. oblique | cardinal |
| `b3` | radial vs. tangential | radial |
| `b4` | polar-cardinal vs. polar-oblique | polar-cardinal |

Coefficients are reported as **2·b** — the pro-minus-con difference, in percent signal change,
directly comparable to the difference scores the figures plot.

**Demeaning removes the blank.** The "blank" condition is full-field pink noise, not a baseline
([`local_qc/DATA_QUALITY.md`](local_qc/DATA_QUALITY.md) §1), so removing each vertex's mean over the
four orientations means only orientation *differences* enter. All four predictors sum to zero across
the four conditions in both experiments, so no intercept is needed and the demeaning does not
distort the design.

**Angle convention.** `θ` is the orientation of the **bars** in conventional visual-field degrees
(0° = right horizontal meridian, CCW positive), so **0° is a horizontal grating and 90° a vertical
one**; `θ_V` is the polar angle of the vertex's pRF centre. This is *not* the "direction of
luminance variation" convention, which is rotated 90° and would flip the sign of both first-harmonic
terms. Full chain: [`STIMULUS_CONVENTIONS.md`](STIMULUS_CONVENTIONS.md).

**At the eight wedge centres the four predictors reduce exactly to the ±1/0 asymmetry codes the ROI
analysis uses** (`cleanroom/lme_codes.m`), so the harmonic route and the ROI route estimate the same
quantities. `cleanroom/test_harmonic_model.m` asserts this, and must pass before any result is
interpreted.

**What the model cannot represent.** Every term depends on θ_V only through cos(2θ_V), sin(2θ_V) and
cos(4θ_V), all of period 180°, so the model is invariant to θ_V → θ_V + 180° and cannot express
upper/lower or left/right field asymmetries. Such effects are averaged, not fitted. The full account
of the model — its identification structure, adequacy tests, and the geometry-versus-context result
it was built for — is [`supplement/SUPPLEMENT_harmonic_model.md`](supplement/SUPPLEMENT_harmonic_model.md).

---

## 3. Why the ROI-average route was replaced

The earlier V1 route took each observer's asymmetry as the mean over their eight polar-angle ROIs.
That is only meaningful if every observer contributes the same eight. In V1 they do — the sparsest
cell holds 22 vertices. Outside V1 they do not, and two facts make the gaps dangerous:

- **The asymmetry varies strongly with polar angle.** For `dg` polar-cardinal−polar-oblique the
  fitted profile alternates ±0.13–0.26 with cardinal/oblique ROI, against a mean effect of 0.04.
- **The loss is not random with respect to polar angle.** Map boundaries lie on the meridians, and
  vertices near a boundary have the poorest pRF fits, so cells empty along exactly the axis the
  asymmetry varies along.

Deleting meridian ROIs from V1 in simulation biases the ROI route by **+0.071** on `dg` polc−polo,
whose true value is 0.040 — cell loss nearly triples it.

**Measured in full 2026-08-25, and the larger effect is on radial−tangential.** Deleting MT's own 30
empty cells from V1 moves the ROI route's `dg` rad−tang from 0.103 to **−0.070** — a sign flip —
against 0.120 → 0.186 for the specification. Deleting the vertical meridian and deleting the
horizontal meridian push the ROI estimate in *opposite* directions (−0.048 and +0.273), which is the
signature of an unbalanced mean rather than of lost precision. The bias is systematic: it comes
almost entirely from the loss every observer shares, and randomising the observer-specific part
leaves the mean shift unchanged. → [`MISSING_DATA.md`](MISSING_DATA.md) §2–4

Note what is *not* broken. The inclusion mask is orientation-independent, so all four orientations of
a cell vanish together and the within-ROI subtraction is never partially missing. Only the average
over ROIs breaks.

---

## 4. Gain, per observer × map

The pRF-gain rescaling divides each observer's BOLD by their own mean gain and multiplies the group
gain back in. It is a per-observer scalar applied before observers are combined, so it touches no
within-observer quantity — but it does **not** commute with precision weighting (scaling `y_i` and
`σ_i` by `c_i` changes τ̂² and hence the weights), so it is applied at the observer boundary, never to
a finished group estimate.

Gain is computed per observer **per map** (`server_extract/collect_gain_areas.m`). Using a V1-derived
scalar everywhere would special-case V1 as the source; per-map scale factors correlate with the
V1 one at r = 0.88 (V2), 0.73 (V3), 0.66 (hV4), so V1 is a fair proxy for V2 and a poor one for hV4.

**The rescaling is amplitude-neutral per map** — the geometric mean of the scale factors is exactly 1
(verified to 12 decimals in every map), so it equalises observers *within* a map without altering
that map's response level. The V1 > V2 > V3 decline in gain, and the asymmetry attenuation up the
hierarchy, are therefore not artefacts of the normalisation.

| | V1 | V2 | V3 | hV4 (2–10°) |
|---|---|---|---|---|
| group geometric mean, % BOLD | 4.369 | 3.843 | 3.188 | 3.654 |
| across-observer range | 1.67× | 1.75× | 1.80× | 1.77× |

Gain is the smallest of the four decisions — worth at most 0.0127 in the group estimate (gain versus
no gain). The mov/stat protocols are combined by **geometric** mean, superseding `gainSummary.csv`'s
arithmetic average. Validation: V1 4–8° reproduces `gainSummary.csv` for all eight observers to
**5.3e-15** via the independent `dg_computeGain` path.
See [`../AnalysisCode/01_calculate_observer_gain/README.md`](../AnalysisCode/01_calculate_observer_gain/README.md)
for what gain means and how it is computed.

---

## 5. Equal-coverage weighting at 45°

Vertex density across the eight ROIs spans 5.4× in V1, 3.4× V2, 2.4× V3, 95× hV4. Natural weighting
therefore estimates an average over **cortex**; equal coverage estimates an average over the **visual
field**. When the question is polar-angle dependence, the cortical average confounds the asymmetry
with cortical magnification.

**This is an estimand choice, not numerical hygiene.** The proof is V2, where the design correlation
`r(b1,b3)` is 0.001 under natural weighting — no collinearity to fix — and the estimate still moves
by 0.0275.

That equal coverage also orthogonalises the design is a second, weaker, map-specific benefit: mean
`r(b1,b3)` under natural weighting is 0.349 in V1, 0.182 in V3, 0.375 in hV4, but 0.001 in V2. What
matters is not the density range but its alignment with cos(2θ_V): V1's excess sits at 0° and 180°
where cos(2θ_V) = +1 and does not cancel; V2's sits at obliques and does.

**45°, not finer.** The bin width controls leverage, since `w = 1/count` gives a bin holding one
vertex the same total weight as a bin holding two hundred. The ratio below is between **individual
vertices** — do not confuse it with the *precision*-weight ratio across observers that the §7
coverage criterion thresholds at 25; they are different quantities on different units:

| map | width | min count | max **vertex**-weight ratio | eff N | r(b1,b3) |
|---|---|---|---|---|---|
| V1 | 45° | 33 | **11×** | 1005 | 0.075 |
| V1 | 15° | 2 | 101× | 742 | 0.016 |
| V3 | 45° | 6 | 27× | 433 | 0.086 |
| V3 | 15° | 1 | 89× | 249 | 0.035 |

A residual `r(b1,b3)` of 0.075 is a variance inflation of 1.006 — nothing. A 100× leverage on a
single vertex is a real hazard. Finer binning does not even help asymptotically: with continuous θ_V
in V1 the residual runs natural 0.349, 8 bins 0.075, 24 bins 0.016, 48 bins 0.009, **96 bins 0.023,
360 bins 0.080** — inverse-count weighting on fine bins is a noisy density estimate, so it turns back
up. A von Mises kernel gives 0.026, no better. 45° also keeps **one** binning in the pipeline and
makes the estimand exactly "equal weight per polar-angle ROI", which is what the ROI route computes.

**Empty weighting bins are a non-event.** They are not analysis units; only vertices that exist are
indexed, so an empty bin never divides and never creates a missing value (verified: 0 non-finite
weights across all maps and observers). Bin width has nothing to do with the missing-data problem.

---

## 6. Continuous θ_V, equal weighting of observers, and identifiability

**θ_V is continuous.** Binned θ_V puts cos(4θ_V) at only two values, making `b4 = b2 × cos(4θ_V)` a
two-point design — exactly orthogonal at full coverage, exactly degenerate (VIF = ∞) once one ROI
class is lost. Continuous θ_V never degenerates: VIF 4.12 even with all four cardinal ROIs removed.

**Nothing is near degenerate at real coverage** — max |r| = 0.467, max VIF 1.28 across all maps,
weightings and θ_V choices. But the second-harmonic pair deserves care. Per cell, polc−polo is
exactly ± card−obl with the sign alternating cardinal/oblique (maximum difference over all cells:
0.000, both experiments), so the two are one measurement separated only by how cells combine across
ROIs — which is what ROI loss damages. **In a map with class-structured ROI loss they are not two
findings**, under either route.

**Equal weighting of observers is primary.** With 8 runs, σ̂ carries ~7 df and scatters even when
every observer has identical true precision: the null max/min over 8 observers has median 2.20× and
95th percentile 3.55×, against an observed median of 3.37×. Only 14 of 32 rows exceed the null, so
most of what you would weight by is estimation noise. Simulation at this dataset's operating point
puts precision weighting **4% worse** than equal weighting; it only pays when τ is small *and* the
spread is large. Empirically it changes little — max |equal − precision| is 0.029 (V1), 0.016 (V2),
0.017 (V3), 0.024 (hV4) — and where it does most work is the τ̂² = 0 rows, which are the least
trustworthy. Both columns are reported in every table; this is only about which is quoted. The
mechanics of the precision weighting are in
[`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) §4.

A real and predictable structure sits underneath: **matched-frame contrasts are ~3× noisier** (mean
σ 0.081 versus 0.027) because the contrast projects onto correlated trial-wise noise, and they carry
more genuine between-observer spread. But that variation is *between asymmetries*, and precision
weighting only ever uses the between-observer part.

---

## 7. Coverage — which maps can be reported by polar angle

All eight maps were extracted and all have per-map gain. Whether each can support a
*polar-angle-resolved* analysis is a separate, measurable question, decided on a fixed criterion
rather than case by case:

> at most **2 empty** (observer × ROI) cells of 64, a **median of ≥ 20 vertices** per cell, a
> **maximum precision-weight ratio < 25**, and a **per-map gain that exists**.

The measured coverage table and which maps pass are in [`RESULTS.md`](RESULTS.md) §5. In short:
**V1, V2 and V3 at 4–8°, plus V3a at 2–10°.**

Failing the criterion says **this design cannot resolve that map by polar angle** — not that the map
has no asymmetries. Whole-ROI analyses of hV4 and MT remain available and are unaffected; that is
what `local_qc/group_addv4mt.m` already did.

**What the criterion is really screening for, measured 2026-08-25.** The three thresholds are
proxies; the quantity behind them is now measured directly. Reproducing each map's whole coverage
profile in V1 — its empty cells *and* its vertex counts, `subsample_cells.m` — gives a 90% band on
`dg` rad−tang of 0.018 wide under V3's coverage and **0.131 wide under MT's**, against an effect of
0.119. A map whose band is wider than the effect cannot resolve it however the cells are combined.
Note that this is **sparsity, not holes**: V1 with MT's 30 cells deleted still holds 7045 vertices,
MT holds 545, a median of one per cell, and the specification largely absorbs the holes on their own.
→ [`MISSING_DATA.md`](MISSING_DATA.md) §6

Turning that measurement into a calibration of the thresholds themselves has not been done; nothing
reported depends on it, since V3 and MT are separated by a factor of seven in band width.

The gain clause matters only in principle: MT and MST lack a per-map gain for some observers, and
scoring them with a V1 scalar would reintroduce the V1 special-casing §4 removed. They fail on
coverage anyway, so it changes no map's status.

---

## 8. Running it

```bash
matlab -batch "cd('Reproduction/cleanroom'); run_spec_outputs"
```

One command, two to three minutes. Needs `~/dg_collect/` (run-wise betas and per-map gain — see
[`server_extract/README.md`](server_extract/README.md)); it does not touch
`Support/allsubjectsTable.csv`. Everything reported is regenerated by this; nothing is transcribed
by hand.

**Three variants** are produced, because route and weighting are orthogonal knobs and the comparison
is only readable if they change one at a time (`cleanroom/spec_variants.m`):

| tag | route | weighting | |
|---|---|---|---|
| **`spec`** | harmonic model, continuous θ_V | equal | **primary** — the settled specification |
| `roi` | eight polar-angle wedges | equal | `spec` vs `roi` isolates the **route** |
| `roipw` | eight polar-angle wedges | precision | `roi` vs `roipw` isolates the **weighting** |

There is deliberately no harmonic-plus-precision variant: adding it would leave no pair differing in
exactly one thing. All twelve resulting figures are laid out side by side, figure by figure, in
[`supplement/FIGURE_VARIANTS.md`](supplement/FIGURE_VARIANTS.md).

**Outputs.** Figures in `supplement/figures/` (PNG + PDF), tables in `supplement/` (CSV), both
tracked. `<tag>` ∈ {`spec`, `roi`, `roipw`}:

| file | what it is |
|---|---|
| `Figure_5_<tag>_dg`, `Figure_6_<tag>_da` | Figures 5 and 6, V1 4–8° |
| `Figure_5_6_<tag>_profile` | the polar-angle profile |
| `Figure_S5_<tag>_hierarchy` | the extrastriate supplement figure |
| `spec_asymmetries_<tag>_v1_4-8.csv` | the eight V1 asymmetries |
| `spec_context_<tag>_v1_4-8.csv` | the four V1 context effects |
| `spec_perobserver_<tag>_v1_4-8.csv` | every observer's value behind those |
| `spec_areas_coverage_<tag>.csv` | eight maps × two bands, coverage and reportability |
| `spec_areas_{asymmetries,context,perobserver}_<tag>.csv` | the same quantities in every map × band, flagged |
| `spec_areas_trend_<tag>.csv` | the hierarchy trend |

Every table carries `variant`, `route` and `weighting` columns and **both** the equal-weighted
(`eq_*`) and precision-weighted (`pw_*`) estimates in every row. Reporting only one would hide that
the choice is a choice. `Reproduction/figures/` is the git-ignored scratch destination for older
diagnostic plots and holds no copy of these — one image, one location.

Individual pieces are callable without the full sweep: `spec_profiles('area','V2')`,
`spec_tables('area','V2','variant','roi')`, `spec_areas_summary('variant','roi')`.

**Where the time goes**, measured: the 500-resample bootstrap over runs inside
`diagnose_within_observer_error`, which measures each observer's within-observer σ — ~75% of every
call (4.7 s per map × band in V1, against 1.2 s at `nBoot = 50`). File I/O is negligible; the
harmonic fit itself is ~2.7 s. `nBoot` is the only real lever, and lowering it degrades σ, which is
what the precision weighting rests on, so it stays at 500.

### Code

| file | role |
|---|---|
| `cleanroom/spec_profiles.m` | fits the specification per observer; returns the model and observed profiles |
| `cleanroom/spec_tables.m` | the V1 tables, four interval methods |
| `cleanroom/spec_areas_summary.m` | the eight-map sweep, coverage, reportability, hierarchy trend |
| `cleanroom/spec_group.m` | the group estimator, equal or precision — **one** definition |
| `cleanroom/spec_axis_limits.m` | the shared axis limits — **one** definition, spanning all variants |
| `cleanroom/spec_variants.m` | the three variants and what each isolates |
| `cleanroom/load_runbetas_area.m` | the analysed vertex set — **one** definition |
| `cleanroom/plot_fig5_6_spec.m`, `plot_spec_profile.m`, `plot_spec_hierarchy.m` | the figures |
| `cleanroom/run_spec_outputs.m` | the driver |
| `cleanroom/harmonic_*.m`, `fit_harmonic_vertex.m` | the model itself |
| `cleanroom/test_harmonic_model.m` | assertions — must pass first |
| `cleanroom/diagnose_within_observer_error.m` | per-observer σ from bootstrapping runs |
| `cleanroom/precision_weighted_table.m`, `precision_weighted_cells.m` | the precision-weighted estimates |
| `server_extract/collect_gain_areas.m`, `collect_runwise_betas_areas.m` | the extractions this needs |

**Missing-data diagnostics**, run on demand and not part of `run_spec_outputs` —
[`MISSING_DATA.md`](MISSING_DATA.md) §8:

| file | role |
|---|---|
| `cleanroom/cell_occupancy.m` | vertices per (observer × ROI) for any map, from labels alone |
| `cleanroom/diagnose_cell_loss.m` | delete a map's empty cells from V1; both routes, null, class test |
| `cleanroom/diagnose_loss_structure.m` | is the shift systematic, and how much data does each map have |
| `cleanroom/diagnose_pooled_fit.m`, `spec_pooled.m` | one group fit versus averaging per-observer fits |
| `cleanroom/subsample_cells.m` | thin a loaded dataset to a target vertex count per cell |

`spec_profiles`, `diagnose_within_observer_error` and `load_runbetas_area` carry an additive
`dropCells` argument for these. It defaults to empty, and empty is the specification unchanged: the
V1 baseline reproduces `supplement/spec_asymmetries_spec_v1_4-8.csv` exactly with the option in place.

**Two guards keep the figures and the tables from drifting apart.** `spec_profiles` asserts that its
per-observer asymmetries equal the ones `diagnose_within_observer_error` computes (agreement 3e-16 in
V1), and `spec_tables` asserts that its equal-weighted means equal `precision_weighted_table`'s. Both
fire on mismatch rather than reporting two versions of the same number.

---

## 9. Traps in this code

Current guidance, not history. Each of these produced wrong numbers silently.

- **Bootstrap draws must be generated one row at a time.** Pre-generating them as a matrix fills
  column-major and consumes the random stream in a different order than the original per-iteration
  call, **silently changing every bootstrap SE** — no error, no warning, different numbers. Draws do
  have to be generated up front, since the fitting path reseeds the global stream; generate them a
  row at a time.
- **Do not mask a full-surface vector with a mask defined over `vertIndex`.** An indexing bug in the
  gain summary did exactly that and read the wrong elements. It was caught only because an
  independent path (`dg_computeGain`) existed to check against.
- **Read significance off the equal-weighted column**, which is primary. Reading it off the
  precision-weighted column once produced a spurious report that adopting continuous θ_V moved
  `dg` polc−polo across *p* = .05.
- **σᵢ must be measured on the rescaled data.** Rescaling observer *i* by *sᵢ* multiplies their
  measurement SE by *sᵢ* too; taking σᵢ from unrescaled data and applying the weights to rescaled
  data would systematically over-weight high-gain observers. The clean-room scales `asymSE` by the
  same factor as `asym`, so this is handled — but it is the easy mistake.
- **Do not "align" the wedge dimension of `medianBOLDpa` with `[0 45 90 …]`.** It is stored in
  Benson order and the downstream code expects that; "fixing" it introduces a bug into working code.
  → [`STIMULUS_CONVENTIONS.md`](STIMULUS_CONVENTIONS.md) §3.
