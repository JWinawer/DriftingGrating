# Missing data — what cell loss does, and why some maps are not reported by polar angle

**Measured 2026-08-25.** Everything here is a simulation *in V1*, where the answer on complete data
is known, so any change is attributable to the manipulation and nothing else.

The analysis divides each map into 64 **cells** — 8 observers × 8 polar-angle ROIs. In V1 every cell
is populated. Outside V1 some hold no vertex that passes the inclusion filter, and
[`SPECIFICATION.md`](SPECIFICATION.md) §7 refuses to report a map by polar angle when too many are
empty. This document measures what that loss actually does, and separates two things the coverage
criterion had been treating as one.

**The short version.** Cell loss has two distinct consequences, and they need separating.

1. **Holes bias the estimate, systematically and in a direction fixed by which ROIs are lost.**
   The bias comes almost entirely from the loss that *every observer shares*. It is modest under the
   settled specification (worst case 0.066, on `dg` rad−tang, whose value is 0.119) and severe under
   the ROI-average route (0.173, a sign flip).
2. **Sparsity is a much bigger problem than holes, and it is the reason MT is excluded.** V1 with
   MT's holes still holds 7045 vertices. MT holds 545, a median of one per cell. Reproduce MT's
   whole coverage profile in V1 and the 90% band on `dg` rad−tang is 0.131 wide — wider than the
   effect being measured.

Fitting one group model to all observers' vertices instead of averaging per-observer fits does
**not** help with either, and costs the per-observer estimates the within-subject tests need.

---

## 1. What the missing cells look like

MT at 4–8° is the worst case among the eight maps: **30 empty cells of 64**. The pattern is not
scattered.

| ROI | 0° | 45° | **90°** | 135° | 180° | 225° | **270°** | 315° |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| observers with an empty cell | 2 | 3 | **8** | 2 | 1 | 3 | **8** | 3 |

**MT loses the entire vertical meridian — both wedges, all eight observers.** The other 14 empty
cells are scattered over particular observers, and one observer (sub-wlsubj123) loses all eight
wedges and so leaves the analysis altogether. Deleting the same 30 cells from V1 destroys 4030 of
its 11075 vertices.

`cell_occupancy.m` computes this per-cell table for any map from `labels_*.mat` and `ret_*.mat`
alone — occupancy depends on the label, eccentricity and pRF R², never on the betas — so it needs no
`runbetas_areas_*` extraction and reproduces the empty-cell and median columns of
`supplement/spec_areas_coverage_spec.csv` exactly.

---

## 2. Deleting those cells from V1

`diagnose_cell_loss.m`. Group estimates over the seven observers finite in **both** conditions, so
the delta is the cell loss and not the observer that disappeared.

| exp | asymmetry | spec: base → deleted | Δ | ROI route: base → deleted | Δ |
|---|---|--:|--:|--:|--:|
| dg | horiz−vert | −0.557 → −0.623 | −0.066 | −0.558 → −0.575 | −0.017 |
| dg | card−obl | −0.217 → −0.198 | +0.020 | −0.220 → −0.215 | +0.005 |
| dg | rad−tang | 0.120 → 0.186 | +0.066 | 0.103 → **−0.070** | **−0.173** |
| dg | polc−polo | 0.071 → 0.053 | −0.018 | 0.041 → 0.075 | +0.034 |
| dg−da | horiz−vert | −0.313 → −0.361 | −0.048 | −0.321 → −0.362 | −0.041 |
| dg−da | card−obl | −0.181 → −0.161 | +0.020 | −0.195 → −0.166 | +0.029 |

Two quantities cross *p* = .05 under the specification: `dg` polc−polo (.013 → .111) and the
card−obl context effect (.029 → .052). The horiz−vert context effect, the strongest claim in
[`RESULTS.md`](RESULTS.md), survives in both routes.

At **2–10°**, where MT's pattern (n = 19) leaves all eight observers standing, the split is stark:

- **Specification: every one of the twelve quantities moves by ≤ 0.015 and no *p* crosses .05.**
- **ROI route: `dg` rad−tang goes 0.117 → −0.024**, a sign flip, *p* .003 → .663. `da` horiz−vert
  loses significance (.008 → .058).

---

## 3. The bias is systematic, not noise

The shifts in §2 go up for some asymmetries and down for others, which looks like noise. It is not.
An extrastriate map's loss has a **common** part and an **idiosyncratic** part, and only the common
part can bias the group. `diagnose_loss_structure.m` Part A holds the vertical meridian fixed —
empty for all eight observers — and randomises the 14 scattered cells, 200 draws each way.

| quantity | MT actual Δ | vertical meridian only | permute observers | sd | permute obs + ROIs |
|---|--:|--:|--:|--:|--:|
| dg horiz−vert | −0.066 | **−0.044** | −0.043 | 0.024 | −0.048 |
| dg card−obl | +0.020 | **+0.007** | +0.005 | 0.012 | +0.006 |
| dg rad−tang | +0.066 | **+0.040** | +0.038 | 0.019 | +0.040 |
| da rad−tang | +0.019 | **+0.027** | +0.027 | 0.014 | +0.028 |
| dg−da horiz−vert | −0.048 | **−0.030** | −0.019 | 0.022 | −0.030 |
| dg−da rad−tang | +0.047 | **+0.013** | +0.011 | 0.022 | +0.012 |

**The randomised means land on the vertical-meridian-only column in every row.** The scattered cells
contribute nothing on average; they only add scatter, sd ≈ 0.01–0.025. The whole mean shift is the
loss every observer shares, and its sign per asymmetry is fixed by which ROIs are gone. "Some up,
some down" is four asymmetries responding differently to one deterministic cause.

MT's actual configuration is a somewhat unlucky draw: `dg` rad−tang +0.066 sits at the top of the
randomisation band [0.007, 0.067], as does the rad−tang context effect.

**Do not read this as "loss is random, so it averages out."** It does not average out. It is
reproducible: delete the same cells and you get the same answer.

---

## 4. Why — the wedges stop balancing

`diagnose_cell_loss.m`'s class test deletes a whole polar-angle class from *every* observer.

| deleted | route | horiz−vert | card−obl | rad−tang | polc−polo |
|---|---|--:|--:|--:|--:|
| (none) | spec | −0.548 | −0.221 | 0.119 | 0.072 |
| vertical meridian | spec | −0.592 | −0.214 | 0.159 | 0.077 |
| horizontal meridian | spec | −0.592 | −0.228 | 0.104 | 0.061 |
| obliques | spec | −0.496 | −0.241 | 0.088 | 0.093 |
| (none) | ROI | −0.548 | −0.221 | 0.104 | 0.039 |
| vertical meridian | ROI | −0.545 | −0.231 | **−0.048** | 0.115 |
| horizontal meridian | ROI | −0.597 | −0.236 | **0.273** | 0.110 |
| obliques | ROI | −0.481 | −0.181 | 0.078 | **−0.181** |

Under the ROI route, deleting the vertical meridian and deleting the horizontal meridian push
rad−tang in **opposite** directions. That is the signature of an unbalanced mean, not of lost
precision. A radial grating is horizontal on the horizontal meridian and vertical on the vertical
meridian, so the ROI average over eight wedges is only the intended quantity while the two classes
balance; break the balance and V1's large horiz−vert asymmetry (−0.55) leaks straight into the
polar-frame estimate. The same mechanism moves polc−polo under oblique deletion, since polc−polo is
exactly ± card−obl with the sign alternating cardinal/oblique ([`SPECIFICATION.md`](SPECIFICATION.md)
§6).

The specification fits continuous θ_V and never forms that mean, so it stays within ~0.04 in all
three cases. This is the measured version of the argument in
[`SPECIFICATION.md`](SPECIFICATION.md) §3; it confirms the +0.071 figure quoted there for polc−polo
and adds the larger and more dangerous effect, which is on rad−tang.

---

## 5. A single pooled group fit does not help

`spec_pooled.m` and `diagnose_pooled_fit.m`. The specification fits each observer separately and
averages the eight fits. The alternative stacks every observer's vertices into one weighted least
squares, with weights on two levels: equal coverage over the eight ROIs within an observer, and
equal total weight across observers so that map size does not tilt the fit.

**The cross-observer normalisation is necessary and it works.** Without it, pooling weights observers
by vertex count — `harmonic_weights` rescales to mean 1, so `sum(w)` = vertex count, a 1.59× range
across these observers.

| dg, no deletion | horiz−vert | card−obl | rad−tang | polc−polo |
|---|--:|--:|--:|--:|
| pooled, count-weighted | −0.565 | −0.232 | 0.127 | 0.069 |
| pooled, observers equalised | −0.556 | −0.217 | 0.120 | 0.063 |
| specification average | −0.557 | −0.217 | 0.120 | 0.071 |

Equalising is worth up to 0.016 and lands the pooled fit on the specification's answer. On complete
data the two estimators agree to ≤ 0.008 across all eight asymmetries.

**But pooling buys no robustness.** Under MT's 4–8° pattern:

| quantity | avg (spec) Δ | pooled, obs equal Δ | pooled, cell equal Δ |
|---|--:|--:|--:|
| dg horiz−vert | −0.066 | −0.068 | −0.050 |
| dg rad−tang | +0.066 | +0.059 | +0.060 |
| dg−da horiz−vert | −0.048 | −0.052 | −0.024 |
| dg−da card−obl | +0.020 | +0.027 | +0.027 |

Weighted least squares is linear in *y*, so with M_i = X_i′W_iX_i the pooled fit is exactly

> **b**_pool = (Σ M_i)⁻¹ (Σ M_i **b**_i)

a **matrix-weighted average of the same per-observer fits** the specification averages equally. Any
robustness it could have must come from reweighting observers, and there is almost none available:
under MT's pattern the observer shares run 0.137–0.148 against 0.143 uniform, even for an observer
down to five of eight ROIs. Equal-coverage weighting within an observer has already removed the
vertex-count dependence, so a holed observer still carries nearly full weight.

The deeper reason is §3: MT's damaging loss is **common to all observers**. Pooling combines *across*
observers and cannot recover information missing from every one of them. When every observer loses
the same class the two pooled weightings become algebraically identical — they differ by one common
constant, which weighted least squares ignores — and both track the averaged fit.

**Tested where it should win.** Idiosyncratic loss is the case the pooling intuition rests on. RMSE
of `dg` rad−tang against the no-deletion estimate, 100 draws per row:

| loss | k of 8 | avg | pooled, obs | pooled, cell |
|---|--:|--:|--:|--:|
| idiosyncratic | 2 | 0.0115 | 0.0117 | 0.0119 |
| idiosyncratic | 4 | 0.0233 | 0.0240 | 0.0235 |
| idiosyncratic | 6 | 0.0464 | 0.0431 | 0.0416 |
| common | 4 | 0.0283 | 0.0277 | 0.0277 |
| common | 6 | 0.0537 | 0.0534 | 0.0534 |

Pooling is ~7% better only at the most extreme severity, six of eight ROIs gone from every observer.
At matched k, idiosyncratic loss is *less* damaging than common loss for both estimators — random
per-observer loss partly averages out, shared loss does not.

**And pooling has a cost.** A single group fit has no per-observer estimate, so the within-observer
*t* test that every context claim rests on (standing fact 6) has nothing to run on. `spec_pooled`
uses a delete-one-observer jackknife instead, and the context effect becomes a difference of two
group fits rather than the average of within-observer differences. Comparable in spirit, not
identical. **The specification's per-observer route stays primary.**

One incidental result worth keeping: `equalcell` — *not* compensating an observer for their missing
ROIs — beats `equalobserver` on the Cartesian terms (`dg` horiz−vert Δ −0.050 vs −0.068; the
horiz−vert context effect −0.024 vs −0.052). Renormalising an observer's weights over the ROIs they
happen to have inflates exactly the surviving wedges that are unbalanced. It does not help on
rad−tang, so it is a lead rather than a recommendation.

---

## 6. Why MT is not reported by polar angle

**Everything above deletes cells and leaves V1's density intact everywhere else.** That tests the
holes only, and the holes are not MT's main problem. V1 with MT's 30 cells deleted still holds 7045
vertices; **MT holds 545 in total, a median of one per cell.**

`diagnose_loss_structure.m` Part B reproduces each map's *whole* coverage profile in V1 —
`subsample_cells.m` thins every cell to that map's actual vertex count, holes and sparsity together
— and reports the spread of the estimate over 200 draws. At 4–8°:

| map | empty | median/cell | vertices | observers | `dg` rad−tang: mean [5th, 95th] | band width |
|---|--:|--:|--:|--:|---|--:|
| V3 | 1 | 80 | 5160 | 8 | 0.124 [0.114, 0.132] | 0.018 |
| V3a | 4 | 18.5 | 2005 | 8 | 0.126 [0.103, 0.147] | 0.044 |
| hV4 | 24 | 4 | 1031 | 8 | 0.138 [0.084, 0.190] | 0.106 |
| **MT** | 30 | 1 | 545 | **7** | **0.187 [0.124, 0.255]** | **0.131** |

V1's own value is 0.119, from 11075 vertices. **With MT's coverage the 90% band is wider than the
effect being measured**, and the mean is biased up by more than half. V3's band is 0.018. That is
where the §7 criterion separates them, and it is separating them on the right thing.

Two reasons this is a **lower bound** on MT's real problem:

- The subsample keeps V1's own responses, and V1's vertices have better pRF fits and larger BOLD
  than an extrastriate map's. It reproduces the count, not the reliability.
- MT at 4–8° has one observer with **zero** surviving vertices, so under the NaN convention the
  group estimate does not exist at all. That is why MT's 4–8° rows in
  `supplement/spec_areas_asymmetries_spec.csv` are NaN rather than numbers.

At 2–10°, where MT does produce numbers, the band is still [0.097, 0.199] against V1's 0.124 — and
MT still fails on all three of empty cells (19 > 2), median vertices (9 < 20) and weight ratio
(48 > 25).

**A small caveat on reading this table.** The subsample reproduces each map's polar-angle *density
profile*, not just its total, so the mean can shift a little even where the band is tight — V3
coverage gives horiz−vert −0.563 against V1's −0.548. The band width is the quantity to read.

### What the exclusion does and does not say

It says **this design cannot resolve MT by polar angle**. It is not a claim that MT has no
orientation asymmetries, and it says nothing about whole-ROI analyses of MT, which remain available
and unaffected — `local_qc/group_addv4mt.m` already does one, and both once-flagged observers were
cleared partly on MT motion selectivity ([`local_qc/DATA_QUALITY.md`](local_qc/DATA_QUALITY.md) §2).
The excluded numbers stay in `supplement/spec_areas_*.csv` behind a `reportable` flag;
`plot_spec_hierarchy.m` filters on that flag, which is why MT is absent from
`Figure_S5_spec_hierarchy`.

---

## 7. What is still open

Two things, both listed with everything else in [`../AGENTS.md`](../AGENTS.md) §5: calibrating the
[`SPECIFICATION.md`](SPECIFICATION.md) §7 thresholds against the band width §6 measures directly, and
testing whether `equalcell` weighting (§5, last paragraph) helps in the maps that are actually
reported. Nothing currently reported depends on either.

---

## 8. How to reproduce

From `cleanroom/` in MATLAB. None of this is part of `run_spec_outputs`; it is diagnostic and is run
on demand.

| command | what it does |
|---|---|
| `cell_occupancy('area','MT')` | the 8 × 8 per-cell vertex table for any map |
| `diagnose_cell_loss('donor','MT','eccRange',[4 8])` | §2 and §4 — delete a map's empty cells from V1, both routes, with a random-deletion null and the class test |
| `diagnose_loss_structure('donor','MT')` | §3 and §6 — the randomisation, and the coverage-density test across maps |
| `diagnose_pooled_fit('donor','MT')` | §5 — pooled versus per-observer, the leverage table, the severity sweep |
| `spec_pooled` | the pooled estimator itself; returns both estimators side by side |
| `subsample_cells` | thins a loaded dataset to a target count per cell |

All of them need `~/dg_collect/`, including `gain_areas_summary.csv` — a copy is tracked at
`supplement/gain_areas_summary.csv` and `observer_gain_weights` looks for it next to
`gainSummary.csv`. Without it the run falls back to the V1-derived scalar and every number sits
about 1% high.

`spec_profiles`, `diagnose_within_observer_error` and `load_runbetas_area` take an additive
`dropCells` argument for this work. It defaults to empty, empty means the settled specification
unchanged, and the V1 baseline reproduces `supplement/spec_asymmetries_spec_v1_4-8.csv` exactly with
the option in place.

**Guards.** The per-observer fits recovered from the pooled path match `spec_profiles` to 6e-15; the
two pooled weightings are bit-identical on complete data; the leverage columns sum to 1; and
`diagnose_pooled_fit`'s `avg` column reproduces `diagnose_cell_loss`'s harmonic arm exactly from an
independent in-memory code path.
