# Reproduction findings

Running log of what the reproduction has surfaced. Numbers are group-mean differences
(z-scored, sigma units), averaged over the 8 polar-angle wedges then over 8 subjects, from the
clean-room (Path A) starting from `Support/allsubjectsTable.csv`.

## 1. Seven of eight asymmetries reproduce the manuscript almost exactly

Filter: V1, `4<=ecc<=8`, `pRF_r2>0.1` → **11,075 vertices** (matches the CSV's own strict
inclusion count). Z-scored contrast = `(orientation - blank)/beta_std` per vertex.

| asymmetry | dg computed | dg manuscript | da computed | da manuscript |
|---|---|---|---|---|
| horizontal − vertical | **−1.155** | −1.155 | −0.041 ⚠ | −0.45 |
| cardinal − oblique | −0.398 | −0.40 | −0.064 | −0.06 |
| radial − tangential | 0.226 | 0.23 | 0.603 | 0.60 |
| polar-card − polar-obl | 0.056 | 0.06 | 0.173 | 0.17 |

All direct asymmetries and **all dg derived (PA-dependent) asymmetries** match. This validates
the whole data path: loading, the inclusion filter, per-vertex z-scoring, wedge binning by
`pRF_angle_bin`, median aggregation, and the local-orientation / frame mapping.

## 2. Open discrepancy: `da` horizontal − vertical (Fig 6A)

The one mismatch is the **derived** (polar-stimulus → Cartesian-frame) H-V asymmetry for the
polar experiment: clean-room **−0.04** (asymmetry essentially absent) vs manuscript **−0.45**.

**Mechanism (established, not yet definitively attributed):** the derived asymmetries depend on
the per-wedge reference polar angle. The repo's `compute_derivativeDirections.m` /
`deriveLocalMotionfromUVM.m` use the PA order `[90 45 0 315 270 225 180 135]`, whereas
`meanWithinLabel.m` stores the wedge dimension in the order `[0 45 90 135 180 225 270 315]`.
These disagree at the cardinal meridians.

Evidence from the clean-room, toggling only the reference-angle ordering:

| quantity | correct ref `[0..315]` | repo ref `[90 45 0 315 …]` | manuscript |
|---|---|---|---|
| dg radial − tangential | **0.226** | 0.003 | 0.23 |
| da horizontal − vertical | −0.041 | **0.446** (|.| matches) | −0.45 |

The per-θ profiles show why: `dg radTan` is strongly **+** at the vertical meridians and **−** at
the horizontal meridians; averaging preserves +0.23 under the correct order but cancels under
the repo order. `da HV` is dominated by the radial-tangential effect projected onto the H/V
axis, which flips sign between horizontal and vertical meridians, so it **cancels to ~0** under
the correct order but does **not** cancel under the repo order.

So the manuscript's `dg radTan` (0.23) tracks the **correct** convention while its `da HV`
(−0.45) tracks the **repo PA-ordering** convention — an internal inconsistency between two
derived quantities. (The clean-room hand-replication of the repo order gives +0.45, matching the
manuscript's magnitude but not its sign, so the exact repo indexing is not yet fully mirrored.)

**Implication if confirmed:** the reported `da` horizontal-vs-vertical asymmetry (Fig 6A, and its
weight in the Fig 7 LME) may be an artifact of the polar-angle-ordering in the derived-direction
code rather than a genuine effect; the physically-correct value is ~0 (Cartesian asymmetry fully
absent for polar gratings, which would *strengthen* the paper's thesis).

**To resolve definitively (Path B):** run the actual `compute_derivativeDirections.m` on the
bridged group arrays and compare **per-θ** against the clean-room. That will show exactly which
wedge-to-reference pairing the real pipeline uses and settle whether −0.45 is real or an ordering
artifact.

## 3. Corroboration from the manuscript's own Fig 7 (LME)

The joint LME (Fig 7) and the independent analysis (Figs 5/6) use the *same* 256 data points, so
they should agree in direction. They do for every asymmetry **except da H-V**:

| da asymmetry | manuscript Fig 6 (independent) | manuscript Fig 7 (LME) | clean-room (both) |
|---|---|---|---|
| horizontal − vertical | **−0.45** | **0.02** | −0.04 |
| cardinal − oblique | −0.06 | 0.06 | −0.06 |
| radial − tangential | 0.60 | 0.60 | 0.60 |
| polar-card − polar-obl | 0.17 | 0.18 | 0.17 |

The manuscript's **own LME reports da H-V ≈ 0.02** (absent) — agreeing with the clean-room (~0)
and contradicting its own independent Fig 6A value of −0.45. So the two manuscript analyses are
internally inconsistent for exactly this quantity, and the clean-room sides with the LME. This
strengthens the read that the independent Fig 6A −0.45 is an artifact of the derived-direction
PA-ordering, not a real effect.

(Aside: the clean-room LME Δ reproduces the independent means *exactly* because the balanced
design makes the four predictors orthogonal. The manuscript's LME weights are ~5–15% attenuated
vs its own independent estimates — the signature of slightly correlated predictors, again
consistent with a misaligned derived-direction coding.)

## Status
- [x] Path A data path validated (7/8 exact)
- [x] Discrepancy localized + corroborated by manuscript's own Fig 6-vs-Fig 7 inconsistency
- [ ] `da` H-V resolved via Path B per-θ comparison (run actual repo code)
