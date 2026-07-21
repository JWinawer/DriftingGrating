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

## 4. RESOLVED (Path B): the −0.45 is a cardinal-meridian label swap

Running the **actual repo function** `compute_derivativeDirections.m` on the bridged da
`medianBOLDpa` (Path B, `bridge/resolve_da_HV.m`) gives, per polar angle:

```
              0°     45°    90°    135°   180°   225°   270°   315°
repo  code : -0.169  0.236  0.876  0.522  0.011  0.321  1.23   0.54
clean-room:   0.169  0.236 -0.876  0.522 -0.011  0.321 -1.23   0.54
```

**Identical at the oblique meridians (45/135/225/315); exactly sign-flipped at the four
cardinal meridians (0/90/180/270).** That is the signature of the repo routine **swapping the
horizontal and vertical labels at the cardinal locations**, caused by its polar-angle reference
order `[90 45 0 315 270 225 180 135]` not matching the data's `[0 45 90 135 180 225 270 315]`.

Group means: repo code **+0.446**, clean-room **−0.041**, manuscript Fig 6A **−0.45**, manuscript
Fig 7B LME **0.02**. The manuscript's |0.45| is reproduced *only* by the repo's derived-direction
routine; it is an **ordering artifact**. The physically-correct value is ~0 (Cartesian H-V
asymmetry absent for polar gratings), which the manuscript's own Fig 7 LME already reports.
(The repo yields +0.45 vs Fig 6A's −0.45 — a pro/con plotting-sign convention, immaterial to the
artifact.)

**Scope of the artifact:** it affects the *derived* (reference-frame-mismatched) asymmetries —
da H-V/cardinal-oblique and dg radial-tangential/polar-cardinal — because those go through
`compute_derivativeDirections.m`. It does **not** touch the direct asymmetries (dg H-V/card-obl,
da rad-tang/polcard). Whether the reported dg radial-tangential (0.23) and da cardinal-oblique
are also affected should be checked the same way; the clean-room values are the correct reference.

## 5. Full Path B sweep — all 8 asymmetries through the real repo code

`bridge/run_pathB_values.m` runs the actual repo functions (`retrieveProConIdx.m` for direct
asymmetries, `compute_derivativeDirections.m` for derived) on the bridged group array;
`bridge/run_pathB_figures.m` regenerates Figs 5/6 through the real `plot1/plot2` plotting
functions. **Cross-path check:** the bridged `medianBOLDpa` orientation rows equal the
clean-room array to `0.000e+00` over every orientation/wedge/subject/variant — so all
differences below come purely from the existing code, not the bridge.

z-scored group means (equal-weight over wedges):

| experiment | asymmetry | kind | existing code | clean-room (correct) | manuscript |
|---|---|---|---|---|---|
| dg | horiz − vert | direct  | −1.155 | −1.155 | −1.155 |
| dg | card − obl   | direct  | −0.398 | −0.398 | −0.40 |
| dg | **rad − tang** | **derived** | **0.003** ⚠ | **0.226** | **0.23** |
| dg | polc − polo  | derived | 0.056  | 0.056  | 0.06 |
| da | **horiz − vert** | **derived** | **0.446** ⚠ | **−0.041** | **−0.45** |
| da | card − obl   | derived | −0.064 | −0.064 | −0.06 |
| da | rad − tang   | direct  | 0.603  | 0.603  | 0.60 |
| da | polc − polo  | direct  | 0.173  | 0.173  | 0.17 |

**What this shows.** The current `compute_derivativeDirections.m` is wrong for **both
first-harmonic derived asymmetries** — it collapses `dg` radial-tangential to ~0 (should be
0.23) and inflates `da` horizontal-vertical to |0.45| (should be ~0). The second-harmonic
derived asymmetries (`dg` polar-card, `da` card-obl) are invariant to the cardinal-meridian
swap and come out correct.

**The manuscript is internally inconsistent on the two first-harmonic derived values:** its
Fig 5C `dg` radial-tangential (0.23) matches the **correct** convention (so Fig 5C was *not*
produced by the current code), while its Fig 6A `da` horizontal-vertical (−0.45) matches the
**current buggy** code. The two figures were generated with different derivations.

**Corrected conclusions (trust the clean-room):**
- `da` horizontal-vertical ≈ 0 → Cartesian asymmetry **absent** for polar gratings (strengthens
  the paper's thesis; the reported −0.45 is spurious).
- `dg` radial-tangential ≈ 0.23 (the manuscript's Fig 5C value is already correct); note that
  re-running the *current* code would wrongly report ~0.
- Both are 1st-harmonic derived quantities that pass through `compute_derivativeDirections.m`;
  fixing the polar-angle reference order there would make the code reproduce the clean-room
  values for every asymmetry.

## Status
- [x] Path A data path validated (7/8 exact)
- [x] Discrepancy localized + corroborated by manuscript's own Fig 6-vs-Fig 7 inconsistency
- [x] **RESOLVED via Path B**: the first-harmonic *derived* asymmetries (`da` H-V, `dg` rad-tang)
      carry a cardinal-meridian label swap (PA-ordering artifact) in
      `compute_derivativeDirections.m`; correct values are ~0 and 0.23 respectively
- [x] Cross-path validated (bridged vs clean-room group arrays identical); Figs 5/6 regenerated
      through the real `plot1/plot2`
- [ ] Suggested fix: align the polar-angle reference order in `compute_derivativeDirections.m`
      (and `deriveLocalMotionfromUVM.m`) with `meanWithinLabel.m`'s `[0 45 90 …]` ordering
