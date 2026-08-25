# Stimulus and angle conventions

**What each stimulus was, what each condition index means, and which polar-angle frame each part of
the pipeline uses.** Verified end to end — experiment code → design matrix → GLM → CSV → figures —
against sub-0255's raw files at machine precision.

**Verdict: the `AnalysisCode` pipeline is correct on all eight asymmetries.** A polar-angle-ordering
bug once reported against `compute_derivativeDirections.m` was **not a bug**; it came from two errors
inside this reproduction, both since fixed (§5).

Six conventions change along the chain (§4). Every one is handled correctly in `AnalysisCode`, and
every one is a place a reimplementation can go wrong.

---

## 1. What was on the screen

`ExperimentCode/Config/trialtypes.csv` codes each trial as `mainCondition` (0 = blank, 1 = static,
2 = motion) plus an `OrientationDir` id ∈ {0, 45, 90, 135}.

**Cartesian experiment (`dg`).** The base texture is a **vertical** grating
(`Config/constConfig.m:50-55`), rotated by `mod(90 − id, 180)` in PsychToolbox screen coordinates:

| id | PTB rotation | rendered bars | visual-field orientation |
|---|---|---|---|
| 0 | 90° | horizontal | 0° |
| 45 | 45° | `/` | 45° |
| 90 | 0° | vertical | 90° |
| 135 | 135° | `\` | 135° |

So **the id *is* the grating orientation in conventional visual-field degrees**; the `mod(90−id,180)`
map is exactly the screen y-flip.

**Polar experiment (`da`).** `Stim/my_stim.m:89-113` switches the same ids onto polar gratings:

| id | shader setting | stimulus | local orientation |
|---|---|---|---|
| 90 | `circularFrequency = 0` | pinwheel | **radial** |
| 0 | `radialFrequency = 0` | annulus | **tangential** |
| 45 | both /√2, +circ | spiral | oblique |
| 135 | both /√2, −circ | spiral | oblique |

At the **upper vertical meridian** radial = vertical = 90° and tangential = horizontal = 0°, so the
`dg` and `da` ids coincide exactly there. **The shared id is the stimulus's local orientation at the
upper vertical meridian.** That is the organising principle of the whole design, and why the
derivation routine is named `deriveLocalMotionfromUVM.m`.

## 2. Condition index ↔ stimulus

`createStimMap.m` gives canonical names; `general/jsons/CONTRASTS.json` fixes the index order.
Indices 26–29 are the four stationary orientations:

| idx | `dg_contrast_name` | `da_contrast_name` | UVM-local orientation |
|---|---|---|---|
| 26 | `s0_v_b` | `sannulus_v_b` | 0° |
| 27 | `s90_v_b` | `spinwheel_v_b` | 90° |
| 28 | `s45_v_b` | `scspiral_v_b` | **45°** |
| 29 | `s135_v_b` | `sccspiral_v_b` | **135°** |

The `dg` and `da` names share an index because they share the UVM-local orientation — which is why
`meanWithinLabel.m:36` deliberately uses the `dg` names for both experiments. Note in particular
**cspiral = 45°, ccspiral = 135°**, settled directly by the data: sub-0255's
`polexp_cspiral_stationary` column matches `betamaps` id 45 to 5.3e-15 and id 135 not at all.

## 3. The polar-angle frame — the crux

Two conventions are in play:

- **Benson**: 0° = upper vertical meridian, increasing clockwise, signed ±180 (positive = right
  visual field).
- **Conventional**: 0° = right horizontal meridian, counter-clockwise.

`general/map_theta.m` converts between them — `conventional = mod(90 − benson, 360)` — and is its own
inverse. The two consumers of the pRF angle differ, **deliberately**:

| consumer | source | `map_theta` applied? | resulting frame |
|---|---|---|---|
| `createTables.m:75` → `allsubjectsTable.csv` | `*.angle_adj.mgz` | **yes** | conventional |
| `meanWithinLabel.m:118` → `meanBOLDpa` / `medianBOLDpa` | `*.angle_adj.mgz` | **no** | Benson |

`meanWithinLabel.m:73` bins with `polarAngles = [0, 45, 90, 135, 180, −135, −90, −45]` in **Benson**
degrees. Mapping that index order to conventional degrees gives

```
[90, 45, 0, 315, 270, 225, 180, 135]
```

which is *verbatim* what three downstream files independently assume —
`compute_derivativeDirections.m:48`, `plot1_experimentalCond.m:121` (commented *"manually converted
based on the order of polarAngles above (Noah's convention)"*) and `lme1_fit.m:88`. **The original
pipeline is self-consistent across four files.** `deriveLocalMotionfromUVM.m` then rotates the
UVM-referenced orientation by `(θ − 90)`, which is the correct physics, and its `ismember`
reindexing is an involutive permutation.

> ⚠️ **Do not "align" the wedge dimension of `medianBOLDpa`/`meanBOLDpa` with `[0 45 90 …]`.** That
> would introduce the cardinal-meridian swap into working code and corrupt both first-harmonic
> derived asymmetries.

**Three independent confirmations that the CSV is conventional and the wedge arrays Benson:**

1. **Cortical magnification.** Vertices per CSV bin (V1, 4–8°, R² > 0.1): 0→2882, 45→1412, 90→**536**,
   135→1077, 180→2195, 225→1112, 270→841, 315→1020. Maxima at 0/180, global minimum at 90, and
   90 < 270 — the established V1 surface-area profile. So bins 0/180 are the horizontal meridian.
2. **Radial bias.** The 2nd harmonic of `dg` (horizontal − vertical) across bins peaks at phase
   ≈ **−10°**, i.e. at the horizontal meridian, which is only possible if bin 0 is the HM.
3. **Anatomy.** In left V2d — which represents the *right lower* field — 95.6% of vertices have CSV
   `pRF_angle` ∈ 270–360° (circular mean 309°); only 0.1% fall in 90–180°, where they would sit under
   Benson.

## 4. The six places a convention changes

1. **Orientation ids and motion ids measure different things.** `OrientationDir` id *is* the grating
   orientation (0 = horizontal); `MotionDir` id is the *direction of motion* (0 = rightward), whose
   grating is perpendicular. Hence design column 1 is `M_0` = *vertical* grating moving rightward,
   while column 9 is `S_0` = *horizontal* grating. Handled by `format_desmats.m` reading
   `trialMat(:,3)` for static and `(:,4)` for motion. Benign downstream: "distance 0 from the polar
   angle = radial" is true for a radial *orientation* and for radial *motion* alike.
2. **Screen vs visual-field angles.** The PsychToolbox y-down flip, confined to `ExperimentCode`;
   never propagates.
3. **The same id is a different stimulus in `dg` vs `da`.** id 0 = horizontal grating but *annulus*;
   id 90 = vertical but *pinwheel*. Both `results.mat` files use the *same* field names (`s0Vb`,
   `s90Vb`, …), i.e. the GLM names live in **id space**, not semantic space.
4. **id space → semantic names** happens once, at `createTables.m:126-151`, and differently per
   experiment. This is where `s0Vb` becomes `cartexp_horizontal_stationary` *or*
   `polexp_annulus_grating_stationary`.
5. **Benson → conventional polar angle, on one branch only** (§3). **This is the trap.**
6. **z-scoring.** CSV betas are non-z-scored but ship `beta_std`, so the z-scored contrast is
   recoverable; the saved `results.mat` contrasts (Nov 2024) are already z-scored — verified for
   sub-0255 to 8.9e-16. `main_singlesub.m` now has `normalize = 0`, so re-running the pipeline today
   would not reproduce those saved files. That is correct rather than a trap, since z-scoring was
   decided against ([`METHOD_DECISIONS.md`](METHOD_DECISIONS.md) §1) — but the flag is still a
   hand-edited constant and would be better as a recorded parameter.

## 5. The retracted bug report

The `FINDINGS.md` report of a polar-angle-ordering bug in `compute_derivativeDirections.m` is **fully
retracted**. Two errors *inside this reproduction*, each flipping a different half of the wedges,
produced it:

- **The spirals were swapped** in `cleanroom/config_repro.m` (ccspiral→45°, cspiral→135°; §2 has it
  the other way). At the four oblique wedges the spirals supply the locally-horizontal and
  locally-vertical stimuli, so this flipped the sign of `da` (H−V) there.
- **The bridge fed the wrong frame** — `bridge/build_group_matrices_fromCSV.m` built `medianBOLDpa`
  in *conventional* order and handed it to the repo's Benson-expecting functions. That relation is a
  reflection about 45°, so it flipped exactly the four cardinal wedges.

Combined they flip all eight wedges, turning `da` H−V from −0.446 into +0.446 — the unexplained sign
flip the original report noticed and set aside as a "plotting-sign convention". Both are fixed;
neither ever existed in `AnalysisCode`. Recomputing all eight asymmetries directly from the CSV with
the §2 semantics reproduces the draft's values.

## 6. End-to-end verification (sub-0255)

**Experiment → design matrix.** `S09_const_file_Run1.mat` stores the `containers.Map` objects
actually used at scan time, for both experiments, identical to `constConfig.m`. `expDes.trialMat` has
52 trials = 13 conditions × 4 repeats.

**Design matrix → GLM.** Rebuilding the trial→condition sequence from the design files alone, using
`format_desmats.m`'s labelling rule, reproduces GLMsingle's `stimorder` exactly: 468 trials (9 runs)
for `dg`, 416 (8 runs) for `da`, **36 / 32 trials per condition, perfectly balanced**. Averaging
`results.allevents.modelmd` over those trials reconstructs `betamaps`.

**GLM → CSV.** Each reconstructed `betamaps` column equals its `allsubjectsTable.csv` column to
≤ 5.3e-15 (CSV text round-trip precision), for all 13 conditions of both experiments, plus
`dg_beta_std` / `da_beta_std`.

## 7. The one defensive change worth making

`meanWithinLabel.m` is the only consumer of `angle_adj.mgz` that does not call `map_theta`, and the
Benson ordering it produces is re-derived by hand in three separate downstream files. Storing the
wedge centres alongside the data — or converting once at the source and updating the three
`[90,45,0,…]` literals together — would remove the trap. That is a refactor; the current code's
*output* is correct, so it is optional.
