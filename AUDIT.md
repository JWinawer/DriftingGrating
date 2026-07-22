# Stimulus-design audit: experiment code → GLM → CSV → figures

**Date:** 2026-07-22 · **Scope:** Figures 5–8 · **Verdict: the manuscript is correct on all
eight asymmetries. The original `AnalysisCode` pipeline is correct. The
polar-angle-ordering "bug" reported in `Reproduction/FINDINGS.md` is not a bug — it was
produced by two errors inside the reproduction itself.**

This audit traces the orientation/polar-angle conventions along the whole chain. §§1–5 reach
the verdict from repo code plus `Support/allsubjectsTable.csv` alone. **§6 then confirms every
link at machine precision against sub-0255's raw experiment files and GLM output**, including
the two facts §§1–5 could only infer: the spiral identities and the Benson polar-angle frame.
§7 lists every convention change along the chain; §8 the remaining open items.

---

## 1. Ground truth: what was actually put on the screen

`ExperimentCode/Config/trialtypes.csv` codes each trial as `mainCondition`
(0=blank, 1=static, 2=motion) plus an `OrientationDir` id ∈ {0,45,90,135}.

**Cartesian experiment (`dg`).** The base texture is built at
`Config/constConfig.m:50-55` as `repmat(signal,[N,1])` with `signal` varying along *x* —
a **vertical** grating. `constConfig.m:269-270` maps id → PsychToolbox rotation:

```matlab
orientationids = 0:45:315; ptborientation = {90, 45, 0, 135, ...};   % PTB = mod(90-id,180)
```

`Screen('DrawTexture', ..., angle)` rotates clockwise in screen coordinates (y down).
Rotating a vertical grating by the mapped angle gives:

| id | PTB rotation | rendered bars | visual-field orientation |
|---|---|---|---|
| 0   | 90° | horizontal | 0° |
| 45  | 45° | `/`        | 45° |
| 90  | 0°  | vertical   | 90° |
| 135 | 135°| `\`        | 135° |

So **the id *is* the grating orientation in conventional visual-field degrees**; the
`mod(90-id,180)` map is exactly the screen y-flip.

**Polar experiment (`da`).** `Stim/my_stim.m:89-113` switches the same ids onto polar gratings:

| id | shader setting | stimulus | local orientation |
|---|---|---|---|
| 90  | `circularFrequency = 0` | pinwheel | **radial** |
| 0   | `radialFrequency = 0`   | annulus  | **tangential** |
| 45  | both /√2, +circ         | spiral   | oblique |
| 135 | both /√2, −circ         | spiral   | oblique |

At the **upper vertical meridian** radial = vertical = 90° and tangential = horizontal = 0°.
So the `dg` and `da` ids coincide precisely at the UVM: **the shared id is the stimulus's
local orientation at the upper vertical meridian.** This is the organizing principle of the
whole design, and it is why the derivation routine is named `deriveLocalMotionfromUVM.m`.

## 2. GLM stage: condition index ↔ stimulus

`AnalysisCode/01_process_singlesubjectGLM/createStimMap.m` gives canonical names and
`ori` / `oriOffset` (offset from radial) values; `general/jsons/CONTRASTS.json` fixes the
index order. Indices 26–29 are the four stationary orientations:

| idx | `dg_contrast_name` | `da_contrast_name` | UVM-local orientation |
|---|---|---|---|
| 26 | `s0_v_b`   | `sannulus_v_b`  | 0° |
| 27 | `s90_v_b`  | `spinwheel_v_b` | 90° |
| 28 | `s45_v_b`  | `scspiral_v_b`  | **45°** |
| 29 | `s135_v_b` | `sccspiral_v_b` | **135°** |

The dg and da names occupy a *shared* index because they share the UVM-local orientation —
which is exactly why `meanWithinLabel.m:36` deliberately uses the dg names for both
experiments ("they apply to both DG and DA"). Note in particular **cspiral = 45°,
ccspiral = 135°**.

## 3. The polar-angle frame — the crux

Two conventions are in play:

- **Benson** (`Noah's convention`): 0° = upper vertical meridian, increasing clockwise,
  signed ±180 (positive = right visual field).
- **Conventional**: 0° = right horizontal meridian, counter-clockwise.

`general/map_theta.m` converts between them: `conventional = mod(90 − benson, 360)`. It is
an **involution** (its own inverse).

The two consumers of the pRF angle differ, deliberately:

| consumer | source | `map_theta` applied? | resulting frame |
|---|---|---|---|
| `01_.../createTables.m:75` → `allsubjectsTable.csv` | `*.angle_adj.mgz` | **yes** | conventional |
| `03_.../meanWithinLabel.m:118` → `medianBOLDpa` | `*.angle_adj.mgz` | **no** | Benson |

`meanWithinLabel.m:73` bins with `polarAngles = [0, 45, 90, 135, 180, -135, -90, -45]` in
**Benson** degrees. Mapping that index order to conventional degrees gives

```
[90, 45, 0, 315, 270, 225, 180, 135]
```

which is *verbatim* what three downstream files independently assume:

- `04_.../compute_derivativeDirections.m:48` — `polarAngles = [90, 45, 0, 315, ...]`
- `04_.../plot1_experimentalCond.m:121` — `anglevals = [90, 45, 0, 315, ...]`, commented
  *"these were manually converted based on the order of polarAngles above (Noah's convention)"*
- `04_.../lme1_fit.m:88` — same array

**The original pipeline is self-consistent across four files.** `deriveLocalMotionfromUVM.m`
then rotates the UVM-referenced orientation by `(θ − 90)`, which is the correct physics, and
its `ismember` reindexing is an involutive permutation, so it is applied correctly.

### Empirical confirmation that the CSV is in the conventional frame

Two independent checks on `Support/allsubjectsTable.csv` (V1, 4≤ecc≤8, R²>0.1; 11,075 vertices):

1. **Cortical magnification.** Vertices per bin: 0→2882, 45→1412, 90→**536**, 135→1077,
   180→2195, 225→1112, 270→841, 315→1020. Maxima at bins 0/180, global minimum at 90, and
   bin 90 < bin 270. That is the established V1 surface-area profile — most area at the
   horizontal meridian, least at the vertical, less at the upper than lower vertical meridian.
   So bins 0/180 = horizontal meridian.
2. **Radial bias.** Fitting the 2nd harmonic of dg (horizontal − vertical) across bins gives
   amplitude 0.23 at phase ≈ **−10°** — i.e. (H−V) peaks at the horizontal meridian, as the
   radial bias requires, and only possible if bin 0 is the HM.

Both confirm: **CSV `pRF_angle` is conventional**, hence `angle_adj.mgz` is Benson, hence the
`medianBOLDpa` wedge dimension is Benson-ordered.

## 4. Independent recomputation reproduces the manuscript exactly

Recomputing all eight asymmetries directly from the CSV — per subject × wedge median of the
z-scored `(orientation − blank)/beta_std`, equal-weighted over the 8 wedges, then averaged
over the 8 subjects — using the §2 stimulus semantics:

| asymmetry | recomputed | manuscript |
|---|---|---|
| dg horizontal − vertical | −1.155 | −1.16 |
| dg cardinal − oblique | −0.398 | −0.40 |
| dg radial − tangential | +0.226 | +0.23 |
| dg polar-card − polar-obl | +0.056 | +0.06 |
| **da horizontal − vertical** | **−0.446** | **−0.45** |
| da cardinal − oblique | −0.064 | −0.06 |
| da radial − tangential | +0.603 | +0.60 |
| da polar-card − polar-obl | +0.173 | +0.17 |

All eight match, **including the disputed `da` horizontal−vertical**. This recomputation
lives entirely in the CSV/conventional frame, so it is independent of how one resolves the
Benson question — the manuscript's numbers are the physically correct ones either way.

## 5. What went wrong in `Reproduction/`

Two independent errors, each flipping a different half of the wedges:

**Error 1 — the spirals are swapped.** `cleanroom/config_repro.m:74-80` declares
`oriNames = {'pinwheel','annulus','ccspiral','cspiral'}` with `oriAngle = [90 0 45 135]`,
i.e. ccspiral→45°, cspiral→135°. Per `CONTRASTS.json` (§2) it is the reverse. At the four
**oblique** wedges the spirals are what supply the locally-horizontal and locally-vertical
stimuli, so this flips the sign of `da` (H−V) at 45°/135°/225°/315°:

```
theta          0      45      90     135     180     225     270     315     mean
correct    +0.169  -0.236  -0.876  -0.522  -0.011  -0.321  -1.230  -0.540  -0.446
FINDINGS   +0.169  +0.236  -0.876  +0.522  -0.011  +0.321  -1.230  +0.540  -0.041
```

Re-introducing only this swap reproduces the reproduction's −0.041 exactly.

**Error 2 — the bridge feeds the wrong frame.** `cleanroom/config_repro.m:31` sets
`cfg.paBins = [0 45 90 135 180 225 270 315]` (conventional) and comments *"Matches the data
order used by the original meanWithinLabel.m"* — it does not; `meanWithinLabel` is
Benson-ordered (§3). `bridge/build_group_matrices_fromCSV.m:16` builds `medianBOLDpa` in
that conventional order and hands it to the repo's Benson-expecting functions. The
conventional↔Benson relation is a reflection about 45°, identity at 45°/225° and swapping
0↔90 and 180↔270, so it flips exactly the four **cardinal** wedges. This is what collapsed
`dg` radial−tangential from 0.226 to 0.003, and it is the "identical at the obliques,
sign-flipped at the cardinals" signature reported in `FINDINGS.md` §4.

**Combined**, errors 1 and 2 flip all eight wedges, turning −0.446 into **+0.446** — the
reproduction's "repo code" value, whose unexplained sign flip `FINDINGS.md` §2 noticed and
set aside as a "pro/con plotting-sign convention."

Both errors originate in the reproduction. Neither exists in `AnalysisCode`.

## 6. End-to-end verification against sub-0255 raw files

Raw files for **sub-0255** (`Support/sub-0255/{dg,da}/`) confirm every link at machine
precision. Scripts in the session scratchpad; summary:

**Experiment → design matrix.** `S09_const_file_Run1.mat` stores the `containers.Map` objects
actually used at scan time, for both experiments: `maporientation` = `{0→90, 45→45, 90→0,
135→135}` and `mapdirection` = `{0→180, 45→135, 90→90, …}` — identical to `constConfig.m`.
`expDes.trialMat` has 52 trials = 13 conditions × 4 repeats, `OrientationDir` ∈ {0,45,90,135}
on static trials and `MotionDir` ∈ {0…315} on motion trials.

**Design matrix → GLM.** Rebuilding the trial→condition sequence from the design files alone,
using `format_desmats.m`'s labelling rule, reproduces GLMsingle's `stimorder` exactly:
468 trials (9 runs) for `dg`, 416 (8 runs) for `da`, **36 / 32 trials per condition, perfectly
balanced**. Averaging `results.allevents.modelmd` over those trials reconstructs `betamaps`.

**GLM → CSV.** Each reconstructed `betamaps` column equals its `allsubjectsTable.csv` column
to ≤5.3e-15 (CSV text round-trip precision), for all 13 conditions of both experiments, plus
`dg_beta_std` / `da_beta_std`. The **spiral assignment is settled directly by the data**:

```
CSV polexp_cspiral_stationary  vs betamaps S_45  (id 45)  : 5.3e-15   MATCH
CSV polexp_cspiral_stationary  vs betamaps S_135 (id 135) : 5.3e+00
CSV polexp_ccspiral_stationary vs betamaps S_135 (id 135) : 5.1e-15   MATCH
CSV polexp_ccspiral_stationary vs betamaps S_45  (id 45)  : 5.3e+00
```

So **cspiral = id 45, ccspiral = id 135**, confirming `CONTRASTS.json` and `createTables.m`
and refuting `config_repro.m` (§5, Error 1) from the raw data.

**Polar-angle frame, confirmed anatomically.** In `summaryTables_wleftV2d/`, left V2d — which
represents the *right lower* visual field — has 95.6% of its vertices at `pRF_angle` ∈
270–360° (circular mean 309°); only 0.1% fall in 90–180°, where they would sit under Benson.
Together with `createTables.m:74` (*"angle converted from Benson coordinates to 0 to 360
(cc from right horizontal)"*) this **confirms the one assumption §3 could not check in-repo**:
`angle_adj.mgz` is Benson, the CSV is conventional, and `medianBOLDpa` is Benson-ordered.
(Independently reproduced and confirmed by JW, 2026-07-22.)

**The saved GLM results are z-scored.** For sub-0255 the stored contrasts satisfy
`results.contrasts.s0Vb == (betamaps(:,9) − betamaps(:,13)) ./ std(betamaps,0,2)` to 8.9e-16,
and are *not* equal to the raw difference (off by ~4). So these `results.mat` files
(Nov 2024) were produced with `normalize = 1`. `main_singlesub.m:131` now reads
`normalize = 0; % changed to 0 on 7-14-2026`, so **re-running the pipeline today would not
reproduce the saved files** — see §8. It also confirms the CSV route is exactly equivalent:
`(orientation − blank)/beta_std` from the CSV *is* the pipeline's z-scored contrast.

## 7. Convention changes along the chain

Six places where a convention changes. All are handled correctly by `AnalysisCode`, but each
is a place a reimplementation can go wrong — and two of them are where the reproduction did.

1. **Orientation ids and motion ids measure different things.** `OrientationDir` id *is* the
   grating orientation (0 = horizontal); `MotionDir` id is the *direction of motion*
   (0 = rightward), whose grating is perpendicular. Hence design column 1 is `M_0` =
   *vertical* grating moving rightward, while column 9 is `S_0` = *horizontal* grating.
   Handled by `format_desmats.m` reading `trialMat(:,3)` for static and `(:,4)` for motion.
   Downstream this is benign: "distance 0 from the polar angle = radial" is true both for a
   radial *orientation* and for radial *motion*.
2. **Screen vs visual-field angles.** `maporientation`/`mapdirection` apply
   `mod(90−id,180)` / `mod(180−id,360)`, the PsychToolbox y-down flip. Confined to
   `ExperimentCode`; never propagates.
3. **The same id is a different stimulus in `dg` vs `da`.** id 0 = horizontal grating but
   *annulus*; id 90 = vertical but *pinwheel*. They coincide only at the **upper vertical
   meridian**. Both `results.mat` files use the *same* field names (`s0Vb`, `s90Vb`, …), i.e.
   the GLM names live in **id space**, not semantic space — which is exactly why
   `meanWithinLabel.m:36` uses the dg names for both experiments.
4. **id space → semantic names** happens once, at `createTables.m:126-151`, and differently
   per experiment. This is where `s0Vb` becomes `cartexp_horizontal_stationary` *or*
   `polexp_annulus_grating_stationary`.
5. **Benson → conventional polar angle, on one branch only.** `createTables.m:75` applies
   `map_theta`; `meanWithinLabel.m` does not. Three downstream files re-apply the conversion
   by hand as `[90 45 0 315 270 225 180 135]`. **This is the trap the reproduction fell into.**
6. **z-scoring.** CSV betas are non-z-scored (`betas_nonzscored.mat`) but ship `beta_std` so
   the z-scored contrast is recoverable; `results.mat` contrasts are already z-scored.

## 8. Still open

- **`Fig 7` LME.** `FINDINGS.md` §3 reports the manuscript's Fig 7B `da` H−V as ≈0.02 against
  the independent −0.45. `lme1_fit.m:88` uses the *correct* angle convention, so this is not
  the ordering issue. Note that the 0.02 is a transcription from the PDF by the earlier agent
  that I have not verified (the eight §4 values all did verify). Worth confirming it is a real
  inconsistency before chasing it.
- **z-scoring provenance.** The saved results are z-scored but `main_singlesub.m:131` now has
  `normalize = 0` (changed 2026-07-14). Anyone re-running stage 01 will silently produce
  non-z-scored betas and different Fig 7 statistics. This flag should be made an explicit,
  recorded parameter rather than a hand-edited constant — and the manuscript Methods, which
  commit to z-scoring, match the *saved* files, not the current default.
- **The z-scoring decision itself** (`Reproduction/NEXT_STEPS.md`) remains open; this audit
  only establishes which variant the existing figures came from.

## 9. Recommendation

**Do not "fix" `compute_derivativeDirections.m`.** Aligning its `polarAngles` with
`[0 45 90 …]`, as `FINDINGS.md` currently suggests, would *introduce* the cardinal-meridian
swap into working code and corrupt both first-harmonic derived asymmetries.

The one substantive improvement worth making is defensive, not corrective: `meanWithinLabel.m`
is the only consumer of `angle_adj.mgz` that does not call `map_theta`, and the Benson
ordering it produces is re-derived by hand in three separate downstream files. Storing the
wedge centres alongside the data (or converting once at the source and updating the three
`[90,45,0,…]` literals together) would remove the trap that the reproduction fell into —
but that is a refactor, and the current code's *output* is correct.
