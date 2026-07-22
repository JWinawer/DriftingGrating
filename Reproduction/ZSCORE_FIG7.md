# Why z-scoring reverses radTan vs H−V in Fig 7B (polar gratings)

**Question.** In the z-scored Fig 7, the largest polar-grating asymmetry is radial−tangential
(0.603) with horizontal−vertical second (−0.446). Without z-scoring the order reverses
(0.150 vs −0.211). What accounts for it?

**Answer in one line.** The reversal is a **between-subject reweighting**, not a vertex-level
effect: z-scoring weights each subject by 1/`beta_std`, and in the polar experiment `beta_std`
is essentially the subject's overall response amplitude (r = +0.94). The subjects with the
weakest polar-grating responses happen to be the ones with the largest radial−tangential and
the least negative horizontal−vertical, so a **single positive covariance inflates radTan and
shrinks H−V at the same time** — purely because the two effects have opposite sign. The two most heavily up-weighted observers turn out to be the two
whose polar-grating data are anomalous (§3a), and rebuilding the divisor from conditions that
exclude the blank puts H−V back on top (§3b).

Reproduce §1–2 and §4–5 with `cleanroom/diagnose_zscore_fig7.m`, and §3a–3b with
`cleanroom/diagnose_response_signs.m` (which needs `cleanroom/load_allconditions.m`, a variant
of `load_and_filter` that caches all 13 conditions per experiment).

---

## 1. It is a subject-level effect, not a vertex-level one

Replacing the per-vertex divisor with progressively coarser versions (`da`):

| divisor | HV | cardObl | radTan | polcard |
|---|---|---|---|---|
| D0 raw | −0.211 | −0.030 | 0.150 | 0.034 |
| D1 one global scalar | ×3.08 | ×3.08 | ×3.08 | ×3.08 |
| **D2 one scalar per subject** | **×2.69** | ×2.76 | **×4.54** | ×6.29 |
| D3 per subject × wedge | ×2.01 | ×2.08 | ×3.85 | ×4.51 |
| D4 true per vertex | ×2.11 | ×2.17 | ×4.02 | ×5.16 |

A **single number per subject (D2) already produces the whole differential** — in fact
slightly more of it than the real per-vertex z-scoring does. Within-wedge reweighting of
individual vertices contributes nothing to the reversal. The same ladder for `dg` shows no
comparable spread (all four asymmetries land at ×1.9–2.4), so this is specific to `da`.

## 2. The exact algebra

For a per-subject divisor, amplification decomposes exactly as

```
amp = mean_i(1/s_i)  +  cov(a_i, 1/s_i) / mean_i(a_i)
```

| | mean(a) | cov(a, 1/s) | cov/mean(a) | amp |
|---|---|---|---|---|
| radTan | **+0.150** | +0.188 | +1.25 | 4.54 |
| H−V | **−0.211** | +0.127 | −0.60 | 2.69 |

with `mean(1/s) = 3.29` as the generic factor. The covariance is **positive for both**:
low-`beta_std` subjects have more positive values of *both* asymmetries. Because radTan is a
positive effect it gets inflated above 3.29; because H−V is a negative effect the identical
covariance pushes its amplification below 3.29. One covariance, two opposite consequences.

## 3. Which subjects, and what `beta_std` actually is

| subject | da `beta_std` | mean resp | med pRF R² | raw radTan | raw H−V | weight under z |
|---|---|---|---|---|---|---|
| sub-0037 | **0.186** | **−0.02** | 0.66 | **+0.354** | **+0.026** | **20.4%** |
| sub-0201 | 0.206 | −0.28 | 0.60 | +0.162 | −0.192 | 18.4% |
| sub-0250 | 0.288 | +0.27 | 0.65 | +0.383 | −0.266 | 13.2% |
| sub-wlsubj123 | 0.284 | +0.39 | 0.55 | +0.237 | −0.051 | 13.4% |
| sub-0255 | 0.328 | +0.20 | 0.62 | +0.263 | −0.329 | 11.6% |
| sub-wlsubj124 | 0.374 | +0.59 | 0.61 | +0.047 | −0.137 | 10.2% |
| sub-0426 | 0.501 | +0.53 | 0.36 | +0.001 | −0.423 | 7.6% |
| sub-0395 | **0.739** | **+1.35** | 0.62 | **−0.244** | −0.319 | **5.1%** |

(uniform weighting would be 12.5% each. **"mean resp"** = median over vertices of
`mean(4 stationary orientations) − blank`, i.e. the mean of the four quantities the figures are
built from — the blank is already subtracted, and the 8 motion conditions do not enter it.
See §3a for why two subjects are negative.)

At the subject level `beta_std` is **overall responsiveness**, not data quality:

- corr(`beta_std`, median mean stimulus−blank response) = **+0.936**
- corr(`beta_std`, median pRF R²) = **−0.281** (if anything, *negative*)
- corr(`da beta_std`, `dg beta_std`) = +0.55 Pearson but only **+0.17 Spearman** — it is largely
  experiment-specific, not a stable trait of the observer.

The consequence is uncomfortable. **sub-0395 has by far the strongest polar-grating drive**
(median blank-subtracted response +1.35, vs ≈0 for the top two) and is the one subject whose
radial−tangential goes the *wrong* way (−0.244); z-scoring cuts its weight to 5.1%.
**sub-0037 has essentially no net response to the polar stimuli at all** (−0.02) yet carries the
largest radTan (+0.354) and the only positive (wrong-signed) H−V; z-scoring raises its weight
to 20.4%. sub-0201's net response is actually negative (−0.28) and it gets 18.4%.
So the effective ordering of subject weights runs roughly opposite to how well the polar
stimuli drove V1.

## 3a. The two negative subjects are real, and they are the two most up-weighted

V1 showing *less* response to a grating in its receptive field than to blank is surprising
enough to check. It is not an averaging artefact — the quantity is already a blank-subtracted
contrast — and it survives pulling in all 13 GLMsingle conditions
(`diagnose_response_signs.m`). First, a definition check: `beta_mean` and `beta_std` are exactly
the mean and std of the 13 conditions (agreement to 6e-15), so the published divisor **includes
the blank and the 4 analysed stationary conditions**.

Median raw betas per subject, polar experiment:

| subject | stat×4 | motion×8 | blank | stat−blank | motion−blank | blank's rank among the 13 |
|---|---|---|---|---|---|---|
| **sub-0037** | −0.104 | −0.102 | −0.091 | **−0.023** | **−0.030** | 8 of 13 (mid-pack) |
| **sub-0201** | −0.393 | −0.525 | −0.122 | **−0.278** | **−0.402** | **13 of 13 (the largest)** |
| sub-0255 | −0.025 | +0.345 | −0.279 | +0.201 | +0.631 | 2 |
| sub-wlsubj123 | −0.069 | +0.057 | −0.544 | +0.391 | +0.645 | 1 |
| sub-wlsubj124 | −0.089 | +0.253 | −0.726 | +0.590 | +0.960 | 1 |
| sub-0395 | +0.179 | +1.009 | −1.103 | +1.352 | +2.196 | 1 |
| sub-0426 | −0.172 | +0.607 | −0.735 | +0.527 | +1.288 | 1 |
| sub-0250 | −0.019 | +0.174 | −0.344 | +0.270 | +0.544 | 2 |

(A negative beta is not itself alarming: ~47% of all betas are negative and their median is
+0.04, so the beta zero-point sits near the middle of the conditions. Only stimulus-vs-blank
is interpretable.)

The two subjects fail in **different** ways, and neither is fixed by dropping the blank from an
average:

- **sub-0201** — the blank is the *largest* of all 13 betas. All 12 stimulus conditions,
  motion included, fall below blank. The same holds in the Cartesian experiment
  (blank rank 12 of 13; stat−blank −0.203, motion−blank −0.225), so this is a subject-level
  problem present in both datasets, not a one-off blank estimate.
- **sub-0037** — a flat profile: stat (−0.104), motion (−0.102) and blank (−0.091) are all
  effectively equal and the blank sits mid-pack, i.e. V1 did not differentiate anything in the
  polar experiment. Yet the *same subject* responds strongly to Cartesian gratings
  (motion−blank = **+0.997**). So this is specific to the polar session.

These are exactly the two subjects the published normalization weights most heavily (20.4% and
18.4% against 12.5% uniform) — and necessarily so: their `beta_std` is small *because* they did
not respond, and `beta_std` is the divisor.

## 3b. Removing the blank from the divisor reverses the reversal

Because `beta_std` spans all 13 conditions, the blank — usually the most extreme of the 13 — is
a large part of what it measures. A divisor built from the **8 motion conditions only** drops
the blank *and* is independent of the 4 stationary conditions being analysed, so it does
per-vertex gain control without being contaminated by the effect or by the blank estimate:

| divisor (polar experiment) | H−V | cardObl | radTan | polcard | ordering |
|---|---|---|---|---|---|
| none (raw) | −0.211 | −0.030 | 0.150 | 0.034 | H−V larger |
| **`beta_std`, all 13 (published)** | −0.446 | −0.064 | **0.603** | 0.173 | **radTan larger** |
| **std of the 8 motion conditions** | **−0.793** | −0.133 | 0.725 | 0.209 | **H−V larger** |

Subject weights flatten markedly under the motion-only divisor (7.8–15.2%, versus 5.1–20.4%
under `beta_std`); sub-0037 in particular drops from 20.4% to 13.0%, because its motion-only std
(0.187) is no smaller than its all-13 std (0.186) — every other subject's ratio is 1.3–2.6.

So **the published ordering is not reproduced by an equally principled but uncontaminated
normalizer.** Under the motion-only divisor 6 of 8 leave-one-out subsets put H−V first and the
subject bootstrap gives P(radTan > |H−V|) = **0.43**. Across all three divisors that
probability is 0.28 (raw), 0.69 (published), 0.43 (motion-only) — all straddling chance, and the
published choice is the only one of the three that puts radTan on top.

The Cartesian experiment is unaffected: H−V is largest under all three divisors.

(A fourth divisor, `mean(8 motion) − blank`, is *not* usable — it passes through zero for many
vertices and flips signs; it is shown in the script only to document that.)

## 4. Why H−V is the natural loser: the two are harmonics of one pattern

For polar gratings, pinwheel is radial and annulus tangential at *every* polar angle, so:

- **radTan is the polar-angle MEAN of (pinwheel − annulus).**
- **H−V is the cos(2θ) MODULATION of that same quantity** (+ the matching spiral term), because
  pinwheel is locally horizontal at PA 0/180 and locally vertical at PA 90/270.

Verified exactly:

```
raw     : radTan = +0.150 ;  H−V = cos2θ(rad−tang) −0.101 + cos2θ(spirals) −0.110 = −0.211
zscored : radTan = +0.603 ;  H−V = cos2θ(rad−tang) −0.244 + cos2θ(spirals) −0.202 = −0.446
```

These are not two independent findings — they are the DC and second-harmonic components of a
single polar-angle profile, so anything that changes the *shape* of that profile trades one
against the other. Per-wedge rad−tang:

| PA | 0 | 45 | 90 | 135 | 180 | 225 | 270 | 315 |
|---|---|---|---|---|---|---|---|---|
| raw | −0.046 | 0.078 | 0.318 | 0.034 | −0.086 | 0.247 | 0.357 | 0.298 |
| z-scored | **+0.169** | 0.373 | 0.876 | 0.281 | **−0.011** | 0.929 | 1.230 | 0.978 |

Raw rad−tang is *negative* on the horizontal meridian (PA 0 and 180) and strongly positive on
the vertical meridian. The reweighting largely removes that negative lobe, which simultaneously
raises the mean (radTan ↑) and flattens the cos2θ modulation (|H−V| ↓). That is the reversal.

## 5. Neither ordering is statistically established

- **Leave-one-subject-out (z-scored):** dropping **sub-0037** flips the order back
  (radTan 0.447 vs |H−V| 0.552). All seven other drops preserve it. The z-scored result rests
  on one observer.
- **Bootstrap over the 8 subjects:** P(radTan > |H−V|) = **0.69** z-scored, **0.28** raw.

So the rank order is not a reliable feature of the data in *either* variant. It is within the
sampling noise of n=8.

---

## Recommendation

1. **Do not frame the paper around which polar-grating asymmetry is "largest."** The bootstrap
   (0.69 vs 0.28) says the dataset cannot support that ranking either way, and the z-scored
   ordering hinges on a single observer.
2. The claim the data *do* support, in both variants and with the sign preserved throughout, is
   the one in AGENTS.md: **the polar-frame asymmetries appear for polar gratings and the
   Cartesian-frame asymmetries weaken** relative to the Cartesian experiment. Compare each
   asymmetry *across experiments*, not against each other within a panel.
3. **If z-scoring is retained** (the Methods currently commit to it), the manuscript should say
   plainly that it weights observers by the inverse of their overall response amplitude, and the
   per-subject weights should be reported. As it stands the normalization silently down-weights
   the best-driven observer 4× relative to the least-driven one.
4. **A per-vertex divisor is not what is doing the work here** — a per-subject scalar reproduces
   the entire effect. If the intent of z-scoring was vertex-level gain control, that intent is
   not what is changing Fig 7; the observer-level rescaling is.
5. **If a per-vertex divisor is wanted, use one built from the 8 motion conditions** rather than
   all 13 (§3b). It excludes the blank, is independent of the conditions being analysed, and
   yields far more uniform observer weights. Note that it does *not* rescue the z-scored
   ordering — it restores H−V as the largest polar-grating asymmetry.
6. **Consider whether sub-0201 and sub-0037 belong in the polar analysis at all** (§3a).
   sub-0201's blank beta exceeds all 12 stimulus conditions in *both* experiments; sub-0037 shows
   no differentiation whatsoever in the polar experiment despite a strong Cartesian response.
   Whatever is decided, it should be decided on data-quality grounds and stated — not left to be
   applied implicitly, and in reverse, by the choice of normalizer.

## Caveats

- The `included` flag is 1 for all 11,075 vertices surviving the V1 / 4–8° / R² > 0.1 filter, so
  it plays no role.
- Amplification factors for individual subjects are unstable where the raw effect is near zero
  (sub-0426's raw radTan is +0.001), which is why the group-level covariance identity in §2 is
  the right summary rather than per-subject ratios.
