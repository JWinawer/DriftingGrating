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

Underneath, this is a question about **what units make observers commensurable** (§6) — a real
question, since percent BOLD is scanner-dependent, and one that is separate from data quality.
Normalising is the right call, but it has a precondition: an observer whose gain cannot be
measured cannot be normalised. **Apply both and the reversal disappears** — at n=6 all five
normalisers, `beta_std` included, put H−V first (§8). Still provisional on the GLM fits, which
nothing in the pipeline currently checks (§7).

| section | script |
|---|---|
| §1–2, §4–5 | `cleanroom/diagnose_zscore_fig7.m` |
| §3a–3b | `cleanroom/diagnose_response_signs.m` (via `cleanroom/load_allconditions.m`, a variant of `load_and_filter` caching all 13 conditions) |
| §6 precision | `cleanroom/compare_subject_weighting.m` |
| §6 units / gain | `cleanroom/diagnose_gain_normalization.m` |
| §8 | `cleanroom/diagnose_exclusion_x_normalization.m` |

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

## 6. Units and precision are two separate decisions

An earlier version of this section framed the choice as "how should observers be weighted" and
concluded that down-weighting sub-0395 "discards signal, not noise." **That framing was wrong**,
and the correction matters, so it is set out here rather than quietly edited away.

### Why the weighting framing is slippery

Percent BOLD is not a well-defined biological quantity. Observers differ in how large a response
any stimulus produces — 1% versus 3% is ordinary — and changing the pulse sequence or field
strength would change those numbers with no change in the underlying neural response. Averaging
raw values therefore weights observers by an arbitrary gain, which is not obviously better than
any other choice. Normalizing each observer before averaging is the standard response, and has a
direct precedent in single-unit work, where each neuron's responses are scaled so its peak is 1
before averaging across neurons.

Crucially, **normalizing and reweighting are the same operation seen from two sides**. Averaging
normalized per-observer values with *equal* weight is arithmetically identical, up to one global
constant, to a 1/gain-weighted average of the raw values:

```
equal-weighted mean of normalised radTan = 0.6807
1/gain-weighted mean of raw radTan       = 0.2072
ratio 3.2859 = mean(1/gain) exactly
```

They rank the four asymmetries identically. So "z-scoring gives sub-0037 20.4% of the weight
against a uniform 12.5%" is **not by itself an objection** — that is simply what any
normalization looks like when viewed in un-normalized units. The objection has to be to the
*divisor*, not to the reweighting. The §3 weight table should be read as a description of what
the normalization does, not as evidence that it is unfair.

There are really two orthogonal decisions:

- **Units** — in what units are observers commensurable? (raw % change, or normalized by gain)
- **Precision** — should better-measured observers count more? (equal vs inverse-variance)

`beta_std` conflates them: it is a units choice that also happens to reweight, and it was never
stated as either.

### On precision, the answer is clear

Within-observer bootstrap SEs are 0.015–0.033 in raw units, while observers disagree by far more:

| | between-subject SD | mean within-subject SE | variance ratio |
|---|---|---|---|
| da radTan | 0.209 | 0.020 | **104×** |
| da H−V | 0.152 | 0.021 | 52× |
| dg H−V | 0.200 | 0.025 | 64× |
| dg radTan | 0.085 | 0.023 | 14× |

Between-observer variance exceeds within-observer measurement variance by one to two orders of
magnitude, so inverse-variance weights are nearly uniform and precision weighting barely moves
the answer (da radTan 0.154 vs 0.150 equal-weighted). **Precision is not what is in dispute** —
no observer here is meaningfully noisier than the others at the level of the group estimate.

### On units, the premise fails empirically

Gain normalization requires that each observer *have* an estimable gain. Test it three ways
(`diagnose_gain_normalization.m`):

**1. Is gain a property of the observer, or of the session?** If it were a stable individual
trait it should transfer between the two experiments. Correlations across the 8 observers:

| gain estimate | Pearson | Spearman |
|---|---|---|
| `beta_std` (all 13) | 0.55 | 0.17 |
| std of 8 motion conditions | 0.45 | 0.38 |
| mean motion drive over blank | 0.67 | 0.36 |
| peak stimulus response (ephys-style) | 0.68 | 0.36 |

Weak, and at n=8 none is distinguishable from zero. Whatever these divisors measure is
substantially **session-specific**, not a stable individual trait. That is the assumption the
normalization rests on, and it is not supported here.

**2. Is it estimable at all?** For the polar experiment, blank-referenced gain estimates for
**sub-0037** (mean motion drive −0.03, peak +0.29) and **sub-0201** (−0.40 and −0.10) are zero or
negative. Their V1 did not respond to the stimuli, so they have no measurable gain. Normalizing
by a non-positive divisor sign-flips the observer, which is why both the mean-motion and
ephys-style rows below are uninterpretable rather than merely different:

| normalizer (polar, per-observer scalar) | H−V | radTan | ordering | |
|---|---|---|---|---|
| raw | −0.211 | 0.150 | H−V larger | |
| `beta_std`, all 13 (published) | −0.568 | **0.681** | **radTan larger** | |
| std of 8 motion conditions | −1.043 | 0.907 | H−V larger | |
| mean motion drive | −0.262 | −1.368 | — | **invalid: divisor ≤ 0** |
| peak response (ephys-style) | 0.093 | 0.071 | — | **invalid: divisor ≤ 0** |

Note *which* divisors stay positive. Standard-deviation-based ones always do, because a standard
deviation is positive by construction — including for an observer whose V1 did not respond, where
it returns the **noise level** in place of a gain and then divides by it. `beta_std` does not
avoid the failure that breaks the other two; it hides it.

**3. Would normalization even reconcile the dissenting observer?** No, and this is where the
original "signal, not noise" claim was closest to right for the wrong reason. Gain is positive, so
no normalization can flip a sign:

| normalizer | sub-0395 radTan | group mean | ratio |
|---|---|---|---|
| raw | −0.244 | 0.150 | −1.63 |
| `beta_std` (published) | −0.331 | 0.681 | −0.49 |
| std of 8 motion | −0.789 | 0.907 | −0.87 |

sub-0395 dissents in *every* set of units, and proportionally no less after normalization.
Normalization does not bring that observer into line with the group — it reduces the observer's
leverage on the mean. Those are different things, and only the first would justify the change.

### Where that leaves it

The conceptual case for normalizing is sound: raw percent BOLD is scanner-dependent and averaging
it weights observers arbitrarily. What fails is the **execution available in this dataset**. A
defensible gain normalization needs a divisor that is (a) independent of the effect being
measured, (b) positive and stably estimable for every observer, and (c) a property of the observer
rather than the session. No divisor derivable from these 13 conditions satisfies all three, and
`beta_std` satisfies none of them cleanly — it contains the four analyzed conditions and the
blank, it substitutes noise for gain in two observers, and it transfers across experiments at
r = 0.55/0.17.

**The constructive route** is to estimate gain from something genuinely independent: the
retinotopy scan is a separate session with its own stimulus, and `prfvista_mov` response
amplitude could supply a per-observer gain that is not contaminated by the orientation conditions
and does not collapse when V1 fails to respond to *these* stimuli. That is testable and is worth
doing before the z-scoring question is settled either way. Failing that, report raw units and say
that observers were weighted equally as a deliberate choice.

## 7. Open: nobody has checked the GLM fits

This analysis has been screening the wrong quality metric. `allsubjectsTable.csv` carries exactly
one quality column, **`pRF_r2` — the retinotopy model fit, not the GLM fit** — and the inclusion
filter (`pRF_r2 > 0.1`) is built on it. **No GLMsingle fit-quality metric enters the pipeline at
any stage.** Given §3a, that is a gap: a vertex can pass the filter on the strength of a good pRF
fit while its 13 condition betas are essentially unconstrained.

The metrics do exist, but **only one subject's are checkable from this repository.**
`Support/sub-0255/{dg,da}/results.mat` are the only GLM outputs present; a full sweep of the repo
turns up no other `results.mat`, and neither `allsubjectsTable.csv` (nor the second copy under
`Support/summaryTables_wleftV2d/`, which has an identical column set) carries a GLM metric. The
claim that the *other* seven subjects have the same fields is an **inference from
`main_singlesub.m`**, which saves `results.allevents = modelOut{1,4}` unconditionally under the
`glmsingle` HRF setting — not something that has been observed. It needs `/Volumes/Vision`
mounted to confirm.

What sub-0255's `results.mat` does contain (both experiments, 270,291 vertices):

| field | what it gives you |
|---|---|
| `R2` | variance explained per vertex — the direct GLM quality measure |
| `R2run` | per-run R², so run-to-run consistency and bad runs are visible |
| `FRACvalue` | ridge fraction chosen per vertex; low = heavy shrinkage = poorly constrained |
| `noisepool` | which vertices GLMsingle itself classified as noise |
| `HRFindex` / `xvaltrend` | which HRF was selected, and the cross-validation trend |
| `meanvol` | mean EPI intensity — flags dropout |

Recommended next step, in priority order:

1. Extract, for all 8 subjects × both experiments, the distribution of `R2`, `FRACvalue`,
   `noisepool` membership and `meanvol` **within the analysed V1 4–8° patch** (this needs the
   V1 label to map `results.mat` vertices onto CSV rows — the CSV has no vertex index column,
   which is itself worth fixing).
2. Ask specifically whether **sub-0201 and sub-0037** are outliers on those metrics in the polar
   experiment. sub-0201's blank beta exceeding all 12 stimulus conditions in both experiments is
   the kind of thing a bad-run or motion problem produces; sub-0037's total lack of condition
   differentiation in `da` but strong response in `dg` points at a session-specific failure.
   `R2run` would show a single bad run.
3. Decide inclusion/exclusion on those grounds, **before** settling the z-scoring question — the
   two are entangled, since z-scoring is currently acting as an implicit *inverse* quality
   weighting.
4. Consider adding a GLM-`R2` column to `allsubjectsTable.csv` so the vertex filter can screen on
   it alongside `pRF_r2`.

## 8. Normalisation and exclusion are one decision, and together they resolve it

Accepting that observers should be brought into commensurate units — which is the right call
(§6) — carries a precondition: **an observer whose gain cannot be measured cannot be normalised.**
In the polar experiment sub-0037 (−0.03) and sub-0201 (−0.40) have non-positive blank-referenced
gain, so for them the normalisation is undefined. `beta_std` conceals this by returning a positive
number regardless, and then divides by it.

So the choice is not "normalise or not" but "normalise, having first removed the observers for
whom normalisation is undefined". Crossing the two decisions
(`diagnose_exclusion_x_normalization.m`); `P` is a subject bootstrap of P(radTan > |H−V|):

**All 8 observers**

| normaliser | H−V | radTan | P | ordering |
|---|---|---|---|---|
| raw | −0.211 | 0.150 | 0.28 | H−V larger |
| **`beta_std`, all 13 (published)** | −0.568 | **0.681** | **0.63** | **radTan larger** |
| std of 8 motion | −1.043 | 0.907 | 0.39 | H−V larger |
| mean motion drive | — | — | — | *invalid: divisor ≤ 0* |
| peak response | — | — | — | *invalid: divisor ≤ 0* |

**n = 6 (gain estimable for every observer)**

| normaliser | H−V | radTan | P | ordering |
|---|---|---|---|---|
| raw | −0.254 | 0.114 | 0.10 | H−V larger |
| `beta_std`, all 13 (published) | −0.625 | 0.460 | 0.23 | **H−V larger** |
| std of 8 motion | −1.213 | 0.725 | 0.17 | H−V larger |
| mean motion drive | −0.284 | 0.237 | 0.31 | H−V larger |
| peak response (ephys-style) | −0.204 | 0.158 | 0.27 | H−V larger |

Two things happen at once. Removing those two observers makes **every** normaliser well-defined —
the blank-referenced ones become positive throughout, so five schemes can be compared instead of
three. And all five then **agree**: horizontal−vertical is the largest polar-grating asymmetry,
including under the manuscript's own `beta_std`, which flips from 0.681 to 0.460 against an H−V of
−0.625.

**The radTan-largest result exists only in the configuration that normalises two observers whose
gain is unmeasurable.** It does not survive its own methodology applied consistently.

Strength of evidence, stated honestly: P runs 0.10–0.31 across the five, so the *direction* is
unanimous but the margin is moderate — this is good evidence against radTan-largest, not proof of
H−V-largest. The recommendation in §5 stands: do not build the framing on the ranking.

**This is provisional on the GLM audit (§7).** The exclusion is currently justified on
response-amplitude grounds alone. If the GLMsingle metrics show sub-0037 and sub-0201 are fine,
the reasoning above needs revisiting — though it is hard to see how an observer with no
stimulus-driven response can be normalised by its response amplitude whatever the fit statistics
say.

## Recommendation

0. **Normalise, but exclude first** (§8). Bringing observers into commensurate units is
   justified — percent BOLD is scanner-dependent. Its precondition is a measurable gain, which
   sub-0037 and sub-0201 lack in the polar experiment. Apply both and all five normalisers agree
   that H−V is the largest polar asymmetry, the manuscript's own `beta_std` included. The
   radTan-largest result survives only when two unnormalisable observers are normalised anyway.
   *Provisional on the GLM audit (§7).*
1. **Even so, do not frame the paper around which polar-grating asymmetry is "largest."** The
   bootstraps run 0.10–0.31 at n=6 and 0.28–0.63 at n=8: unanimous in direction once the
   precondition is met, but never decisive. The ranking is not a load-bearing result.
2. The claim the data *do* support, in both variants and with the sign preserved throughout, is
   the one in AGENTS.md: **the polar-frame asymmetries appear for polar gratings and the
   Cartesian-frame asymmetries weaken** relative to the Cartesian experiment. Compare each
   asymmetry *across experiments*, not against each other within a panel.
3. **If z-scoring is retained** (the Methods currently commit to it), the manuscript should say
   plainly that observers are normalized by their own overall response amplitude, and should
   justify the divisor. Note this is a **units** claim, not a fairness one (§6): normalizing per
   observer is a legitimate response to percent BOLD being scanner-dependent. What needs
   defending is `beta_std` specifically.
4. **A per-vertex divisor is not what is doing the work here** — a per-subject scalar reproduces
   the entire effect. If the intent of z-scoring was vertex-level gain control, that intent is
   not what is changing Fig 7; the observer-level rescaling is.
5. **If a per-vertex divisor is wanted, use one built from the 8 motion conditions** rather than
   all 13 (§3b). It excludes the blank, is independent of the conditions being analysed, and
   yields far more uniform observer weights. Note that it does *not* rescue the z-scored
   ordering — it restores H−V as the largest polar-grating asymmetry.
6. **Separate the units question from the precision question, and settle precision first** (§6).
   Precision is not in dispute: between-observer variance exceeds within-observer measurement
   variance by 14–104×, so inverse-variance weighting converges on equal weighting and changes
   nothing. Units are the live question, and the case for normalizing per observer is sound in
   principle — percent BOLD is scanner-dependent. What fails is every divisor available in these
   13 conditions: none is simultaneously effect-independent, positive for all 8 observers, and
   stable across sessions.
6a. **Try estimating observer gain from the retinotopy scan** (§6, "Where that leaves it"). It is
   an independent session with its own stimulus, so `prfvista_mov` response amplitude would not be
   contaminated by the orientation conditions and would not collapse when V1 fails to respond to
   *these* stimuli. This is the one route that could make a principled normalization possible;
   it is worth testing before the z-scoring question is settled either way.
7. **Consider whether sub-0201 and sub-0037 belong in the polar analysis at all** (§3a).
   sub-0201's blank beta exceeds all 12 stimulus conditions in *both* experiments; sub-0037 shows
   no differentiation whatsoever in the polar experiment despite a strong Cartesian response.
   Whatever is decided, it should be decided on data-quality grounds and stated — not left to be
   applied implicitly, and in reverse, by the choice of normalizer.
8. **Check the GLM fits before settling any of this** (§7). No GLMsingle quality metric currently
   enters the pipeline — the only quality filter is on the *pRF* fit. The metrics are sitting
   unused in each subject's `results.mat`.

## Caveats

- The `included` flag is 1 for all 11,075 vertices surviving the V1 / 4–8° / R² > 0.1 filter, so
  it plays no role.
- Amplification factors for individual subjects are unstable where the raw effect is near zero
  (sub-0426's raw radTan is +0.001), which is why the group-level covariance identity in §2 is
  the right summary rather than per-subject ratios.
