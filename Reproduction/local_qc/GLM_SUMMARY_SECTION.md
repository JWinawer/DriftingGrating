# The GLM-summary section — reliability, route, and the Figure 4 controls

*What quantifies the first-stage GLM estimates for the Results section "GLM summary", plus draft
manuscript text. Computed by
[`../cleanroom/splithalf_reliability.m`](../cleanroom/splithalf_reliability.m); per-observer values
in `splithalf_reliability.csv` / `.mat`. Written 2026-08-19.*

All numbers here are on the settled route: **across-vertex mean, with the observer pRF-gain
rescaling** (JW, 2026-08-19). Result 3 below is what established that route.

| section | what it settles |
|---|---|
| Why reliability rather than R² | why the section's current three R² numbers are the wrong currency |
| Method | the split-half over runs |
| Which route the manuscript uses | mean vs. median, gain vs. none — **resolved** |
| Results 1–2 | split-half reliability; measurement error against the effect |
| Results 3–4 | the asymmetry tables on the current route; gain × precision weighting |
| Result 5 | the Fig-4A run-mismatch control |
| Draft text | prose for the manuscript |
| Figure 4C | why it carries no quantitative claim |

## Why reliability rather than R²

The section currently reports three different quantities all called R² — variance explained in a
group-mean ROI time series (0.42 / 0.28), fit to a trial-triggered average (0.37–0.21), and
vertexwise GLM R² (Fig 4C) — from single examples, and none of them is the quantity the paper's
claims rest on.

R² is also the wrong currency for this dataset. The reference condition is full-field pink noise
present throughout every run ([`REPORT.md`](REPORT.md) §1), so GLM R² measures differentiation
*among* conditions, not responsiveness. Median V1 R² in the analysed 4–8° band is 4.7% (dg) and
4.4% (da) against 3.3% for the whole surface — a reader who takes that at face value concludes the
data are bad. What Figures 5–8 require is that the analysed profile be **reproducible**, which is
measured by resampling the measurement, i.e. splitting the runs.

## Method

The design is balanced — 8 runs, 4 trials per condition per run — so a 4-vs-4 run split is balanced
by construction. For each observer and experiment: recompute the analysed profile (4 orientations ×
8 polar-angle wedges) from each half, correlate the two halves, Fisher-z average over **all 35**
balanced splits, and Spearman–Brown correct to the full run count. Two observers depart from 8 runs
(sub-0255 dg has 9, sub-0395 da has 6); for the odd count each run is dropped in turn and the even
procedure averaged over the 9 drops, which corrects to 8 rather than 9 and is therefore
conservative.

Reliability is reported at two spatial levels × two contrast levels:

| | includes overall response level | orientation-differential only |
|---|---|---|
| **per vertex** (nVert 1126–1786) | `vertexFull` | `vertexDiff` |
| **ROI profile** (32 values — what Figs 5–8 analyse) | `roiFull` | `roiDiff` |

**What the correlation is computed over.** For `roiFull`, the 32 numbers are the response to each
of the 4 orientations in each of the 8 polar-angle ROIs, for one observer in one experiment — the
same 32 numbers that Figures 5/6 plot (8 observers × 32 = 256 values per experiment). Each half of the runs yields its own 4 × 8 matrix; *r* is the Pearson correlation of those
two matrices flattened to 32-vectors. For `roiDiff` it is the same 32 values after subtracting,
within each ROI, that ROI's mean over the 4 orientations — so each group of 4 sums to zero and the
32 values carry 24 independent degrees of freedom. For the two `vertex` rows the same thing is done
one level down: nVert × 4 orientation contrasts, flattened (≈4500–7100 values per observer).

The `Diff` variants remove each vertex's / each wedge's mean across the four orientations. That
strips response-level differences which are reliable but carry no orientation information, leaving
exactly the subspace that carries all four asymmetries. **`roiDiff` is the strict test and the one
to report.**

**Validation.** This route (run-averaged GLMsingle single-trial betas) reproduces the manuscript ROI
pipeline — see below — and the split-half within-observer SEs agree with the independent
bootstrap-over-runs estimator in [`REPORT.md`](REPORT.md) §2.7 to ±0.005. The vertex counts also
match `gainSummary.csv`'s independently-computed `nVertices` exactly for all 8 observers,
confirming the same V1 / 4–8° / pRF R² > 0.1 selection.

## Which route the manuscript uses — resolved

Four candidate routes were compared against the eight manuscript asymmetries, summing the absolute
deviation over all eight values:

| route | Σ&#124;route − manuscript&#124; |
|---|--:|
| median across vertices, unweighted | 0.128 |
| median across vertices, gain-weighted | 0.142 |
| mean across vertices, unweighted | 0.040 |
| **mean across vertices, gain-weighted** | **0.015** |

| exp | asymmetry | manuscript | medUnwt | medWtd | meanUnwt | **meanWtd** |
|---|---|--:|--:|--:|--:|--:|
| dg | horiz−vert | −0.55 | −0.480 | −0.476 | −0.554 | **−0.553** |
| dg | card−obl | −0.22 | −0.204 | −0.192 | −0.237 | **−0.223** |
| dg | rad−tang | 0.10 | 0.116 | 0.112 | 0.109 | **0.105** |
| dg | polc−polo | 0.04 | 0.029 | 0.027 | 0.042 | **0.040** |
| da | horiz−vert | −0.22 | −0.211 | −0.218 | −0.215 | **−0.221** |
| da | card−obl | −0.03 | −0.030 | −0.027 | −0.032 | **−0.029** |
| da | rad−tang | 0.15 | 0.150 | 0.148 | 0.151 | **0.151** |
| da | polc−polo | 0.04 | 0.034 | 0.033 | 0.040 | **0.040** |

**Mean + gain weighting matches every manuscript value to ±0.003.** The aggregator is what does the
work (the median route misses dg horiz−vert by 0.07); the gain weighting is a small refinement that
closes the residual on card−obl and polc−polo.

Notes on the two components:

- **Mean vs. median. Applied 2026-08-19.** `meanWithinLabel.m` saves both `meanBOLDpa` and
  `medianBOLDpa`. `lme1_fit.m` was already reading the mean; `plot_NeuralAsymmetries.m` passed the
  median and has been switched, so Figs 5/6 and Fig 7 are no longer computed from different
  aggregates. The clean-room switches via `cfg.aggregator`. `plot1_/plot2_experimentalCond.m` and
  the `compute_derivative*` functions still *name* their first parameter `medianBOLDpa` while
  receiving mean data — deliberately left, because upstream also edits those files and the rename
  would create merge conflicts for no functional gain.
- **Gain weighting.** Upstream commits a498f08 (2026-07-31) and e2acbf7 (2026-08-07) divide each
  observer's BOLD by their own mean pRF gain and multiply the group-mean gain back in
  (`retrieveObserverGainWeights.m`), applied to Figures 5, 6 and 7. Gains are read from
  `gainSummary.csv` (now at `~/dg_collect/`). **The Methods say the group statistics are multiplied
  by the *geometric* mean gain; the code uses the *arithmetic* mean** (`groupGain = mean(gainWeights)`
  in `lme1_fit.m`). For these 8 observers that is 4.4417 vs 4.3952 — a 1% scale difference, so it
  changes no conclusion, but the sentence should be corrected to match the code.

---

## Result 1 — split-half reliability

Spearman–Brown corrected, per observer. **These values are invariant to the gain weighting** (a
per-observer scalar cancels in a within-observer correlation), so they do not depend on that choice
at all:

| observer | nRun dg/da | vertexFull dg | da | vertexDiff dg | da | roiFull dg | da | **roiDiff dg** | **da** |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| sub-0037 | 8 / 8 | 0.734 | 0.121 | 0.797 | 0.301 | 0.899 | 0.328 | **0.916** | **0.500** |
| sub-0201 | 8 / 8 | 0.485 | 0.693 | 0.404 | 0.202 | 0.675 | 0.578 | **0.713** | **0.546** |
| sub-0255 | 9 / 8 | 0.631 | 0.714 | 0.748 | 0.707 | 0.823 | 0.894 | **0.854** | **0.877** |
| sub-wlsubj123 | 8 / 8 | 0.863 | 0.753 | 0.811 | 0.395 | 0.916 | 0.833 | **0.944** | **0.691** |
| sub-wlsubj124 | 8 / 8 | 0.924 | 0.602 | 0.710 | 0.144 | 0.931 | 0.649 | **0.880** | **0.556** |
| sub-0395 | 8 / 6 | 0.944 | 0.903 | 0.894 | 0.668 | 0.976 | 0.893 | **0.976** | **0.800** |
| sub-0426 | 8 / 8 | 0.707 | 0.734 | 0.713 | 0.482 | 0.752 | 0.756 | **0.845** | **0.507** |
| sub-0250 | 8 / 8 | 0.615 | 0.718 | 0.690 | 0.651 | 0.663 | 0.775 | **0.797** | **0.823** |
| **median** | | 0.721 | 0.716 | 0.731 | 0.438 | 0.861 | 0.766 | **0.867** | **0.623** |

**Headline: `roiDiff` median 0.87 (range 0.71–0.98) for Cartesian and 0.62 (range 0.50–0.88) for
polar; ≥ 0.50 in all 16 sessions.**

**The ROI level is the right level, and the data say so.** Per-vertex orientation-differential
reliability is much lower (median 0.73 dg, **0.44** da) than the ROI profile. Each analysed value is
a mean over ~1100–1800 vertices, so per-vertex reliability is not what any claim depends on;
reporting it would understate the reliability of the actual measurements. Worth stating the
reasoning explicitly in the text, to pre-empt "why not per vertex".

## Result 2 — measurement error against the effect being tested

Across-vertex mean, gain-rescaled (geometric). `SNR` is |group mean| / within-observer SE;
`within/total` is the fraction of the **variance across observers** attributable to measurement
error; `p disatt` removes that variance and re-tests.

| | asymmetry | group mean | SD across obs | within-obs SE | SNR | within/total | p observed | p disatt |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| **dg** | horiz−vert | −0.547 | 0.158 | 0.127 | 4.3 | 65% | <0.0001 | <0.0001 |
| | card−obl | −0.221 | 0.158 | 0.073 | 3.0 | 21% | 0.0056 | 0.0030 |
| | rad−tang | 0.104 | 0.080 | 0.037 | 2.8 | 21% | 0.0081 | 0.0044 |
| | polc−polo | 0.039 | 0.044 | 0.018 | 2.1 | 17% | 0.0396 | 0.0274 |
| **da** | horiz−vert | −0.218 | 0.176 | 0.041 | 5.3 | 6% | 0.0098 | 0.0086 |
| | card−obl | −0.029 | 0.045 | 0.020 | 1.4 | 21% | 0.1129 | 0.0811 |
| | rad−tang | 0.150 | 0.170 | 0.136 | 1.1 | 64% | 0.0416 | 0.0043 |
| | polc−polo | 0.040 | 0.105 | 0.073 | 0.5 | 48% | 0.3185 | 0.1795 |

> **`within/total` is the wrong statistic to headline.** It is a ratio of *between-observer*
> variance, so it is large either when
> measurement error is large **or when observers genuinely agree with one another** — and those
> are opposite situations. The two cells above 50% are one of each:
>
> - **dg horiz−vert (65%) is the good case.** The effect is large and near-identical in every
>   observer, so the between-observer SD is small (0.158) and there is little true variation for
>   measurement error to be a small fraction *of*. The true between-observer SD is τ = 0.093. Per
>   observer, |effect| / own SE runs 1.9 to 7.9 (median 5.0), so each observer individually
>   resolves the effect several times over. That is why it is negative in **all 8 observers**
>   (−0.275 to −0.720) and in **all 8 ROIs** (−0.373 to −0.692), and why the group *t* is −9.8.
>   Nothing here is noise-limited.
> - **da rad−tang (64%) is the real limitation.** There the within-observer SE (0.136) is
>   genuinely comparable to the effect (0.150), so SNR is 1.1 — a single observer barely resolves
>   it, and disattenuation moves *p* from 0.042 to 0.004.
>
> **Report SNR and the consistency counts, not `within/total`.** SNR separates the two cases;
> the ratio conflates them.

**What survives as a claim:** the effects are large relative to their own measurement error
(SNR ≥ 2.8 for all four Cartesian-experiment asymmetries and for da horiz−vert), disattenuation
changes no conclusion — every asymmetry significant at *p* < 0.05 stays significant and every null
stays null — and the only genuinely noise-limited cells are da rad−tang and da polc−polo, both of
which a better measurement would *strengthen*.

---

## Result 3 — the asymmetry and context-effect tables on the current route

Regenerated by [`../cleanroom/asymmetry_tables.m`](../cleanroom/asymmetry_tables.m); CSVs at
`../supplement/asymmetry_tables_{asymmetries,context}.csv`. Bootstrap CI = 1000 resamples of
observers (the Figs 5/6 method); *t* CI on 7 df. `obs`/`ROI` count how many of the 8 observers and
8 ROIs share the group sign.

**Four asymmetries, per experiment** (% signal change):

| exp | asymmetry | mean | bootstrap 95% CI | t 95% CI | t | p | obs | ROI |
|---|---|--:|---|---|--:|--:|--:|--:|
| dg | horiz−vert | −0.547 | [−0.646, −0.442] | [−0.679, −0.416] | −9.82 | <0.0001 | 8/8 | 8/8 |
| dg | card−obl | −0.221 | [−0.321, −0.119] | [−0.354, −0.089] | −3.95 | 0.0056 | 7/8 | 8/8 |
| dg | rad−tang | 0.104 | [0.054, 0.157] | [0.037, 0.171] | 3.66 | 0.0081 | 8/8 | 6/8 |
| dg | polc−polo | 0.039 | [0.013, 0.070] | [0.002, 0.076] | 2.52 | 0.0396 | 6/8 | 4/8 |
| da | horiz−vert | −0.218 | [−0.327, −0.110] | [−0.365, −0.071] | −3.51 | 0.0098 | 7/8 | 8/8 |
| da | card−obl | −0.029 | [−0.059, −0.001] | [−0.066, 0.009] | −1.81 | 0.1129 | 6/8 | 6/8 |
| da | rad−tang | 0.150 | [0.033, 0.245] | [0.007, 0.292] | 2.49 | 0.0416 | 7/8 | 6/8 |
| da | polc−polo | 0.040 | [−0.028, 0.109] | [−0.048, 0.127] | 1.07 | 0.3185 | 4/8 | 6/8 |

**Four context effects** (dg − da), paired across observers:

| asymmetry | mean | bootstrap 95% CI | t 95% CI | t | p | obs |
|---|--:|---|---|--:|--:|--:|
| horiz−vert | −0.329 | [−0.493, −0.194] | [−0.528, −0.130] | −3.92 | 0.0058 | 8/8 |
| card−obl | −0.193 | [−0.286, −0.098] | [−0.314, −0.071] | −3.74 | 0.0073 | 7/8 |
| rad−tang | −0.046 | [−0.154, 0.112] | [−0.224, 0.132] | −0.61 | 0.5604 | 6/8 |
| polc−polo | −0.001 | [−0.074, 0.086] | [−0.103, 0.102] | −0.01 | 0.9898 | 3/8 |

Two things moved relative to the median-route numbers in [`../LME.md`](../LME.md) §5 (now carrying
a dated banner with these values): the **horiz−vert context effect grew from −0.268 to −0.329**,
and **da card−obl's bootstrap CI now excludes zero** ([−0.059, −0.001]) while its *t* interval does
not ([−0.066, 0.009], p = 0.11) — the CI-method disagreement LME.md already flags. Report the *t*
interval for that cell, or say explicitly which method is used. The paper's core claims — both
Cartesian context effects significant, both polar ones null — are unchanged.

---

## Result 4 — gain rescaling *and* precision weighting together

Regenerated in [`../LME.md`](../LME.md) §5 / `../supplement/precision_weighted.csv`.

The two are different kinds of operation and stacking them is defensible: the **gain rescaling is a
data transformation** that removes a nuisance multiplicative factor (overall responsiveness),
measured from an independent protocol, so that observers' effects are on a common scale; the
**precision weighting is an estimator choice** about how to combine already-comparable estimates.
Scale and precision are separate problems, so this is not weighting twice for the same reason.

Three things to be careful about:

1. **σᵢ must be measured on the rescaled data.** Rescaling observer *i* by *sᵢ* multiplies their
   measurement SE by *sᵢ* too. Taking σᵢ from unrescaled data and applying the weights to rescaled
   data would systematically over-weight high-gain observers. The clean-room code scales `asymSE`
   by the same factor as `asym`, so this is handled — but it is the easy mistake.
2. **The two partially cancel, by construction.** A low-gain observer is scaled *up*, which
   inflates their effect *and* their noise; precision weighting then pulls them back *down*. Net
   reweighting is milder than either alone. Visible in the numbers: for da rad−tang the weight
   spread rises from 2.95× (median route, no gain rescaling) to 9.89× now, yet the estimate moves
   by about the same amount (0.150 → 0.176).
3. **It still changes nothing that matters.** Precision weighting moves no context effect and no
   Cartesian-frame asymmetry. The single cell whose status changes is da rad−tang
   (0.150 [0.007, 0.292] → 0.176 [0.065, 0.286]), which is marginal by every route and should not
   be reported as though the weighting settled it.

So: keep the gain rescaling in the main analysis, keep precision weighting as the supplementary
robustness check it already is, and state in the supplement that σᵢ is computed post-rescaling.

---

## Result 5 — the run-mismatch control, quantified without a server pull

[`../cleanroom/run_mismatch_local.m`](../cleanroom/run_mismatch_local.m) and
[`../cleanroom/check_fig4a.m`](../cleanroom/check_fig4a.m).

Fig 4A's mismatched-run comparison was a single anecdote (two runs, one ROI). It can be
generalised from local data. For a **fixed** (not refit) prediction — which is why R² can go
negative — writing ρ = corr(ŷ_c, ŷ_m) and *k* = ‖ŷ_m‖/‖ŷ_c‖:

> R²_m = R²_c · (2ρ*k* − *k*²)

This ratio depends only on the predicted time series, so it is fully local and valid at **any**
averaging level, provided ρ and *k* are computed at that level.

**The sign result needs nothing from the server.** R²_m < 0 whenever *k* > 2ρ. The fixed-seed
designs are near-orthogonal — median ρ = 0.006 (dg), 0.020 (da), with *k* ≈ 1.04 — so the ratio is
**negative for every run pair in all 16 sessions** (range −1.81 to −0.58). The one session worth
checking before claiming universality is **sub-0201/da**, where ρ = 0.36 and only 66% of run pairs
are negative: some of its run-pair designs are genuinely correlated.

**Validated against the two runs already in Fig 4A**, using the figure's own R²_c:

| panel | ρ | *k* | R²_c (fig) | R²_m predicted | R²_m (fig) | gap |
|---|--:|--:|--:|--:|--:|--:|
| run 2 observed, run 5 design | −0.060 | 0.817 | 0.42 | −0.322 | −0.25 | +0.072 |
| run 5 observed, run 2 design | −0.060 | 1.223 | 0.28 | −0.460 | −0.53 | −0.070 |

Agreement to ±0.07 with the two errors **opposite in sign** (mean +0.001), so the neglected term
⟨e, ŷ_m⟩ behaves like noise rather than a systematic bias. (An earlier prediction here — that shared
event timing would make it reliably positive — is not supported by these two points.) The check
also *explains* the figure: the striking difference between −0.25 and −0.53 is almost entirely *k*,
the ratio of predicted amplitudes (0.817 vs its reciprocal 1.223), not a difference in data quality.

**What still needs measured data: absolute ROI-level magnitude.** Converting the ratio to an
absolute R² requires R²_c *at the same averaging level*. GLMsingle's `R2run` is **per vertex**,
while Fig 4A's 0.42/0.28 are for a series averaged over vertices *and* observers — far higher, by
an amount set by the **cross-vertex noise correlation**, which cannot be determined from predictions
alone. An earlier version of `run_mismatch_local.m` scaled its ROI estimate by the median
per-vertex `R2run` (≈0.05 where ≈0.42 was needed, an ~8× error); that column has been removed and
only the scale-free ratio is reported at ROI level.

### Calibration result (sub-0037, dg, 8 runs — 56 run pairs)

Fetched and analysed 2026-08-19 ([`../cleanroom/calibrate_run_mismatch.m`](../cleanroom/calibrate_run_mismatch.m)).
Predictions use the **overall** condition betas (mean across runs), not a run's own, to avoid
overfitting the run being predicted.

**The sign claim is confirmed: measured R²_m is negative in 56 of 56 run pairs**, at both the
per-vertex and ROI-mean level. That is the claim the manuscript needs, and it is now measured
rather than inferred.

**The failing assumption is the other one.** Decomposing the gap (per vertex, normalised by ‖y‖²):

| term | measured | the formula assumes |
|---|--:|---|
| 2⟨e, ŷ_m⟩ / ‖y‖² | **−0.010** (SD 0.037) | 0 |
| 2⟨e, ŷ_c⟩ / ‖y‖² | **−0.079** | 0 |

⟨e, ŷ_m⟩ — the term flagged as the risk — is genuinely negligible. The whole bias comes from
⟨e, ŷ_c⟩, dismissed above as a minor "exact for OLS, approximate for ridge" caveat. **That was
backwards.**

**The fix is to rewrite the formula.** Using *Q* = ‖ŷ_c‖²/‖y‖² in place of R²_c,

> R²_m = *Q* · (2ρ*k* − *k*²)

needs only ⟨e, ŷ_m⟩ = 0, with no assumption about ⟨e, ŷ_c⟩. Here *Q* = 0.131 against
R²_c = 0.051 (ratio 2.54), and substituting it cuts the mean gap from −0.087 to **−0.014** — which
matches the measured ⟨e, ŷ_m⟩ term of −0.010. The accounting closes.

**But the estimator predicts the average, not individual pairs.** Correlation between measured and
predicted R²_m across the 56 pairs is only 0.15 (per vertex) / 0.12 (ROI), and the *Q* rescaling
increases per-pair scatter (gap SD 0.081 → 0.177) even as it fixes the mean. Use it for a group
statement, never to predict a specific run pair. *Q* is also not locally available — it needs
‖y‖² — so the 2.54 ratio is itself a calibration constant that may not transfer across observers.

### What averaging over vertices buys — the blocked quantity, now measured

| quantity | value |
|---|--:|
| measured R²_c, per-vertex median | 0.051 |
| measured R²_c, ROI mean (single observer) | 0.181 |
| **ratio** | **3.52×** |
| mean cross-vertex residual correlation | **0.273** |

Averaging 1232 vertices buys only **3.5×** in R², nowhere near what independent noise would give,
because the residuals are strongly correlated across vertices. This is the quantity that blocks
converting the scale-free ratio into an absolute ROI-level R², and it cannot be derived from
predictions — it had to be measured.

> **Note on comparability.** Measured per-vertex R²_c (0.051) is well below GLMsingle's own
> `R2run` (0.174) for two reasons, neither a bug: GLMsingle fits 52 **single-trial** betas per run
> where this uses 13 condition betas fit across runs, and the extraction projects out only
> polynomials (following `createTTaveTable.m`) while GLMsingle also projects noise-pool PCs. The
> calibration is internally consistent because it uses the measured R²_c throughout, but these
> values are not interchangeable with `R2run`.

**Bottom line for the manuscript.** State the run-mismatch control as the qualitative claim — the
mismatched-design prediction fits worse than the mean in every run pair — which is now measured on
56 pairs plus the 2 in Fig 4A. Quote exact magnitudes only from measured data, for whatever subset
is reported, not from the local estimator.

### Original pull note

**Calibration pull.** `../server_extract/collect_timeseries.m` fetches one session
(sub-0037, dg; 8 runs, 2.35 GB). Completed in ~10 min; per-run times 52–172 s, i.e. an effective
1.7–3.4 MB/s. Measured mount throughput has ranged 0.8–3.4 MB/s
within one session of testing (Abu Dhabi → New York, and a recent router change in the building), so
treat any single throughput measurement as unreliable for planning.

---

## Draft text for the section

*Keep the existing Fig 4A–C paragraphs as illustration; these follow them.*

> Because the reference condition was full-field pink noise rather than a mean-luminance screen, and
> this background was present continuously throughout each run, the variance explained by the GLM
> reflects differentiation among the stimulus conditions rather than overall responsiveness: a
> region that responds strongly but similarly to the gratings and to the pink-noise background
> yields low variance explained. We therefore assessed the quality of the first-stage estimates by
> their reliability — the reproducibility of the response estimates that enter the asymmetry
> analyses — rather than by variance explained.
>
> Each observer's runs were split in half. Because the design was balanced (four trials per
> condition in every run), a 4-vs-4 split of the eight runs is balanced by construction. We
> recomputed the analysed profile — the response to each of the four orientations in each of the
> eight polar-angle ROIs — separately from each half, correlated the two halves, averaged over all
> 35 balanced splits, and applied the Spearman–Brown correction to the full run count. Because
> differences in overall responsiveness among ROIs are reliable but carry no orientation
> information, we report the reliability of the profile after removing each ROI's mean response
> across the four orientations, which leaves exactly the component that carries the four
> asymmetries. This orientation-differential profile had a split-half reliability of r = 0.87
> (range across observers 0.71–0.98) for the Cartesian experiment and r = 0.62 (range 0.50–0.88)
> for the polar experiment, exceeding 0.5 in all 16 sessions (Table SX). We report reliability at
> the level of the ROI profile rather than of individual vertices because each analysed value is a
> summary over more than a thousand vertices; the reliability of a single vertex is correspondingly
> lower (median r = 0.73 and 0.44 for the two experiments) but is not the quantity on which any of
> the following analyses depends.
>
> The same run splits give the measurement error of each asymmetry directly, so each effect can be
> expressed relative to the noise in its own measurement. In the Cartesian experiment the
> horizontal–vertical asymmetry was 4.3 times its within-observer standard error, and the
> cardinal–oblique, radial–tangential and polar cardinal–polar oblique asymmetries were 3.0, 2.8
> and 2.1 times theirs; in the polar experiment the horizontal–vertical asymmetry was 5.3 times its
> standard error and the radial–tangential asymmetry 1.1 times. Individual observers therefore
> resolve the larger asymmetries on their own, which is why the horizontal–vertical difference in
> the Cartesian experiment was negative in all eight observers and at all eight locations. Removing
> the measurement variance entirely and re-testing changes no conclusion: every asymmetry
> significant at p < 0.05 remains so and every null remains null. The exception in the other
> direction is the radial–tangential asymmetry in the polar experiment, whose measurement error is
> comparable to the effect itself (p = 0.042 observed, 0.004 disattenuated); that effect is limited
> by measurement precision rather than by the number of observers, and we treat it as marginal
> throughout.
>
> As an independent check that the GLM recovered known response properties, we examined two regions
> about which this study makes no claims. Area MT responded more strongly to moving than to
> stationary gratings in all 16 sessions (+0.25 to +1.10% signal change), and area V4 responded more
> strongly to gratings than to the pink-noise background in all 16 sessions (+0.06 to +0.82%)
> (Table SY). The pattern is physiologically coherent — MT prefers motion, V4 prefers coherent
> gratings over noise — and confirms that the sessions contain well-differentiated responses,
> including in the two observers whose V1 responses were weakest.

**Fig 4D suggestion.** A split-half scatter — the 32 orientation-differential profile values from
one half against the other, one colour per observer, pooled over the 35 splits or shown for a
representative split — puts the reliability claim in the same "show the data" register as 4A–C.

**Supplementary tables.** Table SX = `splithalf_reliability.csv` (per-observer reliability).
Table SY = the MT / V4 columns of `glm_summary.csv` ([`REPORT.md`](REPORT.md) §2.6).

**Left for Rania: the geometric vs. arithmetic mean gain.** The manuscript should keep the
**geometric** mean (JW's recommendation, which she has agreed to adopt); `lme1_fit.m`,
`plot1_experimentalCond.m` and `plot2_experimentalCond.m` currently use `mean(gainWeights)` and
should be changed to `exp(mean(log(gainWeights)))`. The clean-room already defaults to geometric
(`cfg.gainMean`). The difference is a single scalar common to all observers — 4.3952 vs 4.4417 for
these 8, a factor of 0.9895 — so it shifts every reported effect size by ~1% and leaves every
correlation, variance ratio, *t* and *p* value **exactly** unchanged.

## Figure 4C — no quantitative claim attached

**4–8° is the window where spatial frequency is most closely matched between the two experiments,
not the stimulated extent.** The stimulus spanned roughly **0.5–12°**, so the model should fit
across that whole range and fall off only beyond it — 4–8° is not expected to be diagnostic of
where the model fits better, and no R²-by-eccentricity claim should be built on it. *(A "R²
concentrated in the stimulated 4–8° band" reading of this map was withdrawn 2026-08-19, JW.)*

Tested directly (V1 median GLM R² by eccentricity, dg, all 8 observers):

| ecc (°) | 0–0.5 | 0.5–1 | 1–2 | 2–3 | 3–4 | 4–6 | 6–8 | 8–10 | 10–12 | 12–14 | 14–18 | 18–25 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| median R² | 4.1 | 3.9 | 4.0 | 4.6 | 5.4 | 5.3 | 4.5 | 4.5 | 4.2 | 3.6 | 6.5 | 5.8 |

**The group profile is essentially flat** (3.6–6.5) with no drop beyond 12°, so the map does not
support a spatial-specificity claim either — and on reflection **a flat profile is what should be
expected here**, for three reasons (JW, 2026-08-19):

1. **R² is sign-blind.** A vertex *suppressed* by the stimulus has as much cross-condition variance
   as a driven one and earns the same R². Variance explained therefore cannot localise the
   stimulus even in principle.
2. **pRF extent, not just centre.** A pRF whose centre lies beyond the stimulus can still overlap
   it. Measured: at 14–18°, median σ = 2.17 and **59%** of V1 vertices have (ecc − 2σ) < 12°; at
   18–25°, σ = 2.98 but only **18%** do. So extent accounts for the near periphery, not the far.
3. **Surround suppression, with no way to see it.** Vertices well outside the stimulus would
   normally show a *negative* response, and with a sparser design carrying occasional long blanks
   one would see a positive→negative transition with eccentricity. The rapid event-related design
   plus the absence of any true blank makes that transition essentially invisible. Measured:
   median stationary−blank stays **positive at every eccentricity** and rises rather than falls
   (+0.17 at 4–8°, +0.32 at 14–18°, +0.42 at 18–25°) — no crossing at all.

So the flat profile is a property of the measure and the design, not evidence about the model. Only the two high-R² observers show the expected
shape — sub-0037 and sub-0395 rise to 16–32% across ~3–12° — while the other six sit near 4%
everywhere, which is what the pink-noise reference predicts (little cross-condition variance to
explain). **A further reason not to read anything into the outer bins** (JW, 2026-08-19): the pRF *mapping*
stimulus itself only extended to ~12°. So a vertex with an estimated pRF centre beyond 12° must
nonetheless have responded to the mapping stimulus, which filled the same aperture as the gratings —
it is an extrapolated fit from a vertex that *did* respond, not evidence of a receptive field
outside the stimulated region. Vertices with pRF centres > 12° are therefore **not an unstimulated
control**, and no eccentricity contrast built from them can be. A valid comparison would be against
the rest of the brain with no pRF screener applied at all.

The far-periphery bins also rest on few vertices (median n ≈ 100–120 beyond 14°, against
≈1350 at 4–8°), so their elevation is not reliable. pRF *fit quality* there is adequate (median pRF
R² 0.52–0.56), but median σ is non-monotonic across the outer bins (0.91 → 0.20 → 2.17 → 2.98),
which points to poorly-constrained fits at the V1 label edge rather than a real size progression —
another reason not to lean on 14–25°.

### Would a split-half reliability map do better than the R² map?

Tested per vertex over all V1 vertices, 35 balanced splits, Fisher-z averaged, Spearman–Brown
corrected:

| observer | GLM R² (4–8° / >12°) | split-half reliability (4–8° / >12°) |
|---|---|---|
| sub-0037 | 16.6 / 6.9 — 2.41× | 0.862 / 0.643 — 1.34× |
| sub-0395 | 26.4 / 23.3 — 1.13× | 0.918 / 0.785 — 1.17× |
| sub-0201 | 3.5 / 3.4 — 1.02× | 0.426 / 0.424 — 1.00× |

> **These ratios are not a valid test and should not be quoted** — see the pRF-mapping caveat
> above. The >12° group is not an unstimulated control. They are recorded only to show that
> nothing here motivates building the map.

**As a localiser, reliability is no better than R².** That is expected once the earlier reasoning is
applied to it: reliability is *also* sign-blind, so a reliably *suppressed* vertex scores as high as
a reliably driven one. Restricting to the 4 stimuli of interest does not help, because the problem
was never which conditions were included.

**Where reliability does win is interpretability of magnitude.** In the responsive patch it reads
0.86 and 0.92 — numbers a reader can evaluate directly — where R² reads 16.6% and 26.4%, which look
poor to anyone unaware of the pink-noise reference; and it correctly flags sub-0201 at 0.43. It
answers "are these estimates trustworthy?" rather than "where is the stimulus?", and
trustworthiness is what this section needs to establish. But that is an argument for the reliability
**statistic** (Result 1), not for rendering it on a surface.

**Two practical reasons not to build the map.** `collect_runwise_betas.m` saved **V1 only**
(`v1Index`, 5158–7546 vertices per observer), so a whole-surface reliability map needs `modelmd`
re-pulled — 425 MB × 16 ≈ 6.8 GB. And per-vertex reliability over just 4 orientation values is
intrinsically noisy (raw half-vs-half *r* ranges 0.27–0.85 across observers); it behaves far better
as the ROI-level statistic already reported.

### Recommendation

Either **drop Fig 4C**, or present it as an **illustration for an example observer** with no
quantitative claim attached about how well or where the model fits. Do not spend a pull on a
reliability map. The section's claims rest on the ROI-level reliability statistics (Results 1–2).

A moving-vs-stationary map (`allmValls`, already local for all 16 sessions) *would* be a genuinely
diagnostic surface figure — MT+ positive in 8/8 observers, hV4 near zero, whole-surface median
0.000, and the pink-noise reference cancels exactly since both conditions share it. It is set aside
because **this paper analyses only the stationary conditions**, so a motion figure would be off-topic
here. Worth remembering for a paper where it is not.

---

## What is still open here

Carried in [`../../AGENTS.md`](../../AGENTS.md) §7 along with everything else open in the project:

- **Finish the run-mismatch calibration** (Result 5). The sign result is settled and measured on 56
  run pairs; what a further pull adds is the cross-vertex noise correlation needed for absolute
  ROI-level magnitudes. Check **sub-0201/da** (ρ = 0.36, only 66% of pairs negative) before stating
  the claim over all run pairs.
- **Geometric vs. arithmetic mean gain** — left for Rania, above.
- **Figure 4C** — drop it, or present it as an illustration with no quantitative claim.

*The "merge or triage the 4 upstream commits" item that stood here is closed: a498f08, e2acbf7 and
the rest are in this history, including the `lme1_fit.m` ROI-order fix, the ROI set expanded to 8,
`01_calculate_observer_gain/`, and `02_ttave/createTTaveTable.m`.*
