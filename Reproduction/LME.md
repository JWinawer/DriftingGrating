# The linear mixed-effects model (Fig 7) — what it adds, and why it should come out

**Resolution: the LME estimates nothing the beta-subtraction route does not, and its standard error
is wrong in a known direction. Recommend omitting it from the manuscript, reporting the
within-observer error summary instead (see [`_archive/GLM_QUALITY.md`](_archive/GLM_QUALITY.md)), and pointing to this
repository for the precision-weighted check.**

## 1. The estimates are identical, not merely similar

Both main-text routes reduce to the same number. Running the ROI route and `fit_lme_fig7` over the
same `M` from `bin_and_aggregate`, **all four asymmetries in both experiments agree to < 2e-16**:

| exp | asymmetry | wedge-average route | LME fixed effect |
|---|---|---|---|
| dg | H−V | −0.479829 | −0.479829 |
| dg | card−obl | −0.204375 | −0.204375 |
| dg | rad−tang | 0.115887 | 0.115887 |
| dg | polc−polo | 0.029248 | 0.029248 |
| da | H−V | −0.211370 | −0.211370 |
| da | card−obl | −0.029590 | −0.029590 |
| da | rad−tang | 0.150079 | 0.150079 |
| da | polc−polo | 0.033570 | 0.033570 |

Two conditions make this exact, and both hold:

- the four codes from `lme_codes` are **exactly orthogonal**, to each other and to the intercept
  (cross-product matrix diagonal, 16/32/16/32; every column sums to 0, both experiments);
- the design is **completely balanced** — no empty wedge/observer cells in either experiment,
  minimum 22 vertices per wedge.

Under balance plus orthogonality, GLS with compound symmetry, OLS, and the per-observer mean all
coincide for within-observer contrasts. Note the dependence on balance: deleting a single wedge for
a single observer makes the routes diverge by ~1.2e-2. Not a concern at the published inclusion
criteria; it would become one if a stricter filter emptied a wedge.

(The two routes as *published* also differ in taking the median vs the mean within each wedge. That
is a separate, accidental difference and should be made consistent; in this repository both routes
are fed the same wedge median from `bin_and_aggregate`.)

## 2. The random intercept cannot affect the four asymmetries

With `(1|subject)` the only thing shrunk is the observer intercept, and each asymmetry contrast sums
to zero within an observer, so the intercept cancels out of it exactly — shrunk or not. There is
nothing for the Gaussian prior to act on.

Consequently **treating observer as a fixed effect changes nothing**: fixed-observer OLS returns the
same four estimates to ~1e-16, spending 7 more DF on intercepts (DFE 244 vs 251) while still fitting
one asymmetry slope for everyone. Random-vs-fixed is not the lever.

## 3. Where the routes do differ: the standard error

| route | DF | SE(2b) | p |
|---|---|---|---|
| LME `(1\|subject)`, residual DF | 502 | — | 0.0009 |
| paired *t* across the 8 observers | 7 | — | 0.024 |

The LME tests the effect against wedge-level observations rather than against observers, with no
Satterthwaite or Kenward–Roger correction, giving *p* smaller by 5–25×. The lever is **whether the
asymmetry slope varies by observer, and how DF is computed**: fitting random slopes and querying
`fixedEffects(lme,'DFMethod','satterthwaite')` gives DF ≈ n−1 and very nearly reproduces the paired
test. `fitlme`'s default *residual* DF reports the wedge-level number whatever the random structure.

## 4. Would trial-wise input make it compute something new?

GLMsingle retains a beta per vertex per trial; these are averaged to condition-wise betas. Feeding
the trial-level values (averaged across vertices within a wedge) to the LME was tested by simulation
at the real design size (8 observers × 8 wedges × 4 orientations × 32 trials = 8,192 rows):

| specification | DF | SE(2b) | p | estimate |
|---|---|---|---|---|
| trial LME `(1\|subject)` | 8187 | 0.0199 | 3e−126 | −0.483248 |
| trial LME `(1\|subject)+(1\|run)` | 8187 | 0.0191 | 1e−135 | −0.483248 |
| trial LME random slopes + `(1\|run)`, Satterthwaite | 9.2 | 0.0296 | 4e−08 | −0.483248 |
| paired *t* across 8 observers | 7 | 0.0295 | 8e−07 | −0.483248 |

**The estimate does not move.** Balance and orthogonality do the same work at trial level. Feeding
trials to the current specification makes the anti-conservatism far worse (DF 502 → 8,187), and
adding `(1|run)` does not rescue the fixed-effect test. Only random slopes fix it, and then the
answer lands on the paired *t*.

The one genuinely new thing is a variance decomposition, since condition-wise input has no
replication at all — one number per cell, so the residual is not measurement error but everything
at once. In simulation the trial-level model recovers trial SD 0.606 (true 0.600) and run SD 0.184
(true 0.200), but the between-observer SD of the slope 0.032 (true 0.050): **the quantity it would
uniquely add is the one it estimates worst**, because it rests on 8 groups.

**Caution on trial-level independence.** GLMsingle single-trial betas are not independent
replicates: trials within a run share the HRF fit, nuisance regressors, drift and motion, and the
fractional-ridge shrinkage couples them at each vertex. Averaging across vertices collapses the
spatial dimension and leaves the trial coupling untouched. This is the same class of error as
resampling vertices, one level up — which is why `diagnose_within_observer_error.m` resamples
**runs**, which are much closer to exchangeable.

## 5. Should noisier observers count less? (precision weighting)

They should in principle, and **no version of the LME as fed can do it**. With one number per cell
there is no replication, so σ² is a single shared parameter and a reliable observer is
algebraically indistinguishable from a noisy one. Trial-level input does not fix it either while
the residual is homoscedastic and the design balanced — every observer's cell mean then has the same
nominal precision by construction. It needs a per-observer residual variance (`varIdent` in R's
`nlme`; `fitlme` will not do it directly).

Done properly, with per-observer reliabilities measured across runs and weights
*wᵢ* = 1/(τ² + σᵢ²). **This is the table to report** — it is generated by
`cleanroom/precision_weighted_table.m` and written to
`supplement/precision_weighted.csv`. Intervals are 95% on 7 df.

> **Regenerated 2026-08-19 on the published route (across-vertex MEAN + observer pRF-gain
> rescaling).** The two tables below were computed with the across-vertex *median* and no gain
> rescaling. Both inputs have since changed (see
> [`local_qc/GLM_SUMMARY_SECTION.md`](local_qc/GLM_SUMMARY_SECTION.md)), so the current numbers
> are those in `supplement/precision_weighted.csv`, reproduced here; the superseded tables are
> kept below. **The conclusion is unchanged**: precision weighting moves no context effect and no
> Cartesian-frame asymmetry, and da rad−tang remains the one cell whose *status* it changes.
>
> | exp | asymmetry | equal-weighted [95% CI] | precision-weighted [95% CI] | τ | mean σᵢ | *w* max/min |
> |---|---|---|---|---|---|---|
> | dg | horiz−vert | −0.547 [−0.679, −0.416] | −0.549 [−0.669, −0.428] | 0.086 | 0.132 | 3.63× |
> | dg | card−obl | −0.221 [−0.354, −0.089] | −0.213 [−0.345, −0.082] | 0.138 | 0.078 | 1.39× |
> | dg | rad−tang | 0.104 [0.037, 0.171] | 0.112 [0.047, 0.178] | 0.070 | 0.040 | 1.89× |
> | dg | polc−polo | 0.039 [0.002, 0.076] | 0.040 [0.004, 0.077] | 0.040 | 0.019 | 1.38× |
> | da | horiz−vert | −0.218 [−0.365, −0.071] | −0.214 [−0.361, −0.068] | 0.170 | 0.046 | 1.23× |
> | da | card−obl | −0.029 [−0.066, 0.009] | −0.032 [−0.068, 0.005] | 0.039 | 0.021 | 1.80× |
> | da | rad−tang | 0.150 [0.007, 0.292] | **0.176 [0.065, 0.286]** | 0.072 | 0.154 | 9.89× |
> | da | polc−polo | 0.040 [−0.048, 0.127] | 0.040 [−0.044, 0.124] | 0.077 | 0.071 | 2.69× |
>
> | asymmetry | equal-weighted [95% CI] | *p* | precision-weighted [95% CI] | *p* | *w* max/min |
> |---|---|---|---|---|---|
> | horiz−vert | −0.329 [−0.528, −0.130] | 0.0058 | −0.335 [−0.530, −0.140] | 0.0048 | 1.71× |
> | card−obl | −0.193 [−0.314, −0.071] | 0.0073 | −0.183 [−0.303, −0.062] | 0.0089 | 1.49× |
> | rad−tang | −0.046 [−0.224, 0.132] | 0.56 | −0.061 [−0.219, 0.097] | 0.39 | 4.84× |
> | polc−polo | −0.001 [−0.103, 0.102] | 0.99 | −0.001 [−0.102, 0.099] | 0.97 | 2.14× |
>
> Note `w` max/min for da rad−tang rises from 2.95× to 9.89×: the gain rescaling and the
> precision weighting interact (see the "double weighting" note in
> [`local_qc/GLM_SUMMARY_SECTION.md`](local_qc/GLM_SUMMARY_SECTION.md)).

**The four asymmetries, per experiment** (% signal change) — *superseded, median route:*

| exp | asymmetry | equal-weighted [95% CI] | precision-weighted [95% CI] | τ | mean σᵢ | *w* max/min |
|---|---|---|---|---|---|---|
| dg | horiz−vert | −0.480 [−0.647, −0.313] | −0.471 [−0.634, −0.308] | 0.155 | 0.126 | 2.22× |
| dg | card−obl | −0.204 [−0.334, −0.074] | −0.194 [−0.323, −0.065] | 0.134 | 0.078 | 1.40× |
| dg | rad−tang | 0.116 [0.045, 0.187] | 0.122 [0.052, 0.191] | 0.072 | 0.045 | 1.90× |
| dg | polc−polo | 0.029 [−0.017, 0.075] | 0.028 [−0.017, 0.074] | 0.049 | 0.025 | 1.40× |
| da | horiz−vert | −0.211 [−0.338, −0.085] | −0.209 [−0.336, −0.082] | 0.145 | 0.044 | 1.15× |
| da | card−obl | −0.030 [−0.068, 0.009] | −0.028 [−0.067, 0.010] | 0.040 | 0.024 | 1.41× |
| da | rad−tang | 0.150 [−0.024, 0.324] | **0.176 [0.011, 0.340]** | 0.160 | 0.134 | 2.95× |
| da | polc−polo | 0.034 [−0.062, 0.129] | 0.031 [−0.062, 0.125] | 0.090 | 0.070 | 1.83× |

**The four context effects** (dg − da) — *superseded, median route:*

| asymmetry | equal-weighted [95% CI] | *p* | precision-weighted [95% CI] | *p* | *w* max/min |
|---|---|---|---|---|---|
| horiz−vert | −0.268 [−0.490, −0.047] | 0.024 | −0.267 [−0.487, −0.048] | 0.024 | 1.62× |
| card−obl | −0.175 [−0.280, −0.069] | 0.006 | −0.159 [−0.263, −0.055] | 0.009 | 1.67× |
| rad−tang | −0.034 [−0.239, 0.170] | 0.70 | −0.055 [−0.252, 0.143] | 0.54 | 2.32× |
| polc−polo | −0.004 [−0.119, 0.110] | 0.93 | −0.007 [−0.119, 0.106] | 0.90 | 1.56× |

**Every context effect is unchanged**, which is what the manuscript's claims rest on, and seven of
the eight asymmetries are unchanged. The exception is **da rad−tang**, which the weighting moves
from 0.150 [−0.024, 0.324] to 0.176 [0.011, 0.340] — i.e. from including zero to excluding it. That
cell is marginal by every route (see §7) and should not be reported as though weighting settled it;
it is also the cell with the most heterogeneous reliability (2.95×), and the reweighting acts by
down-weighting sub-0426 (σ = 0.251) and sub-0395 (σ = 0.158, 6 runs rather than 8) — sub-0395 being
the observer the leave-one-out analysis already singles out. So precision weighting is a principled
version of that probe. Note it does **not** make the rad−tang *context* effect significant
(p = 0.54 vs 0.70), which is the claim that matters.
**Why the weights compress even though the reliabilities do not.** Notation, as in random-effects
meta-analysis: **σᵢ²** is observer *i*'s within-observer measurement variance (how much their own
estimate would move on re-measurement, here from resampling runs), and **τ²** is the
between-observer variance of the *true* effects — the heterogeneity, i.e. the random-effect variance
component, a single number shared by the whole sample. It is obtained by subtraction,
τ̂² = var(yᵢ) − mean(σᵢ²).

An observer's measured value is not a noisy reading of the group mean; it is a noisy reading of
*their own* true effect, which is itself a draw from a population that genuinely varies. So its
variance as an estimator of the group mean is τ² + σᵢ², which is where the weight comes from.

Worked, for da rad−tang (the most heterogeneous case): τ̂² = 0.0256 (τ = 0.160), mean σᵢ² = 0.0179.

```
sigma_i^2       ranges 0.0045 -> 0.0630     ratio 14.1x
tau^2 + sigma_i^2  ranges 0.0300 -> 0.0886     ratio  2.95x
```

Adding a constant to every term shrinks the ratios toward 1, and here that constant is larger than
the mean σᵢ². So a 14× spread in reliability becomes a 3× spread in weight. Were τ² zero — every
observer sharing one true effect, all disagreement being noise — the weights would spread the full
14.1×, and the noisiest observer would drop from 5.5% to 1.8%.

Two consequences. A perfectly measured observer still has a weight ceiling of 1/τ², because one
exact draw from a varying population still does not pin down its mean. And a noisy observer is not
worthless: sub-0426 falls from 12.5% under equal weighting to 5.5%, real down-weighting, but bounded.
**You can only discount the part of an observer's deviation that is noise, and here most of it is
not.** Note this is compression, not a fixed cap — the ratio is (τ² + σ²ₘₐₓ)/(τ² + σ²ₘᵢₙ), always
smaller than σ²ₘₐₓ/σ²ₘᵢₙ when τ² > 0, but an arbitrarily noisy observer would still approach zero
weight.

Two caveats. **τ̂² is estimated from 7 df**, so the entire weighting scheme rests on a quantity this
design cannot measure well; each σᵢ² comes from only 8 runs. At n = 8, estimated-weight GLS can add
variance rather than remove it. And the disattenuation ceiling (p 0.024 → 0.014 for H−V, see
`_archive/GLM_QUALITY.md`) bounds what any of this can achieve.

Worth noting: precision weighting down-weights sub-0426 (SE 0.251) and sub-0395 (SE 0.158, and only
6 runs), and sub-0395 is the observer the leave-one-out analysis singles out. So it is a principled
version of that probe — it moves in the same direction and still does not approach significance
(p = 0.54 vs 0.70), which strengthens the "uninformative, not null" reading of the polar-frame
result rather than undermining it.
### How it is computed, and why not with `fitlme`

The precision-weighted estimate is computed **in closed form** in
`cleanroom/precision_weighted_table.m`, not through MATLAB's mixed-model functions. This is forced,
not a shortcut: `fitlme` has no mechanism for **known, per-observer** measurement variances. It
estimates a single shared residual variance, and its `'Weights'` argument supplies *relative*
precisions that are then multiplied by that estimated variance — so the σᵢ² measured from runs
cannot enter as fixed quantities. On da rad−tang:

| route | estimate | |
|---|---|---|
| closed form (random-effects, σᵢ² known) | **0.1756** | the correct estimator |
| `fitlme y~1+(1\|obs)`, `Weights` = 1/σᵢ² | 0.1995 | the τ² = 0 case — over-weights clean observers |
| `fitlme y~1+(1\|obs)`, no weights | 0.1501 | equal weighting |

`fitlme` returns either extreme but not the correct intermediate. (With one row per observer the
random intercept is also not separable from the residual, so that formulation is degenerate
regardless — it simply returns the weighted least-squares answer, and reports an estimated residual
variance of 2.06 where the meta-analytic model fixes the scale at 1.)

What the closed form implements is the standard **random-effects meta-analysis** estimator: observed
effects with known sampling variances, plus a between-study variance. R has it as `metafor::rma`;
MATLAB has no built-in equivalent. This matters for how the analysis is described — it is not
"we declined to fit a mixed model", it is "the mixed model that can use measured within-observer
error is a random-effects meta-analysis, which MATLAB's LME functions cannot express."

**The τ² estimator barely matters.** The table uses a simple method-of-moments estimate,
τ̂² = max(0, var(yᵢ) − mean(σᵢ²)). The classic DerSimonian–Laird form uses a weighted *Q* statistic
and gives noticeably different τ̂² in places (da rad−tang: 0.160 vs 0.116) but almost identical
weighted means (0.1756 vs 0.1825; largest discrepancy across all eight asymmetries is 0.007). The
conclusion does not depend on the choice.


## 6. Recommendation

Omit the LME from the manuscript. As specified it reports numbers already obtained by subtraction,
with a standard error that is wrong in a known direction; respecified correctly it reports the same
numbers with a standard error that converges on the paired *t*. Report instead:

- the within-observer error summary, in the data-quality section — it shows both that the GLM data
  are good and that precision weighting is not needed (`_archive/GLM_QUALITY.md`);
- a one-line statement that a precision-weighted analysis changes no conclusion, pointing to this
  repository for the numbers above.
## 7. ⚠️ Bootstrap vs *t* intervals disagree on two polar asymmetries

Surfaced while building the §5 table, and it needs a decision before submission. The harmonic model
reports **bootstrap percentile** intervals throughout (`HARMONIC_MODEL.md`, supplement §S3: "95%
confidence intervals from 1,000 bootstrap resamples of the eight observers"). On the same
per-observer values, a *t* interval on 7 df disagrees for two of the eight asymmetries — both in the
polar experiment:

| exp | asymmetry | mean | bootstrap pct 95% | *t* interval (7 df) | *t* p |
|---|---|---|---|---|---|
| dg | horiz−vert | −0.480 | [−0.613, −0.353] | [−0.647, −0.313] | <0.001 |
| dg | card−obl | −0.204 | [−0.305, −0.105] | [−0.334, −0.074] | 0.007 |
| dg | rad−tang | 0.116 | [0.064, 0.174] | [0.045, 0.187] | 0.006 |
| dg | polc−polo | 0.029 | [−0.005, 0.066] | [−0.017, 0.075] | 0.177 |
| da | horiz−vert | −0.211 | [−0.308, −0.112] | [−0.338, −0.085] | 0.006 |
| da | **card−obl** | −0.030 | **[−0.062, −0.002]** | **[−0.068, 0.009]** | **0.114** |
| da | **rad−tang** | 0.150 | **[0.004, 0.279]** | **[−0.024, 0.324]** | **0.081** |
| da | polc−polo | 0.034 | [−0.041, 0.105] | [−0.062, 0.129] | 0.432

The Cartesian experiment agrees on all four. The two polar cells that flip are exactly the ones the
supplement currently reports as non-zero:

- §S5.1 gives polar rad−tang as 0.162 [0.025, 0.264] and polar card−obl as −0.048 [−0.102, −0.004];
- §S5.3 says of radial−tangential "0.119 ... for Cartesian gratings and 0.162 ... for polar
  gratings — **both clearly non-zero**".

"Clearly" is not supported for the polar case: *p* = 0.081 by a *t*-test. The per-observer values are
0.354, 0.162, 0.263, 0.237, 0.047, **−0.244**, 0.001, 0.383 — one strong negative (sub-0395) in a
sample of eight. Percentile bootstrap intervals are known to have poor coverage at n = 8, and the
percentile method does not account for uncertainty in the spread; the *t* interval is the more
conservative and more standard choice at this sample size.

**What this does and does not touch.** It does *not* touch the context effects — those were already
tested with a paired *t* (horiz−vert p = 0.024, card−obl p = 0.006) and are unaffected. It does not
touch any Cartesian-frame result. What it touches is the secondary statement that the polar-frame
asymmetries are individually non-zero in the polar experiment, which supports the reading that
radial/tangential is a robust local property present in both experiments. That statement should
either be softened or restated against a *t* interval.

Reproduce with the disagreement check in `cleanroom/precision_weighted_table.m` and the per-observer
values in `diagnose_within_observer_error.m`.


## Code

`cleanroom/fit_lme_fig7.m` (the model as published), `cleanroom/lme_codes.m` (the asymmetry codes),
`cleanroom/compute_asymmetries.m` (the subtraction route), `cleanroom/diagnose_context_asymmetry.m`
(the within-subject tests and the DF comparison), `cleanroom/diagnose_within_observer_error.m`
(per-observer reliabilities from runs, feeding the precision weights).
