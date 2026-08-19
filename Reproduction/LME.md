# The linear mixed-effects model (Fig 7) — what it adds, and why it should come out

**Resolution: the LME estimates nothing the beta-subtraction route does not, and its standard error
is wrong in a known direction. Recommend omitting it from the manuscript, reporting the
within-observer error summary instead (see [`GLM_QUALITY.md`](GLM_QUALITY.md)), and pointing to this
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
*wᵢ* = 1/(τ² + σᵢ²):

| exp | asymmetry | equal-w | precision-w | weight spread |
|---|---|---|---|---|
| dg | H−V | −0.4798 | −0.4712 | 2.22× |
| dg | card−obl | −0.2044 | −0.1940 | 1.40× |
| dg | rad−tang | 0.1159 | 0.1216 | 1.90× |
| dg | polc−polo | 0.0292 | 0.0285 | 1.40× |
| da | H−V | −0.2114 | −0.2091 | 1.15× |
| da | card−obl | −0.0296 | −0.0283 | 1.41× |
| da | rad−tang | 0.1501 | 0.1756 | 2.95× |
| da | polc−polo | 0.0336 | 0.0314 | 1.83× |

And on the context differences the claims rest on:

| asymmetry | equal-w | *p* | precision-w | *p* |
|---|---|---|---|---|
| H−V | −0.2685 | 0.024 | −0.2673 | 0.024 |
| card−obl | −0.1748 | 0.006 | −0.1590 | 0.009 |
| rad−tang | −0.0342 | 0.70 | −0.0546 | 0.54 |
| polc−polo | −0.0043 | 0.93 | −0.0065 | 0.90 |

**No conclusion changes.**

**Why the weights compress even though the reliabilities do not.** The weight is 1/(τ² + σᵢ²) and τ²
is common to every observer. The per-observer SEs genuinely differ a great deal — da rad−tang runs
0.067 to 0.251, a 14× spread in variance — but τ² sits underneath all of them as a floor. Since
within-observer error is only 23–39% of the total, τ² is the larger term and the weights can spread
by at most ~3×, not 14×. The general point: **precision weighting is bounded by the
between-to-within variance ratio, not by how much the reliabilities differ.** An observer whose
effect genuinely differs from the group cannot be down-weighted for being noisy, because most of why
they differ is not noise.

Two caveats. τ² is estimated from 7 df and each σᵢ² from 8 runs, so the weights are themselves
noisy; at n = 8, estimated-weight GLS can add variance rather than remove it. And the disattenuation
ceiling (p 0.024 → 0.014 for H−V, see `GLM_QUALITY.md`) bounds what any of this can achieve.

Worth noting: precision weighting down-weights sub-0426 (SE 0.251) and sub-0395 (SE 0.158, and only
6 runs), and sub-0395 is the observer the leave-one-out analysis singles out. So it is a principled
version of that probe — it moves in the same direction and still does not approach significance
(p = 0.54 vs 0.70), which strengthens the "uninformative, not null" reading of the polar-frame
result rather than undermining it.

## 6. Recommendation

Omit the LME from the manuscript. As specified it reports numbers already obtained by subtraction,
with a standard error that is wrong in a known direction; respecified correctly it reports the same
numbers with a standard error that converges on the paired *t*. Report instead:

- the within-observer error summary, in the data-quality section — it shows both that the GLM data
  are good and that precision weighting is not needed (`GLM_QUALITY.md`);
- a one-line statement that a precision-weighted analysis changes no conclusion, pointing to this
  repository for the numbers above.

## Code

`cleanroom/fit_lme_fig7.m` (the model as published), `cleanroom/lme_codes.m` (the asymmetry codes),
`cleanroom/compute_asymmetries.m` (the subtraction route), `cleanroom/diagnose_context_asymmetry.m`
(the within-subject tests and the DF comparison), `cleanroom/diagnose_within_observer_error.m`
(per-observer reliabilities from runs, feeding the precision weights).
