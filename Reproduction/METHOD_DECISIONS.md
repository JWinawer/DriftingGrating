# Closed methodological decisions

Five choices that were argued out and are now settled. Each is recorded because reversing it, or
re-deriving it, would waste time — and because two of them (§1, §5) change what can be claimed.

The analysis these feed into is [`SPECIFICATION.md`](SPECIFICATION.md); the numbers are in
[`RESULTS.md`](RESULTS.md).

---

## 1. Raw percent signal change, not z-scored *(decided 2026-07-24)*

Z-scoring here meant dividing each vertex's 13 GLM betas by `beta_std`, the standard deviation across
those 13 conditions, before any analysis (the `normalize` flag in
`01_process_singlesubjectGLM/main_singlesub.m`, now `0`). The intent was to remove per-observer BOLD
gain.

**Why it is wrong for this dataset.**

- **The "blank" is full-field pink noise, not a mean-luminance baseline**
  ([`local_qc/DATA_QUALITY.md`](local_qc/DATA_QUALITY.md) §1). A divisor built from the spread across
  13 conditions is a gain estimate only if those responses are anchored to a true baseline. They are
  not, so `beta_std` conflates overall BOLD gain with orientation-tuning strength and motion
  sensitivity. Two observers with equal gain but different tuning get different `beta_std`, and
  z-scoring forces them to match — erasing a real difference as though it were nuisance gain.
- **It reintroduces the blank into quantities that were free of it.** The orientation asymmetries are
  differences between stimulus conditions, so the blank cancels exactly. Dividing by `beta_std` puts
  the pink-noise blank and the motion conditions back into the denominator.
- **No valid substitute normaliser exists in this data.** The retinotopy model stores no gain map,
  and no divisor available from the 13 conditions is simultaneously effect-independent, positive for
  every observer, and stable across sessions — blank-referenced gain estimates are ≤ 0 for sub-0037
  and sub-0201 in the polar experiment. Std-based divisors stay positive only because a standard
  deviation always is, which substitutes the noise level for a gain.
- **Precision was never the argument for it.** Between-observer variance exceeds within-observer
  measurement variance for every asymmetry (roughly 3–17×, measured across runs), so
  inverse-variance weighting converges on near-equal weighting and changes nothing.

**The one conclusion that turns on this.** The radial/tangential **context** effect reverses under
z-scoring: in the z-scored variant the difference between experiments reaches significance, in the
raw variant it does not. The raw analysis is adopted for the reasons above, which are independent of
that result — but the dependence should be stated, not buried. The z-scored variant still runs
(`run_harmonic_model(true)`, plus `cleanroom/diagnose_zscore_fig7.m`,
`compare_subject_weighting.m`, `diagnose_exclusion_x_normalization.m`); the mechanism of the reversal
is a between-observer reweighting, described in [`_archive/ZSCORE_FIG7.md`](_archive/ZSCORE_FIG7.md)
§1–§5.

**Consequence for the draft.** Methods drop the "beta weights for each vertex were standardized"
statement and the σ-unit in-text statistics; figures use the raw variants. The pink-noise caveat
applies to raw and z-scored analyses equally, so it is not a mark against the raw one.

---

## 2. Mean across vertices, with per-observer gain rescaling *(decided 2026-08-19)*

`meanWithinLabel.m` saves both `meanBOLDpa` and `medianBOLDpa`. Four candidate routes were compared
against the eight values in the draft, summing absolute deviation over all eight:

| route | Σ&#124;route − draft&#124; |
|---|--:|
| median across vertices, unweighted | 0.128 |
| median across vertices, gain-weighted | 0.142 |
| mean across vertices, unweighted | 0.040 |
| **mean across vertices, gain-weighted** | **0.015** |

**Mean + gain weighting matches every value in the draft to ±0.003.** The aggregator does the work —
the median route misses `dg` horiz−vert by 0.07 — and the gain weighting is a small refinement that
closes the residual on card−obl and polc−polo.

`lme1_fit.m` was already reading the mean; `plot_NeuralAsymmetries.m` passed the median and has been
switched, so Figures 5/6 and Figure 7 are no longer computed from different aggregates. That mismatch
was the whole of the apparent disagreement between the two routes.

**Known cosmetic wart, left alone deliberately.** `plot1_/plot2_experimentalCond.m` and the
`compute_derivative*` functions still name their first parameter `medianBOLDpa` while receiving mean
data. Renaming was skipped because upstream also edits those files.

**Still to change upstream: geometric, not arithmetic, mean gain.** The Methods say the group
statistics are multiplied by the *geometric* mean gain; `lme1_fit.m`, `plot1_experimentalCond.m` and
`plot2_experimentalCond.m` use `mean(gainWeights)`. For these 8 observers that is 4.4417 against
4.3952 — a scalar common to every observer, so it shifts every effect size by ~1% and leaves every
correlation, variance ratio, *t* and *p* **exactly** unchanged. The clean-room already defaults to
geometric.

---

## 3. No mixed model; Figure 7 is removed *(decided 2026-08-19)*

The Figure-7 LME estimates nothing the beta-subtraction route does not, and its standard error is
wrong in a known direction.

**The estimates are identical, not merely similar.** Running the ROI route and `fit_lme_fig7` over
the same input, all four asymmetries in both experiments agree to **< 2e-16**. Two conditions make
this exact, and both hold: the four asymmetry codes are exactly orthogonal to each other and to the
intercept (Gram matrix `diag(16, 32, 16, 32)`), and the design is completely balanced. Under balance
plus orthogonality, GLS with compound symmetry, OLS, and the per-observer mean all coincide for
within-observer contrasts. *(The dependence on balance is real — deleting a single wedge for a single
observer makes the routes diverge by ~1.2e-2. Not a concern at these inclusion criteria.)*

**The random intercept cannot affect the asymmetries.** Each asymmetry contrast sums to zero within
an observer, so the intercept cancels exactly — shrunk or not. Treating observer as a fixed effect
returns the same four estimates to ~1e-16. Random-vs-fixed is not the lever.

**Where the routes differ is the standard error, and the LME's is anti-conservative.**

| route | DF | *p* (dg horiz−vert) |
|---|---|---|
| LME `(1\|subject)`, residual DF | 502 | 0.0009 |
| paired *t* across the 8 observers | 7 | 0.024 |

The LME tests the effect against wedge-level observations rather than against observers, with no
Satterthwaite or Kenward–Roger correction — *p* smaller by 5–25×. Fitting random slopes and querying
`DFMethod','satterthwaite'` gives DF ≈ n−1 and very nearly reproduces the paired test; `fitlme`'s
default *residual* DF reports the wedge-level number whatever the random structure. **The LME
*p*-values must not be quoted.**

**Trial-level input does not rescue it.** Simulated at the real design size (8 observers × 8 wedges ×
4 orientations × 32 trials): the estimate does not move at all, and the anti-conservatism gets far
worse (DF 502 → 8,187). Only random slopes fix it, and then the answer lands on the paired *t*. The
one genuinely new thing trial-level input offers is a variance decomposition — and the component it
would uniquely add (between-observer SD of the slope) is the one it estimates worst, because it rests
on 8 groups. **Caution:** GLMsingle single-trial betas are not independent replicates — trials within
a run share the HRF fit, nuisance regressors, drift and motion, and fractional-ridge shrinkage
couples them at each vertex. This is the same class of error as resampling vertices, one level up,
which is why `diagnose_within_observer_error.m` resamples **runs**.

**What replaces it.** Report the within-observer error summary in the data-quality section
([`local_qc/RELIABILITY.md`](local_qc/RELIABILITY.md)) and a one-line statement that a
precision-weighted analysis changes no conclusion, pointing here for the numbers.

**All context tests are within subject** — form the per-observer difference first, then test across
the 8 observers. That is the correct analysis for a balanced orthogonal within-subject design, and it
is what is reported everywhere.

**One consequence is still open.** Figure 8 plots this same LME's estimate in red against the
observer data in black, so whether Figure 8 keeps the overlay has not been decided
([`../AGENTS.md`](../AGENTS.md) §5).

---

## 4. Precision weighting: reported, not primary

Should noisier observers count less? In principle yes — and **no version of the LME as fed can do
it.** With one number per cell there is no replication, so σ² is a single shared parameter and a
reliable observer is algebraically indistinguishable from a noisy one. Trial-level input does not fix
it either while the residual is homoscedastic and the design balanced. It needs a per-observer
residual variance, which `fitlme` will not express.

**What is done instead.** Per-observer reliabilities σᵢ measured across runs, weights
*wᵢ* = 1/(τ² + σᵢ²), computed **in closed form** in `cleanroom/precision_weighted_table.m`. This is
the standard **random-effects meta-analysis** estimator — observed effects with known sampling
variances plus a between-study variance. R has it as `metafor::rma`; MATLAB has no built-in
equivalent. So the honest description is not "we declined to fit a mixed model", it is "the mixed
model that can use measured within-observer error is a random-effects meta-analysis, which MATLAB's
LME functions cannot express." Confirming this on `da` rad−tang: the closed form gives 0.176;
`fitlme` with `Weights = 1/σᵢ²` gives 0.200 (the τ² = 0 case, over-weighting clean observers) and
without weights 0.150 (equal weighting) — it returns either extreme but not the correct intermediate.

**Why the weights compress even though the reliabilities do not.** σᵢ² is observer *i*'s
within-observer measurement variance; τ² is the between-observer variance of the *true* effects, a
single number shared by the sample, obtained by subtraction τ̂² = var(yᵢ) − mean(σᵢ²). An observer's
measured value is not a noisy reading of the group mean; it is a noisy reading of *their own* true
effect, which is itself a draw from a population that genuinely varies. So its variance as an
estimator of the group mean is τ² + σᵢ². Worked, for `da` rad−tang (the most heterogeneous case):

```
sigma_i^2          ranges 0.0045 -> 0.0630     ratio 14.1x
tau^2 + sigma_i^2  ranges 0.0300 -> 0.0886     ratio  2.95x
```

Adding a constant to every term shrinks the ratios toward 1, and here that constant exceeds the mean
σᵢ². **A 14× spread in reliability becomes a 3× spread in weight.** Two consequences: a perfectly
measured observer still has a weight ceiling of 1/τ², because one exact draw from a varying
population still does not pin down its mean; and a noisy observer is not worthless — sub-0426 falls
from 12.5% to 5.5%, real down-weighting but bounded. **You can only discount the part of an
observer's deviation that is noise, and here most of it is not.**

**Why it is not primary.** Two reasons. τ̂² is estimated from 7 df and each σᵢ² from only 8 runs, so
at n = 8 estimated-weight GLS can add variance rather than remove it — the simulation and the null
distribution behind that judgement are in [`SPECIFICATION.md`](SPECIFICATION.md) §6. And empirically
it moves no context effect and no Cartesian-frame asymmetry ([`RESULTS.md`](RESULTS.md) §8). The one
cell whose *status* it changes is `da` rad−tang, which is marginal by every route and should not be
reported as though weighting settled it.

**The τ² estimator barely matters.** The tables use method-of-moments,
τ̂² = max(0, var(yᵢ) − mean(σᵢ²)). DerSimonian–Laird gives noticeably different τ̂² in places
(`da` rad−tang: 0.160 vs 0.116) but almost identical weighted means (largest discrepancy across all
eight asymmetries: 0.007).

**Stacking gain rescaling and precision weighting is defensible, and mild.** They are different kinds
of operation — the gain rescaling is a *data transformation* removing a nuisance multiplicative
factor measured from an independent protocol; the precision weighting is an *estimator choice* about
how to combine already-comparable estimates. Scale and precision are separate problems, so this is
not weighting twice for the same reason. They also partially cancel by construction: a low-gain
observer is scaled up, which inflates their effect *and* their noise, and precision weighting then
pulls them back down.

---

## 5. *t* intervals, not percentile bootstrap

At n = 8 the percentile bootstrap has poor coverage and does not account for uncertainty in the
spread. ***t* on n−1 df is primary throughout.** Both are computed and both appear in every CSV, with
a `ci_methods_disagree` flag, so the choice stays visible rather than being made silently.

Under the settled specification exactly one V1 cell disagrees — `da` card−obl, bootstrap
[−0.083, −0.002] excluding zero against *t* [−0.092, 0.011], *p* = .105. Report the *t* reading.
(`da` rad−tang used to disagree too; the specification's per-map gain rescaling tightened it and both
intervals now exclude zero, *p* = .028.)

**What this touches.** Not the context effects — those were always tested with a paired *t*. Not any
Cartesian-frame result. It touches only the secondary statement that the polar-frame asymmetries are
individually non-zero in the polar experiment.

**Note for the supplement.** `supplement/SUPPLEMENT_harmonic_model.md` reports percentile bootstrap
intervals throughout, since that is what the harmonic-model code computes. Where a claim there turns
on an interval excluding zero — the polar experiment's cardinal−oblique in particular — the *t*
reading in [`RESULTS.md`](RESULTS.md) §2 is the one to report.

---

## Code

`cleanroom/fit_lme_fig7.m` (the model as the draft specifies it), `lme_codes.m` (the asymmetry
codes), `compute_asymmetries.m` (the subtraction route), `diagnose_context_asymmetry.m` (the
within-subject tests and the DF comparison), `diagnose_within_observer_error.m` (per-observer
reliabilities from runs), `precision_weighted_table.m` / `precision_weighted_cells.m` (the
closed-form estimator), `fit_cell_meta.m`.
