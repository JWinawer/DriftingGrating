# Why we do not z-score

**Decided 2026-07-24. This is the standing summary; it replaces the working documents, which are in
[`_archive/`](_archive/README.md).** All analyses use raw beta weights in percent signal change.

## What z-scoring meant here

Each vertex's 13 GLM betas were divided by `beta_std`, the standard deviation across those 13
conditions, before any analysis (`01_process_singlesubjectGLM/main_singlesub.m`, the `normalize`
flag). The intent was to remove per-observer BOLD gain, which is scanner- and session-dependent and
not of interest.

## Why it is wrong for this dataset

**The "blank" condition is full-field pink noise, not a mean-luminance baseline.** This is the fact
that settles it ([`local_qc/REPORT.md`](local_qc/REPORT.md) §1). A divisor built from the spread
across 13 conditions is a gain estimate only if those responses are anchored to a true baseline.
They are not, so `beta_std` is the spread among *contrast-pattern* responses — it conflates overall
BOLD gain with orientation-tuning strength and motion sensitivity. Two observers with equal gain but
different tuning get different `beta_std`, and z-scoring forces them to match, erasing a real
difference as though it were nuisance gain.

**It reintroduces the blank into quantities that were free of it.** The orientation asymmetries are
differences between stimulus conditions, so the blank cancels exactly. Dividing by `beta_std` puts
the pink-noise blank and the motion conditions back into the denominator. Z-scoring can only add
distortion.

**No valid substitute normalizer exists in this data.** The retinotopy model stores no gain map, and
no divisor available from the 13 conditions is simultaneously effect-independent, positive for every
observer, and stable across sessions — blank-referenced gain estimates are ≤ 0 for sub-0037 and
sub-0201 in the polar experiment. Std-based divisors stay positive only because a standard deviation
always is, which substitutes the noise level for a gain. The appropriate default is not to normalize.

**Precision was never the argument for it.** Between-observer variance exceeds within-observer
measurement variance for every asymmetry (roughly 3–17×, measured across runs), so inverse-variance
weighting converges on near-equal weighting and changes nothing. Normalizing is a *units* choice that
also happens to reweight observers; it was never stated as either.

## What follows from the decision

- Methods drop the "beta weights for each vertex were standardized" statement and the σ-unit in-text
  statistics; figures use the raw variants, which already exist.
- In the polar experiment the largest asymmetry is horizontal−vertical, the raw result. (This was
  the headline difference in Figure 7, which is itself being removed — [`LME.md`](LME.md).)
- The pink-noise-baseline caveat applies to raw and z-scored analyses equally, so it is not a mark
  against the raw one. Draft manuscript text: [`local_qc/manuscript_caveat_paragraph.md`](local_qc/manuscript_caveat_paragraph.md).

## The one conclusion that still turns on this

The **radial/tangential context effect** reverses under z-scoring: in the z-scored variant the
difference between experiments reaches significance, in the raw variant it does not
([`HARMONIC_MODEL.md`](HARMONIC_MODEL.md) Result 3). The raw analysis is the one adopted, for the
reasons above, which are independent of that result — but the dependence should be stated, not
buried.

## If you need the mechanism

Why z-scoring reverses the radTan/H−V rank order in Fig 7B — a between-observer reweighting, and the
algebra of the reversal — is in [`_archive/ZSCORE_FIG7.md`](_archive/ZSCORE_FIG7.md) §1–§5. Its
§7–§8 are superseded and its §6 precision table was corrected. Diagnostics that still run:
`cleanroom/diagnose_zscore_fig7.m`, `compare_subject_weighting.m`,
`diagnose_exclusion_x_normalization.m`; `run_harmonic_model(true)` fits the z-scored variant.
