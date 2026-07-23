# GLMsingle fit quality, all 8 observers × both experiments

Results of the audit set up in [`NEXT_STEPS.md`](NEXT_STEPS.md) ("per-subject GLM fit quality").
Data extracted on the server 2026-07-23 with
[`server_extract/extract_for_transfer.m`](server_extract/extract_for_transfer.m), returned in
`Support/glm_qc_for_transfer/`, audited with `cleanroom/audit_glm_quality.m`.

```matlab
cd Reproduction/cleanroom
A = audit_glm_quality('../../Support/glm_qc_for_transfer');
```

## 1. The delivery

16/16 files, every row `ok` in `summary.csv` — no `missing`, no `bad-structure`, so every
observer was run under `hRF_setting = 'glmsingle'` and has the TYPED output. The **V1-restriction
branch worked** (it was untested when written, since no labels or pRF maps were on the machine
where it was authored): every observer has a real patch, 1126–1786 vertices, and `fieldsMissing`
is empty throughout. Files are ~50 kB rather than the ~10 MB fallback.

One incompatibility, since fixed: `extract_glm_qc` writes the field `scope`, the standalone
server script writes the same information as `patchNote`, and `audit_glm_quality` read only the
former. It now accepts either. The claim at NEXT_STEPS.md that the audit "reads its output
directly" was not true as written.

## 2. Median GLM R² in the V1 patch

| subject | dg (cartesian) | da (polar) | da − dg |
|---|---|---|---|
| sub-0037 | **16.99** | **2.83** | **−14.16** |
| sub-0201 | 3.46 | 4.50 | +1.04 |
| sub-0250 | 6.33 | 3.59 | −2.74 |
| sub-0255 | 3.57 | 4.32 | +0.75 |
| sub-0395 | 27.43 | 20.56 | −6.87 |
| sub-0426 | 5.02 | 10.99 | +5.97 |
| sub-wlsubj123 | 4.94 | 4.46 | −0.48 |
| sub-wlsubj124 | 4.72 | 3.43 | −1.30 |
| *median* | *4.98* | *4.39* | |

Run counts are not uniform: sub-0395 has 6 polar runs, sub-0255 has 9 cartesian, everyone else
has 8.

## 3. Why the raw table is the wrong scale — and the right one

On the raw numbers sub-0037's polar R² of 2.83 is the lowest but unremarkable: sub-wlsubj124 sits
at 3.43 and the group median is 4.39. An IQR test flags **no one** in the bad direction; the only
flagged observer is sub-0395, flagged for being unusually *good*.

R² is not comparable across sessions on an absolute scale — it depends on run count, motion, coil
and SNR. But each file also carries the median R² over the **whole cortical surface**, most of
which is non-visual and should not respond to gratings at all. That is a per-session noise floor.
Dividing the V1 patch by it asks the question that matters: **did V1 fit better than average
cortex in this session?**

| subject | dg (cartesian) | da (polar) |
|---|---|---|
| sub-0037 | 4.96 | **0.92** |
| sub-0201 | **1.15** | **1.16** |
| sub-0250 | 1.58 | 1.11 |
| sub-0255 | 1.08 | 1.29 |
| sub-0395 | 9.10 | 5.24 |
| sub-0426 | 1.55 | 2.56 |
| sub-wlsubj123 | 1.47 | 1.27 |
| sub-wlsubj124 | 1.73 | 1.19 |
| *median* | *1.56* | *1.23* |

1.0 means V1 fit no better than average cortex, i.e. no measurable response to the stimulus.

Note how low this is across the board. Only sub-0395 (both experiments), sub-0037 (`dg`) and
sub-0426 (`da`) are much above ~1.5 — for most observer-sessions, V1 is fit only slightly better
than non-visual cortex. That is a general weakness of this dataset, not a two-observer problem.

## 4. The two anomalous observers fail differently

**sub-0037 — confirmed session-specific failure.** Its polar ratio of 0.92 is the only value below
1.0 anywhere in the dataset: V1 was fit *worse* than the brain-wide baseline. Its cartesian
session is 4.96, second-best in the set. This is exactly the prediction from ZSCORE_FIG7.md §3a
(strong in `dg`, no differentiation in `da`), and it is independent of the beta analysis that
generated the prediction. GLMsingle's own ridge fraction agrees: `FRACvalue` = 0.25 for the
cartesian session, 0.05 — maximum shrinkage, the algorithm's verdict that there is little signal
to fit — for the polar one.

It is **not** a bad run. All 8 polar runs sit at median R² 2.4–3.8 (against 11.5–19.9 for all 8
cartesian runs). There is nothing to drop; the whole session is uniformly flat.

**sub-0201 — prediction not confirmed.** The expected bad-run/motion signature is absent. Its
`runSpread` is the *smallest* in the entire group (0.27 in `dg`, 0.89 in `da`) — every run is
uniformly mediocre, which is the opposite of a motion artefact. Its ratios are 1.15 and 1.16, the
two weakest non-failures in the set, and they are weak in **both** experiments. In the raw polar
table it is mid-pack, 6th of 8.

So sub-0201 is a quiet responder throughout, not a botched polar session. Its blank beta exceeding
all 12 stimulus betas in both experiments remains unexplained by GLM fit quality.

## 5. What this changes

[`ZSCORE_FIG7.md`](ZSCORE_FIG7.md) §8 excludes sub-0037 and sub-0201 together, on the grounds that
neither has a measurable gain in the polar experiment, and reports that all five normalisers then
agree that H−V is the largest polar asymmetry. **This audit supports half of that.**

- **sub-0037:** independent, non-circular support for exclusion. Its polar session did not produce
  a usable V1 response by a measure that never touches the 13 condition betas.
- **sub-0201:** no support. Excluding it is still defensible on the blank-beta grounds, but that
  is the *same* observation the exclusion is meant to justify — the argument is circular in a way
  it is not for sub-0037. §8 should say so rather than treating the two symmetrically.

Worth checking: does the §8 result survive excluding **sub-0037 only**? If H−V remains the largest
polar asymmetry under all five normalisers with sub-0201 retained, the conclusion holds on
non-circular grounds and the §8 framing can simply be tightened. If it does not, the result depends
on an exclusion that this audit cannot independently justify, and that has to be stated plainly.

## 6. Open: is sub-0037's polar session a processing error?

A processing error — wrong design matrix, misaligned stimulus timing, swapped or mislabelled runs
— produces a poor fit in V1 while the rest of the brain looks normal. That is precisely the
observed pattern, so **this data does not distinguish "the observer did not respond" from "the
session was processed wrong."** What it does rule out is a partial failure: with all 8 runs
uniformly flat there is no bad run to drop.

The discriminating check is to pull that session's design matrix and stimulus timing and compare
against a session that worked. It is worth doing before excluding: sub-0037 has the second-best
cartesian data in the set, so losing the observer costs real power, and if the polar session is
recoverable it should be recovered rather than dropped.

## 7. Also worth noting

`prfvista_mov` saves no gain map — see the retinotopy section of [`NEXT_STEPS.md`](NEXT_STEPS.md),
now closed negative. The inventory is identical for all 8 observers: `angle`, `angle_adj`,
`eccen`, `sigma`, `vexpl`, `x`, `y`.

Step 5 of the original task — adding a GLM-`R2` column to `allsubjectsTable.csv` so the vertex
filter can screen on it alongside `pRF_r2` — remains open and is now better motivated: §3 shows
that for most observer-sessions V1 barely outperforms non-visual cortex, so a vertex can pass
`pRF_r2 > 0.1` on a good retinotopy fit while its 13 condition betas are essentially unconstrained.
