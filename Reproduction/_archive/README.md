# `_archive` — retired working documents

**Nothing here is current.** Two documents are kept, and only because current code and documents cite
them for content that was never restated elsewhere. Do not act on anything here without checking the
superseding document first. Start from [`../../AGENTS.md`](../../AGENTS.md).

| document | what it was | why it is still here |
|---|---|---|
| `FINDINGS.md` | reported a polar-angle-ordering bug in `compute_derivativeDirections.m` | **Fully retracted.** The original code is correct; the discrepancy came from two bugs inside the reproduction itself. [`../STIMULUS_CONVENTIONS.md`](../STIMULUS_CONVENTIONS.md) is the correct account. Kept because `test_harmonic_model.m`, `validate_against_manuscript.m` and `run_pathB_figures.m` refer back to it as the thing they guard against. |
| `ZSCORE_FIG7.md` | why z-scoring reverses the radial/tangential vs. horizontal/vertical rank order | **Closed** — z-scoring was abandoned 2026-07-24 ([`../METHOD_DECISIONS.md`](../METHOD_DECISIONS.md) §1). §1–§5 remain the only account of the *mechanism* of that reversal, which is why the file survives; §7–§8 are superseded, and its §6 precision table was corrected in 2026-08 (the original within-observer SEs came from resampling vertices, which is invalid). |

Other working documents from this phase — `NEXT_STEPS.md`, `GLM_QUALITY.md`, and the reproduction
documents superseded by the 2026-08-25 consolidation (`EXTRASTRIATE.md`, `SPEC_FIGURES.md`,
`HARMONIC_MODEL.md`, `LME.md`, `WHY_NOT_ZSCORE.md`, `AUDIT.md`, `local_qc/REPORT.md`,
`local_qc/GLM_SUMMARY_SECTION.md`) — were deleted rather than kept. Every claim they carried is
stated in current form in the document that replaced it, and git history is the record.

## Two things these two documents say that are wrong repo-wide

- **They call the work "published."** It is not. This is a manuscript in preparation — no journal
  version, no preprint. Read every "published" here as "in the draft at the time."
- **They cite `GLM_QUALITY.md` for the within-observer error summary.** That material is now in
  [`../local_qc/RELIABILITY.md`](../local_qc/RELIABILITY.md). Cite that.
