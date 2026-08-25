# `_archive` — retired reproduction documents

Working documents from earlier phases, kept as a record. **Nothing here is current.** Each was
superseded; the pointers below say by what. Do not act on anything in this folder without checking
the superseding document first. Start from [`../../AGENTS.md`](../../AGENTS.md).

| document | what it was | status |
|---|---|---|
| `FINDINGS.md` | reported a polar-angle-ordering bug in `compute_derivativeDirections.m` | **fully retracted.** The original code is correct; the discrepancy came from two bugs in the reproduction itself. [`../AUDIT.md`](../AUDIT.md) is the correct account. Kept because `AUDIT.md`, `test_harmonic_model.m` and `validate_against_manuscript.m` refer back to it. |
| `NEXT_STEPS.md` | the task setups behind the z-scoring and GLM-quality audits | **closed.** z-scoring decided in [`../WHY_NOT_ZSCORE.md`](../WHY_NOT_ZSCORE.md); GLM quality answered in [`../local_qc/REPORT.md`](../local_qc/REPORT.md). Its one live item — a GLM-`R²` column in `allsubjectsTable.csv` — is carried in `AGENTS.md` §7. |
| `ZSCORE_FIG7.md` | why z-scoring reverses the radTan/H−V rank order in Fig 7B | **closed.** Z-scoring was abandoned 2026-07-24. §1–§5 remain the account of the *mechanism*, which [`../HARMONIC_MODEL.md`](../HARMONIC_MODEL.md) still cites for the one conclusion that turns on the choice; §7–§8 were already superseded. Its §6 precision table was corrected 2026-08-18 — the original within-observer SEs came from resampling vertices, which is invalid. |
| `GLM_QUALITY.md` | the **first** GLM fit-quality audit, on the 4–8°-filtered extraction | **superseded** by [`../local_qc/REPORT.md`](../local_qc/REPORT.md), the second pass on the unfiltered re-extraction this document's §6a called for, whose §2.1 table is strictly richer. Unique to it: the `fracAbove` scale argument (§3) and the median-understates-sessions point (§4a), both about the first extraction. Its §7 measurement-reliability section was **moved** to `REPORT.md` §2.7 — cite that, not this. |

Retired 2026-08-18, per the leading-underscore convention for one's own set-aside material.

---

## The correction convention changed, 2026-08-24

Superseded text is now **deleted or moved here**, not kept in place below a dated banner. A short
dated note stays where the claim was, saying what changed; git history is the record. Warnings
about what *not* to do are current guidance and are kept. See [`../../CLAUDE.md`](../../CLAUDE.md).

The documents in this folder predate that change and still carry layered revision blocks. That is
fine — they are archived, and nothing here is current.

What moved under the new convention, all recoverable from git history rather than copied here,
because every claim they carried is stated in its current form in the document that replaced them:

| removed from | what it was |
|---|---|
| `LME.md` §5 | the median-route, un-rescaled precision-weighting tables, superseded 2026-08-19 by the mean + gain route |
| `HARMONIC_MODEL.md` | the revision banners 2026-08-17b/c, and the withdrawn recommendation to read the `ref` (median) column |
| `EXTRASTRIATE.md` §9 | three withdrawn claims from earlier passes; the section now keeps only the code traps that could recur |
| `local_qc/REPORT.md` §2.1 | the quoted superseded table caption |
| `local_qc/GLM_SUMMARY_SECTION.md` | the "merge the 4 upstream commits" item (those commits are merged) and the Fig-4C correction's quoted earlier wording |
| `server_extract/RUNME.md` | folded whole into that folder's `README.md` |

The same pass made `AGENTS.md` the single entry point and the single home for open items, and
collapsed material duplicated across `AGENTS.md`, `Reproduction/README.md` and individual document
banners to one statement each.

## Two things these documents say that are wrong repo-wide

- **They call the work "published."** It is not. This is a manuscript in preparation — no journal
  version, no preprint. Across the current documents "published route / values / asymmetries" was
  replaced with "manuscript route / values / asymmetries" on 2026-08-24; the files in this folder
  were left as they were, since nothing here is current. Read every "published" in this folder as
  "in the draft at the time."
- **`GLM_QUALITY.md` is cited elsewhere for the within-observer error summary.** That section moved
  to `../local_qc/REPORT.md` §2.7. Cite that.
