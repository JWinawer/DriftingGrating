# `_archive` — retired reproduction documents

Working documents from earlier phases of the reproduction, kept as a record. **Nothing here is
current.** Each was superseded by a later document; the pointers below say by what. Do not act on
anything in this folder without checking the superseding document first.

| document | what it was | status |
|---|---|---|
| `FINDINGS.md` | reported a polar-angle-ordering bug in `compute_derivativeDirections.m` | **fully retracted.** The original code is correct; the discrepancy came from two bugs in the reproduction itself. `../../AUDIT.md` is the correct account. Kept only because `AUDIT.md`, `test_harmonic_model.m` and `validate_against_manuscript.m` refer back to it. |
| `NEXT_STEPS.md` | the task setups behind the z-scoring and GLM-quality audits | **closed.** z-scoring decided in `../local_qc/REPORT.md` §4 (use raw); GLM quality answered there and in `../GLM_QUALITY.md`. Its one live item — adding a GLM-`R²` column to `allsubjectsTable.csv` — is carried in `../README.md`'s open list. |
| `ZSCORE_FIG7.md` | why z-scoring reverses the radTan/H−V rank order in Fig 7B | **closed.** Z-scoring was abandoned 2026-07-24 (`../local_qc/REPORT.md` §4). §1–§5 remain the account of the *mechanism*, which `../HARMONIC_MODEL.md` still cites for the one conclusion that turns on the choice; §7–§8 were already superseded. Its §6 precision table was corrected 2026-08-18 — the original within-observer SEs came from resampling vertices, which is invalid. |

Retired 2026-08-18, per the leading-underscore convention for one's own set-aside material.
