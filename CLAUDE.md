# CLAUDE.md

**Read [`AGENTS.md`](AGENTS.md) first.** It is the entry point for this repository: what the target
analyses are, the code map, the current state of every open and closed question, and an annotated
index of every other document.

This file is a stub so that agents looking for `CLAUDE.md` by convention find their way there.
`AGENTS.md` holds the content and is referenced by name throughout the repo and its scripts — keep
it as the single source of truth rather than duplicating any of it here.

Two conventions worth knowing before editing anything:

- **`_archive/` folders hold retired material.** `Reproduction/_archive/` is not current; check the
  superseding document named in its `README.md` before acting on anything there.
- **Corrections replace the superseded text — they are not stacked on top of it.** When a claim is
  withdrawn or re-measured, rewrite the passage to state the current position, and leave a short
  dated note saying what changed. **Delete the superseded prose, tables and numbers rather than
  keeping them below the note**; git history is the record, and stacking them makes the documents
  progressively harder to read. When a whole document is superseded, move it to `_archive/` and say
  in that folder's `README.md` what replaced it.

  Two things this does *not* mean. Do not rewrite silently — the dated note stays, so a reader can
  see that a number moved and go find the old one. And do not delete a **warning**: "this statistic
  conflates two situations", "this comparison is not a valid test", "generating bootstrap draws as
  a matrix silently changes every SE" are current guidance about what not to do, not history, even
  though they are phrased around a past mistake. Keep those; drop the superseded version of the
  claim itself.

  (Changed 2026-08-24. The earlier convention kept superseded versions in place below a dated
  banner, which is why some documents carry layered revision blocks.)
