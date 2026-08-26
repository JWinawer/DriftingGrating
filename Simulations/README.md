# `Simulations` — optic-flow illustrations

Two standalone MATLAB scripts that estimate local motion vectors from a movie and render them as a
quiver field over the frames, writing an `.avi`. They exist to **illustrate** that natural viewing
produces radially structured optic flow — self-motion, object motion and eye movements each giving a
different field — as motivation for the radial/tangential asymmetries.

**Nothing here feeds the analyses.** No number in
[`../Reproduction/RESULTS.md`](../Reproduction/RESULTS.md) comes from these scripts, and no
document depends on them. They are figure/illustration material only.

| file | what it does |
|---|---|
| `movie_motionVectors.m` | block-matched motion vectors over a central ROI, averaged over a sliding window, drawn as a quiver field |
| `movie_motionVectors_alternate.m` | the same with a different ROI/parameter choice |

**To run one you must edit it first.** Each hard-codes an absolute `videoFile` path on a machine
that is not this one (`/Users/rje257/...`), and the source movies are not in the repo. Set
`videoFile` to a movie you have; everything else is self-contained base MATLAB plus
`VideoReader`/`VideoWriter`.
