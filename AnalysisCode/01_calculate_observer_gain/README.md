# `01_calculate_observer_gain` — per-observer pRF gain

Computes the per-vertex pRF **gain**, the divisor behind the per-observer rescaling applied to
Figures 5–7. The rescaling divides each observer's BOLD by their own mean gain and multiplies the
group mean back in (`04_plot_betaAsymmetries/retrieveObserverGainWeights.m`).

**Status: done.** Gain was computed at the server for all 8 observers and the results are in
`gainSummary.csv` (at `~/dg_collect/`). The V1-only restriction in `dg_computeGain` has since been
widened — `Reproduction/server_extract/collect_gain_areas.m` computes gain over all eight visual
maps, so a per-observer × **map** factor is available, which is what the settled specification uses
(`Reproduction/SPECIFICATION.md` §4). That path reproduces this one in V1 4–8° to **5.3e-15**, which
is the argument for keeping both.

One difference worth knowing: the mov/stat protocols are combined by **geometric** mean in the
current specification, superseding `gainSummary.csv`'s arithmetic average.

## Why this module exists at all

The `ret_sub-*.mat` files copied to New York contain only the prfvista summary maps
(`x, y, eccen, angle, angle_adj, sigma, vexpl`). They do **not** contain the fitted pRF amplitude
(`model.beta`) or the stimulus, and there is no way to recover the amplitude from what is in them.
Gain needs both, so it has to be computed at the server, from

```
<derivatives>/prfvista_mov/<subject>/<ses-XX>/*_task-prf_acq-normal_run-99_results.mat
```

which is the same file `vistaPRF2MAP.m` reads to write the `.mgz` maps. Its vertex order is left
hemisphere then right, identical to `ret_<subject>.mat`, so the gain vector drops straight into the
existing analysis.

## What gain means here

vistasoft's `rfGaussian2d` makes Gaussians of **unit height**, not unit area, so the fitted `beta`
shrinks as roughly 1/σ² and is not comparable across vertices with different pRF sizes. Instead we
rebuild the model's predicted time series — already in percent BOLD, already HRF-convolved — and
take its largest excursion away from the DC baseline:

```
gain = max_t | beta(1) * (allstimimages * RF)(t) |
```

That is size-independent for pRFs inside the stimulus aperture (verified: over a 16× range of
`beta`, gain is constant to 6 decimal places).

**Caveats worth knowing.**

- **Gain is stimulus-dependent.** It is the peak predicted BOLD *given the stimulus that was run*.
  For a sweeping bar, a large pRF is never fully stimulated, so gain still falls with pRF size —
  that is real, not an artifact, but it is not a stimulus-free property of the vertex.
- If the stimulus never leaves the pRF unstimulated, there is no time point where the prediction
  returns to baseline. Use `'metric', 'range'` (max minus min) in that case.
- `vexpl` in these files is a **fraction**, not a percent (max ~0.86), and is NaN at unfit vertices.

## Files

| file | what it is |
|---|---|
| `dg_inspectPrfResults.m` | read-only inspection — confirms the results files carry `model.beta` and `params.analysis.allstimimages` |
| `dg_computeGain.m` | computes and saves one `gain_<subject>.mat` per subject (~1 MB each) |
| `rmModelGain.m` | the gain calculation itself; also lives in vistasoft at `mrBOLD/Analysis/retinotopyModel/` |
| `dg_compareGainROI.m` | mean/median gain in a V1 mask, moving (`prfvista_mov`) vs. stationary (`prfvista`) carrier |
| `dg_plotGainScatter.m` | scatters the two carriers against each other |
| `computeMotminusStaticDotsPerSubject.m` | a separate per-subject scalar from the `mt+2` localiser, unrelated to pRF gain |

`rmModelGain.m` is self-contained. It reads the model fields directly and builds its own Gaussians,
via local `modelGet` and `gauss2d` subfunctions that reproduce vistasoft's `rmGet` and
`rfGaussian2d` exactly for the parameters used here — checked against the real vistasoft functions
on a machine that had them, agreeing to 5e-15. (That check, `test_rmModelGain`, is not in this
repo; it needs vistasoft on the path.) It handles both the linear and the CSS/exponent pRF
models, chunks over vertices to bound memory, and returns 0 for unfit vertices.

## Requirements

**None beyond base MATLAB** for the gain path itself — no vistasoft, no MathWorks toolboxes. Copy
this folder, `addpath` it, and go. Verified by stripping the path down to base MATLAB and running
the pipeline: `requiredFilesAndProducts` reports only `dg_inspectPrfResults`, `dg_computeGain` and
`rmModelGain`, and the single product "MATLAB".

## Re-running it at the server

```matlab
addpath('/path/to/this/folder');

% 1. Inspect (30 s, reads variable headers only)
diary ~/dg_inspect.txt
dg_inspectPrfResults('root', '/path/to/derivatives/prfvista_mov');
diary off

% 2. Compute
dg_computeGain('root', '/path/to/derivatives/prfvista_mov', 'outdir', '~/dg_gain');
```

Step 2 writes one `gain_<subject>.mat` per subject — one value per vertex, percent BOLD, whole
surface, lh then rh — plus an `info` struct recording source file, metric, date, vertex count and
units. Roughly 8 MB total, against tens of GB for the source files.
