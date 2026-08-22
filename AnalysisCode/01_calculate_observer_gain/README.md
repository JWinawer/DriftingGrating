# pRF gain — handoff to whoever is at the server

## The problem

The `ret_sub-*.mat` files copied to New York contain only the prfvista summary
maps: `x, y, eccen, angle, angle_adj, sigma, vexpl`. They do **not** contain the
fitted pRF amplitude (`model.beta`) or the stimulus, and there is no way to
recover the amplitude from what is in them.

Gain needs both. So it has to be computed at the server, from

```
<derivatives>/prfvista_mov/<subject>/<ses-XX>/*_task-prf_acq-normal_run-99_results.mat
```

which is the same file `vistaPRF2MAP.m` reads to write the `.mgz` maps. Its
vertex order is left hemisphere then right, identical to `ret_<subject>.mat`, so
the gain vector drops straight into the existing analysis.

## What gain means here

vistasoft's `rfGaussian2d` makes Gaussians of **unit height**, not unit area, so
the fitted `beta` shrinks as roughly 1/sigma² and is not comparable across
vertices with different pRF sizes. Instead we rebuild the model's predicted time
series — already in percent BOLD, already HRF-convolved — and take its largest
excursion away from the DC baseline:

```
gain = max_t | beta(1) * (allstimimages * RF)(t) |
```

That is size-independent for pRFs inside the stimulus aperture (verified: over a
16x range of `beta`, gain is constant to 6 decimal places).

## Requirements

**None beyond base MATLAB.** No vistasoft, no MathWorks toolboxes. Copy this
folder, `addpath` it, and go. Verified by stripping the path down to base MATLAB
and running the whole pipeline: `requiredFilesAndProducts` reports the three
files in this folder and the single product "MATLAB".

## Do this

### Step 1 — inspect (30 seconds, do this first)

```matlab
addpath('/path/to/this/handoff/folder');
diary ~/dg_inspect.txt
dg_inspectPrfResults('root', '/path/to/derivatives/prfvista_mov');
diary off
```

Send back `~/dg_inspect.txt`. It is a few KB. It reads variable headers only.

This settles the one open question: whether the results files actually carry
`model.beta` and `params.analysis.allstimimages`. Nobody could check from New
York because the volume is not mounted there. **Do not skip this** — step 2 is
pointless if those fields are absent, and the inspect output tells us what to do
instead.

### Step 2 — compute and send back

```matlab
addpath('/path/to/this/handoff/folder');
dg_computeGain('root', '/path/to/derivatives/prfvista_mov', ...
               'outdir', '~/dg_gain');
```

This writes one `gain_<subject>.mat` per subject, each about 1 MB — roughly 8 MB
in total, versus tens of GB for the source files. Send back the whole `~/dg_gain`
folder.

Each file contains:

- `gain` — one value per vertex, percent BOLD, whole surface, lh then rh
- `info` — source file, metric, date, vertex count, units

## Files here

| file | what it is |
|---|---|
| `dg_inspectPrfResults.m` | step 1, read-only inspection |
| `dg_computeGain.m` | step 2, computes and saves the gain vectors |
| `rmModelGain.m` | the gain calculation itself; also lives in vistasoft at `mrBOLD/Analysis/retinotopyModel/` |

`rmModelGain.m` is self-contained. It reads the model fields directly and builds
its own Gaussians, via local `modelGet` and `gauss2d` subfunctions that reproduce
vistasoft's `rmGet` and `rfGaussian2d` exactly for the parameters used here —
checked against the real vistasoft functions in `test_rmModelGain`, agreeing to
5e-15. It handles both the linear and the CSS/exponent pRF models, chunks over
vertices to bound memory, and returns 0 for unfit vertices.

## Caveats worth knowing

- **Gain is stimulus-dependent.** It is the peak predicted BOLD *given the
  stimulus that was run*. For a sweeping bar, a large pRF is never fully
  stimulated, so gain still falls with pRF size — that is real, not an artifact,
  but it is not a stimulus-free property of the vertex.
- If the stimulus never leaves the pRF unstimulated, there is no time point where
  the prediction returns to baseline. Use `'metric', 'range'` (max minus min) in
  that case.
- `vexpl` in these files is a **fraction**, not a percent (max ~0.86), and is NaN
  at unfit vertices.
