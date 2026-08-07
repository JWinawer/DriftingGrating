# AGENTS.md

Orientation notes for anyone (human or AI agent) picking up this repository. This
file summarizes what the project is, how the pieces fit together, and the
non-obvious gotchas that have already cost real debugging time.

## What this project is

This repo supports a manuscript in progress:

> **"Local orientation asymmetries in V1 depend on global stimulus properties"**
> Ezzo, Carrasco, Rokers, Winawer — draft at `Support/draft.pdf`

The study measures fMRI (V1) BOLD responses to two families of grating stimuli:
- **"dg" = drifting grating** → the **Cartesian** experiment (horizontal, vertical,
  and two diagonal gratings)
- **"da" = drifting annulus** → the **Polar** experiment (pinwheel, annulus,
  clockwise/counterclockwise spiral gratings)

These two protocol names (`dg`, `da`) are used throughout the MATLAB code and
folder/derivative naming — whenever you see `dg`/`da` as a `projectName`
variable, that's Cartesian vs. Polar.

The central finding: four orientation asymmetries (vertical>horizontal,
oblique>cardinal, radial>tangential, polar-cardinal>polar-oblique) are strongest
when the asymmetry's reference frame matches the global stimulus's reference
frame (Cartesian asymmetries largest for Cartesian gratings, polar asymmetries
largest for polar gratings), even though local orientation/spatial-frequency
content is matched across the two stimulus types at 6° eccentricity.

Read `Support/draft.pdf` for full scientific context before making analysis
decisions — it explains *why* the pipeline is structured the way it is (e.g.
why there's a subject-wise gain correction, why V1 is binned into 8 polar-angle
sectors, why four specific asymmetries were chosen).

## Repository layout

```
AnalysisCode/     MATLAB pipeline: GLM fitting -> ROI/polar-angle asymmetry
                  stats -> linear mixed-effects models -> plotting.
Models/           Python/Jupyter: image-computable normalization models
                  (steerable-pyramid energy + divisive normalization) that
                  attempt to explain the *sign inversion* of Cartesian
                  asymmetries. Implements the paper's Methods subsection
                  "Normalization as a description of contextual suppression".
ExperimentCode/   Psychtoolbox stimulus presentation code (not analysis).
Simulations/      Standalone motion-vector demo scripts.
Reproduction/     Scripts for pulling data from a remote server.
Support/          draft.pdf (the manuscript) and a summary CSV table.
```

### AnalysisCode/ subfolder pipeline order

Numeric prefixes indicate pipeline stage, but they are not strictly run
top-to-bottom for every analysis — the two chains that matter most this year:

1. **`01_process_singlesubjectGLM/`** — per-subject GLMsingle fit
   (`main_singlesub.m`), produces per-subject `results.mat` contrast files.
2. **`01_calculate_observer_gain/`** — computes each observer's mean pRF gain
   from the retinotopy protocols (`dg_computeGain.m`, `dg_compareGainROI.m`),
   used later to correct for cross-subject differences in raw BOLD scale (see
   *Observer gain correction* below). Output: `gainSummary.mat`/`.csv` under
   `derivatives/summaryTables/`.
3. **`03_process_groupBetas/meanWithinLabel.m`** — reads every subject's
   per-ROI FreeSurfer labels + GLM betas, bins by polar angle, and writes the
   group-level matrices `meanBOLD(pa)`/`medianBOLD(pa)` that everything
   downstream consumes. **This is the one script that actually defines the ROI
   dimension's order** (see ROI ordering gotcha below).
4. **`04_plot_betaAsymmetries/`** — the main analysis + figures:
   - `plot_NeuralAsymmetries.m` → `plot1_experimentalCond.m` (polar plots) /
     `plot2_experimentalCond.m` (pairwise, equal-PA-weighted plots): model-free
     per-asymmetry summaries.
   - `lme1_fit.m`: the linear mixed-effects model (all 4 asymmetries fit
     jointly per ROI), plus a subject-bootstrap for CIs, plus three sets of
     figures (per-asymmetry-across-ROIs, the "master" ROI-1 summary, and a
     per-ROI version of the master figure).
   - `lme2_ploteachDirLoc.m`: polar plots of the fitted model vs. raw data,
     per absolute direction/location — reads `lme1_fit.m`'s saved
     `modeldata.mat`/`LME_bold.mat`, has no logic of its own.

Data (`meanBOLD*.mat`, GLM results, FreeSurfer labels, etc.) lives **outside
this git repo**, on a network volume at
`/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/`. The repo only holds
code and small support tables.

## Known gotchas (read before editing the MATLAB pipeline)

### ROI ordering is positional, not declarative
`AnalysisCode/general/jsons/ROIS_ALL.json` defines the ROI list (`filename`,
`plotname`, `index`). **The `index` field is only meaningful if
`meanWithinLabel.m` was run with that exact ROI list, in that exact order** —
it does not derive from an existing canonical mapping. `meanWithinLabel.m`
writes `meanBOLD(ci, ri, si)` etc. using the ROI list's *position* (`ri`), so
the Nth entry in the JSON becomes the Nth column of every saved
`meanBOLD*`/`medianBOLD*` matrix. Downstream, `plot1_experimentalCond.m` /
`plot2_experimentalCond.m` correctly re-derive the right column via
`projectSettings.roi_idx{ri}`; `lme1_fit.m` was fixed to do the same (previously
it used the raw loop counter, which only ever worked by coincidence). **If you
add/remove/reorder ROIs in `ROIS_ALL.json`, you must re-run
`meanWithinLabel.m`** for every affected project — there is no way to
reconcile a changed ROI list against already-saved matrices. The current
`ROIS_ALL.json` (8 ROIs: V1, V2, V3, V3A, V3B, hV4, pMT, pMST) was restored
from a truncated version (see `ROIS_ALL_2026-08-05.json` backup in the same
folder) — if you need the original 9-ROI canonical order (which included
`hMTcomplex`), it's recoverable from git history (`ROIS_ALL copy.json`,
commit `e37e1be` and earlier).

### Observer gain correction
Each subject's BOLD values get multiplied by `subjectScale(i) = groupGain /
gain_i`, where `gain_i` is that subject's mean pRF gain (from
`retrieveObserverGainWeights.m` / `gainSummary.mat`) and `groupGain =
mean(gainWeights)` (arithmetic mean of the whole subject set). This:
- **down-weights high-gain subjects, up-weights low-gain subjects** — the
  substantive correction.
- multiplying back by `groupGain` afterward is a **pure unit-restoration
  step** — mathematically inert for every inferential statistic (t/F-stats,
  p-values); it only rescales the reported effect sizes back toward original
  BOLD-like units. Swapping arithmetic mean for geometric mean (or any other
  positive constant) would rescale every number by the same tiny constant
  factor and change nothing about significance.
- For `dg` (13 subjects total), only the **first 8** are used in the main
  analysis scripts (`lme1_fit.m`, `plot_NeuralAsymmetries.m`) — hardcoded
  index `[1:8]` — to match the 8 subjects who completed the `da` experiment.
  `lme2_ploteachDirLoc.m` has no subject-selection logic of its own; it
  inherits whatever subset `lme1_fit.m` used when it saved `modeldata.mat`.

### MATLAB working-directory requirement
`AnalysisCode/general/setup_user.m` requires MATLAB's current folder to be
`AnalysisCode/` itself (or a path containing it) — it explicitly checks
`pwd` and errors otherwise (`"AnalysisCode folder containing setup.json must
be the parent folder."`). Scripts that call `setup_user()`
(`meanWithinLabel.m`, GLM-fitting scripts) must be run with that working
directory. Scripts that don't call `setup_user()` (e.g. `lme1_fit.m`) still
rely on `addpath(genpath(pwd))`, so make sure `AnalysisCode/` (and all
subfolders) are on the MATLAB path before running anything, regardless of
which tool/method is used to launch the script.

### Figure-export loops need `drawnow` before `print()`
`lme1_fit.m` previously had a real bug where two figures generated back-to-back
in the same loop (`mainSubset`/`derivedSubset` asymmetries) were saved with
*identical* content — colors, legend, and data all from the wrong figure —
because `print()` grabbed a stale render before MATLAB finished drawing the
legend/boxcharts. Fixed by calling `drawnow;` immediately before every
`print(...)` inside a loop that creates multiple figures. Any new
figure-generation loop added to this pipeline should do the same.

## Models/ (Python) — normalization simulations

`model_DriftingGratings.ipynb` implements the paper's Methods subsection
*"Normalization as a description of contextual suppression"* and its
(currently incomplete) Results counterpart *"Normalization model describes
context-dependent suppression"* (draft.pdf, Methods p.13, Results p.25 — the
results section still has a placeholder figure, "Fig 9. Image computable
model? ... schematic?").

Pipeline: build Cartesian/polar stimulus images → run through a Steerable
Pyramid (4 orientations × 6 spatial frequencies, via the `plenoptic` library)
→ artificially double the cardinal-orientation channels' energy (simulating
V1's overrepresentation of cardinal-tuned neurons) → test whether various
divisive-normalization formulations can turn that overrepresentation into a
*suppressed* net response (matching the "inverted" BOLD asymmetry seen for
cardinal-vs-oblique in the Cartesian experiment, but not the polar one).

Models implemented in `helper_functions/utils_image.py`:
1. `div_normalization(..., tuned=False)` — untuned divisive normalization.
2. `div_normalization(..., tuned=True)` — feature/orientation-tuned normalization.
3. `div_normalization(..., tuned=False, q_exp=2)` — untuned + exponentiated
   (superlinear) suppression.
4. `normalization_byAnisotropy(...)` — suppression scales with the standard
   deviation of energy across orientation channels (adapted from Fang et al.).
5. `normalization_byStimHomogeneity(...)` — **not yet in the manuscript text**;
   a newer, exploratory model based on cosine similarity between center/surround
   feature vectors. Treat as unpublished/in-progress.

Per the draft: models 1–2 fail to reproduce the suppression; models 3–4
succeed, and only for Cartesian (not polar) stimuli — mirroring the fMRI data.
Beyond the manuscript's 4 named models, the notebook's own internal numbering
(see below) also includes an anisotropy variant that skips spatial collapse
and the not-yet-published `normalization_byStimHomogeneity` model.

**The notebook's title cell (cell 0) is a running DONE / open-questions / TO-DO
log kept by the author — read it first, it is more current than any static
summary of the code.** Highlights as of this writing: radial-orientation
overrepresentation has also been added (in addition to cardinal) with
behaviorally-derived weights; there's an open question about how to choose the
4D Gaussian convolution bandwidth for the untuned/tuned models; and the TO-DO
list includes adding a polar-cardinal asymmetry version, a literature review of
suppressive-exponent normalization, and testing predictions for contrast/
adaptation/natural-scene manipulations.

The notebook otherwise has no markdown documentation of the code itself and
ends in empty cells — it is a working draft, not a finished analysis. See the
installation instructions at the top of the notebook and in `README.md`
before running it.

### Performance: the normalization models are GPU-accelerated, the pyramid step is not

Profiling (2026-08-07) found the notebook's ~12-minute runtime was ~99% spent
in `div_normalization` (Models 1–3) and `normalization_byStimHomogeneity`
(Model 6) — specifically in their `scipy.ndimage.gaussian_filter` calls with
a ~90px-radius kernel — not in the steerable-pyramid step, which is <1% of
total time. `utils_image.py` now has `gaussian_filter_gpu()` (backed by
`get_fast_device()`: CUDA > MPS > CPU), a drop-in torch-based replacement for
`scipy.ndimage.gaussian_filter` used inside those two functions, verified to
match scipy's output to ~1e-7 relative error. This cut total runtime from
~715s to ~44s on this Intel Mac (MPS backend) — a ~16x speedup — with no
change to model behavior. `div_normalization` and
`normalization_byStimHomogeneity` both take an optional `device=` kwarg
(defaults to auto-detect).

Implementation note: `gaussian_filter_gpu` cannot just move the SciPy call's
axis-length-based reasoning to torch — SciPy's per-axis cost scales with
*kernel radius*, not axis length, so even the tiny 4–6-element SF/orientation
axes are expensive to filter at this sigma. It instead builds a dense `(L,
L)` linear operator per axis (matching scipy's default half-sample-symmetric
`mode='reflect'` boundary behavior, including the periodic folding needed
when the kernel radius exceeds the axis length) and applies it via matmul —
padding-based conv1d was tried first and blew MPS's buffer limit (~9GB) when
padding a small axis while other axes stayed huge.

The steerable pyramid step is intentionally left on CPU regardless of
`DEVICE`: `plenoptic`'s `SteerablePyramidFreq` calls `torch.fft.fft2` on
complex tensors internally, and PyTorch's MPS backend does not support
complex dtypes for FFT ops as of `torch==2.2.2` (confirmed by direct test:
`torch.fft.fft2` on a real MPS tensor raises "Trying to convert ComplexFloat
to the MPS backend but it does not have support for that dtype"). This may
be fixed in newer torch releases (worth re-checking on Apple Silicon setups
pinned to `torch==2.9.0`, e.g. via `Models/requirements.txt`), but isn't
worth chasing further given the step is <1% of runtime either way.

### Dependency pins are platform-specific — do not `pip install` unpinned

`Models/requirements.txt` (Apple Silicon / arm64) and
`Models/requirements-intel-mac.txt` (Intel / x86_64) exist because a plain
`pip install plenoptic torch ...` silently grabs versions that break the
notebook, in two independent ways that were both hit and diagnosed in this
repo's history:

1. **`plenoptic>=2.0` removed the `plenoptic.simulate` module** that this
   notebook imports (`from plenoptic.simulate import SteerablePyramidFreq`) —
   this fails with `ModuleNotFoundError: No module named 'plenoptic.simulate'`,
   not an obvious version-mismatch error. Both requirements files pin
   `plenoptic==1.3.1`, matching a real working environment (confirmed via
   `pip show`/`pkgutil.iter_modules` that `2.1.0` genuinely lacks this
   submodule; `1.3.1` has it).
2. **PyTorch dropped Intel-macOS (x86_64) wheels after the `2.2.x` series.**
   On an Intel Mac, no `torch` version at all is installable for Python
   ≥3.12, and the newest installable version for Python 3.11 is `2.2.2`.
   `torch==2.2.2` was compiled against NumPy 1.x's ABI, so installing NumPy
   2.x alongside it produces `UserWarning: Failed to initialize NumPy:
   _ARRAY_API not found` and real crash risk (NumPy's own message: "may
   crash"). This further constrains Intel Macs to `numpy<2` (verified pin:
   `1.26.4`) and, transitively, to an older `opencv-python` release
   (`4.9.0.80`) since `opencv-python>=4.10ish` requires `numpy>=2`, which
   would directly conflict with the `numpy<2` pin. None of this applies on
   Apple Silicon, which has current wheels for everything.

If you see either failure mode again (missing `plenoptic.simulate`, or the
NumPy ABI warning cascading into an unrelated `ModuleNotFoundError`), the fix
is almost always "an unpinned `pip install` grabbed a too-new package" — check
`pip show <package>` for the installed version against the pinned
`requirements*.txt`, don't just reinstall the same unpinned command again.

When updating these pins (e.g. to allow a newer `plenoptic` once its API
stabilizes), verify on **both** architectures before merging — `python3 -c
"import platform; print(platform.machine())"` distinguishes them
(`arm64` vs `x86_64`), and a clean import check looks like:
```python
import torch, numpy
from plenoptic.simulate import SteerablePyramidFreq  # must not raise
```
with no NumPy ABI warning printed.

## When picking up new work here

- Check `Support/draft.pdf` for the current state of the manuscript text
  before assuming what's "done" — sections can be explicitly marked
  incomplete (placeholder figures, missing quantitative comparisons).
- Anything under `AnalysisCode/04_plot_betaAsymmetries/` that reads
  `meanBOLD*`/`medianBOLD*` assumes those matrices already exist for the
  current ROI list; regenerate via `meanWithinLabel.m` first if in doubt.
- Prefer extending the existing per-project (`dg`/`da`) branching pattern
  (`if strcmp(projectName, 'dg') ... elseif strcmp(projectName, 'da')`) rather
  than introducing a new convention — it's used consistently across every
  script in `04_plot_betaAsymmetries/`.
