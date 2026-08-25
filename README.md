# DriftingGrating

See `AGENTS.md` for a fuller summary of the project (the analysis pipeline,
what is settled, what is open, and known gotchas). The manuscript is **in
preparation** — not published, not preprinted; it is a Google Doc reachable
through the shortcut in `Manuscript/`.

## ExperimentCode

DESCRIPTION
ExperimentCode folder contains the scripts to run the experiment.

INSTRUCTIONS BEFORE RUNNING
It is advised to open the scrConfig.m script (inside config folder) and add the screen specs for correct presentation.

INSTRUCTIONS TO RUN
run expLauncher.m from inside the ExperimentCode folder

TESTED ON
MATLAB R2020a

## AnalysisCode

MATLAB pipeline that turns the preprocessed fMRI data into the orientation
asymmetry statistics and figures reported in the manuscript. See
`AnalysisCode/README.rtf` for a per-script description, and `AGENTS.md` for
pipeline order and known gotchas (ROI ordering, observer gain correction,
required working directory).

## Models

`Models/model_DriftingGratings.ipynb` is a **Jupyter notebook** implementing
the image-computable normalization simulations described in the manuscript's
Methods subsection "Normalization as a description of contextual suppression"
(see `AGENTS.md` for details on what it does and its current, unfinished
state).

### Installation

Running this notebook requires Jupyter itself plus several scientific-Python
and deep-learning packages (including
[`plenoptic`](https://github.com/plenoptic-org/plenoptic) for the
steerable-pyramid model). Recommended setup, using a dedicated conda
environment named `viz`.

**First, check which Mac architecture you have** (this matters — see
"Apple Silicon vs. Intel" below):

```bash
python3 -c "import platform; print(platform.machine())"
# arm64  -> Apple Silicon (M1/M2/M3/...)
# x86_64 -> Intel
```

```bash
# 1. Create and activate a new conda environment named "viz"
#    (Apple Silicon: python 3.11-3.13 all work; Intel: use 3.11 -- see below)
conda create -n viz -c conda-forge python=3.11 -y
conda activate viz

# 2. Install Jupyter itself (needed to open/run the notebook)
pip install jupyter ipykernel

# 3. Install the notebook's pinned dependencies -- pick the file matching
#    your Mac's architecture (from the `platform.machine()` check above):
pip install -r Models/requirements.txt            # Apple Silicon (arm64)
pip install -r Models/requirements-intel-mac.txt   # Intel (x86_64)

# 4. Register "viz" as a selectable Jupyter kernel
python -m ipykernel install --user --name viz --display-name "Python (viz)"
```

#### Apple Silicon vs. Intel: why there are two requirements files

`Models/requirements.txt` (Apple Silicon / arm64) and
`Models/requirements-intel-mac.txt` (Intel / x86_64) pin **different**
`numpy`/`torch`/`torchvision`/`opencv-python` versions, verified working on
each platform. This isn't arbitrary — PyTorch stopped publishing macOS
x86_64 wheels after the `2.2.x` series, and `torch==2.2.2` was compiled
against NumPy 1.x's ABI (installing NumPy 2.x alongside it produces a
"compiled using NumPy 1.x" warning and real crash risk, and can cascade into
`ModuleNotFoundError`s in packages that depend on `torch`+`numpy` both
importing cleanly). Intel Macs are therefore capped at Python 3.11 +
`torch==2.2.2` + `numpy<2` + an older `opencv-python` release (`4.12.0.88`+
requires `numpy>=2`, which conflicts). Apple Silicon has no such ceiling and
can use current releases of everything. See `AGENTS.md` for the full story.

Also important: **`plenoptic` must be `<2.0`** (both requirements files pin
`1.3.1`) — `plenoptic>=2.0` removed the `plenoptic.simulate` module this
notebook imports, so `pip install plenoptic` without a version pin will
silently grab a version that's missing `SteerablePyramidFreq` entirely.

### Running the notebook

```bash
conda activate viz
jupyter notebook   # or: jupyter lab
```

Then open `Models/model_DriftingGratings.ipynb`. Its saved kernel metadata
already points at the `viz`/"Python (viz)" kernel (from step 4 above), so it
should launch with the right environment automatically — no need to
manually switch kernels via the Kernel menu, though it's worth glancing at
the kernel name shown in the top-right of the notebook UI to confirm.

You technically don't have to `conda activate viz` before launching
`jupyter` — kernels registered via `ipykernel install --user` are visible to
*any* Jupyter server on the machine, not just the one launched from that
environment. But activating `viz` first is the simplest way to guarantee
`jupyter` itself resolves to a Jupyter installation that can see the kernel,
especially if you have other Jupyter installs (e.g. a base conda
environment) that might otherwise shadow it on your `PATH`.

A `ModuleNotFoundError` after following these steps (e.g. for `ipdb`) almost
always means Jupyter is running a different kernel/environment than `viz` —
check the kernel name shown in the notebook UI, and confirm with
`jupyter kernelspec list` that `viz` points at the right environment
(`.../envs/viz/bin/python`).
