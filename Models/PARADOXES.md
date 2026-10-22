# Three paradoxes, read through the steerable-pyramid model

**What this is.** This uses the models in this folder *qualitatively*. The aim is to show which
model component can produce each puzzling pattern in the V1 data, and which component rules it out.
It is **not** a fit to the fMRI data, and no model parameter was tuned to the data. Written
2026-10-21 against `main` at 547ca09 (Rania's Models/ update). The repo's model code is unchanged.
All new code is in [`paradoxes/`](paradoxes/).

## The three paradoxes

V1, matched 7 observers, raw % signal change (from [`RESULTS.md`](../Reproduction/RESULTS.md)).
`dg` = Cartesian gratings, `da` = polar gratings.

| | dg | da | context (dg − da) |
|---|--:|--:|--:|
| cardinal − oblique | **−0.193** | −0.044 (n.s.) | **−0.149** |
| radial − tangential | **+0.091** | **+0.198** | **−0.107** |

1. **Cardinals are overrepresented in V1, but they give a weaker BOLD response than obliques.**
   Fang, Cai and Lu (2022) found both cardinal and radial overrepresentation in macaque V1
   ([doi:10.1073/pnas.2113407119](https://doi.org/10.1073/pnas.2113407119)).
2. **Radial orientations are also overrepresented, and they give a *stronger* response.**
3. **Both effects depend on context.** The cardinal deficit is large for Cartesian gratings and
   near zero for polar gratings. The radial advantage is about twice as large for polar gratings.

One point about the data matters for paradox 1. The table's two Cartesian numbers can be split
apart. Horizontal − vertical in dg is −0.514. So, relative to the obliques, vertical is about
**+0.06** and horizontal about **−0.45**. (This is arithmetic on the group means:
V − obl = (card−obl) − (H−V)/2.) The "cardinal deficit" is therefore really a **horizontal**
deficit. The models below have one gain shared by both cardinal channels. They cannot express this.

## How the model was read out

The readout is the same as the notebook's. Take the pyramid energy, apply the imbalance and the
normalization, then average over SF and sum over orientation channels. Then average within each
45° polar-angle wedge at 4–8° eccentricity. The four asymmetries are formed from those wedge values
exactly as in the fMRI analysis, with equal weight on every wedge.
[`paradox_lib.py`](paradoxes/paradox_lib.py) re-implements each normalization pool so that its
parameters can be varied. Each re-implementation matches the repo function to better than 1e-5 at
the repo's settings ([`00_checks.py`](paradoxes/00_checks.py)).

Model numbers are given as **% of that model's mean response**. Only their **signs and ordering**
can be compared with the data. Absolute size cannot be compared: the data have no true baseline,
because the "blank" is pink noise (AGENTS.md, standing fact 1).

## Result table

At the repo's parameters: cardinal gain 2, radial gain 1.5, σ = 0.01 (0.1 for M5), 100% contrast.
✓ means the model's sign matches the data.

| model variant | P1: dg card < obl | P2: rad > tang, both sets | P3: card context (dg < da) | P3: rad context (da > dg) | all |
|---|:-:|:-:|:-:|:-:|:-:|
| linear, no normalization | ✗ +40 | ✓ | ✗ | ✗ (equal) | ✗ |
| M1 untuned, p = q = 1 | ✗ +9 | ✓ | ✓ | ✗ | ✗ |
| M2 untuned, q = 2 | ✓ −39 | ✗ da −12 | ✓ | ✗ | ✗ |
| M3 tuned, p = q = 1 | ✗ +55 | ✓ | ✗ | ✗ | ✗ |
| M4 anisotropy, global | ✓ −10 | ✓ | ✓ | ✓ | ✓ but da card +34 |
| M5 homogeneity (repo) | ✓ −24 | ✗ dg −3 | ✓ | ✓ | ✗ |
| M2 + radial gain on the output | ✓ −28 | ✓ | ✓ | ✓ | ✓ but da card +20 |
| **M4, 3° window + radial gain on the output** | **✓ −17** | **✓** | **✓ (da −1)** | **✓ (+18 vs +28)** | **✓** |

Full numbers: [`fig1_sign_table.png`](paradoxes/figures/fig1_sign_table.png). "Radial gain on the
output" means the cardinal gain enters before normalization and the radial gain multiplies the
normalized response. In other words, the radial bias is kept out of the suppressive pool. This is
the split the notebook's own open question suggests: "cardinal = more neurons, polar = feedback or
lateral connections".

## The mechanism, in plain terms

**Paradox 1 needs suppression that grows faster than excitation.** In divisive normalization a
response is drive ÷ (σ + pooled drive). Extra cardinal neurons raise both the top and the bottom of
that fraction. When both rise in step (p = q, as in M1), the ratio can only flatten toward equality.
It never reverses ([`fig4_mechanism.png`](paradoxes/figures/fig4_mechanism.png) A). That is the
ordinary saturating contrast response (Heeger 1992,
[doi:10.1017/s0952523800009640](https://doi.org/10.1017/s0952523800009640)). A reversal needs the
pool to grow faster than the drive. That is supersaturation, which Peirce (2007) found in a minority
of macaque V1 neurons ([doi:10.1167/7.6.13](https://doi.org/10.1167/7.6.13)). Three repo components
do this:

- **q > p with the gain applied before the exponent** (M2). Overrepresentation then acts like extra
  contrast in those channels. If it is modelled instead as extra neurons of unchanged gain, M2
  inverts only when normalization is very strong (σ = 0.001; at σ = 0.01 it gives +5.5).
- **a homogeneity-gated exponent** (M5), which is a form of q > p.
- **an anisotropy pool** (M4). Here the suppression is the spread of energy across orientations,
  not the total. An oblique grating leaks some energy into the boosted cardinal channels. That leak
  raises the summed response more than it raises the spread. **This depends on how broad the
  channels' orientation tuning is.** With narrow tuning the leak, and the inversion, would shrink.

Tuned normalization (M3) makes paradox 1 *worse* at every setting tested (+55). An oblique grating's
cross-orientation pool contains the boosted cardinal leak, so obliques are suppressed more.

**Paradox 3 (cardinal) comes almost for free from any spatially broad pool.** A Cartesian grating
has one orientation everywhere. So the pool around any point "sees" the boosted cardinal channel,
and suppresses it. A polar grating's local orientation turns with polar angle. At 6°, a pool wider
than about 2° already mixes many orientations. At 12° the mixture is nearly balanced
([`fig4_mechanism.png`](paradoxes/figures/fig4_mechanism.png) B). That pool cannot tell a cardinal
location from an oblique one, so the extra cardinal drive passes through unsuppressed. This is why
M1, M2, M4 and M5 all give a more negative cardinal effect in dg than in da. The size of the
context effect follows the pool size
([`fig3_pool_size.png`](paradoxes/figures/fig3_pool_size.png)). With a 1° pool the context effect is
small, and da cardinal is still clearly negative. da cardinal reaches ≈ 0 at about 3° for M4, and
between 1° and 3° for M2. It is clearly *positive* by 6°. **The data (da ≈ 0) therefore point to a
pool of a few degrees at 6° eccentricity.** That is the scale of the suppressive surround:
Cavanaugh, Bair and Movshon (2002) report surrounds about 3× the classical RF diameter.

**Paradox 2 conflicts with paradox 1 if both biases enter the same pool.** The two biases differ
only in where they sit. The cardinal bias is the same everywhere. The radial bias turns with polar
angle. The two stimulus sets swap these roles. A Cartesian grating hits the cardinal bias
everywhere. A pinwheel hits the radial bias everywhere. So any pool that suppresses the cardinal
bias in dg also suppresses the radial bias in da. This is what M2 and M5 do: in da they weaken or
reverse the radial effect. The data are not symmetric under this swap. Cardinal is suppressed in
dg, while radial is *enhanced* in da. **So the two biases must enter the model differently.** Two
ways work:

- Keep the radial bias out of the suppressive pool (the "radial on output" rows). The radial effect
  is then positive in both sets for every model, σ and gain set tested.
- Use a pool that ignores a gain shared equally by all orientations (M4). Averaged over the whole
  field, a pinwheel's radial boost falls on every orientation channel equally. A std-across-
  orientation pool then does not see it. This also lets the polar gratings escape normalization
  almost entirely, which is why global M4 leaves da cardinal at +34.

## What is robust and what is not

**Robust** (holds at every σ from 0.001 to 1, at all three gain sets tested, and at every pool or
window size tested):

- With no normalization, or with p = q (M1), overrepresentation always gives a *larger* response.
  These models cannot produce paradox 1.
- Tuned normalization (M3) always makes cardinals *stronger*.
- Whenever paradox 1 is produced, the cardinal effect is less negative in da than in dg. The gap
  is small with a 1° pool and grows with pool size. That is, paradox 3 follows from paradox 1 plus
  a broad pool.
- With the radial gain outside the pool, the radial effect is positive in both sets.
- Gain sets: the repo's (2, 1.5), behavioural values (1.3, 1.1), and macaque V1 values (1.11, 1.15)
  give the same signs. Smaller gains give smaller effects.

**Depends on parameters:**

- **Every inversion needs strong normalization**, with σ well below the pooled drive. All models
  return to the linear answer (cardinal stronger) once σ ≥ 0.1
  ([`fig2_sigma.png`](paradoxes/figures/fig2_sigma.png)). **The model stimulus is at 100% contrast,
  but the experiment used 50%.** Energy scales as contrast², so at 50% the repo's σ acts like
  a σ that is 4× larger (p = q) to 16× larger (q = 2). At that equivalent value, the best variant's
  dg cardinal effect is ≈ 0 or slightly positive, and M2 and M5 no longer invert. So with the repo's σ values
  **at the real contrast, none of the models reproduces paradox 1.** It needs σ ≲ 0.005 in
  100%-contrast units. σ is a free parameter, so this is not fatal, but it should be stated.
- **The da cardinal value goes from negative to ≈ 0 to strongly positive with pool size** (1° →
  3° → 6°). Matching the data's ≈ 0 needs a pool of about 3°.
- **The radial effect being larger in da is the weakest result.** It appears in the best variant
  (+28 vs +18) only because the polar gratings are less normalized overall, so the output gain
  multiplies a larger response. That rests on the overall response level, which the data cannot
  check (pink-noise baseline). The harmonic supplement also reports that polar gratings give about
  half the overall asymmetry amplitude. Treat this one as **not explained**.
- **M4's inversion depends on orientation bandwidth**, as above. Untested beyond the pyramid's own
  bandwidth (4 orientations, order 3).

**Not addressed:** horizontal vs vertical. No variant has separate H and V gains. A version with
g_H > g_V under q > p would be expected to give H weaker than V, but **that is a guess, not
tested.**

## Problems found in the repo model (not fixed; flagged for Rania)

1. **M5 creates a horizontal–vertical asymmetry from nothing**: +7% of mean, H > V, with a balanced
   population. Its pool is a Gaussian blur along the orientation-channel index with reflecting
   ends, so channels 0 and 3 (90°, 135°) are treated differently from 1 and 2 (45°, 0°). The channel
   order already steps by −45°, so blurring that axis with `mode='wrap'` is the correct circular
   pool. With wrap the artifact is exactly zero, and the other signs barely change (fig1, "circular
   ori" row).
2. **The model's 100% contrast vs the experiment's 50%** (above).
3. Polar gratings match 1 cpd at ~7.2° (the code's "half radius"), not at 6° as in the experiment.
   This is minor for signs, but worth fixing if the model is ever compared quantitatively.

## Reproducing

Python 3.11 with the repo's `requirements-intel-mac.txt` (system Python 3.9 cannot install
plenoptic 1.3.1). From `Models/paradoxes/`:

```bash
python 00_checks.py     # ~20 s: conventions, energy scale, match to repo functions
python 01_variants.py   # ~2 min: the 6 repo models x 4 imbalances
python 02_sweeps.py     # ~21 min on an Intel i9 CPU: all sweeps; logs progress to _cache/02_sweeps.log
python 03_figures.py    # figures/ from _cache/sweeps.json
```

`_cache/` (energies, ~0.6 GB) is git-ignored.

## References

- Cavanaugh JR, Bair W, Movshon JA (2002). Nature and interaction of signals from the receptive
  field center and surround in macaque V1 neurons. *J Neurophysiol* 88:2530–2546. (Cited from
  `SUMMARY.md`; not re-checked here.)
- Fang C, Cai X, Lu HD (2022). Orientation anisotropies in macaque visual areas. *PNAS*
  119:e2113407119. [doi](https://doi.org/10.1073/pnas.2113407119)
- Heeger DJ (1992). Normalization of cell responses in cat striate cortex. *Vis Neurosci*
  9:181–197. [doi](https://doi.org/10.1017/s0952523800009640)
- Peirce JW (2007). The potential importance of saturating and supersaturating contrast response
  functions in visual cortex. *J Vis* 7(6):13. [doi](https://doi.org/10.1167/7.6.13)

Related, not relied on above:

- Furmanski & Engel (2000) found the opposite fMRI pattern for isolated gratings: cardinal *larger*
  than oblique in V1 ([doi](https://doi.org/10.1038/75702)).
- Essock et al. (2003) and Hansen & Essock (2004) describe a behavioural "horizontal effect":
  horizontal content is seen *worst* in broadband images, and they explain it by
  overrepresentation-driven suppression ([doi](https://doi.org/10.1016/s0042-6989(03)00142-1),
  [doi](https://doi.org/10.1167/4.12.5)).
- Sasaki et al. (2006) report the fMRI radial bias ([doi](https://doi.org/10.1016/j.neuron.2006.07.021)).
- Roth, Heeger and Merriam (2018) show that aperture edges ("vignetting") can create apparent
  orientation biases ([doi](https://doi.org/10.7554/eLife.37241)). In this model the aperture edges
  are far from 4–8°: with a balanced population, the radial − tangential effect for Cartesian
  gratings is 0.0%, so vignetting does not contribute here.
- Coen-Cagli, Kohn and Schwartz (2015) on homogeneity-gated surround suppression, which motivates
  M5 ([doi](https://doi.org/10.1038/nn.4128)).

Citations checked on PubMed on 2026-10-21, except Cavanaugh et al.
