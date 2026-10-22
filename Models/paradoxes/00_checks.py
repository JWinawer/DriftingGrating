"""Step 0: does the model run, are the conventions what the code assumes, and do
the re-implemented pools match the originals? Prints energy magnitudes too,
because every exponent in the normalization models acts on absolute energy."""
import time, json, numpy as np
import paradox_lib as L

t0 = time.time()
print("channel bar orientations:", L.check_channel_orientations())
E = L.stimulus_energy()
print(f"energies built in {time.time()-t0:.1f}s")
masks = L.wedge_masks()
bar = L.local_bar_orientation(E, masks)
print("dominant local bar orientation per wedge (rows = stimuli)")
print(" wedges   ", L.WEDGES)
for i, names in enumerate([L.CART_NAMES, L.POLAR_NAMES]):
    for j, n in enumerate(names):
        print(f" {n:8s}", bar[i, j].astype(int))
_, ecc = L.angle_maps(); band = (ecc >= 4) & (ecc <= 8)
for (i, j), e in E.items():
    tot = e.sum(axis=(0, 1))[band]
    print(f"set {i} stim {j}: summed energy in 4-8 deg: mean {tot.mean():.3e}  "
          f"max single channel {e[..., band].max():.3e}")
t1 = time.time()
print("re-implementation vs original (max relative error):", L.check_against_originals(E))
print(f"checks took {time.time()-t1:.1f}s; total {time.time()-t0:.1f}s")
L.os.makedirs(L.CACHE, exist_ok=True); np.save(L.os.path.join(L.CACHE, "bar_ori.npy"), bar)
