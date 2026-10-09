"""Independent numpy ground surfaces, implemented from scratch.

  open50 / open90 : 0.5 m grid-minimum + grey-erosion (morphological opening)
                    with a 5x5 (2.5 m) or 9x9 (4.5 m) flat element, then
                    slope-corrected back up by +win/2 * tan(theta) using the
                    slope field derived from the SMRF DEM.
  cloth05         : grid-level Cloth Simulation Filter (Zhang et al. 2016),
                    cloth_dir = -1 (upward looking). The cloth starts above the
                    highest terrain and is lowered at most `step` per iteration
                    toward  max(terrain in 3x3 neighbourhood) - dmin; friction
                    stops it where it first touches. Grid-level: 4 orders of
                    magnitude cheaper than PDAL's per-point CSF, which needs
                    >9 GB on this tile and was OOM-killed.
  gmn05           : plain 0.5 m grid minimum + 3x3 median (control).
"""
import numpy as np, os, time
from scipy import ndimage

G = "/tmp/opencode/cloud-research/out/grids"; O = "/tmp/opencode/cloud-research/out/dem"
m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5, NX5, NY5 = (int(m5["res"]), int(m5["x0"]), int(m5["y0"]),
                        int(m5["nx"]), int(m5["ny"]))
mn5 = np.load(f"{G}/minz_5000.npy").astype(np.float64) * 1e-4
cnt5 = np.load(f"{G}/cnt_5000.npy")
valid = cnt5 > 0
slope5 = np.load(f"{O}/slope05.npy")
t0 = time.time()

d = ndimage.distance_transform_edt(~valid)
_, ii = ndimage.distance_transform_edt(~valid, return_indices=True)
terr = np.where(np.isfinite(mn5), mn5, np.nan)[tuple(ii)]
terr = np.where(d <= 4, terr, np.nan)
print("terrain pct", np.nanpercentile(terr[valid], [5, 50, 95]).round(1))
print("slope pct (from SMRF)", np.nanpercentile(slope5[valid], [5, 50, 95]).round(1))

gmn = ndimage.generic_filter(terr, np.nanmedian, size=3, mode="nearest")
np.save(f"{O}/dem_gmn05.npy", np.where(valid, gmn, np.nan))

for K in (5, 9):
    er = ndimage.grey_erosion(gmn, size=K, mode="nearest")
    lift = (K - 1) / 2.0 * 0.5 * np.tan(np.deg2rad(slope5))
    op = ndimage.generic_filter(er + lift, np.nanmedian, size=3, mode="nearest")
    np.save(f"{O}/dem_open{K}0.npy", np.where(valid, op, np.nan))
    print(f"open{K}0 median change vs gmn {np.nanmedian((op - gmn)[valid]):+.2f} m  "
          f"vs grid-min {np.nanmedian((op - terr)[valid]):+.2f} m  ({time.time()-t0:.0f}s)")


def cloth_sim(terrain, step=0.65, dmin=0.5, iterations=600):
    tmax = np.nanmax(terrain)
    c = np.full(terrain.shape, tmax + 1.0)
    loc = ndimage.maximum_filter(terrain, size=3, mode="nearest")   # highest nearby terrain
    target = loc - dmin
    for it in range(iterations):
        c = np.maximum(target, c - step)
        if it % 100 == 0:
            print(f"  cloth it {it} mean {c.mean():.1f}", flush=True)
    return c, target

c, target = cloth_sim(terr)
np.save(f"{O}/dem_cloth05.npy", np.where(valid, np.minimum(c, terr), np.nan))
print("cloth05 median vs grid-min", np.nanmedian((np.minimum(c, terr) - terr)[valid]).round(2),
      f"({time.time()-t0:.0f}s)")
for nm in ["gmn05", "open50", "open90", "cloth05"]:
    a = np.load(f"{O}/dem_{nm}.npy"); v = np.isfinite(a) & valid
    print(f"{nm:8s} median {np.nanmedian(a[v]):8.2f}  p5 {np.nanpercentile(a[v],5):8.2f}  "
          f"p95 {np.nanpercentile(a[v],95):8.2f}")
print("total", time.time() - t0)