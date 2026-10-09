"""Build 0.25 m and 0.5 m per-cell min-z / count grids for the whole tile, cached.

Also caches the flat cell index per point (as int32) so later passes are cheap.
"""
import numpy as np, time, os, sys
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, N

OUT = "/tmp/opencode/cloud-research/out/grids"
os.makedirs(OUT, exist_ok=True)

X = load_xyz()
x, y, z = X["x"].astype(np.int64), X["y"].astype(np.int64), X["z"].astype(np.int64)
A = load_attr()
fm = np.asarray(A["fm"]).copy()

print("x range", x.min(), x.max(), "y", y.min(), y.max(), "z", z.min(), z.max(), flush=True)

t0 = time.time()
for res in (2500, 5000):  # 0.25 m, 0.50 m  (raw units of 1e-4 m)
    x0 = (x.min() // res) * res
    y0 = (y.min() // res) * res
    nx = (x.max() // res) - (x0 // res) + 1
    ny = (y.max() // res) - (y0 // res) + 1
    ix = ((x - x0) // res).astype(np.int64)
    iy = ((y - y0) // res).astype(np.int64)
    flat = (iy * nx + ix).astype(np.int64)
    nc = nx * ny
    print(f"res {res/10000:g}m grid {nx}x{ny}={nc}  sorting...", flush=True)

    order = np.argsort(flat, kind="stable")
    fs = flat[order]
    starts = np.searchsorted(fs, np.arange(nc + 1))
    # min z = z of first element of each group
    zmin = np.full(nc, np.nan)
    cnt = (starts[1:] - starts[:-1]).astype(np.int32)
    nonempty = starts[:-1] < starts[1:]
    zmin[nonempty] = z[order[starts[:-1][nonempty]]]
    np.save(f"{OUT}/minz_{res}.npy", zmin.reshape(ny, nx))
    np.save(f"{OUT}/cnt_{res}.npy", cnt.reshape(ny, nx))
    print("  done", time.time() - t0, "empty cells", int((cnt == 0).sum()), flush=True)
    with open(f"{OUT}/grid_{res}.txt", "w") as fh:
        fh.write(f"res={res}\nx0={x0}\ny0={y0}\nnx={nx}\nny={ny}\n")
    if res == 25:
        np.save(f"{OUT}/flat25.npy", flat.astype(np.int32))
        np.save(f"{OUT}/fm.npy", fm)
print("total", time.time() - t0)