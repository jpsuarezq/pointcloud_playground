import numpy as np, sys, time, os
from scipy import ndimage
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, N, SCALE

G = "/tmp/opencode/cloud-research/out/grids"
O = "/tmp/opencode/cloud-research/out/dem"
X = load_xyz()
x = np.asarray(X["x"], dtype=np.int64); y = np.asarray(X["y"], dtype=np.int64)
z = np.asarray(X["z"], dtype=np.float64) * SCALE
t0 = time.time()
for res, tag in ((10000, "100"), (20000, "200")):
    x0 = (x.min() // res) * res; y0 = (y.min() // res) * res
    nx = (x.max() // res) - (x0 // res) + 1; ny = (y.max() // res) - (y0 // res) + 1
    ix = ((x - x0) // res).astype(np.int64); iy = ((y - y0) // res).astype(np.int64)
    flat = iy * nx + ix; nc = nx * ny
    o = np.argsort(flat, kind="stable"); fs = flat[o]
    st = np.searchsorted(fs, np.arange(nc + 1))
    g = np.full(nc, np.nan); ne = st[:-1] < st[1:]
    g[ne] = z[o[st[:-1][ne]]]
    g = g.reshape(ny, nx); cnt = (st[1:] - st[:-1]).astype(np.int32).reshape(ny, nx)
    np.save(f"{G}/minz_{res}.npy", g * 1e4)
    np.save(f"{G}/cnt_{res}.npy", cnt)
    open(f"{G}/grid_{res}.txt", "w").write(f"res={res}\nx0={x0}\ny0={y0}\nnx={nx}\nny={ny}\n")
    v = cnt > 0
    d = ndimage.distance_transform_edt(~v); _, ii = ndimage.distance_transform_edt(~v, return_indices=True)
    ff = np.where(d <= 4, g[tuple(ii)], np.nan)
    ff = ndimage.generic_filter(ff, np.nanmedian, size=3, mode="nearest")
    np.save(f"{O}/dem_gmin{tag}0.npy", ff)
    print(tag, g.shape, "empty", int((~v).sum()), time.time() - t0, flush=True)