"""Numpy ground surfaces on the 0.25 m grid.

Built (all metres, y-major (ny,nx), NaN = no-data):
  gmin025  - grid minimum per 0.25 m cell, holes nearest-filled, 3x3 median
  pmf025   - Zhang et al. (2003) progressive morphological filter applied to
             the grid-minimum surface (canopy-penetrating, purely geometric)
  fm0025   - grid minimum using ONLY points labelled PredSemantic_FM==0
             (a classifier-driven surface, independent of morphology)
  gmin100 / gmin200 - coarse 1 m / 2 m grid minimum
  pdal-smrf025 / pdal-pmf025 / pdal-csf025 - written by scripts/02_pdal_dems.sh
"""
import numpy as np, os, sys, time
from scipy import ndimage

sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, N, SCALE, OX, OY

G = "/tmp/opencode/cloud-research/out/grids"
O = "/tmp/opencode/cloud-research/out/dem"
os.makedirs(O, exist_ok=True)

meta = dict(l.split("=") for l in open(f"{G}/grid_2500.txt"))
RES = int(meta["res"]); X0 = int(meta["x0"]); Y0 = int(meta["y0"])
NX, NY = int(meta["nx"]), int(meta["ny"])
RESM = RES * SCALE   # 0.25 m
print("grid", NY, NX, "res", RESM, flush=True)

cnt = np.load(f"{G}/cnt_2500.npy")
valid = cnt > 0
dist = ndimage.distance_transform_edt(~valid)
_, inds = ndimage.distance_transform_edt(~valid, return_indices=True)
interp_mask = valid | (dist <= 8.0)      # allow fill up to 2 m from data
np.save(f"{O}/interp_mask025.npy", interp_mask)
np.save(f"{O}/valid025.npy", valid)
print("interp coverage", round(float(interp_mask.mean()), 4), flush=True)


def cellmin(x, y, z, pmask=None):
    """pmask: optional per-POINT boolean keep mask."""
    ix = ((x - X0) // RES).astype(np.int64)
    iy = ((y - Y0) // RES).astype(np.int64)
    np.clip(ix, 0, NX - 1, out=ix); np.clip(iy, 0, NY - 1, out=iy)
    if pmask is not None:
        keep = pmask
        ix, iy, z = ix[keep], iy[keep], z[keep]
    flat = iy * NX + ix
    nc = NX * NY
    o = np.argsort(flat, kind="stable")
    fs = flat[o]
    st = np.searchsorted(fs, np.arange(nc + 1))
    g = np.full(nc, np.nan)
    ne = st[:-1] < st[1:]
    g[ne] = z[o[st[:-1][ne]]]
    return g.reshape(NY, NX)


def fill_holes(g, vmask, iters=8):
    if vmask.all():
        return g.copy(), vmask.copy()
    d = ndimage.distance_transform_edt(~vmask)
    _, inds = ndimage.distance_transform_edt(~vmask, return_indices=True)
    filled = g[tuple(inds)]
    return np.where(d <= iters, filled, np.nan), (vmask | (d <= iters))


X = load_xyz(); A = load_attr()
t0 = time.time()
x = np.asarray(X["x"], dtype=np.int64)
y = np.asarray(X["y"], dtype=np.int64)
z = np.asarray(X["z"], dtype=np.float64) * SCALE
fm = np.asarray(A["fm"])
print("loaded", time.time() - t0, flush=True)

# ---------- grid minimum ----------
gm = cellmin(x, y, z)
np.save(f"{O}/dem_gmin025_raw.npy", gm)
gmf, gmv = fill_holes(gm, valid)
gmf = ndimage.generic_filter(gmf, np.nanmedian, size=3, mode="nearest")
np.save(f"{O}/dem_gmin025.npy", gmf)
np.save(f"{O}/valid_gmin025.npy", gmv)
print("gmin done", time.time() - t0, flush=True)

# ---------- progressive morphological filter ----------
def pmf_open(grid, slope_deg, max_r=12, h0=0.10, hmax=1.5):
    """Zhang et al. 2003 progressive morphological filter on a 2-D grid.

    For k = 1..K, opening with a (2k+1)^2 flat SE; a cell is raised to the
    opened value when (opened - original) > h_k, h_k = h0 + (hmax-h0)*(k-1)/(K-1).
    """
    K = max_r
    out = grid.copy()
    raised = np.zeros_like(grid, dtype=bool)
    for k in range(1, K + 1):
        size = 2 * k + 1
        er = ndimage.grey_erosion(out, size=size, mode="nearest")
        op = ndimage.grey_dilation(er, size=size, mode="nearest")
        h = h0 + (hmax - h0) * (k - 1) / (K - 1)
        flag = op - out > h
        raised |= flag & np.isfinite(out)
        out = np.where(flag, op, out)
    return out, raised

dem_pmf, raised = pmf_open(gmf, 30.0, max_r=12)
# a real opening also removes isolated spikes downward; final clean-up
dem_pmf = ndimage.generic_filter(dem_pmf, np.nanmedian, size=3, mode="nearest")
np.save(f"{O}/dem_pmf025.npy", dem_pmf)
np.save(f"{O}/raised_pmf025.npy", raised)
print("pmf done", time.time() - t0, "cells raised", int(raised.sum()), flush=True)

# ---------- FM0-only surface ----------
gm0 = cellmin(x, y, z, pmask=(fm == 0))
d0 = ndimage.distance_transform_edt(np.isnan(gm0))
np.save(f"{O}/dem_fm0025_raw.npy", gm0)
gm0f, _ = fill_holes(gm0, ~np.isnan(gm0), iters=8)
gm0f = ndimage.generic_filter(gm0f, np.nanmedian, size=3, mode="nearest")
np.save(f"{O}/dem_fm0025.npy", gm0f)
np.save(f"{O}/hasfm0_025.npy", ~np.isnan(gm0))
print("fm0 done", time.time() - t0, "cells with FM0", int((~np.isnan(gm0)).sum()), flush=True)

# ---------- coarse grids ----------
for res, tag in ((10000, "100"), (20000, "200")):
    m = np.load(f"{G}/minz_{res}.npy").astype(np.float64) * SCALE
    c = np.load(f"{G}/cnt_{res}.npy")
    v = c > 0
    dd = ndimage.distance_transform_edt(~v)
    _, ii = ndimage.distance_transform_edt(~v, return_indices=True)
    ff = np.where(dd <= 4, m[tuple(ii)], np.nan)
    ff = ndimage.generic_filter(ff, np.nanmedian, size=3, mode="nearest")
    np.save(f"{O}/dem_gmin{tag}0.npy", ff)
    print(tag, "done", time.time() - t0, flush=True)
print("total", time.time() - t0)