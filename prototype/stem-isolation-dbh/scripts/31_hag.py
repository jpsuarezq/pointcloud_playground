"""Build a 0.25 m ground DEM from the SMRF ground points and compute height
above ground (HAG) for every point of the tile.

DEM = grid-minimum of SMRF ground points per 0.25 m cell, holes nearest-filled
up to 2 m, then a 3x3 median. Bilinear-sampled per point. ok = cell within the
filled footprint. Writes out/hag/hag.f4 and out/hag/ok.b1.
"""
import os, json, time, numpy as np, laspy
from scipy import ndimage

OUT = "/tmp/opencode/dbh/out"
C = f"{OUT}/cache"
H = f"{OUT}/hag"
os.makedirs(H, exist_ok=True)
RES = 0.25

meta = json.load(open(f"{C}/meta.json"))
n = meta["n"]
xyz = np.memmap(f"{C}/xyz.f4", dtype=np.float32, mode="r", shape=(n, 3))
x = np.asarray(xyz[:, 0], dtype=np.float64) + meta["ox"]
y = np.asarray(xyz[:, 1], dtype=np.float64) + meta["oy"]
z = np.asarray(xyz[:, 2], dtype=np.float64) + meta["oz"]

x0 = np.floor(x.min() / RES) * RES
y0 = np.floor(y.min() / RES) * RES
nx = int(np.ceil((x.max() - x0) / RES)) + 1
ny = int(np.ceil((y.max() - y0) / RES)) + 1
print("grid", ny, nx, "res", RES, flush=True)

# ---- ground points
gx, gy, gz = [], [], []
with laspy.open(f"{OUT}/ground_smrf025.laz") as f:
    for pts in f.chunk_iterator(8_000_000):
        gx.append(np.asarray(pts.x, dtype=np.float64))
        gy.append(np.asarray(pts.y, dtype=np.float64))
        gz.append(np.asarray(pts.z, dtype=np.float64))
gx = np.concatenate(gx); gy = np.concatenate(gy); gz = np.concatenate(gz)
print("ground points", len(gx), flush=True)

ix = np.clip(((gx - x0) // RES).astype(np.int64), 0, nx - 1)
iy = np.clip(((gy - y0) // RES).astype(np.int64), 0, ny - 1)
flat = iy * nx + ix
order = np.argsort(flat, kind="stable")
fs = flat[order]
st = np.searchsorted(fs, np.arange(nx * ny + 1))
dem = np.full(nx * ny, np.nan)
ne = st[:-1] < st[1:]
dem[ne] = gz[order[st[:-1][ne]]]
dem = dem.reshape(ny, nx)
valid = np.isfinite(dem)
print("valid cells", int(valid.sum()), "of", dem.size, flush=True)

d = ndimage.distance_transform_edt(~valid)
_, inds = ndimage.distance_transform_edt(~valid, return_indices=True)
filled = dem[tuple(inds)]
dem = np.where(d <= 8, filled, np.nan)          # fill holes up to 2 m
dem = ndimage.generic_filter(dem, np.nanmedian, size=3, mode="nearest")
np.save(f"{H}/dem_smrf025.npy", dem)
np.save(f"{H}/grid.npy", np.array([x0, y0, nx, ny]))

# ---- bilinear sample at every point
def bilin(a, rowf, colf):
    i0 = np.floor(rowf).astype(np.int64); j0 = np.floor(colf).astype(np.int64)
    fx = rowf - i0; fy = colf - j0
    i0c = np.clip(i0, 0, a.shape[0] - 1); i1c = np.clip(i0 + 1, 0, a.shape[0] - 1)
    j0c = np.clip(j0, 0, a.shape[1] - 1); j1c = np.clip(j0 + 1, 0, a.shape[1] - 1)
    return (a[i0c, j0c] * (1 - fx) * (1 - fy) + a[i0c, j1c] * fx * (1 - fy) +
            a[i1c, j0c] * (1 - fx) * fy + a[i1c, j1c] * fx * fy)

rowf = (y - y0) / RES
colf = (x - x0) / RES
inb = (rowf >= 0) & (rowf < ny - 1.001) & (colf >= 0) & (colf < nx - 1.001)
gz_s = bilin(dem, rowf, colf)
hag = z - gz_s
ok = inb & np.isfinite(hag) & (hag > -50) & (hag < 200)
hag = np.where(ok, hag, 0.0)
np.save(f"{H}/hag.npy", hag.astype(np.float32))
np.save(f"{H}/ok.npy", ok)
print("ok frac", round(float(ok.mean()), 4), "hag pct", np.percentile(hag[ok], [5, 50, 95]).round(1))
print("done")
