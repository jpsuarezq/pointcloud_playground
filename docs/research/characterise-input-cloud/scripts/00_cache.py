"""Stream the whole LAS once into flat .npy memmaps in /tmp/opencode/cloud-research/out/cache.

Avoids re-reading 2.16 GB for every experiment. Coordinates stored as int32
raw LAS units (scale 1e-4, offsets from header) -> exact, compact.
"""
import os, time, numpy as np, laspy

SRC = "/home/usuario/workspaces/NubePuntos/pointclouds/segmentado_RGB_color.las"
OUT = "/tmp/opencode/cloud-research/out/cache"
os.makedirs(OUT, exist_ok=True)

f = laspy.open(SRC)
n = f.header.point_count
print("points", n)
sx, sy, sz = f.header.scales
ox, oy, oz = f.header.offsets
print("scales", sx, sy, sz, "offsets", ox, oy, oz)

dt_map = dict(
    xyz=np.dtype([("x", "<i4"), ("y", "<i4"), ("z", "<i4")]),
    rgb=np.dtype([("r", "<u2"), ("g", "<u2"), ("b", "<u2")]),
    attr=np.dtype([("fm", "u1"), ("cls", "u1"), ("intensity", "<u2")]),
)
X = np.memmap(f"{OUT}/xyz.i4", dtype=dt_map["xyz"], mode="w+", shape=(n,))
R = np.memmap(f"{OUT}/rgb.u2", dtype=dt_map["rgb"], mode="w+", shape=(n,))
A = np.memmap(f"{OUT}/attr.u1", dtype=dt_map["attr"], mode="w+", shape=(n,))

i = 0
fm_hist = np.zeros(3, dtype=np.int64)
cls_hist = np.zeros(256, dtype=np.int64)
t0 = time.time()
for pts in f.chunk_iterator(8_000_000):
    m = len(pts)
    sl = slice(i, i + m)
    # NOTE: laspy's .X/.Y/.Z are ALREADY offset-subtracted scaled ints.
    X["x"][sl] = np.asarray(pts.X, dtype=np.int64)
    X["y"][sl] = np.asarray(pts.Y, dtype=np.int64)
    X["z"][sl] = np.asarray(pts.Z, dtype=np.int64)
    fm = np.asarray(pts["PredSemantic_FM"]).astype(np.int64)
    fm_hist += np.bincount(np.clip(fm, 0, 2), minlength=3)
    A["fm"][sl] = fm.astype(np.uint8)
    c = np.asarray(pts.classification).astype(np.uint8)
    cls_hist += np.bincount(c, minlength=256)
    R["r"][sl] = np.asarray(pts.red)
    R["g"][sl] = np.asarray(pts.green)
    R["b"][sl] = np.asarray(pts.blue)
    A["cls"][sl] = c
    A["intensity"][sl] = np.asarray(pts.intensity)
    i += m
    print(f"  {i}/{n}  {time.time()-t0:.0f}s", flush=True)
f.close()
assert i == n, (i, n)
X.flush(); R.flush(); A.flush()
print("done", time.time() - t0)
print("FM histogram", fm_hist, fm_hist / n)
print("classification nonzero:", {k: int(v) for k, v in enumerate(cls_hist) if v})
np.save(f"{OUT}/fm_hist.npy", fm_hist)
np.save(f"{OUT}/cls_hist.npy", cls_hist)
with open(f"{OUT}/meta.txt", "w") as fh:
    fh.write(f"n={n}\nscale={sx}\noffset_x={ox}\noffset_y={oy}\noffset_z={oz}\n"
             f"raw_x = X_local_raw*scale + offset_x  (laspy .X is already offset-subtracted)\n")
print(open(f"{OUT}/meta.txt").read())