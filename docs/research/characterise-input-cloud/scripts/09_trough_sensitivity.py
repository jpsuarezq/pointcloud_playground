"""Sensitivity of the 1.0-1.6 m trough to a vertical offset / alternative of the
ground surface, and to restricting the point set to a single PredSemantic_FM
class. If the trough is a surface artefact it should move or vanish when the
surface is shifted; if it is a property of the capture it should persist.
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE, OX, OY
from scipy import ndimage

G = "/tmp/opencode/cloud-research/out/grids"; O = "/tmp/opencode/cloud-research/out/dem"
F = "/tmp/opencode/cloud-research/out/fig"
m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5, NX5, NY5 = int(m5["res"]), int(m5["x0"]), int(m5["y0"]), int(m5["nx"]), int(m5["ny"])
j = json.load(open(f"{O}/smrf025.json"))
sm = np.where(np.load(f"{O}/smrf025.npy") < -9000, np.nan, np.load(f"{O}/smrf025.npy")).astype(np.float64)
NYs, NXs = sm.shape
cnt5 = np.load(f"{G}/cnt_5000.npy"); slope5 = np.load(f"{O}/slope05.npy")
usable = (cnt5 >= 20) & np.isfinite(slope5) & (slope5 < 60)
ncell = int(usable.sum()); area = ncell * 0.25

X = load_xyz()
xr = np.asarray(X["x"], dtype=np.int64) * SCALE + OX
yr = np.asarray(X["y"], dtype=np.int64) * SCALE + OY
zr = np.asarray(X["z"], dtype=np.float64) * SCALE
fm = np.asarray(load_attr()["fm"]); del X
Xraw = (np.asarray(load_xyz()["x"], dtype=np.int64)); Yraw = np.asarray(load_xyz()["y"], dtype=np.int64)
ix5 = np.clip(((Xraw - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((Yraw - Y5) // R5).astype(np.int32), 0, NY5 - 1)
inb = usable[iy5, ix5]

gx = (xr - j["origin"][0]) / j["res"][0]; gy = (j["origin"][1] - yr) / j["res"][1]
inb &= (gx >= 0) & (gx < NXs - 1.001) & (gy >= 0) & (gy < NYs - 1.001)


def bilin(a, rf, cf):
    i0 = np.floor(rf).astype(np.int64); j0 = np.floor(cf).astype(np.int64)
    fx = rf - i0; fy = cf - j0
    i0c = np.clip(i0, 0, a.shape[0] - 1); i1c = np.clip(i0 + 1, 0, a.shape[0] - 1)
    j0c = np.clip(j0, 0, a.shape[1] - 1); j1c = np.clip(j0 + 1, 0, a.shape[1] - 1)
    return (a[i0c, j0c] * (1 - fx) * (1 - fy) + a[i0c, j1c] * fx * (1 - fy) +
            a[i1c, j0c] * (1 - fx) * fy + a[i1c, j1c] * fx * fy)


BIN = 0.05; NB = 100; edges = np.round(np.arange(NB + 1) * BIN, 3)


def trough(dens):
    seg = dens[10:80]                      # 0.5 - 4.0 m
    k = int(np.argmin(seg)) + 10
    pk = dens[36:70]
    kp = int(np.argmax(pk)) + 36
    return float(edges[k]), float(dens[k]), float(edges[kp]), float(dens[kp]), float(dens[kp] / dens[k])


rows = []
for off in (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0):
    a = sm - off
    hag = zr - bilin(a, gy, gx)
    sel = inb & np.isfinite(hag)
    b = np.clip((hag[sel] / BIN).astype(np.int32), 0, NB - 1)
    d = np.bincount(b, minlength=NB) / (area * BIN)
    hm, dm, hp, dp, r = trough(d)
    print(f"surface SMRF{off:+.2f}m : trough {hm:.2f}m d={dm:6.1f}  peak {hp:.2f}m d={dp:6.1f}  ratio {r:5.2f}  d1.25={d[25]:.1f}", flush=True)
    rows.append(dict(surface=f"smrf{off:+.2f}", trough_h=hm, d_min=dm, peak_h=hp, d_max=dp, ratio=r,
                     d1_25=float(d[25]), prof=[round(float(v), 1) for v in d]))

# same but restricted to canopy points only (FM2) and to FM1 only
for cls in (0, 1, 2):
    sel = inb & (fm == cls)
    b = np.clip(((zr - bilin(sm, gy, gx))[sel] / BIN).astype(np.int32), 0, NB - 1)
    d = np.bincount(b, minlength=NB) / (area * BIN)
    hm, dm, hp, dp, r = trough(d)
    print(f"FM{cls} only      : trough {hm:.2f}m d={dm:6.1f}  peak {hp:.2f}m d={dp:6.1f}  ratio {r:5.2f}  total {d.sum()*BIN:.0f} pts/m2", flush=True)
    rows.append(dict(surface=f"fm{cls}only", trough_h=hm, d_min=dm, peak_h=hp, d_max=dp, ratio=r,
                     d1_25=float(d[25]), prof=[round(float(v), 1) for v in d]))

json.dump(dict(edges=edges.tolist(), rows=rows), open(f"{F}/trough_sensitivity.json", "w"), indent=1)
print("wrote", f"{F}/trough_sensitivity.json")