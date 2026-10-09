"""Vertical point-density profile vs height above ground, per ground surface,
slope-stratified. 05 m bins 0-6 m.
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE
from scipy import ndimage

G = "/tmp/opencode/cloud-research/out/grids"
O = "/tmp/opencode/cloud-research/out/dem"
HC = "/tmp/opencode/cloud-research/out/hag"
F = "/tmp/opencode/cloud-research/out/fig"
os.makedirs(F, exist_ok=True)

m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5 = int(m5["res"]), int(m5["x0"]), int(m5["y0"])
NX5, NY5 = int(m5["nx"]), int(m5["ny"])
m25 = dict(l.split("=") for l in open(f"{G}/grid_2500.txt"))
X25, Y25 = int(m25["x0"]), int(m25["y0"])
cnt5 = np.load(f"{G}/cnt_5000.npy")

# ---- slope from SMRF 0.25 m surface
sm = np.where(np.load(f"{O}/smrf025.npy") < -9000, np.nan, np.load(f"{O}/smrf025.npy"))[::-1, :]
_, ii = ndimage.distance_transform_edt(~np.isfinite(sm), return_indices=True)
smf = sm[tuple(ii)]
smb = ndimage.uniform_filter(smf, size=5, mode="nearest")
gyy, gxx = np.gradient(smb, 0.25, 0.25)
slope025 = np.degrees(np.arctan(np.hypot(gxx, gyy)))
np.save(f"{O}/slope025.npy", slope025)
print("slope pct (0.25m grid)", np.nanpercentile(slope025, [5, 25, 50, 75, 95, 99]).round(1))

rr = (Y5 + (np.arange(NY5) + 0.5) * R5 - Y25) / 2500.0
cc = (X5 + (np.arange(NX5) + 0.5) * R5 - X25) / 2500.0
RR, CC = np.meshgrid(rr, cc, indexing="ij")
slope5 = ndimage.map_coordinates(slope025, [RR, CC], order=1, mode="nearest")
np.save(f"{O}/slope05.npy", slope5)
usable = (cnt5 >= 20) & np.isfinite(slope5) & (slope5 < 60)
ncell = int(usable.sum()); area = ncell * 0.25
print("usable 0.5 m cells", ncell, "area", round(area, 1),
      "median slope", round(float(np.nanmedian(slope5[usable])), 1),
      "mean cos", round(float(np.nanmean(np.cos(np.deg2rad(slope5[usable])))), 4))

X = load_xyz()
x = np.asarray(X["x"], dtype=np.int64); y = np.asarray(X["y"], dtype=np.int64)
fm = np.asarray(load_attr()["fm"]); del X
ix5 = np.clip(((x - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((y - Y5) // R5).astype(np.int32), 0, NY5 - 1)
inmask = usable[iy5, ix5]

BIN = 0.05; NB = 120; edges = np.round(np.arange(NB + 1) * BIN, 3)
SLOPEBINS = [(0, 15), (15, 25), (25, 35), (35, 45), (45, 90)]
slope_lab = np.digitize(slope5, [15, 25, 35, 45])

res = {"area_m2": area, "n_cells": ncell, "bin": BIN, "edges": edges.tolist(),
       "mean_cos_theta": float(np.nanmean(np.cos(np.deg2rad(slope5[usable])))),
       "median_slope": float(np.nanmedian(slope5[usable])), "profiles": {}}

for surf in sys.argv[1:]:
    hag = np.load(f"{HC}/hag_{surf}.npy"); ok = np.load(f"{HC}/ok_{surf}.npy")
    sel = inmask & ok
    b = np.clip((hag[sel] / BIN).astype(np.int32), 0, NB - 1)
    sl = slope_lab[iy5[sel], ix5[sel]]
    f = fm[sel]
    P = {"all": np.bincount(b, minlength=NB).tolist()}
    for i, (a, c) in enumerate(SLOPEBINS):
        m = sl == i
        P[f"slope_{a}_{c}"] = np.bincount(b[m], minlength=NB).tolist()
        P[f"n_cells_slope_{a}_{c}"] = int((slope_lab[usable] == i).sum())
    for c in (0, 1, 2):
        P[f"fm{c}"] = np.bincount(b[f == c], minlength=NB).tolist()
    res["profiles"][surf] = P
    d = np.array(P["all"]) / (area * BIN)
    print(f"{surf:10s} total {d.sum()*BIN:.0f} pts/m2  min at {edges[int(np.argmin(d[3:80]))+3]:.2f}-"
          f"{edges[int(np.argmin(d[3:80]))+4]:.2f} m  dmin {d[3:80].min():.0f}  "
          f"peak2-4m {d[40:80].max():.0f}  ratio {d[40:80].max()/d[3:80].min():.2f}", flush=True)

json.dump(res, open(f"{F}/density_profile2.json", "w"))

surfs = list(res["profiles"])
L = ["bin_m  " + "  ".join(f"{s:>11s}" for s in surfs) + "   min/max"]
for i in range(0, 80):
    row = np.array([res["profiles"][s]["all"][i] / (area * BIN) for s in surfs])
    L.append(f"{edges[i]:.2f}  " + "  ".join(f"{v:11.1f}" for v in row) +
             f"   {row.max()/max(row.min(),1e-9):6.2f}")
open(f"{F}/density_profile2.txt", "w").write("\n".join(L))
print("\n".join(L[::2][:42]))