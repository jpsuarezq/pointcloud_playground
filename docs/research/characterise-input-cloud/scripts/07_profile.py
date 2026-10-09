"""Vertical point-density profile vs height above ground, computed against
several independent ground surfaces.

Density is normalised by (number of usable 0.5 m cells) x (cell area) x (bin
thickness), and additionally slope-corrected by 1/cos(theta).
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, N, SCALE, OX, OY
from scipy import ndimage

G = "/tmp/opencode/cloud-research/out/grids"
O = "/tmp/opencode/cloud-research/out/dem"
HC = "/tmp/opencode/cloud-research/out/hag"
F = "/tmp/opencode/cloud-research/out/fig"
os.makedirs(F, exist_ok=True)

meta = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5 = int(meta["res"]); X5 = int(meta["x0"]); Y5 = int(meta["y0"])
NX5, NY5 = int(meta["nx"]), int(meta["ny"])
cnt5 = np.load(f"{G}/cnt_5000.npy")
print("0.5m grid", NY5, NX5, "occupied", (cnt5 > 0).mean())

dem5raw = np.load(f"{O}/dem_gmin1000.npy")   # 1 m grid-min for slope only
# ---- slope from the SMRF 0.25 m surface, nan-aware, resampled to the 0.5 m grid
sm = np.load(f"{O}/smrf025.npy")
sm = np.where(sm < -9000, np.nan, sm)[::-1, :]          # row 0 = south, shape (1631,800)
_, ii = ndimage.distance_transform_edt(~np.isfinite(sm), return_indices=True)
smf = sm[tuple(ii)]
smb = ndimage.uniform_filter(smf, size=5, mode="nearest")
gyy, gxx = np.gradient(smb * SCALE, 0.25, 0.25)
slope = np.degrees(np.arctan(np.hypot(gxx, gyy)))
np.save(f"{O}/slope025.npy", slope)
print("slope pct", np.nanpercentile(slope, [5, 25, 50, 75, 95, 99]).round(1))

# resample to 0.5 m cell centres of the cnt5 grid
X0250, Y0250 = int(dict(l.split("=") for l in open(f"{G}/grid_2500.txt"))["x0"]), \
               int(dict(l.split("=") for l in open(f"{G}/grid_2500.txt"))["y0"])
NY25 = slope.shape[0]
rr = (Y5 + (np.arange(NY5) + 0.5) * R5 - Y0250) / 2500.0
cc = (X5 + (np.arange(NX5) + 0.5) * R5 - X0250) / 2500.0
RR, CC = np.meshgrid(rr, cc, indexing="ij")
slope5 = ndimage.map_coordinates(slope, [RR, CC], order=1, mode="nearest")
np.save(f"{O}/slope05.npy", slope5)
usable = (cnt5 >= 20) & np.isfinite(slope5) & (slope5 < 60)
ncell = int(usable.sum())
area = ncell * 0.25
print("usable 0.5m cells", ncell, "area m2", round(area, 1),
      "median slope", np.nanmedian(slope5[usable]).round(1),
      "mean cos(theta)", np.nanmean(np.cos(np.deg2rad(slope5[usable]))).round(4))

X = load_xyz()
x = np.asarray(X["x"], dtype=np.int64); y = np.asarray(X["y"], dtype=np.int64)
fm = np.asarray(load_attr()["fm"])
del X
ix5 = np.clip(((x - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((y - Y5) // R5).astype(np.int32), 0, NY5 - 1)
inmask = usable[iy5, ix5]
print("points in usable cells", inmask.mean().round(4), "->", int(inmask.sum()))

edges = np.round(np.arange(0, 20.05, 0.1), 2)
BIN = 0.1
out = {"n_usable_cells": ncell, "area_m2": area,
       "mean_cos_theta": float(np.nanmean(np.cos(np.deg2rad(slope5[usable])))),
       "median_slope_deg": float(np.nanmedian(slope5[usable])),
       "bins": edges.tolist(), "profiles": {}}

for surf in sys.argv[1:]:
    hag = np.load(f"{HC}/hag_{surf}.npy")
    ok = np.load(f"{HC}/ok_{surf}.npy")
    sel = inmask & ok
    h = hag[sel]
    f = fm[sel]
    b = np.clip(((h) / BIN).astype(np.int32), 0, len(edges) - 2)
    tot = np.bincount(b, minlength=len(edges) - 1)
    dens = tot / (area * BIN)
    dens_slope = dens / out["mean_cos_theta"]
    prof = {"total": tot.tolist(),
            "density": dens.tolist(),
            "density_slopecorr": dens_slope.tolist()}
    for c in (0, 1, 2):
        totc = np.bincount(b[f == c], minlength=len(edges) - 1)
        prof[f"fm{c}_count"] = totc.tolist()
        prof[f"fm{c}_density"] = (totc / (area * BIN)).tolist()
    out["profiles"][surf] = prof
    print(f"{surf}: " + "  ".join(
        f"{edges[i]:.1f}-{edges[i+1]:.1f}m:{dens[i]:.0f}" for i in range(0, 60, 5)), flush=True)

json.dump(out, open(f"{F}/density_profile.json", "w"))

# ---------------- readable table ----------------
surfs = list(out["profiles"])
lines = ["bin_m  " + "  ".join(f"{s:>12s}" for s in surfs) + "   ratio_smrf_over_gmin"]
for i in range(0, 60):
    row = [out["profiles"][s]["density"][i] for s in surfs]
    r = row[0] / row[1] if "smrf025" in surfs and len(surfs) > 1 and row[1] > 0 else float("nan")
    lines.append(f"{edges[i]:.1f}-{edges[i+1]:.1f}  " + "  ".join(f"{v:12.1f}" for v in row) + f"   {r:6.2f}")
open(f"{F}/density_profile.txt", "w").write("\n".join(lines))
print("\n".join(lines[:22]))
print("wrote", f"{F}/density_profile.txt")