"""Is the 1.0-1.6 m density trough physical? Stratify the vertical profile by
(a) terrain slope and (b) canopy closure, using the SMRF ground surface.
Closure = (local DSM - local DTM) per 0.5 m cell; low-closure cells have
directly-visible ground, so every ground-surface method must agree there.
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE

G = "/tmp/opencode/cloud-research/out/grids"
O = "/tmp/opencode/cloud-research/out/dem"
HC = "/tmp/opencode/cloud-research/out/hag"
F = "/tmp/opencode/cloud-research/out/fig"

m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5, NX5, NY5 = int(m5["res"]), int(m5["x0"]), int(m5["y0"]), int(m5["nx"]), int(m5["ny"])
m25 = dict(l.split("=") for l in open(f"{G}/grid_2500.txt"))
RES, X25, Y25 = int(m25["res"]), int(m25["x0"]), int(m25["y0"])
NY25 = int(m25["ny"])

# ---- DSM at 0.25 m, then aggregate to 0.5 m
X = load_xyz()
x = np.asarray(X["x"], dtype=np.int64); y = np.asarray(X["y"], dtype=np.int64)
z = np.asarray(X["z"], dtype=np.float64) * SCALE
fm = np.asarray(load_attr()["fm"]); del X
ix25 = np.clip(((x - X25) // RES).astype(np.int64), 0, int(m25["nx"]) - 1)
iy25 = np.clip(((y - Y25) // RES).astype(np.int64), 0, NY25 - 1)
dsm = np.full((NY25, int(m25["nx"])), -np.inf)
np.maximum.at(dsm, (iy25, ix25), z)
dsm[~np.isfinite(dsm)] = np.nan
np.save(f"{O}/dsm025.npy", dsm)

# DTM = SMRF, resampled (row0 = south)
sm = np.where(np.load(f"{O}/smrf025.npy") < -9000, np.nan, np.load(f"{O}/smrf025.npy"))[::-1, :]
rr = (Y5 + (np.arange(NY5) + 0.5) * R5 - Y25) / RES
cc = (X5 + (np.arange(NX5) + 0.5) * R5 - X25) / RES
RR, CC = np.meshgrid(rr, cc, indexing="ij")
from scipy import ndimage
dtm5 = ndimage.map_coordinates(sm, [RR, CC], order=1, mode="nearest")
dsm5 = ndimage.map_coordinates(dsm, [RR, CC], order=1, mode="nearest")
closure = dsm5 - dtm5
cnt5 = np.load(f"{G}/cnt_5000.npy")
slope5 = np.load(f"{O}/slope05.npy")
okc = np.isfinite(closure) & (cnt5 > 0)
print("closure pct", np.nanpercentile(closure[okc], [5, 25, 50, 75, 95]).round(1))

ix5 = np.clip(((x - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((y - Y5) // R5).astype(np.int32), 0, NY5 - 1)
inb = okc[iy5, ix5]

BIN = 0.05; NB = 100; edges = np.round(np.arange(NB + 1) * BIN, 3)
SLO = np.digitize(slope5, [15, 25, 35, 45])
CLO = np.digitize(np.clip(closure, -5, 60), [2, 8, 15])
STRATA = {"slope": [("0-15", 0), ("15-25", 1), ("25-35", 2), ("35-45", 3), ("45+", 4)],
          "closure": [("<2m", 0), ("2-8m", 1), ("8-15m", 2), (">15m", 3)]}
LAB = {"slope": SLO, "closure": CLO}

out = {}
for kind, strata in STRATA.items():
    lab = LAB[kind]
    ncell = {nm: int((okc & (lab == i)).sum()) for nm, i in strata}
    out[kind] = {nm: ncell[nm] for nm, _ in strata}
    print(kind, "cells", out[kind], "area", {k: round(v * 0.25) for k, v in out[kind].items()})

for surf in sys.argv[1:]:
    hag = np.load(f"{HC}/hag_{surf}.npy"); ok = np.load(f"{HC}/ok_{surf}.npy")
    for kind, strata in STRATA.items():
        lab = LAB[kind]
        for nm, i in strata:
            nc = int((okc & (lab == i)).sum())
            if nc < 500:
                continue
            sel = inb & ok & (lab[iy5, ix5] == i)
            b = np.clip((hag[sel] / BIN).astype(np.int32), 0, NB - 1)
            cnt = np.bincount(b, minlength=NB)
            dens = cnt / (nc * 0.25 * BIN)
            seg = dens[3:80]                       # 0.15 - 4.0 m
            kmin = int(np.argmin(seg)) + 3
            # local peak in 1.8-3.2 m
            p = slice(36, 64)
            kmax = int(np.argmax(seg[p])) + 36
            out.setdefault(surf, {}).setdefault(kind, {})[nm] = dict(
                ncells=nc, area=nc * 0.25, total=float(dens.sum() * BIN),
                dmin=float(dens[kmin]), hag_min=float(edges[kmin]),
                dmax=float(dens[kmax]), hag_max=float(edges[kmax]),
                ratio=float(dens[kmax] / dens[kmin]),
                d_0_05=float(dens[1]), d_0_55=float(dens[11]),
                d_1_25=float(dens[25]), d_2_55=float(dens[55]),
                prof=[round(float(v), 1) for v in dens])
            print(f"{surf:9s} {kind:8s} {nm:7s} cells {nc:6d} area {nc*0.25:7.0f}  "
                  f"dmin {dens[kmin]:6.1f} @{edges[kmin]:.2f}  dmax {dens[kmax]:6.1f} @{edges[kmax]:.2f}  "
                  f"ratio {dens[kmax]/dens[kmin]:5.2f}  d0.05 {dens[1]:7.1f} d0.55 {dens[11]:6.1f} "
                  f"d1.25 {dens[25]:6.1f} d2.55 {dens[55]:6.1f}", flush=True)
json.dump(out, open(f"{F}/strata_profile.json", "w"), indent=1)