"""One-pass comparison of every candidate ground surface:
coverage, fraction of points below the surface, HAG distribution, and the
depth/position of the 1.0-1.6 m density trough. Also slope-stratified.
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE, OX, OY
from scipy import ndimage

G = "/tmp/opencode/cloud-research/out/grids"; O = "/tmp/opencode/cloud-research/out/dem"
F = "/tmp/opencode/cloud-research/out/fig"
m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5, NX5, NY5 = int(m5["res"]), int(m5["x0"]), int(m5["y0"]), int(m5["nx"]), int(m5["ny"])
cnt5 = np.load(f"{G}/cnt_5000.npy"); slope5 = np.load(f"{O}/slope05.npy")
usable = (cnt5 >= 20) & np.isfinite(slope5) & (slope5 < 60)
ncell = int(usable.sum()); area = ncell * 0.25
slopelab = np.digitize(slope5, [15, 25, 35, 45])

X = load_xyz()
xr = np.asarray(X["x"], dtype=np.int64) * SCALE + OX
yr = np.asarray(X["y"], dtype=np.int64) * SCALE + OY
zr = np.asarray(X["z"], dtype=np.float64) * SCALE
xraw = np.asarray(X["x"], dtype=np.int64); yraw = np.asarray(X["y"], dtype=np.int64)
fm = np.asarray(load_attr()["fm"]); del X
ix5 = np.clip(((xraw - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((yraw - Y5) // R5).astype(np.int32), 0, NY5 - 1)
base = usable[iy5, ix5] & (slopelab[iy5, ix5] < 4)
print("base mask frac", round(float(base.mean()), 4), "area", area)


def bilin(a, rf, cf):
    i0 = np.floor(rf).astype(np.int64); j0 = np.floor(cf).astype(np.int64)
    fx = rf - i0; fy = cf - j0
    i0c = np.clip(i0, 0, a.shape[0] - 1); i1c = np.clip(i0 + 1, 0, a.shape[0] - 1)
    j0c = np.clip(j0, 0, a.shape[1] - 1); j1c = np.clip(j0 + 1, 0, a.shape[1] - 1)
    return (a[i0c, j0c] * (1 - fx) * (1 - fy) + a[i0c, j1c] * fx * (1 - fy) +
            a[i1c, j0c] * (1 - fx) * fy + a[i1c, j1c] * fx * fy)


BIN = 0.05; NB = 100; edges = np.round(np.arange(NB + 1) * BIN, 3)
out = {"area_m2": area, "n_cells": ncell, "edges": edges.tolist(), "surfaces": {}}
REF_P50 = float(os.environ.get("REF_P50", "nan"))
print(f"{'surface':12s} {'cov%':>6s} {'<0%':>6s} {'p1':>7s} {'p25':>7s} {'p50':>7s} "
      f"{'trough':>7s} {'dmin':>7s} {'peak':>6s} {'dmax':>7s} {'ratio':>6s} {'d1.25':>7s}")
for spec in sys.argv[1:]:
    parts = spec.split(":")
    kind, name = parts[0], parts[1]
    delta = float(parts[2]) if len(parts) > 2 else 0.0
    if kind == "tif":
        j = json.load(open(f"{O}/{name}.json")); a = np.load(f"{O}/{name}.npy").astype(np.float64)
        a = np.where(a < -9000, np.nan, a)
        gx = (xr - j["origin"][0]) / j["res"][0]; gy = (j["origin"][1] - yr) / j["res"][1]
        okg = (gx >= 0) & (gx < a.shape[1] - 1.001) & (gy >= 0) & (gy < a.shape[0] - 1.001)
    else:
        res = int(kind)
        gm = dict(l.split("=") for l in open(f"{G}/grid_{res}.txt"))
        gx0, gy0 = int(gm["x0"]), int(gm["y0"])
        fn = name if name.endswith(".npy") else name + ".npy"
        a = np.load(f"{O}/{fn}")
        a = np.where(np.isfinite(a), a, np.nan)
        gy = (yr - gy0 * 1e-4 - OY) * 1e4 / res
        gx = (xr - gx0 * 1e-4 - OX) * 1e4 / res
        okg = True
    hag = zr - bilin(a, gy, gx)
    ok = base & okg & np.isfinite(hag) & (hag > -50) & (hag < 200)
    sel_ok = np.where(ok)[0]
    if np.isfinite(REF_P50) and delta == 0.0:
        delta = REF_P50 - float(np.median(hag[sel_ok]))
        print(f"   [aligned {name} by {delta:+.2f} m to reference p50]")
    hag = hag - delta
    h_all = hag[sel_ok]
    # drop out-of-range HAG instead of clipping, otherwise the +/-50 m guard
    # values pile up in the first/last bin
    keep = ok & (hag >= 0) & (hag < NB * BIN)
    sel = np.where(keep)[0]
    b = (hag[sel] / BIN).astype(np.int32)
    d = np.bincount(b, minlength=NB) / (area * BIN)
    seg = d[10:80]; k = int(np.argmin(seg)) + 10
    pk = d[36:70]; kp = int(np.argmax(pk)) + 36
    rec = dict(coverage=float(ok.mean()), delta=delta, frac_below0=float(np.mean(h_all < 0)),
               pct={q: float(np.percentile(h_all, q)) for q in (1, 5, 25, 50, 75, 95, 99)},
               trough_h=float(edges[k]), d_min=float(d[k]), peak_h=float(edges[kp]),
               d_max=float(d[kp]), ratio=float(d[kp] / d[k]), d_1_25=float(d[25]),
               d_0_55=float(d[11]), prof=[round(float(v), 2) for v in d])
    # slope stratified trough
    sl = slopelab[iy5[sel], ix5[sel]]
    rs = {}
    for i, (lo, hi) in enumerate([(0, 15), (15, 25), (25, 35), (35, 45), (45, 90)]):
        nc = int((usable & (slopelab == i)).sum())
        m = sl == i
        dd = np.bincount(b[m], minlength=NB) / (nc * 0.25 * BIN)
        s2 = dd[10:80]; k2 = int(np.argmin(s2)) + 10
        p2 = dd[36:70]; kp2 = int(np.argmax(p2)) + 36
        rs[f"{lo}-{hi}"] = dict(ncells=nc, trough_h=float(edges[k2]), d_min=float(dd[k2]),
                                peak_h=float(edges[kp2]), d_max=float(dd[kp2]),
                                ratio=float(dd[kp2] / dd[k2]))
    rec["by_slope"] = rs
    out["surfaces"][name] = rec
    print(f"{name:12s} {ok.mean()*100:6.2f} {np.mean(h_all<0)*100:6.2f} {np.percentile(h_all,1):7.2f} "
          f"{np.percentile(h_all,25):7.2f} {np.percentile(h_all,50):7.2f} {edges[k]:7.2f} {d[k]:7.1f} "
          f"{edges[kp]:6.2f} {d[kp]:7.1f} {d[kp]/d[k]:6.2f} {d[25]:7.1f}", flush=True)
    np.save(f"{F}/prof_{name}.npy", d)
json.dump(out, open(f"{F}/surface_comparison.json", "w"), indent=1)
print("wrote", f"{F}/surface_comparison.json")