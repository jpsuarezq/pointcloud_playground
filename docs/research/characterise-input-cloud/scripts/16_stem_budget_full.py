"""Full-tile stem point budget: KNN geometry on all 1.15-1.45 m band points
(533k -- cheap), cluster the woody subset, calibrate points-on-stem vs diameter.
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE, OX, OY
from scipy.spatial import cKDTree

HC = "/tmp/opencode/cloud-research/out/hag"; F = "/tmp/opencode/cloud-research/out/fig"
G = "/tmp/opencode/cloud-research/out/grids"; O = "/tmp/opencode/cloud-research/out/dem"
m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5, NX5, NY5 = int(m5["res"]), int(m5["x0"]), int(m5["y0"]), int(m5["nx"]), int(m5["ny"])
cnt5 = np.load(f"{G}/cnt_5000.npy"); slope5 = np.load(f"{O}/slope05.npy")
usable = (cnt5 >= 20) & np.isfinite(slope5) & (slope5 < 60)
area = int(usable.sum()) * 0.25

X = load_xyz()
xr = np.asarray(X["x"], dtype=np.int64) * SCALE + OX
yr = np.asarray(X["y"], dtype=np.int64) * SCALE + OY
zr = np.asarray(X["z"], dtype=np.float64) * SCALE
xraw = np.asarray(X["x"], dtype=np.int64); yraw = np.asarray(X["y"], dtype=np.int64); del X
hag = np.load(f"{HC}/hag_smrf025.npy"); okh = np.load(f"{HC}/ok_smrf025.npy")
fm = np.asarray(load_attr()["fm"])
ix5 = np.clip(((xraw - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((yraw - Y5) // R5).astype(np.int32), 0, NY5 - 1)
valid = usable[iy5, ix5] & okh

sel = np.where(valid & (hag >= 1.15) & (hag <= 1.45))[0]
P = np.column_stack([xr[sel], yr[sel], zr[sel]])
print("band points", len(sel), " per m2 of plot:", round(len(sel) / area, 2))
t = cKDTree(P); K = 16
dd, nb = t.query(P, k=K + 1, workers=-1)
Q = P[nb]; Q -= Q.mean(1, keepdims=True)
cov = np.einsum("nki,nkj->nij", Q, Q) / K
ev, evec = np.linalg.eigh(cov); ev = np.clip(ev[:, ::-1], 1e-12, None)
lin = (ev[:, 0] - ev[:, 1]) / ev[:, 0]
pla = (ev[:, 1] - ev[:, 2]) / ev[:, 0]
vert = np.abs(evec[:, :, 2][:, 2])
rho_all = K / (np.pi * dd[:, K] ** 2)
np.save(f"{F}/band_pca.npy", np.column_stack(
    [sel.astype(np.int64), lin, pla, vert, rho_all]).astype(np.float32))
print("band median r_k %.3f m -> surface density %.0f pts/m2 (all band surfaces)"
      % (np.median(dd[:, K]), np.median(rho_all)))
woody = (lin > .6) & (vert > .8)
print("band geometry: linear>0.6 %.2f%%  vert>0.8 %.2f%%  woody %.2f%%  leafy(lin<0.25) %.2f%%"
      % (100 * np.mean(lin > .6), 100 * np.mean(vert > .8),
         100 * woody.mean(), 100 * np.mean(lin < .25)))
print("woody fraction per FM class in the band:")
for c in (0, 1, 2):
    m = fm[sel] == c
    if m.sum() > 50:
        print("  FM%d n=%8d woody %.2f%%" % (c, m.sum(), 100 * float(woody[m].mean())))

m = woody
idx = np.where(m)[0]
W = P[idx]
print("woody points in band:", len(W))
pairs = cKDTree(W).query_pairs(0.10, output_type="ndarray")
par = np.arange(len(W))
def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for a, b in pairs:
    ra, rb = find(a), find(b)
    if ra != rb: par[ra] = rb
uq, lab = np.unique(np.array([find(i) for i in range(len(W))]), return_inverse=True)
sizes = np.bincount(lab)
cent = np.array([W[lab == u].mean(0) for u in range(len(uq))])
print("woody clusters:", len(uq), " >=6 pts:", int((sizes >= 6).sum()))
t3 = cKDTree(P[:, :2])
rows = []
for u in range(len(uq)):
    if sizes[u] < 6: continue
    own = W[lab == u]
    d = 2 * np.percentile(np.hypot(*(own - cent[u]).T), 75)
    if d < 0.03 or d > 1.2: continue
    alld = np.hypot(*(P[:, :2] - cent[u][:2]).T)
    rows.append((d, sizes[u], int((alld <= d / 2 + 0.03).sum()), int((alld <= d / 2 + 0.10).sum())))
R = np.array(rows)
print("usable stem candidates:", len(R), " diameter pct", np.percentile(R[:, 0], [10, 25, 50, 75, 90]).round(3))
out = {"area_m2": area, "band_points": int(len(sel)),
       "band_points_per_m2": len(sel) / area,
       "woody_frac_band": float(m.mean()),
       "n_clusters": int(len(uq)), "n_candidates": int(len(R))}
for j, lbl in [(1, "woody pts <0.10 m"), (2, "band pts r=d/2+3cm"), (3, "band pts r=d/2+10cm")]:
    A = np.polyfit(R[:, 0], R[:, j], 1); r = np.corrcoef(R[:, 0], R[:, j])[0, 1]
    print("  %-22s n = %7.1f*d %+6.1f   R=%.3f  median %d" % (lbl, A[0], A[1], r, np.median(R[:, j])))
    out[lbl] = dict(slope=float(A[0]), intercept=float(A[1]), R=float(r), median=float(np.median(R[:, j])))
print("\n  d(cm)   n(woody)  n(r+3cm)  n(r+10cm)   count")
bins = [0.03, 0.07, 0.10, 0.13, 0.16, 0.20, 0.26, 0.40, 1.2]
out["binned"] = {}
for i in range(len(bins) - 1):
    k = (R[:, 0] >= bins[i]) & (R[:, 0] < bins[i + 1])
    if k.sum() >= 8:
        print("  %5.1f-%5.1f %9.0f %10.0f %11.0f %8d" % (bins[i] * 100, bins[i + 1] * 100,
              np.median(R[k, 1]), np.median(R[k, 2]), np.median(R[k, 3]), k.sum()))
        out["binned"][f"{bins[i]}_{bins[i+1]}"] = dict(n=int(k.sum()), woody=float(np.median(R[k, 1])),
                                                      r3=float(np.median(R[k, 2])), r10=float(np.median(R[k, 3])))
A = np.polyfit(R[:, 0], R[:, 2], 1); A1 = np.polyfit(R[:, 0], R[:, 1], 1)
rho = A[0] / (np.pi * 0.30)
print("\ncalibration (r=d/2+3cm): n = %.1f*d %+.1f  -> rho_s = %.0f pts/m2" % (A[0], A[1], rho))
out["rho_s"] = float(rho)
out["budget"] = {}
print("\n  DBH     points on a 0.30 m segment   [both calibrations]")
for d in (0.10, 0.15, 0.20, 0.30, 0.50, 0.80):
    a = max(A[0] * d + A[1], 0); b = max(A1[0] * d + A1[1], 0)
    out["budget"][str(d)] = dict(all_pts=float(a), woody_pts=float(b))
    print("  %5.0f cm %14.1f %20.1f" % (d * 100, a, b))
json.dump(out, open(f"{F}/stem_budget_tilelevel.json", "w"), indent=1)
print("\nwrote", f"{F}/stem_budget_tilelevel.json")