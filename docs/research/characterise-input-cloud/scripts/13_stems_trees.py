"""Stem point budget + tree count.

(1) Point budget. For each vertical stem found in the 1.15-1.45 m DBH band we
    build the horizontal radial profile of point counts around the stem axis.
    For a cylinder the count profile rises until r = d/2 and then plateaus at
    rho_s * pi * d * L. The knee therefore measures the stem diameter and the
    plateau measures the number of points on a 0.30 m tall segment. We fit
    n_points(d) and report the budget for 30/50/80 cm DBH.
    Independent sanity check: rho_s from a dense-block KNN estimate.
(2) Tree count from CHM local maxima (several window sizes), crown connected
    components, and stem clusters.
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE, OX, OY
from scipy import ndimage
from scipy.spatial import cKDTree

G = "/tmp/opencode/cloud-research/out/grids"; O = "/tmp/opencode/cloud-research/out/dem"
HC = "/tmp/opencode/cloud-research/out/hag"; F = "/tmp/opencode/cloud-research/out/fig"
m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5, NX5, NY5 = int(m5["res"]), int(m5["x0"]), int(m5["y0"]), int(m5["nx"]), int(m5["ny"])
cnt5 = np.load(f"{G}/cnt_5000.npy"); slope5 = np.load(f"{O}/slope05.npy")
usable = (cnt5 >= 20) & np.isfinite(slope5) & (slope5 < 60)
area = int(usable.sum()) * 0.25

X = load_xyz()
xr = np.asarray(X["x"], dtype=np.int64) * SCALE + OX
yr = np.asarray(X["y"], dtype=np.int64) * SCALE + OY
zr = np.asarray(X["z"], dtype=np.float64) * SCALE
xraw = np.asarray(X["x"], dtype=np.int64); yraw = np.asarray(X["y"], dtype=np.int64)
fm = np.asarray(load_attr()["fm"]); del X
hag = np.load(f"{HC}/hag_smrf025.npy"); okh = np.load(f"{HC}/ok_smrf025.npy")
ix5 = np.clip(((xraw - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((yraw - Y5) // R5).astype(np.int32), 0, NY5 - 1)
valid = usable[iy5, ix5] & okh
rep = {"area_m2": area, "points_valid": int(valid.sum())}
print("valid points", int(valid.sum()), "area", area, "density",
      round(int(valid.sum()) / area, 1), "pts/m2")

# ---------------- woody fraction in the band, from a DENSE block ----------------
print("\n=== dense-block geometry check ===")
blk = (xr >= 675370) & (xr < 675400) & (yr >= 758815) & (yr < 758845)
print("block points", blk.sum())
bp = np.column_stack([xr[blk], yr[blk], zr[blk]])
t = cKDTree(bp)
K = 16
dd, nb = t.query(bp, k=K + 1, workers=-1)
Q = bp[nb]; Q -= Q.mean(1, keepdims=True)
ev = np.linalg.eigh(np.einsum("nki,nkj->nij", Q, Q) / K)[0]
ev = np.clip(ev[:, ::-1], 1e-12, None)
lin = (ev[:, 0] - ev[:, 1]) / ev[:, 0]
bh = hag[blk]
bandb = (bh >= 1.15) & (bh <= 1.45)
# verticality of the principal direction
_, evec = np.linalg.eigh(np.einsum("nki,nkj->nij", Q, Q) / K)
vert = np.abs(evec[:, :, 2][:, 2])
print("block: n=%d  band n=%d  median r_k=%.3f m  rho_surface(all)=%.0f pts/m2"
      % (len(bp), bandb.sum(), np.median(dd[:, K]), K / (np.pi * np.median(dd[:, K]) ** 2)))
print("  band: median r_k=%.3f m  rho_surface=%.0f pts/m2  linear frac>0.6 = %.1f%%  "
      "vertical frac>0.8 = %.1f%%  woody = %.1f%%"
      % (np.median(dd[bandb, K]), K / (np.pi * np.median(dd[bandb, K]) ** 2),
         100 * np.mean(lin[bandb] > .6), 100 * np.mean(vert[bandb] > .8),
         100 * np.mean((lin[bandb] > .6) & (vert[bandb] > .8))))
rep["block"] = dict(n=int(len(bp)), band_n=int(bandb.sum()),
                    rk_median=float(np.median(dd[bandb, K])),
                    rho_surface_band=float(K / (np.pi * np.median(dd[bandb, K]) ** 2)),
                    woody_frac=float(np.mean((lin[bandb] > .6) & (vert[bandb] > .8))))

# ---------------- cluster stems in the DBH band ----------------
print("\n=== stem clusters in the 1.15-1.45 m band ===")
sel = np.where(valid & (hag >= 1.15) & (hag <= 1.45))[0]
print("band points:", len(sel))
P = np.column_stack([xr[sel], yr[sel], zr[sel]])
tree = cKDTree(P)
pairs = tree.query_pairs(0.12, output_type="ndarray")
n = len(P); par = np.arange(n)
def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for a, b in pairs:
    ra, rb = find(a), find(b)
    if ra != rb: par[ra] = rb
roots = np.array([find(i) for i in range(n)])
uniq, lab = np.unique(roots, return_inverse=True)
sizes = np.bincount(lab)
cents = np.array([P[lab == u].mean(0) for u in range(len(uniq))])
keep = sizes >= 8
cents, csize = cents[keep], sizes[keep]
print(f"clusters (link r=0.12 m, >=8 pts): {len(cents)}; size pct "
      f"{np.percentile(csize,[10,25,50,75,90]).round(0)}")

# radial profile around each cluster centroid
qsel = np.where(valid & (hag >= 1.15) & (hag <= 1.45))[0]
Q = np.column_stack([xr[qsel], yr[qsel]])
t2 = cKDTree(Q)
RAD = np.arange(0.02, 0.42, 0.02)
prof = np.zeros((len(cents), len(RAD)))
for i, c in enumerate(cents):
    nb = t2.query_ball_point(c[:2], RAD[-1])
    if not nb: continue
    nb = np.asarray(nb)
    d = np.hypot(Q[nb, 0] - c[0], Q[nb, 1] - c[1])
    h, _ = np.histogram(d, bins=np.concatenate([[0], RAD]))
    prof[i] = np.cumsum(h)
mean_prof = prof.mean(0)
print("mean cumulative radial profile (counts):")
for r, v in zip(RAD, mean_prof): print(f"  r={r:.2f}  {v:8.1f}")

# knee = smallest r where the profile reaches 95% of its plateau (r>=0.26 mean)
knee = RAD[np.argmax(mean_prof >= 0.9 * mean_prof[-1])]
dbar = 2 * knee
nbar = float(mean_prof[-1])
print(f"\nknee r={knee:.2f} m -> representative DBH {dbar*100:.0f} cm")
print(f"plateau count (points on a {RAD[-1]:.2f} m radius / 0.30 m tall cylinder) = {nbar:.1f}")
# linear fit of plateau count vs knee diameter across clusters
knees = np.array([RAD[np.argmax(prof[i] >= 0.9 * prof[i][-1])] for i in range(len(cents))])
plats = prof[:, -1]
ok = (plats >= 3) & (knees > 0)
A = np.polyfit(knees[ok] * 2, plats[ok], 1)
print(f"per-stem fit: n_points = {A[0]:.2f} * DBH(m) {A[1]:+.2f}   "
      f"(median DBH {2*np.median(knees)*100:.0f} cm, median plateau {np.median(plats):.0f} pts)")
rho_s = A[0] / (np.pi * 0.30)
print(f"implied woody surface density rho_s = {rho_s:.0f} pts/m2")

budget = {}
print("\n  DBH     n points on a 0.30 m segment")
for d in (0.10, 0.15, 0.20, 0.30, 0.50, 0.80):
    npts = A[0] * d + A[1]
    budget[str(d)] = float(npts)
    print(f"  {d*100:5.0f} cm {npts:10.1f}")
rep["point_budget"] = budget
rep["stem_fit"] = dict(slope=float(A[0]), intercept=float(A[1]), rho_s=float(rho_s),
                       median_db_cm=float(2 * np.median(knees) * 100),
                       n_clusters=int(len(cents)), knee_mean_m=float(knee))
rep["radial_profile"] = dict(r=RAD.tolist(), mean=mean_prof.tolist())

# ---------------- tree count ----------------
print("\n=== tree count ===")
chm = np.full((NY5, NX5), -np.inf)
ixc = np.clip(((xraw - X5) // R5).astype(np.int64), 0, NX5 - 1)
iyc = np.clip(((yraw - Y5) // R5).astype(np.int64), 0, NY5 - 1)
m = valid & (hag < 60)
np.maximum.at(chm, (iyc[m], ixc[m]), hag[m])
chm[~np.isfinite(chm)] = np.nan
np.save(f"{O}/chm05.npy", chm)
print("CHM max %.1f m, p50 %.1f, p90 %.1f" % (np.nanmax(chm), np.nanpercentile(chm, 50),
                                             np.nanpercentile(chm, 90)))
for thr in (3, 5, 8):
    for rad in (3, 4, 5, 6, 8):
        chs = ndimage.generic_filter(chm, np.nanmax, size=2 * rad + 1, mode="constant", cval=np.nan)
        ismax = np.isfinite(chm) & (chm >= chs - 1e-9) & (chm > thr)
        _, nl = ndimage.label(ismax)
        print(f"  CHM> {thr} m, window {rad} m: {nl:5d} crowns -> {nl/(area/1e4):6.0f} trees/ha")
        rep.setdefault("chm", {}).setdefault(str(thr), {})[str(rad)] = int(nl)
for thr in (1.6, 2.0):
    selc = valid & (hag >= thr)
    occ = np.zeros((NY5, NX5), bool); occ[iyc[selc], ixc[selc]] = True
    occ = ndimage.binary_closing(occ, np.ones((3, 3)))
    lab2, nl = ndimage.label(occ)
    sizes2 = np.bincount(lab2.ravel())[1:]
    for minsz in (4, 6, 10):
        print(f"  crown CC above {thr} m, >= {minsz} cells: {int((sizes2>=minsz).sum()):5d} "
              f"-> {(sizes2>=minsz).sum()/(area/1e4):6.0f} trees/ha")
        rep.setdefault("crown_cc", {}).setdefault(str(thr), {})[str(minsz)] = int((sizes2 >= minsz).sum())
json.dump(rep, open(f"{F}/stem_and_trees.json", "w"), indent=1)
print("\nwrote", f"{F}/stem_and_trees.json")