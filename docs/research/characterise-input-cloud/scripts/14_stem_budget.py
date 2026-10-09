"""Point budget on stems: measure n_points vs stem diameter directly.

Method: in a dense 60x60 m block, take the DBH band (1.15-1.45 m), label local
geometry by KNN PCA, keep the woody/vertical subset, single-linkage cluster it
at 10 cm; each cluster's horizontal spread gives its diameter; then count how
many band points fall inside a cylinder of that diameter (plus a 5 cm skin).
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE, OX, OY
from scipy.spatial import cKDTree

HC = "/tmp/opencode/cloud-research/out/hag"; F = "/tmp/opencode/cloud-research/out/fig"
G = "/tmp/opencode/cloud-research/out/grids"
m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5, NX5, NY5 = int(m5["res"]), int(m5["x0"]), int(m5["y0"]), int(m5["nx"]), int(m5["ny"])
cnt5 = np.load(f"{G}/cnt_5000.npy"); O = "/tmp/opencode/cloud-research/out/dem"
slope5 = np.load(f"{O}/slope05.npy")
usable = (cnt5 >= 20) & np.isfinite(slope5) & (slope5 < 60)

X = load_xyz()
xr = np.asarray(X["x"], dtype=np.int64) * SCALE + OX
yr = np.asarray(X["y"], dtype=np.int64) * SCALE + OY
zr = np.asarray(X["z"], dtype=np.float64) * SCALE
xraw = np.asarray(X["x"], dtype=np.int64); yraw = np.asarray(X["y"], dtype=np.int64)
del X
hag = np.load(f"{HC}/hag_smrf025.npy"); okh = np.load(f"{HC}/ok_smrf025.npy")
ix5 = np.clip(((xraw - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((yraw - Y5) // R5).astype(np.int32), 0, NY5 - 1)
valid = usable[iy5, ix5] & okh

blk = (xr >= 675360) & (xr < 675420) & (yr >= 758800) & (yr < 758860)
bi = np.where(blk)[0]
print("block points", len(bi))
P = np.column_stack([xr[bi], yr[bi], zr[bi]])
h = hag[bi]
tree = cKDTree(P)
K = 16
dd, nb = tree.query(P, k=K + 1, workers=-1)
Q = P[nb]; Q -= Q.mean(1, keepdims=True)
cov = np.einsum("nki,nkj->nij", Q, Q) / K
ev, evec = np.linalg.eigh(cov)
ev = np.clip(ev[:, ::-1], 1e-12, None)
lin = (ev[:, 0] - ev[:, 1]) / ev[:, 0]
pla = (ev[:, 1] - ev[:, 2]) / ev[:, 0]
vert = np.abs(evec[:, :, 2][:, 2])
rho_all = K / (np.pi * dd[:, K] ** 2)
print("block surface density: median %.0f pts/m2 (all points)" % np.median(rho_all))

band = (h >= 1.15) & (h <= 1.45) & valid[bi]
print("band points in block", band.sum())
for lname, m in [("lin>.6 & vert>.8", (lin > .6) & (vert > .8)),
                 ("lin>.75 & vert>.9", (lin > .75) & (vert > .9))]:
    mm = band & m
    # woody-only neighbourhood density
    idx = np.where(mm)[0]
    print(f"  {lname:20s} {mm.sum():6d} pts ({100*mm.sum()/band.sum():5.2f}% of band)  "
          f"median r_k {np.median(dd[mm,K]):.3f}")
    if len(idx) > 40:
        t2 = cKDTree(P[idx])
        d2, _ = t2.query(P[idx], k=9, workers=-1)
        print(f"      woody-only surface density (k=8): {8/(np.pi*np.median(d2[:,8])**2):.0f} pts/m2 "
              f"(median woody r8 {np.median(d2[:,8]):.3f} m)")

m = band & (lin > .6) & (vert > .8)
idx = np.where(m)[0]
W = np.column_stack([xr[bi][idx], yr[bi][idx]])
pairs = cKDTree(W).query_pairs(0.10, output_type="ndarray")
par = np.arange(len(W))
def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for a, b in pairs:
    ra, rb = find(a), find(b)
    if ra != rb: par[ra] = rb
roots = np.array([find(i) for i in range(len(W))])
uq, lab = np.unique(roots, return_inverse=True)
sizes = np.bincount(lab)
cent = np.array([W[lab == u].mean(0) for u in range(len(uq))])
keep = sizes >= 6
cent, csz = cent[keep], sizes[keep]
print(f"\nwoody clusters (link 0.10 m, >=6 pts): {len(cent)}")

# band points around each cluster
BP = np.column_stack([xr[bi][band], yr[bi][band]])
t3 = cKDTree(BP)
ds, ns = [], []
for c, s, l in zip(cent, csz, np.where(keep)[0]):
    own = W[lab == l]
    dq = np.hypot(*(own - c).T)
    d_est = 2 * np.percentile(dq, 75)          # diameter from the woody points' own spread
    if d_est < 0.03 or d_est > 1.20: continue
    dq = np.hypot(*(BP - c).T)                 # all band points around the same axis
    n = int(np.sum(dq <= d_est / 2 + 0.05))
    ds.append(d_est); ns.append(n)
ds = np.array(ds); ns = np.array(ns)
print(f"usable stem candidates: {len(ds)}  diameter pct "
      f"{np.percentile(ds,[10,25,50,75,90]).round(2)} m")
A = np.polyfit(ds, ns, 1)
print(f"fit  n_points = {A[0]:.1f} * d(m) {A[1]:+.1f}   (R={np.corrcoef(ds,ns)[0,1]:.3f})")
rho_s = A[0] / (np.pi * 0.30)
print(f"implied woody surface density rho_s = {rho_s:.0f} pts/m2")

print("\n  DBH     n points on a 0.30 m segment  (fit)   (from rho_s)")
bud = {}
for d in (0.10, 0.15, 0.20, 0.30, 0.50, 0.80):
    n = A[0] * d + A[1]; n2 = rho_s * np.pi * d * 0.30
    bud[str(d)] = dict(fit=float(n), rho_s=float(n2))
    print(f"  {d*100:5.0f} cm {n:14.1f} {n2:20.1f}")
# binned empirical medians
print("\n  empirical medians by measured diameter bin:")
edges = [0.02, 0.06, 0.10, 0.15, 0.20, 0.30, 0.45, 0.70, 1.20]
for i in range(len(edges) - 1):
    k = (ds >= edges[i]) & (ds < edges[i + 1])
    if k.sum() >= 8:
        print(f"   d {edges[i]*100:5.0f}-{edges[i+1]*100:5.0f} cm  n={k.sum():5d}  "
              f"median points {np.median(ns[k]):6.1f}  (mean {ns[k].mean():6.1f})")
        bud[f"bin_{edges[i]}_{edges[i+1]}"] = dict(n=int(k.sum()), median=float(np.median(ns[k])),
                                                   mean=float(ns[k].mean()))
json.dump(dict(d=ds.tolist(), n=ns.tolist(), fit=list(map(float, A)),
               rho_s=float(rho_s), budget=bud,
               rho_surface_all=float(np.median(rho_all)),
               block_band_points=int(band.sum()),
               block_woody=int(m.sum())),
          open(f"{F}/stem_point_budget.json", "w"), indent=1)
print("\nwrote", f"{F}/stem_point_budget.json")