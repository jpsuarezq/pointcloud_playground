"""What do the 3 values of PredSemantic_FM mean?

Evidence assembled:
  (a) HAG distribution per class (against the SMRF ground)
  (b) RGB / greenness per class
  (c) local PCA geometry (linearity, planarity, verticality) per class and
      per HAG band
  (d) can ANY threshold or combination of FM isolate woody stems at the
      1.3 m DBH band?
"""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, load_rgb, SCALE, OX, OY
from scipy.spatial import cKDTree

G = "/tmp/opencode/cloud-research/out/grids"; O = "/tmp/opencode/cloud-research/out/dem"
HC = "/tmp/opencode/cloud-research/out/hag"; F = "/tmp/opencode/cloud-research/out/fig"
m5 = dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5, X5, Y5, NX5, NY5 = int(m5["res"]), int(m5["x0"]), int(m5["y0"]), int(m5["nx"]), int(m5["ny"])
cnt5 = np.load(f"{G}/cnt_5000.npy"); slope5 = np.load(f"{O}/slope05.npy")
usable = (cnt5 >= 20) & np.isfinite(slope5) & (slope5 < 60)
ncell = int(usable.sum()); area = ncell * 0.25

X = load_xyz()
xr = np.asarray(X["x"], dtype=np.int64) * SCALE + OX
yr = np.asarray(X["y"], dtype=np.int64) * SCALE + OY
zr = np.asarray(X["z"], dtype=np.float64) * SCALE
xraw = np.asarray(X["x"], dtype=np.int64); yraw = np.asarray(X["y"], dtype=np.int64)
del X
fm = np.asarray(load_attr()["fm"])
C = load_rgb()
R = np.asarray(C["r"]).astype(np.float32); G_ = np.asarray(C["g"]).astype(np.float32); B = np.asarray(C["b"]).astype(np.float32)
del C
ix5 = np.clip(((xraw - X5) // R5).astype(np.int32), 0, NX5 - 1)
iy5 = np.clip(((yraw - Y5) // R5).astype(np.int32), 0, NY5 - 1)
base = usable[iy5, ix5]
hag = np.load(f"{HC}/hag_smrf025.npy"); okh = np.load(f"{HC}/ok_smrf025.npy")
sel_all = base & okh
print("valid", sel_all.sum())
rep = {}

# ---------------- (a) HAG distribution ----------------
edges = np.round(np.arange(0, 40.05, 0.25), 3)
print("\n=== (a) HAG distribution per FM class (% of that class, and % of the band) ===")
print("HAG bin        FM0%   FM1%   FM2%   n_FM0    n_FM1    n_FM2   | share of band")
prof = {}
for c in (0, 1, 2):
    h = hag[sel_all & (fm == c)]
    b = np.clip((h / 0.25).astype(np.int32), 0, len(edges) - 2)
    prof[c] = np.bincount(b, minlength=len(edges) - 1)
tot = prof[0] + prof[1] + prof[2]
for i in range(0, 40):
    lo, hi = edges[i], edges[i + 1]
    row = [100 * prof[c][i] / max(prof[c].sum(), 1) for c in (0, 1, 2)]
    sh = [100 * prof[c][i] / max(tot[i], 1) for c in (0, 1, 2)]
    if lo > 16: break
    print(f"{lo:5.2f}-{hi:5.2f} {row[0]:7.2f} {row[1]:7.2f} {row[2]:7.2f} {prof[0][i]:9d} {prof[1][i]:8d} {prof[2][i]:8d}  | {sh[0]:5.1f} {sh[1]:5.1f} {sh[2]:5.1f}")
rep["hag_profile"] = {str(c): prof[c].tolist() for c in (0, 1, 2)}
print("\nquantiles of HAG per class:")
for c in (0, 1, 2):
    h = hag[sel_all & (fm == c)]
    print(f"  FM{c}: n={len(h):9d}  pct " + " ".join(f"{q}%={np.percentile(h,q):7.2f}" for q in (1, 5, 25, 50, 75, 95, 99)))
    rep[f"hag_pct_fm{c}"] = [float(np.percentile(hag[sel_all & (fm == c)], q)) for q in (1, 5, 25, 50, 75, 95, 99)]
print("frac of class below HAG thresholds:")
for c in (0, 1, 2):
    h = hag[sel_all & (fm == c)]
    print(f"  FM{c}: " + "  ".join(f"<{t}m {100*np.mean(h<t):5.2f}%" for t in (0.5, 1.0, 1.3, 2.0, 5.0, 15.0)))
    rep[f"hag_frac_below_fm{c}"] = {str(t): float(np.mean(h < t)) for t in (0.5, 1.0, 1.3, 2.0, 5.0, 15.0)}

# ---------------- (b) RGB ----------------
s = sel_all
mx = R > 0
print("\n=== (b) colour ===")
print("class   n     meanR  meanG  meanB   exG(G-R)/(G+R)  frac  G==R==B")
for c in (0, 1, 2):
    m = s & (fm == c) & mx
    r, g_, b = R[m], G_[m], B[m]
    print(f"FM{c} {m.sum():9d} {r.mean():7.1f} {g_.mean():7.1f} {b.mean():7.1f} "
          f"{np.mean((g_-r)/(g_+r+1e-6)):12.3f}  {np.mean((g_==r)&(r==b))*100:6.1f}%")
    rep[f"rgb_fm{c}"] = dict(n=int(m.sum()), R=float(r.mean()), G=float(g_.mean()), B=float(b.mean()),
                             exg=float(np.mean((g_ - r) / (g_ + r + 1e-6))),
                             grey_frac=float(np.mean((g_ == r) & (r == b))))
band = s & (hag >= 1.15) & (hag <= 1.45)
for c in (0, 1, 2):
    m = band & (fm == c) & mx
    r, g_, b = R[m], G_[m], B[m]
    print(f"1.3m FM{c} {m.sum():9d} {r.mean():7.1f} {g_.mean():7.1f} {b.mean():7.1f} "
          f"{np.mean((g_-r)/(g_+r+1e-6)):12.3f}  {np.mean((g_==r)&(r==b))*100:6.1f}%")
    rep[f"rgb_1m3_fm{c}"] = dict(n=int(m.sum()), R=float(r.mean()), G=float(g_.mean()), B=float(b.mean()),
                                 exg=float(np.mean((g_ - r) / (g_ + r + 1e-6))))

# ---------------- (c)+(d) local PCA geometry ----------------
print("\n=== (c) local PCA geometry (k=16 NN) ===")
rng = np.random.default_rng(0)
idx_pool = np.where(sel_all)[0]
sub = rng.choice(idx_pool, 260_000, replace=False)
P = np.column_stack([xr[sub], yr[sub], zr[sub]])
tree = cKDTree(P)
K = 16
dist, nb = tree.query(P, k=K + 1, workers=-1)
Q = P[nb]                                  # (n,K+1,3)
Q -= Q.mean(1, keepdims=True)
Ccov = np.einsum("nki,nkj->nij", Q, Q) / K
ev, evec = np.linalg.eigh(Ccov)            # ascending
ev = np.clip(ev, 1e-12, None)
l1, l2, l3 = ev[:, 2], ev[:, 1], ev[:, 0]
lin = (l1 - l2) / l1
pla = (l2 - l3) / l1
sca = l3 / l1
vert = np.abs(evec[:, 2, 2])               # |z| component of the principal direction
rk = dist[:, K]                            # radius of the K-NN neighbourhood
np.save(f"{F}/pca_sub.npy", np.column_stack([sub.astype(np.int64), lin, pla, sca, vert, rk]))
print("subsample", len(sub))

def summ(mask, lbl):
    if mask.sum() < 50: return
    print(f"{lbl:22s} n={mask.sum():7d}  lin {np.mean(lin[mask]):.3f}  pla {np.mean(pla[mask]):.3f}  "
          f"vert {np.mean(vert[mask]):.3f}  woody(lin>.6&vert>.8) {100*np.mean((lin[mask]>.6)&(vert[mask]>.8)):5.2f}%  "
          f"leafy(lin<.25) {100*np.mean(lin[mask]<.25):5.2f}%")

fms = fm[sub]; hags = hag[sub]
print("-- whole cloud --")
for c in (0, 1, 2): summ(fms == c, f"FM{c}")
print("-- by HAG band, all classes --")
for lo, hi in [(0, 0.3), (0.3, 1.0), (1.0, 1.6), (1.6, 3.0), (3.0, 8.0), (8.0, 20.0), (20.0, 100)]:
    summ((hags >= lo) & (hags < hi), f"HAG {lo}-{hi}")
print("-- by HAG band x FM, at the DBH band 1.15-1.45 --")
bandm = (hags >= 1.15) & (hags <= 1.45)
for c in (0, 1, 2): summ(bandm & (fms == c), f"1.3m FM{c}")
print("-- by HAG band x FM, 0-0.3 m --")
for c in (0, 1, 2): summ((hags < 0.3) & (fms == c), f"0.3m FM{c}")

# ---------------- (d) separability of stems at 1.3 m ----------------
print("\n=== (d) can FM isolate stems at the 1.3 m DBH band? ===")
wood = (lin > 0.6) & (vert > 0.80)
print(f"in subsample, HAG 1.15-1.45: n={bandm.sum()}  woody fraction {100*wood[bandm].mean():.2f}%")
rep["dbh_band_woody_fraction"] = float(wood[bandm].mean())
for c in (0, 1, 2):
    m = bandm & (fms == c)
    if m.sum() > 20:
        print(f"  FM{c}: n={m.sum():7d}  woody {100*wood[m].mean():6.2f}%  "
              f"(= {100*wood[m].sum()/wood[bandm].sum():5.2f}% of all woody points in the band)")
        rep[f"dbh_fm{c}_woody"] = dict(n=int(m.sum()), woody_frac=float(wood[m].mean()))
# best achievable: use FM=2 as the only rule
m2 = bandm & (fms == 2)
tp = wood[m2].sum()
print(f"rule 'FM==2': keeps {100*wood[m2].sum()/wood[bandm].sum():.2f}% of woody points, "
      f"precision {100*wood[m2].mean():.2f}%  (baseline 'keep all': precision {100*wood[bandm].mean():.2f}%)")
print(f"rule 'FM==0': keeps {100*wood[bandm&(fms==0)].sum()/wood[bandm].sum():.2f}% of woody points, "
      f"precision {100*wood[bandm&(fms==0)].mean() if (bandm&(fms==0)).sum() else float('nan'):.2f}%")
print(f"rule 'FM==1': n={(bandm&(fms==1)).sum()}  woody {100*wood[bandm&(fms==1)].mean() if (bandm&(fms==1)).sum() else float('nan'):.2f}%")
# AUC of FM as a score for woodiness in the band
sc = fms[bandm].astype(np.float64); lb = wood[bandm].astype(int)
order = np.argsort(sc); ranks = np.empty(len(sc)); ranks[order] = np.arange(1, len(sc) + 1)
n1 = lb.sum(); n0 = len(lb) - n1
auc = (ranks[lb == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
print(f"AUC of FM value as a predictor of 'woody' inside the band: {auc:.4f} (0.5 = no information)")
rep["fm_auc_woody_in_band"] = float(auc)

json.dump(rep, open(f"{F}/fm_analysis.json", "w"), indent=1)
print("\nwrote", f"{F}/fm_analysis.json")