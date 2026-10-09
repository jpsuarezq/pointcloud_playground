import json, os, numpy as np, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage

F = "/tmp/opencode/cloud-research/out/fig"
D = "/tmp/opencode/cloud-research/out/dem"
C = json.load(open(f"{F}/surface_comparison.json"))
edges = np.array(C["edges"])
plt.rcParams.update({"font.size": 9, "figure.dpi": 130})

# ---- Fig 1: vertical density profile per ground surface
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
order = ["smrf025", "csf050", "dem_fm0025", "dem_open90", "dem_gmin025"]
lbl = {"smrf025": "PDAL SMRF, 0.25 m", "csf050": "PDAL CSF, 0.5 m",
       "dem_fm0025": "grid-min of PredSemantic_FM=0, 0.25 m",
       "dem_open90": "numpy slope-corrected opening, 0.5 m",
       "dem_gmin025": "numpy grid minimum, 0.25 m (no morphology)"}
col = {"smrf025": "#d62728", "csf050": "#1f77b4", "dem_fm0025": "#2ca02c",
       "dem_open90": "#ff7f0e", "dem_gmin025": "#7f7f7f"}
for k in order:
    if k not in C["surfaces"]: continue
    p = np.array(C["surfaces"][k]["prof"])
    ax[0].plot(edges[:-1], p, lw=1.4, label=lbl[k], color=col[k])
below = {k: C["surfaces"][k]["frac_below0"] for k in order if k in C["surfaces"]}
for k in order:
    if k not in C["surfaces"]: continue
    p = np.array(C["surfaces"][k]["prof"]); p = p / max(p.max(), 1)
    ax[1].plot(edges[:-1], p, lw=1.4, label=lbl[k], color=col[k])
ax[0].set_xlim(0, 6); ax[0].set_xlabel("height above ground (m)")
ax[0].set_ylabel("points / m$^2$ / m"); ax[0].set_title("Vertical density vs ground surface")
ax[0].legend(fontsize=6.5, loc="upper right")
ax[1].set_xlim(0, 6); ax[1].set_ylim(0, 1.05)
ax[1].axvspan(1.0, 1.6, color="orange", alpha=.18)
ax[1].set_xlabel("height above ground (m)"); ax[1].set_ylabel("normalised density")
ax[1].set_title("Normalised; 1.0-1.6 m blind band shaded")
ax[1].text(1.3, 0.15, "blind band", ha="center", fontsize=7, color="#a05000")
for a in ax: a.grid(alpha=.25)
fig.suptitle("Figure 1 - the 1.0-1.6 m trough appears only with lower-envelope ground surfaces", y=1.02)
fig.tight_layout(); fig.savefig(f"{F}/fig1_density_profile.png", bbox_inches="tight")

# ---- Fig 2: PredSemantic_FM composition vs HAG
fm = json.load(open(f"{F}/fm_analysis.json"))
prof = np.array([fm["hag_profile"][str(c)] for c in (0, 1, 2)], dtype=float)
share = 100 * prof / np.maximum(prof.sum(0), 1)
e2 = np.arange(len(share[0])) * 0.25
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
for i, (nm, c_) in enumerate([("FM0", "#2ca02c"), ("FM1", "#ff7f0e"), ("FM2", "#1f77b4")]):
    ax[0].plot(e2, np.maximum(share[i], 1e-3), lw=1.5, color=c_, label=nm)
ax[0].set_xlim(0, 20); ax[0].set_ylim(1e-2, 100); ax[0].set_yscale("log")
ax[0].set_yticks([0.1, 1, 10, 100]); ax[0].set_yticklabels(["0.1", "1", "10", "100"])
ax[0].set_xlabel("height above ground (m)"); ax[0].set_ylabel("% of points in HAG band")
ax[0].set_title("PredSemantic_FM composition vs height (share of band)")
ax[0].legend(loc="center left", fontsize=8)
for c, c_ in zip((0, 1, 2), ("#2ca02c", "#ff7f0e", "#1f77b4")):
    v = fm[f"hag_pct_fm{c}"]
    ax[1].plot([1, 5, 25, 50, 75, 95, 99], v, "o-", color=c_, label=f"FM{c}")
ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set_xticks([1, 5, 25, 50, 95]); ax[1].set_xticklabels(["1", "5", "25", "50", "95"])
ax[1].set_xlabel("percentile of height above ground (m)")
ax[1].set_ylabel("height above ground (m)")
ax[1].set_title("HAG quantiles per FM class"); ax[1].legend(fontsize=8)
for a in ax: a.grid(alpha=.25)
fig.suptitle("Figure 2 - what the three PredSemantic_FM values mean", y=1.03)
fig.tight_layout(); fig.savefig(f"{F}/fig2_fm_semantics.png", bbox_inches="tight")

# ---- Fig 3: trough sensitivity to vertical offset of the ground surface
T = json.load(open(f"{F}/trough_sensitivity.json"))
rows = [r for r in T["rows"] if r["surface"].startswith("smrf")]
fig, ax = plt.subplots(figsize=(7, 3.4))
for r in rows:
    p = np.array(r["prof"]); e = np.array(T["edges"])
    ax.plot(e[:-1], p / max(p.max(), 1), lw=1.2, label=f'{r["surface"]}')
ax.axvspan(0.8, 1.6, color="orange", alpha=.18)
ax.set_xlim(0, 4); ax.set_ylim(0, 1.05)
ax.set_xlabel("height above ground (m)"); ax.set_ylabel("normalised density")
ax.set_title("Figure 3 - trough depth is invariant to a +-0.25 m offset of the ground surface")
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=.25)
fig.tight_layout(); fig.savefig(f"{F}/fig3_trough_offset.png", bbox_inches="tight")

# ---- Fig 4: tree count vs CHM window
T2 = json.load(open(f"{F}/stem_and_trees.json"))
fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
for thr, d in T2.get("chm", {}).items():
    w = sorted(int(k) for k in d)
    ax[0].plot(w, [d[str(k)] for k in w], "o-", label=f"CHM > {thr} m")
ax[0].axhline(500, ls="--", c="k", lw=.8); ax[0].text(6.1, 510, "500", fontsize=7)
ax[0].set_xlabel("CHM local-maximum window (m)"); ax[0].set_ylabel("crowns detected")
ax[0].set_title("Tree count vs window size"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.25)
p = np.array(C["surfaces"]["smrf025"]["prof"]); e = np.array(C["edges"])
ax[1].plot(e[:-1], p, lw=1.4, color="#d62728")
ax[1].axvspan(1.0, 1.6, color="orange", alpha=.25)
ax[1].axvline(1.3, ls=":", c="k", lw=.8)
ax[1].text(1.33, p.max() * .6, "DBH band\n1.3 m", fontsize=7)
ax[1].axvspan(0, 0.5, color="grey", alpha=.15)
ax[1].text(0.05, p.max() * .8, "ground\nsurface\npoints", fontsize=7)
ax[1].set_xlim(0, 6); ax[1].set_xlabel("height above ground (m)")
ax[1].set_ylabel("points / m$^2$ / m")
ax[1].set_title("Where the points are (SMRF ground)"); ax[1].grid(alpha=.25)
fig.tight_layout(); fig.savefig(f"{F}/fig4_trees_and_profile.png", bbox_inches="tight")

# ---- Fig 5: DEM comparison map (hillsshade of surfaces)
names = [("smrf025", "PDAL SMRF 0.25 m"), ("csf050", "PDAL CSF 0.5 m"),
         ("dem_gmin025", "numpy grid-min 0.25 m"), ("dem_fm0025", "numpy grid-min FM0 0.25 m")]
fig, ax = plt.subplots(1, 4, figsize=(13, 4.2))
for a, (k, t) in zip(ax, names):
    if k.startswith("tif"):
        arr = np.load(f"{D}/{k}.npy"); arr = np.where(arr < -9000, np.nan, arr)
        dx, dy = json.load(open(f"{D}/{k}.json"))["res"]
    else:
        arr = np.load(f"{D}/{k}.npy"); dx = dy = 0.25 if "025" in k else 0.5
    A = arr if arr.shape[0] < arr.shape[1] else arr
    gy, gx = np.gradient(np.where(np.isfinite(arr), arr, np.nanmedian(arr)), dy, dx)
    hs = np.arctan(np.hypot(gx, gy))
    im = a.imshow(hs, cmap="terrain", vmin=0, vmax=np.degrees(np.radians(60)))
    a.set_title(f"{t}\n% pts below DEM = {100*C['surfaces'][k]['frac_below0']:.1f}%", fontsize=8)
    a.set_xticks([]); a.set_yticks([])
    fig.colorbar(im, ax=a, fraction=.046, label="slope (deg)")
fig.suptitle("Figure 5 - hillshade of the candidate ground surfaces", y=1.0)
fig.tight_layout(); fig.savefig(f"{F}/fig5_dem_hillshade.png", bbox_inches="tight")
print("figures written")
for f in sorted(os.listdir(F)):
    if f.endswith(".png"): print(" ", f, os.path.getsize(f"{F}/{f}"))