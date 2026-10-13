"""Figures + metrics for the extended (multi-candidate) run over the whole canopy."""
import json, csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

O = "/tmp/opencode/dbh/out"
FIG = "/home/usuario/workspaces/NubePuntos/prototype/stem-isolation-dbh/figures"

rows = json.load(open(f"{O}/results_multi.json"))
tab = list(csv.DictReader(open(f"{O}/tabla_arboles.csv")))
by = {r["tree"]: r for r in rows}
conf = {int(t["id_arbol"]): ("refused" if t["confianza"] not in ("high", "medium", "low") else t["confianza"]) for t in tab}
lat = {int(t["id_arbol"]): float(t["latitud"]) for t in tab}
lon = {int(t["id_arbol"]): float(t["longitud"]) for t in tab}

meas = [r for r in rows if r["status"] == "measured"]
COL = {"high": "#2ca02c", "medium": "#ff7f0e", "low": "#d62728", "refused": "#bbbbbb"}

# ---- fig1: map ----
fig, ax = plt.subplots(figsize=(9, 8))
for c in ("refused", "low", "medium", "high"):
    xs = [lon[r["tree"]] for r in rows if conf[r["tree"]] == c]
    ys = [lat[r["tree"]] for r in rows if conf[r["tree"]] == c]
    s = [ (by[r["tree"]].get("d_ring_cm") or 10) for r in rows if conf[r["tree"]] == c]
    ax.scatter(xs, ys, c=COL[c], s=[max(8, v) for v in s], alpha=0.8, edgecolors="none",
               label=f"{c} (n={len(xs)})")
ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
ax.set_title("Canopy trees by confidence (marker size ∝ DBH)")
ax.legend(); fig.tight_layout(); fig.savefig(f"{FIG}/ext_fig1_map.png", dpi=150); plt.close(fig)

# ---- fig2: DBH distribution by confidence ----
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
bins = np.arange(0, 120, 5)
for c in ("high", "medium", "low"):
    d = [r["d_ring_cm"] for r in meas if r["confidence"] == c]
    axes[0].hist(d, bins=bins, histtype="step", lw=2, color=COL[c], label=f"{c} (n={len(d)})")
axes[0].set_xlabel("DBH (cm)"); axes[0].set_ylabel("trees"); axes[0].legend()
axes[0].set_title("DBH distribution by confidence")
hp = [r["top_m"] for r in meas]; dp = [r["d_ring_cm"] for r in meas]
cc = [COL[r["confidence"]] for r in meas]
axes[1].scatter(dp, hp, c=cc, s=18, alpha=0.8)
axes[1].set_xlabel("DBH (cm)"); axes[1].set_ylabel("total height (m)")
axes[1].set_title("Height vs DBH (coloured by confidence)")
fig.tight_layout(); fig.savefig(f"{FIG}/ext_fig2_dbh.png", dpi=150); plt.close(fig)

# ---- fig3: arc vs DBH ----
fig, ax = plt.subplots(figsize=(8, 6))
for c in ("high", "medium", "low"):
    a = [r["qc"]["arc"] for r in meas if r["confidence"] == c]
    d = [r["d_ring_cm"] for r in meas if r["confidence"] == c]
    ax.scatter(a, d, c=COL[c], s=20, alpha=0.8, label=c)
ax.axvline(0.35, color="k", ls="--", lw=1, label="accept arc 0.35")
ax.axhline(100, color="k", ls=":", lw=1, label="accept d < 1 m")
ax.set_xlabel("arc coverage of the accepted ring"); ax.set_ylabel("DBH (cm)")
ax.set_title("QC: arc coverage vs DBH"); ax.legend()
fig.tight_layout(); fig.savefig(f"{FIG}/ext_fig3_arc_dbh.png", dpi=150); plt.close(fig)

# ---- fig4: outcomes / confidence / seed / slab ----
fig, axes = plt.subplots(1, 4, figsize=(18, 4.2))
st = Counter(r["status"] for r in rows)
axes[0].barh(list(st.keys()), list(st.values()), color="#7f7f7f"); axes[0].set_title(f"Outcome, n={len(rows)}")
cf = Counter(r["confidence"] for r in meas)
axes[1].barh(list(cf.keys()), list(cf.values()), color=[COL[k] for k in cf]); axes[1].set_title("Confidence (measured)")
sd = Counter(r.get("seed_kind") for r in meas)
axes[2].barh(list(sd.keys()), list(sd.values()), color="#1f77b4"); axes[2].set_title("Seed that isolated the axis")
sl = Counter(r.get("slab_half") for r in meas)
axes[3].bar([f"±{k}" for k in sl], list(sl.values()), color="#9467bd"); axes[3].set_title("Slab half-width used (m)")
fig.tight_layout(); fig.savefig(f"{FIG}/ext_fig4_outcomes.png", dpi=150); plt.close(fig)

# ---- fig5: band points vs measurability ----
fig, ax = plt.subplots(figsize=(9, 5))
bins = np.arange(0, 200, 5)
ax.hist([r["n_band"] for r in meas], bins=bins, color="#2ca02c", alpha=0.7, label=f"measured (n={len(meas)})")
ax.hist([r["n_band"] for r in rows if r["status"] == "refused"], bins=bins, color="#bbbbbb",
        alpha=0.7, label=f"refused (n={len(rows)-len(meas)})")
ax.set_yscale("log"); ax.set_xlabel("points in the 1.3 m band"); ax.set_ylabel("trees (log)")
ax.set_title("Band point budget vs measurability"); ax.legend()
fig.tight_layout(); fig.savefig(f"{FIG}/ext_fig5_band_budget.png", dpi=150); plt.close(fig)

# ---- metrics ----
def pct(x): return np.percentile(x, [10, 25, 50, 75, 90]).round(1).tolist()
metrics = dict(
    canopy_trees=len(rows), measured=len(meas), refused=len(rows) - len(meas),
    measured_frac=round(len(meas) / len(rows), 3),
    confidence=dict(Counter(r["confidence"] for r in meas)),
    seed=dict(Counter(r.get("seed_kind") for r in meas)),
    slab=dict(Counter(r.get("slab_half") for r in meas)),
    refused_with_band0=sum(1 for r in rows if r["status"] == "refused" and r.get("n_band", 0) == 0),
    refused_with_band_ge10=sum(1 for r in rows if r["status"] == "refused" and r.get("n_band", 0) >= 10),
    dbh_cm=dict(all=pct([r["d_ring_cm"] for r in meas]),
                high=pct([r["d_ring_cm"] for r in meas if r["confidence"] == "high"]),
                medium=pct([r["d_ring_cm"] for r in meas if r["confidence"] == "medium"]),
                low=pct([r["d_ring_cm"] for r in meas if r["confidence"] == "low"])),
    height_m=pct([r["top_m"] for r in meas]),
    arc=dict(all=[round(float(np.percentile([r["qc"]["arc"] for r in meas], q)), 2) for q in (10, 50, 90)],
             high=[round(float(np.percentile([r["qc"]["arc"] for r in meas if r["confidence"] == "high"], q)), 2) for q in (10, 50, 90)],
             medium=[round(float(np.percentile([r["qc"]["arc"] for r in meas if r["confidence"] == "medium"], q)), 2) for q in (10, 50, 90)]),
    buttress_suspected=sum(1 for r in meas if r.get("buttress_suspected")),
    span_spread_cm_median=float(np.median([r["span_spread_cm"] for r in meas if r.get("span_spread_cm") is not None])),
)
json.dump(metrics, open(f"{O}/metrics.json", "w"), indent=1)
print(json.dumps(metrics, indent=1))
