"""Figures + summary for the DBH prototype (ticket #11)."""
import numpy as np, json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
import dbh_proto as D

F = "/tmp/opencode/cloud-research/out/fig"
O = "/tmp/opencode/cloud-research/out/inst"
res = json.load(open(f"{O}/proto_results.json"))

# ---- fig 1: exemplar 2957 side view with fitted axis ----
z = np.load(f"{F}/tree_2957.npz")
P = z["P"].astype(np.float64); h = z["h"].astype(np.float64); f = z["f"]
c, u, info = D.seed_axis(P, h, f); c, u = D.refine_axis(P, h, c, u)
fig, ax = plt.subplots(figsize=(4.5, 7))
for cls, col, lab in [(2, "#2ca02c", "leaf (FM=2)"), (1, "#d62728", "wood (FM=1)"),
                      (0, "#7f7f7f", "ground/under (FM=0)")]:
    m = f == cls
    if m.sum():
        ax.scatter(P[m, 0], P[m, 2], s=0.1, c=col, label=f"{lab}  n={int(m.sum())}", rasterized=True)
ts = np.array([h.min(), h.max()])
# draw axis over 0.5..30 m of height
zz = np.linspace(c[2], c[2] + 28, 50)
ax.plot(c[0] + (zz - c[2]) * u[0] / u[2], zz, "k-", lw=1.2, label="fitted stem axis")
ax.set_xlabel("UTM X (m)"); ax.set_ylabel("Z (m)")
ax.set_title("Exemplar 2957: isolated stem axis\n(tilt %.1f deg, %d wood pts in 1-5 m)"
             % (D.tilt_deg(u), info["n_wood_low"]))
ax.legend(markerscale=15, fontsize=7, loc="upper left")
ax.set_aspect("equal"); fig.tight_layout(); fig.savefig(f"{F}/proto_fig1_exemplar.png", dpi=150)
plt.close(fig)

# ---- fig 2: plan views at heights for 2957 ----
fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
for ax, h0 in zip(axes, [1.3, 2.0, 3.0, 4.0]):
    m = np.abs(h - h0) <= 0.15
    Q = P[m]; d = D.perp_dist(Q, c, u)
    e1, e2 = D.plane_basis(u); C = Q - c
    xy = np.column_stack([C @ e1, C @ e2])
    for cls, col in [(2, "#2ca02c"), (1, "#d62728"), (0, "#7f7f7f")]:
        s = f[m] == cls
        if s.sum():
            ax.scatter(xy[s, 0], xy[s, 1], s=8, c=col)
    r = D.ring_radius(P, h, c, u, h0, half=0.15, rmax=0.45)
    ex = D.ring_radius(P, h, c, u, h0, half=0.15, rmax=0.60)
    for rr, ls, lab in [(r, "-", "ring"), (ex, "--", "ring 0.6")]:
        if rr.get("ok"):
            th = np.linspace(0, 2 * np.pi, 100)
            ax.plot(rr["r"] * np.cos(th), rr["r"] * np.sin(th), "k" + ls, lw=1,
                    label=f"{lab} d={rr['d']*100:.0f}cm arc={rr['arc']:.2f}")
    ax.set_aspect("equal"); ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2)
    ax.axhline(0, color="k", lw=0.3); ax.axvline(0, color="k", lw=0.3)
    ax.set_title(f"h = {h0} m (n={int(m.sum())})"); ax.legend(fontsize=7)
fig.suptitle("Exemplar 2957: perpendicular plan views around the fitted axis")
fig.tight_layout(); fig.savefig(f"{F}/proto_fig2_plan.png", dpi=150); plt.close(fig)

# ---- fig 3: method agreement ----
m_rows = [r for r in res if r["status"] == "measured" and r.get("d_taper_cm")]
h_rows = [r for r in res if r.get("d_ring_cm") and r.get("d_horiz_cm")]
t_rows = [r for r in res if r.get("d_ring_cm") and r.get("d_taper_cm")]
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
ax = axes[0]
x = [r["d_ring_cm"] for r in h_rows]; y = [r["d_horiz_cm"] for r in h_rows]
ax.scatter(x, y, c="#1f77b4")
lim = [0, max(x + y) * 1.1]
ax.plot(lim, lim, "k--", lw=0.8)
ax.set_xlabel("perpendicular ring d at 1.3 m (cm)"); ax.set_ylabel("naive horizontal slice d (cm)")
ax.set_title("Naive horizontal slice vs axis-perpendicular\n(above the line = elliptical inflation)")
for r in h_rows:
    if r["tree"] == 2957:
        ax.annotate("2957", (r["d_ring_cm"], r["d_horiz_cm"]))
ax = axes[1]
x = [r["d_ring_cm"] for r in t_rows]; y = [r["d_taper_cm"] for r in t_rows]
ax.scatter(x, y, c="#d62728")
lim = [0, max(x + y) * 1.05]
ax.plot(lim, lim, "k--", lw=0.8)
ax.set_xlabel("direct ring d at 1.3 m (cm)"); ax.set_ylabel("taper-extrapolated d at 1.3 m (cm)")
ax.set_title("Taper extrapolation vs direct\n(points far off the line = no agreement)")
fig.tight_layout(); fig.savefig(f"{F}/proto_fig3_agreement.png", dpi=150); plt.close(fig)

# ---- fig 4: QC + measurability ----
measured = [r for r in res if r["status"] == "measured"]
ref = [r for r in res if r["status"] != "measured"]
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
ax = axes[0]
arc = [r["qc"]["arc"] for r in measured]
ax.hist(arc, bins=np.arange(0, 1.05, 0.1), color="#1f77b4")
ax.axvline(0.35, color="k", ls="--", label="accept threshold 0.35")
ax.set_xlabel("arc coverage of accepted ring"); ax.set_ylabel("trees"); ax.legend()
ax.set_title(f"Accepted rings (n={len(measured)})")
ax = axes[1]
sp = [r["span_spread_cm"] for r in measured if r.get("span_spread_cm") is not None]
ax.hist(sp, bins=np.arange(0, 60, 5), color="#ff7f0e")
ax.set_xlabel("d(range) across axis-fitting spans (cm)"); ax.set_ylabel("trees")
ax.set_title("Stability across fitting spans")
ax = axes[2]
from collections import Counter
cc = Counter(r.get("reason", "measured") for r in res)
ax.barh(list(cc.keys()), list(cc.values()), color="#7f7f7f")
ax.set_xlabel("trees"); ax.set_title(f"Outcome, n={len(res)} sampled")
fig.tight_layout(); fig.savefig(f"{F}/proto_fig4_qc.png", dpi=150); plt.close(fig)

# ---- summary ----
print("sampled", len(res), "measured", len(measured), "refused", len(ref))
print("refusal reasons:", dict(Counter(r.get("reason") for r in ref)))
dr = [r["d_ring_cm"] for r in measured]
print("direct ring d cm: pct", np.percentile(dr, [10, 25, 50, 75, 90]).round(1))
arcs = [r["qc"]["arc"] for r in measured]
print("arc: pct", np.percentile(arcs, [10, 50, 90]).round(2))
if h_rows:
    ratio = [r["d_horiz_cm"] / r["d_ring_cm"] for r in h_rows]
    print("horiz/perp ratio: median %.2f  max %.2f" % (np.median(ratio), max(ratio)))
if t_rows:
    tr = [r["d_taper_cm"] / r["d_ring_cm"] for r in t_rows]
    print("taper/direct ratio: median %.2f  range %.2f-%.2f" % (np.median(tr), min(tr), max(tr)))
    print("agreement within 25%:", sum(0.75 < v < 1.25 for v in tr), "/", len(tr))
print("buttress suspected:", sum(1 for r in measured if r.get("buttress_suspected")), "/", len(measured))
