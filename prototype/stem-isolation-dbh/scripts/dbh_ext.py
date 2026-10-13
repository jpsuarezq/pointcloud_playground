"""Extended stem isolation for the whole canopy: seed the axis even when FM
labels no wood.

The v2 prototype measured the ring from ALL points near the axis — wood was only
used to SEED the axis. So a tree that FM left unlabelled is recoverable if the
axis can be seeded another way. This module adds a fallback chain:

  1. FM wood column search (v2 behaviour).
  2. intensity-proxy column search — low intensity (< per-scan threshold derived
     from the FM wood/foliage medians) is a wood proxy.
  3. base anchor — the median XY of the lower-stem points (0.5-2.5 m), vertical
     axis. Needs no labels at all.

Every measurement carries a confidence and the seed used. Nothing is tuned to
ground truth; thresholds are proxies derived from the given FM labels.
"""
import numpy as np
from dbh_proto import (perp_dist, plane_basis, tilt_deg, kasa_circle,
                       ransac_circle, count_bins, ring_radius, refine_axis,
                       radial_ring_profile, fit_taper, _horiz_xy)


def column_search(P, h, mask, vstep=0.5, lo=1.0, hi=5.0, r=0.30, min_bins=4,
                  max_cand=500):
    """Vertical-line Hough search over `mask` candidate points."""
    w = mask & (h >= lo - 0.5) & (h <= hi)
    W, Wh = P[w], h[w]
    info = {"n_cand": int(w.sum())}
    if len(W) < 6:
        info["reason"] = "insufficient_candidates"
        return None, None, info
    xy = W[:, :2]
    cand = xy[:: max(1, len(xy) // max_cand)]
    edges = np.arange(lo, hi, vstep)
    best = None
    for cxy in cand:
        dd = np.hypot(xy[:, 0] - cxy[0], xy[:, 1] - cxy[1])
        near = dd < r
        covered = sum(1 for e in edges if np.any(near & (Wh >= e) & (Wh < e + vstep)))
        if best is None or covered > best[0]:
            best = (covered, int(near.sum()), cxy)
    covered, tot, cxy = best
    info.update(column_bins=int(covered), column_n=int(tot),
                column=[round(float(v), 2) for v in cxy])
    if covered < min_bins:
        info["reason"] = "no_persistent_column"
        return None, None, info
    zc = W[(np.hypot(xy[:, 0] - cxy[0], xy[:, 1] - cxy[1]) < r) & (Wh < 3.0), 2]
    z0 = float(np.median(zc)) if len(zc) else float(np.median(W[:, 2]))
    return np.array([cxy[0], cxy[1], z0]), np.array([0.0, 0.0, 1.0]), info


def base_anchor(P, h, lo=0.5, hi=2.5, min_pts=8):
    m = (h >= lo) & (h <= hi)
    if m.sum() < min_pts:
        return None
    xy = np.median(P[m, :2], axis=0)
    zc = P[m, 2]
    return np.array([xy[0], xy[1], float(np.median(zc))])


def measure_tree_ext(P, h, fm, inten, cloud, tau, tree_id):
    out = {"tree": int(tree_id), "n_pts": int(len(P))}
    if len(P) == 0:
        out.update(status="refused", reason="empty")
        return out

    wood = fm == 1
    t = np.where(cloud == 0, tau[0], tau[1])
    proxy = inten < t
    seed = None
    kind = None

    c, u, info = column_search(P, h, wood)
    if c is not None:
        kind, seed = "wood", (c, u, info)
    else:
        c2, u2, info2 = column_search(P, h, proxy)
        if c2 is not None:
            kind, seed = "intensity", (c2, u2, info2)
    if seed is None:
        c3 = base_anchor(P, h)
        if c3 is not None:
            kind, seed = "base", (c3, np.array([0.0, 0.0, 1.0]), {"seed": "base_anchor"})
    if seed is None:
        out.update(status="refused", reason="no_seed", seed_kind=None)
        return out

    c, u, info = seed
    out["seed_kind"] = kind
    out["seed"] = info
    c, u = refine_axis(P, h, c, u)
    out["axis_tilt_deg"] = round(tilt_deg(u), 1)
    out["axis_persist_bins"] = int(count_bins(P, h, c, u, r=0.45))

    direct = ring_radius(P, h, c, u, 1.3, half=0.15, rmax=0.45)
    out["d_ring_1p3"] = direct
    if not direct["ok"]:
        # widen once: a proxy-seeded axis may sit slightly off the stem
        direct = ring_radius(P, h, c, u, 1.3, half=0.15, rmax=0.65)
        out["d_ring_1p3_wide"] = direct
    if direct["ok"]:
        out["d_ring_cm"] = round(direct["d"] * 100, 1)
        out["qc"] = {"arc": round(direct["arc"], 2),
                     "r_bin_std_cm": round(direct["r_bin_std"] * 100, 2),
                     "n_bins": direct["n_bins"], "n_slab": direct["n_slab"]}
        # what fraction of the accepted ring is FM wood / proxy wood?
        m = np.abs(h - 1.3) <= 0.15
        Q = P[m]; d = perp_dist(Q, c, u); keep = d < 0.45
        if keep.sum():
            out["qc"]["frac_fm_wood"] = round(float((fm[m][keep] == 1).mean()), 2)
            out["qc"]["frac_proxy_wood"] = round(float((inten[m][keep] < np.where(cloud[m][keep] == 0, tau[0], tau[1])).mean()), 2)

    horiz = ransac_circle(_horiz_xy(P, h, 1.3, 0.15, c), r_min=0.015, r_max=0.80)
    out["d_horiz_cm"] = None if horiz is None else round(horiz["d"] * 100, 1)

    prof = radial_ring_profile(P, h, c, u, np.arange(2.0, 5.01, 0.25), half=0.15, rmax=0.45)
    out["ring_profile_2to5"] = [[round(p[0], 2), round(p[1], 3), p[2], round(p[3], 2)] for p in prof]
    tp = fit_taper(prof)
    if tp:
        out["taper"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in tp.items()}
        out["d_taper_cm"] = round(tp["d13_pow"] * 100, 1)

    d13 = direct.get("d") if direct["ok"] else None
    d30 = next((p[1] for p in prof if abs(p[0] - 3.0) < 1e-6), None)
    if d13 and d30:
        out["flare_ratio"] = round(d13 / d30, 2)
        out["buttress_suspected"] = bool(d13 / d30 > 1.35)

    stab = {}
    for lo, hi in [(1.0, 3.0), (1.0, 5.0), (2.0, 5.0), (1.0, 6.0)]:
        cc, uu = refine_axis(P, h, c, u, span=(lo, hi))
        r = ring_radius(P, h, cc, uu, 1.3, half=0.15, rmax=0.45)
        if r["ok"]:
            stab[f"{lo}-{hi}"] = round(r["d"] * 100, 1)
    if stab:
        v = list(stab.values())
        out["span_stability_cm"] = stab
        out["span_spread_cm"] = round(max(v) - min(v), 1)

    # ---- status + confidence
    if direct["ok"]:
        arc = direct["arc"]
        spread = out.get("span_spread_cm")
        if kind == "wood" and arc >= 0.35 and direct["d"] < 1.0:
            out["status"], out["confidence"] = "measured", "high"
        elif arc >= 0.35 and direct["d"] < 1.0:
            out["status"], out["confidence"] = "measured", "medium"
        else:
            out["status"], out["confidence"] = "measured", "low"
            out["low_reasons"] = [x for x, cond in
                                  [("arc<0.35", arc < 0.35),
                                   ("d>=1m", direct["d"] >= 1.0),
                                   ("spread>20cm", (spread or 0) > 20)] if cond]
    else:
        out["status"], out["confidence"] = "refused", "none"
        out["reason"] = direct.get("reason", "qc_reject")
    out["taper_available"] = "taper" in out
    return out
