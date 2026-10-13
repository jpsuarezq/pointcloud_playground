"""Multi-candidate axis seeding.

The dominant recoverable failure is `no_points_near_axis`: the single seed lands
in a gap, so the ring sees nothing. Here we generate several candidate axes,
refine each, measure the 1.3 m ring for each, and pick the best by a QC-driven
score. This keeps the "failure visible" discipline: a candidate only counts as
`valid` when arc >= 0.35 and d < 1 m; otherwise the best-effort ring is emitted
with confidence `low`.
"""
import numpy as np
from sklearn.cluster import DBSCAN
from dbh_proto import (ring_radius, refine_axis, tilt_deg, count_bins,
                       ransac_circle, _horiz_xy, radial_ring_profile, fit_taper,
                       perp_dist)
from dbh_ext import column_search, base_anchor

VALID_ARC = 0.35
VALID_D = 1.0


def band_seed(P, h, base, lo=1.15, hi=1.45):
    """Seed from the band points themselves: densest cluster nearest the base."""
    m = (h >= lo) & (h <= hi)
    if m.sum() < 4:
        return None
    xy = P[m, :2]
    z0 = float(np.median(P[m, 2]))
    if len(xy) >= 8:
        lab = DBSCAN(eps=0.35, min_samples=4).fit_predict(xy)
        best = None
        for c in set(lab.tolist()):
            if c < 0:
                continue
            cx, cy = xy[lab == c].mean(0)
            d = np.hypot(cx - base[0], cy - base[1]) if base is not None else 0.0
            if best is None or d < best[0]:
                best = (d, cx, cy)
        if best is not None:
            return np.array([best[1], best[2], z0])
    cx, cy = xy.mean(0)
    return np.array([cx, cy, z0])


def _candidates(P, h, fm, proxy, tau, cloud):
    cands = []
    c, u, _ = column_search(P, h, fm == 1)
    if c is not None:
        cands.append(("wood", c, u))
    c, u, _ = column_search(P, h, proxy)
    if c is not None:
        cands.append(("intensity", c, u))
    b = base_anchor(P, h)
    if b is not None:
        cands.append(("base", b, np.array([0.0, 0.0, 1.0])))
    bs = band_seed(P, h, b)
    if bs is not None:
        cands.append(("band", bs, np.array([0.0, 0.0, 1.0])))
    return cands


def _ring_best(P, h, c, u):
    """Return the NARROWEST slab that yields a valid ring; else the best-effort
    (most arc) ring for a low-confidence emission."""
    low = None
    for half, rmax in [(0.15, 0.45), (0.15, 0.65), (0.30, 0.45), (0.30, 0.65)]:
        r = ring_radius(P, h, c, u, 1.3, half=half, rmax=rmax)
        if r.get("ok") and r["arc"] >= VALID_ARC and r["d"] < VALID_D:
            return r, half
        if low is None or (r.get("ok") and (not low[0].get("ok") or r["arc"] > low[0]["arc"])):
            low = (r, half)
    return low


def measure_tree_multi(P, h, fm, inten, cloud, tau, tree_id):
    out = {"tree": int(tree_id), "n_pts": int(len(P))}
    if len(P) == 0:
        out.update(status="refused", reason="empty")
        return out
    t = np.where(cloud == 0, tau[0], tau[1])
    proxy = inten < t
    cands = _candidates(P, h, fm, proxy, tau, cloud)
    if not cands:
        out.update(status="refused", reason="no_seed", seed_kind=None)
        return out

    tried = []
    best = None  # (score, valid, kind, c, u, ring, half)
    for kind, c0, u0 in cands:
        c, u = refine_axis(P, h, c0, u0)
        ring, half = _ring_best(P, h, c, u)
        ok = ring.get("ok")
        arc = ring.get("arc", 0.0) if ok else -1.0
        nkeep = ring.get("n_keep", 0) if ok else 0
        valid = ok and arc >= VALID_ARC and ring["d"] < VALID_D
        tried.append({"kind": kind, "ok": bool(ok), "arc": round(arc, 2),
                      "n_keep": int(nkeep), "half": half,
                      "reason": ring.get("reason")})
        # prefer valid, then any ok ring, then the narrowest slab, then more arc
        score = (1 if valid else 0, 1 if ok else 0, -half, arc, nkeep)
        if best is None or score > best[0]:
            best = (score, valid, kind, c, u, ring, half)

    score, valid, kind, c, u, ring, half = best
    out["seed_kind"] = kind
    out["candidates_tried"] = tried
    out["slab_half"] = half
    out["axis_tilt_deg"] = round(tilt_deg(u), 1)
    out["axis_persist_bins"] = int(count_bins(P, h, c, u, r=0.45))
    out["d_ring_1p3"] = ring
    if ring.get("ok"):
        out["d_ring_cm"] = round(ring["d"] * 100, 1)
        out["qc"] = {"arc": round(ring["arc"], 2),
                     "r_bin_std_cm": round(ring["r_bin_std"] * 100, 2),
                     "n_bins": ring["n_bins"], "n_slab": ring["n_slab"]}
        m = np.abs(h - 1.3) <= half
        Q = P[m]; d = perp_dist(Q, c, u); keep = d < 0.45
        if keep.sum():
            out["qc"]["frac_fm_wood"] = round(float((fm[m][keep] == 1).mean()), 2)

    horiz = ransac_circle(_horiz_xy(P, h, 1.3, 0.15, c), r_min=0.015, r_max=0.80)
    out["d_horiz_cm"] = None if horiz is None else round(horiz["d"] * 100, 1)

    prof = radial_ring_profile(P, h, c, u, np.arange(2.0, 5.01, 0.25), half=half, rmax=0.45)
    out["ring_profile_2to5"] = [[round(p[0], 2), round(p[1], 3), p[2], round(p[3], 2)] for p in prof]
    tp = fit_taper(prof)
    if tp:
        out["taper"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in tp.items()}
        out["d_taper_cm"] = round(tp["d13_pow"] * 100, 1)
    d13 = ring.get("d") if ring.get("ok") else None
    d30 = next((p[1] for p in prof if abs(p[0] - 3.0) < 1e-6), None)
    if d13 and d30:
        out["flare_ratio"] = round(d13 / d30, 2)
        out["buttress_suspected"] = bool(d13 / d30 > 1.35)

    stab = {}
    for lo, hi in [(1.0, 3.0), (1.0, 5.0), (2.0, 5.0), (1.0, 6.0)]:
        cc, uu = refine_axis(P, h, c, u, span=(lo, hi))
        r = ring_radius(P, h, cc, uu, 1.3, half=half, rmax=0.45)
        if r["ok"]:
            stab[f"{lo}-{hi}"] = round(r["d"] * 100, 1)
    if stab:
        v = list(stab.values())
        out["span_stability_cm"] = stab
        out["span_spread_cm"] = round(max(v) - min(v), 1)

    if ring.get("ok"):
        spread = out.get("span_spread_cm")
        if valid and kind == "wood":
            out["status"], out["confidence"] = "measured", "high"
        elif valid:
            out["status"], out["confidence"] = "measured", "medium"
        else:
            out["status"], out["confidence"] = "measured", "low"
            out["low_reasons"] = [x for x, cond in
                                  [("arc<0.35", ring["arc"] < VALID_ARC),
                                   ("d>=1m", ring["d"] >= VALID_D),
                                   ("wide_slab", half > 0.15),
                                   ("spread>20cm", (spread or 0) > 20)] if cond]
    else:
        out["status"], out["confidence"] = "refused", "none"
        out["reason"] = ring.get("reason", "qc_reject")
    out["taper_available"] = "taper" in out
    return out
