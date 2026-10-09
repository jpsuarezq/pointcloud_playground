"""Throwaway prototype v2: stem isolation + DBH from per-tree instances.

Question (ticket #11): can a canopy stem be isolated and measured at 1.3 m on
this capture, and does an axis-fitted / taper-extrapolated diameter behave better
than the naive horizontal slice?  The deliverable is a method whose *failure is
visible*: every measurement carries QC and every refusal is explicit.

PROTOTYPE ONLY. No tests, no tuning, no validation. Not for production.

Method v2 (what changed from v1, and why)
- v1's unbounded RANSAC circle fitted near-collinear arcs with arbitrarily large
  circles (diameters of hundreds of metres on real trees). Bounded now.
- v1's wood-only DBSCAN seed locked onto branch clusters (axis tilt 50-87 deg).
  Now the chosen axis must also pass near the trunk's lower-stem horizontal
  centroid and stay near-vertical; wood seeds are a candidate *generator*, the
  base anchor is the verifier.
- radius is estimated as a per-angular-bin lower quantile of axis distance (the
  innermost coherent surface), robust to both far clusters and twigs crossing
  the axis, instead of a global circle fit.
"""
import numpy as np
from sklearn.cluster import DBSCAN

# ---------------------------------------------------------------- geometry
def fit_line(P):
    c = P.mean(0)
    Q = P - c
    _, evec = np.linalg.eigh(Q.T @ Q / len(P))
    u = evec[:, -1]
    return c, u / np.linalg.norm(u)


def perp_dist(P, c, u):
    v = P - c
    return np.linalg.norm(v - np.outer(v @ u, u), axis=1)


def robust_line(P, trim=0.20, iters=5):
    c, u = fit_line(P)
    for _ in range(iters):
        d = perp_dist(P, c, u)
        keep = d <= np.quantile(d, 1 - trim)
        if keep.sum() < 4:
            break
        c, u = fit_line(P[keep])
    return c, u


def plane_basis(u):
    a = np.array([0.0, 0.0, 1.0])
    if abs(u @ a) > 0.9:
        a = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(u, a); e1 /= np.linalg.norm(e1)
    return e1, np.cross(u, e1)


def tilt_deg(u):
    return float(np.degrees(np.arccos(min(1.0, abs(u[2])))))


def kasa_circle(xy):
    x, y = xy[:, 0], xy[:, 1]
    A = np.column_stack([2 * x, 2 * y, np.ones(len(x))])
    (cx, cy, k), *_ = np.linalg.lstsq(A, x * x + y * y, rcond=None)
    return cx, cy, float(np.sqrt(max(k + cx * cx + cy * cy, 0)))


def arc_coverage(xy, cx, cy, nbins=16):
    ang = np.arctan2(xy[:, 1] - cy, xy[:, 0] - cx)
    h, _ = np.histogram(ang, bins=nbins, range=(-np.pi, np.pi))
    return float((h > 0).mean())


def ransac_circle(xy, tol=0.03, r_min=0.015, r_max=0.80, iters=1000, seed=0):
    """Radius-bounded robust circle fit (v1 had no bounds -> runaway radii)."""
    n = len(xy)
    if n < 3:
        return None
    ax, ay = xy[:, 0], xy[:, 1]
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(iters):
        i, j, k = rng.choice(n, 3, replace=False)
        (x1, y1), (x2, y2), (x3, y3) = xy[i], xy[j], xy[k]
        d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
        if abs(d) < 1e-9:
            continue
        s1, s2, s3 = x1 * x1 + y1 * y1, x2 * x2 + y2 * y2, x3 * x3 + y3 * y3
        ux = (s1 * (y2 - y3) + s2 * (y3 - y1) + s3 * (y1 - y2)) / d
        uy = (s1 * (x3 - x2) + s2 * (x1 - x3) + s3 * (x2 - x1)) / d
        r = np.hypot(x1 - ux, y1 - uy)
        if not (r_min <= r <= r_max):
            continue
        cnt = int((np.abs(np.hypot(ax - ux, ay - uy) - r) < tol).sum())
        if best is None or cnt > best[0]:
            best = (cnt, ux, uy, r)
    if best is None:
        return None
    _, ux, uy, r = best
    for _ in range(2):
        inl = np.abs(np.hypot(ax - ux, ay - uy) - r) < tol
        if inl.sum() < 3:
            break
        ux, uy, r = kasa_circle(xy[inl])
        if not (r_min <= r <= r_max):
            break
    inl = np.abs(np.hypot(ax - ux, ay - uy) - r) < tol
    if inl.sum() < 3:
        return None
    rad = np.hypot(ax[inl] - ux, ay[inl] - uy)
    return dict(cx=float(ux), cy=float(uy), d=float(2 * r), r=float(r),
                n_inl=int(inl.sum()), n=int(n),
                residual=float(np.sqrt(np.mean((rad - r) ** 2))),
                arc=arc_coverage(xy[inl], ux, uy))


# ---------------------------------------------------------------- seeding
def count_bins(P, h, c, u, r, lo=0.8, hi=5.5, step=0.5):
    d = perp_dist(P, c, u)
    nb = 0
    for z in np.arange(lo, hi, step):
        if np.any((h >= z) & (h < z + step) & (d < r)):
            nb += 1
    return nb


def base_anchor(P, h):
    lo = (h >= 0.5) & (h <= 2.5)
    if lo.sum() < 8:
        return None
    return np.array([np.median(P[lo, 0]), np.median(P[lo, 1]),
                     np.median(P[lo, 2][h[lo] >= 1.0]) if (h[lo] >= 1.0).any() else np.median(P[lo, 2])])


def seed_axis(P, h, fm, vstep=0.5, lo=1.0, hi=5.0, r=0.30, min_bins=4):
    """Find the trunk as the near-vertical woody column that persists through the
    most height bins (a Hough-like search), not the densest wood cluster.

    Candidate columns are vertical lines through wood points' horizontal
    positions; a column scores by how many 0.5 m bins in [lo, hi] it holds wood
    within r. The winner seeds the axis; a vertical line at the trunk's lower-stem
    centroid is the fallback. This replaces v1's DBSCAN seed, which locked onto
    branch clusters (axis tilt 50-87 deg)."""
    w = (fm == 1) & (h >= lo - 0.5) & (h <= hi)
    W, Wh = P[w], h[w]
    info = {"n_wood_low": int(w.sum())}
    if len(W) < 6:
        info["reason"] = "insufficient_wood"
        return None, None, info
    xy = W[:, :2]
    cand = xy[:: max(1, len(xy) // 500)]
    edges = np.arange(lo, hi, vstep)
    best = None
    for cxy in cand:
        dd = np.hypot(xy[:, 0] - cxy[0], xy[:, 1] - cxy[1])
        near = dd < r
        covered = sum(1 for e in edges if np.any(near & (Wh >= e) & (Wh < e + vstep)))
        if best is None or covered > best[0]:
            best = (covered, float(near.sum()), cxy)
    covered, tot, cxy = best
    info.update(column_bins=int(covered), column_wood=int(tot),
                column=[round(float(v), 2) for v in cxy])
    if covered < min_bins:
        info["reason"] = "no_persistent_wood_column"
        return None, None, info
    zc = W[(np.hypot(xy[:, 0] - cxy[0], xy[:, 1] - cxy[1]) < r) & (Wh < 3.0), 2]
    z0 = float(np.median(zc)) if len(zc) else float(np.median(W[:, 2]))
    info["seed"] = "vertical_column"
    return np.array([cxy[0], cxy[1], z0]), np.array([0.0, 0.0, 1.0]), info


def refine_axis(P, h, c, u, span=(1.0, 5.0), r_max=0.35, iters=5, max_tilt=35.0):
    for _ in range(iters):
        d = perp_dist(P, c, u)
        m = (h >= span[0]) & (h <= span[1]) & (d < r_max)
        if m.sum() < 8:
            break
        cc, uu = robust_line(P[m])
        if tilt_deg(uu) > max_tilt:          # never let the axis roll onto a branch
            break
        c, u = cc, uu
    return c, u


# ---------------------------------------------------------------- diameter
def ring_radius(P, h, c, u, h0, half=0.15, rmax=0.45, nbins=16, q=25):
    """Innermost coherent surface radius, per angular bin.

    Returns radius (m), arc coverage, per-bin spread and representative points.
    Robust to far clusters (excluded by rmax) and to twigs near the axis
    (low quantile, not the minimum)."""
    m = np.abs(h - h0) <= half
    n_slab = int(m.sum())
    if n_slab < 4:
        return {"ok": False, "reason": "too_few_pts", "n": n_slab}
    Q = P[m]
    d = perp_dist(Q, c, u)
    keep = d < rmax
    if keep.sum() < 4:
        return {"ok": False, "reason": "no_points_near_axis", "n": n_slab, "n_keep": int(keep.sum())}
    Qk, dk = Q[keep], d[keep]
    e1, e2 = plane_basis(u)
    C = Qk - c
    xy = np.column_stack([C @ e1, C @ e2])
    ang = np.arctan2(xy[:, 1], xy[:, 0])
    b = np.floor((ang + np.pi) / (2 * np.pi) * nbins).astype(int)
    b = np.clip(b, 0, nbins - 1)
    reps, rep_xy, radii = [], [], []
    for j in range(nbins):
        sel = b == j
        if sel.sum() >= 1:
            rj = np.percentile(dk[sel], q)
            reps.append(rj)
            i = np.argmin(np.abs(dk[sel] - rj))
            rep_xy.append(xy[sel][i])
            radii.append(dk[sel])
    if len(reps) < 4:
        return {"ok": False, "reason": "insufficient_arc", "n": n_slab,
                "n_keep": int(keep.sum()), "n_bins": len(reps)}
    reps = np.array(reps)
    rep_xy = np.array(rep_xy)
    r_est = float(np.median(reps))
    cx, cy, r_fit = kasa_circle(rep_xy)
    fc = ransac_circle(xy, r_min=0.015, r_max=min(rmax, max(r_est * 2.5, 0.10)))
    return {"ok": True, "h0": float(h0), "half": float(half),
            "r": r_est, "d": 2 * r_est, "arc": len(reps) / nbins,
            "n_bins": len(reps), "n_slab": n_slab, "n_keep": int(keep.sum()),
            "r_bin_std": float(np.std(reps)),
            "circle_fit": None if fc is None else {k: round(v, 4) if isinstance(v, float) else v
                                                   for k, v in fc.items()},
            "r_circle_fit": float(r_fit), "d_circle_fit": float(2 * r_fit)}


def radial_ring_profile(P, h, c, u, hs, **kw):
    rows = []
    for h0 in hs:
        r = ring_radius(P, h, c, u, h0, **kw)
        if r["ok"]:
            rows.append((float(h0), r["d"], r["n_bins"], r["arc"], r))
    return rows


def fit_taper(rows):
    """Power-law d = A*h^b over the clean zone, extrapolated to 1.3 m."""
    if len(rows) < 4:
        return None
    hh = np.array([r[0] for r in rows]); dd = np.array([r[1] for r in rows])
    ok = dd > 0.02
    if ok.sum() < 4:
        return None
    hh, dd = hh[ok], dd[ok]
    b, logA = np.polyfit(np.log(hh), np.log(dd), 1)
    A = float(np.exp(logA))
    slope, inter = np.polyfit(hh, dd, 1)
    pred = A * hh ** b
    r2 = float(1 - np.sum((dd - pred) ** 2) / max(np.sum((dd - dd.mean()) ** 2), 1e-12))
    res = np.sqrt(np.mean((dd - pred) ** 2))
    return dict(A=A, b=float(b), d13_pow=float(A * 1.3 ** b),
                d13_lin=float(slope * 1.3 + inter), r2=r2, rmse=float(res), n=len(hh))


def measure_tree(P, h, fm, tree_id):
    out = {"tree": int(tree_id), "n_pts": int(len(P))}
    c, u, info = seed_axis(P, h, fm)
    out["seed"] = info
    if c is None:
        out["status"] = "refused"; out["reason"] = info.get("reason", "no_seed")
        return out
    c, u = refine_axis(P, h, c, u)
    out["axis_tilt_deg"] = round(tilt_deg(u), 1)
    out["axis_persist_bins"] = int(count_bins(P, h, c, u, r=0.45))

    direct = ring_radius(P, h, c, u, 1.3, half=0.15, rmax=0.45)
    out["d_ring_1p3"] = direct
    if direct["ok"]:
        out["d_ring_cm"] = round(direct["d"] * 100, 1)
        out["qc"] = {"arc": round(direct["arc"], 2), "r_bin_std_cm": round(direct["r_bin_std"] * 100, 2),
                     "n_bins": direct["n_bins"], "n_slab": direct["n_slab"]}
    horiz = ransac_circle(_horiz_xy(P, h, 1.3, 0.15, c), r_min=0.015, r_max=0.80)
    out["d_horiz_cm"] = None if horiz is None else round(horiz["d"] * 100, 1)

    prof = radial_ring_profile(P, h, c, u, np.arange(2.0, 5.01, 0.25), half=0.15, rmax=0.45)
    out["ring_profile_2to5"] = [[round(p[0], 2), round(p[1], 3), p[2], round(p[3], 2)] for p in prof]
    tp = fit_taper(prof)
    if tp:
        out["taper"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in tp.items()}
        out["d_taper_cm"] = round(tp["d13_pow"] * 100, 1)
    # buttress: is the 1.3 m ring much wider than the clean 3 m ring?
    d13 = direct.get("d") if direct["ok"] else None
    d30 = next((p[1] for p in prof if abs(p[0] - 3.0) < 1e-6), None)
    if d13 and d30:
        out["flare_ratio"] = round(d13 / d30, 2)
        out["buttress_suspected"] = bool(d13 / d30 > 1.35)
    if out.get("buttress_suspected"):
        dab = ring_radius(P, h, c, u, 2.0, half=0.15, rmax=0.45)
        if dab["ok"]:
            out["d_dab_cm"] = round(dab["d"] * 100, 1)
            out["dab_height_m"] = 2.0

    # stability across axis-fitting spans
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

    ok_direct = direct["ok"] and direct["arc"] >= 0.35 and direct["d"] < 1.0
    if ok_direct:
        out["status"] = "measured"
    else:
        out["status"] = "refused"
        out["reason"] = direct.get("reason", "qc_reject")
    out["taper_available"] = "taper" in out
    return out


def _horiz_xy(P, h, h0, half, c):
    m = np.abs(h - h0) <= half
    return (P[m][:, :2] - c[:2]) if m.sum() else np.empty((0, 2))
