"""Parameterised, automated DBH + height pipeline: one point cloud in, all
outputs out.

    python run_pipeline.py --config pipeline_config.json
    python run_pipeline.py --input /path/cloud.laz --out /path/run

Steps (selectable with --steps): cache, ground, hag, canopy, measure, export,
figures. Everything is driven by a config so a new cloud only needs the input
path (and, if its dimension names differ, the field mapping).

No fitting, no tuning to ground truth: all thresholds are principled proxies and
are locked in the config.
"""
import argparse, csv, json, os, subprocess, sys, time
import numpy as np
import laspy
from laspy import ExtraBytesParams
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dbh_multi import measure_tree_multi

DEFAULTS = dict(
    # field mapping (extra dimensions of the input cloud)
    fields=dict(instance="PredInstance_FM", semantic="PredSemantic_FM",
                score="PredScore_FM", scan_index="Original cloud index"),
    # semantic class ids
    wood_class=1, ground_class=0, foliage_class=2,
    # geometry / measurement
    canopy_top_hag_m=12.0, dbh_height_m=1.3,
    band_half_m=0.15, band_half_wide_m=0.30,
    ring_rmax_m=0.45, ring_rmax_wide_m=0.65,
    valid_arc=0.35, valid_dbh_max_m=1.0,
    # ground / DEM
    smrf_cell_m=0.25, dem_res_m=0.25, dem_fill_m=2.0,
    # export
    export_point_cap=400_000,
    # intensity wood proxy (None -> auto from wood/foliage medians)
    tau=None,
)


def load_config(args):
    cfg = json.loads(json.dumps(DEFAULTS))
    if args.config:
        user = json.load(open(args.config))
        cfg.update({k: v for k, v in user.items() if k != "fields"})
        cfg["fields"].update(user.get("fields", {}))
    if args.input:
        cfg["input"] = args.input
    if args.out:
        cfg["out"] = args.out
    if not cfg.get("input") or not cfg.get("out"):
        sys.exit("need --input and --out (or a config with both)")
    cfg["cache"] = os.path.join(cfg["out"], "cache")
    cfg["hag"] = os.path.join(cfg["out"], "hag")
    cfg["trees"] = os.path.join(cfg["out"], "trees")
    cfg["figures"] = os.path.join(cfg["out"], "figures")
    for d in ("out", "cache", "hag", "trees", "figures"):
        os.makedirs(cfg[d], exist_ok=True)
    json.dump(cfg, open(os.path.join(cfg["out"], "config.json"), "w"), indent=1)
    return cfg


# ------------------------------------------------------------------ cache
def build_cache(cfg):
    C, F = cfg["cache"], cfg["fields"]
    f = laspy.open(cfg["input"]); h = f.header; n = h.point_count
    names = list(h.point_format.dimension_names)

    def has(x): return x in names
    xyz = np.memmap(f"{C}/xyz.f4", np.float32, "w+", shape=(n, 3))
    inst = np.memmap(f"{C}/inst.u4", np.uint32, "w+", shape=(n,))
    sem = np.memmap(f"{C}/sem.u1", np.uint8, "w+", shape=(n,))
    inten = np.memmap(f"{C}/intensity.u2", np.uint16, "w+", shape=(n,))
    scan = np.memmap(f"{C}/scan.u1", np.uint8, "w+", shape=(n,))
    score = np.memmap(f"{C}/score.f4", np.float32, "w+", shape=(n,))
    ox, oy, oz = h.offsets
    i = 0
    for pts in f.chunk_iterator(8_000_000):
        m = len(pts); sl = slice(i, i + m)
        xyz[sl, 0] = np.asarray(pts.x, np.float64) - ox
        xyz[sl, 1] = np.asarray(pts.y, np.float64) - oy
        xyz[sl, 2] = np.asarray(pts.z, np.float64) - oz
        inst[sl] = np.asarray(pts[F["instance"]], np.uint32) if has(F["instance"]) else 0
        sem[sl] = np.asarray(pts[F["semantic"]], np.uint8) if has(F["semantic"]) else 0
        inten[sl] = np.asarray(pts.intensity, np.uint16) if has("intensity") else 0
        scan[sl] = np.asarray(pts[F["scan_index"]], np.uint8) if has(F["scan_index"]) else 0
        score[sl] = np.asarray(pts[F["score"]], np.float32) if has(F["score"]) else 0.0
        i += m
    f.close()
    for a in (xyz, inst, sem, inten, scan, score):
        a.flush()
    json.dump(dict(n=int(n), ox=float(ox), oy=float(oy), oz=float(oz),
                   scales=[float(s) for s in h.scales]), open(f"{C}/meta.json", "w"), indent=1)
    print(f"[cache] {n} points -> {C}")


# ------------------------------------------------------------------ ground + HAG
def build_ground(cfg):
    pipe = os.path.join(cfg["out"], "smrf.json")
    json.dump({"pipeline": [cfg["input"],
                            {"type": "filters.smrf", "cell": cfg["smrf_cell_m"]},
                            {"type": "filters.range", "limits": "Classification[2:2]"},
                            {"type": "writers.las",
                             "filename": os.path.join(cfg["out"], "ground.laz"),
                             "compression": "laszip"}]}, open(pipe, "w"))
    print("[ground] pdal filters.smrf ...")
    subprocess.run(["pdal", "pipeline", pipe], check=True)
    print("[ground] done")


def build_hag(cfg):
    C, H, meta = cfg["cache"], cfg["hag"], None
    meta = json.load(open(f"{C}/meta.json")); n = meta["n"]
    res = cfg["dem_res_m"]
    xyz = np.memmap(f"{C}/xyz.f4", np.float32, "r", shape=(n, 3))
    x = np.asarray(xyz[:, 0], np.float64) + meta["ox"]
    y = np.asarray(xyz[:, 1], np.float64) + meta["oy"]
    z = np.asarray(xyz[:, 2], np.float64) + meta["oz"]
    x0 = np.floor(x.min() / res) * res; y0 = np.floor(y.min() / res) * res
    nx = int(np.ceil((x.max() - x0) / res)) + 1; ny = int(np.ceil((y.max() - y0) / res)) + 1
    gx, gy, gz = [], [], []
    with laspy.open(os.path.join(cfg["out"], "ground.laz")) as f:
        for pts in f.chunk_iterator(8_000_000):
            gx.append(np.asarray(pts.x, np.float64)); gy.append(np.asarray(pts.y, np.float64))
            gz.append(np.asarray(pts.z, np.float64))
    gx = np.concatenate(gx); gy = np.concatenate(gy); gz = np.concatenate(gz)
    ix = np.clip(((gx - x0) // res).astype(np.int64), 0, nx - 1)
    iy = np.clip(((gy - y0) // res).astype(np.int64), 0, ny - 1)
    flat = iy * nx + ix; order = np.argsort(flat, kind="stable"); fs = flat[order]
    st = np.searchsorted(fs, np.arange(nx * ny + 1))
    dem = np.full(nx * ny, np.nan); ne = st[:-1] < st[1:]
    dem[ne] = gz[order[st[:-1][ne]]]; dem = dem.reshape(ny, nx)
    valid = np.isfinite(dem)
    d = ndimage.distance_transform_edt(~valid)
    _, inds = ndimage.distance_transform_edt(~valid, return_indices=True)
    dem = np.where(d <= cfg["dem_fill_m"] / res, dem[tuple(inds)], np.nan)
    dem = ndimage.generic_filter(dem, np.nanmedian, size=3, mode="nearest")
    np.save(f"{H}/dem.npy", dem)

    def bilin(a, rf, cf):
        i0 = np.floor(rf).astype(np.int64); j0 = np.floor(cf).astype(np.int64)
        fx = rf - i0; fy = cf - j0
        i0c = np.clip(i0, 0, a.shape[0] - 1); i1c = np.clip(i0 + 1, 0, a.shape[0] - 1)
        j0c = np.clip(j0, 0, a.shape[1] - 1); j1c = np.clip(j0 + 1, 0, a.shape[1] - 1)
        return (a[i0c, j0c] * (1 - fx) * (1 - fy) + a[i0c, j1c] * fx * (1 - fy) +
                a[i1c, j0c] * (1 - fx) * fy + a[i1c, j1c] * fx * fy)

    rf = (y - y0) / res; cf = (x - x0) / res
    inb = (rf >= 0) & (rf < ny - 1.001) & (cf >= 0) & (cf < nx - 1.001)
    hag = z - bilin(dem, rf, cf)
    ok = inb & np.isfinite(hag) & (hag > -50) & (hag < 200)
    np.save(f"{H}/hag.npy", np.where(ok, hag, 0).astype(np.float32))
    np.save(f"{H}/ok.npy", ok)
    print(f"[hag] grid {ny}x{nx}, ok frac {ok.mean():.3f}")


# ------------------------------------------------------------------ canopy
def canopy_summary(cfg):
    C, H = cfg["cache"], cfg["hag"]
    meta = json.load(open(f"{C}/meta.json")); n = meta["n"]
    inst = np.asarray(np.memmap(f"{C}/inst.u4", np.uint32, "r", shape=(n,)))
    sem = np.asarray(np.memmap(f"{C}/sem.u1", np.uint8, "r", shape=(n,)))
    hag = np.asarray(np.load(f"{H}/hag.npy", mmap_mode="r"))
    NI = int(inst.max()) + 1
    npts = np.bincount(inst, minlength=NI)
    top = np.full(NI, -1e9); np.maximum.at(top, inst, hag)
    lo, hi = cfg["dbh_height_m"] - cfg["band_half_m"], cfg["dbh_height_m"] + cfg["band_half_m"]
    band = np.bincount(inst[(hag >= lo) & (hag <= hi) & (inst > 0)], minlength=NI)
    wood = np.bincount(inst[(hag >= lo) & (hag <= hi) & (inst > 0) & (sem == cfg["wood_class"])], minlength=NI)
    canopy = top > cfg["canopy_top_hag_m"]; canopy[0] = False
    rows = [(int(k), float(top[k]), int(npts[k]), int(band[k]), int(wood[k])) for k in np.where(canopy)[0]]
    rows.sort(key=lambda r: -r[3])
    json.dump(rows, open(os.path.join(cfg["out"], "canopy_summary.json"), "w"))
    print(f"[canopy] {len(rows)} canopy instances (top>{cfg['canopy_top_hag_m']}m), "
          f"{sum(1 for r in rows if r[3] >= 1)} with >=1 band pt")
    return rows


# ------------------------------------------------------------------ measure
def measure_all(cfg, rows):
    C, H = cfg["cache"], cfg["hag"]
    meta = json.load(open(f"{C}/meta.json")); n = meta["n"]
    xyz = np.memmap(f"{C}/xyz.f4", np.float32, "r", shape=(n, 3))
    inst = np.asarray(np.memmap(f"{C}/inst.u4", np.uint32, "r", shape=(n,)))
    sem = np.asarray(np.memmap(f"{C}/sem.u1", np.uint8, "r", shape=(n,)))
    inten = np.asarray(np.memmap(f"{C}/intensity.u2", np.uint16, "r", shape=(n,)))
    scan = np.asarray(np.memmap(f"{C}/scan.u1", np.uint8, "r", shape=(n,)))
    hag = np.asarray(np.load(f"{H}/hag.npy", mmap_mode="r"))

    if cfg.get("tau") is not None:
        tau = tuple(cfg["tau"])
    else:
        tau = {}
        for ci in np.unique(scan):
            s = scan == ci
            tw = np.median(inten[s & (sem == cfg["wood_class"])]) if (s & (sem == cfg["wood_class"])).any() else 0
            tf = np.median(inten[s & (sem == cfg["foliage_class"])]) if (s & (sem == cfg["foliage_class"])).any() else 65535
            tau[int(ci)] = float((tw + tf) / 2) if tf > tw else 65535.0
        tau = (tau.get(0, 65535.0), tau.get(1, tau.get(0, 65535.0)))
    print("[measure] tau =", tau)

    band = {r[0]: r[3] for r in rows}; top = {r[0]: r[1] for r in rows}
    keys = [r[0] for r in rows]
    order = np.argsort(inst, kind="stable"); si = inst[order]
    bounds = np.searchsorted(si, np.arange(int(inst.max()) + 2))
    out = []
    t0 = time.time()
    for j, k in enumerate(keys):
        a, b = bounds[k], bounds[k + 1]; idx = order[a:b]
        P = np.column_stack([np.asarray(xyz[idx, 0], np.float64) + meta["ox"],
                             np.asarray(xyz[idx, 1], np.float64) + meta["oy"],
                             np.asarray(xyz[idx, 2], np.float64) + meta["oz"]])
        r = measure_tree_multi(P, hag[idx], sem[idx], inten[idx], scan[idx], tau, k)
        r["top_m"] = round(top[k], 1); r["n_band"] = int(band[k])
        out.append(r)
        if j % 100 == 0:
            print(f"  {j}/{len(keys)} {time.time()-t0:.0f}s", flush=True)
    json.dump(out, open(os.path.join(cfg["out"], "results.json"), "w"), indent=1, default=str)
    from collections import Counter
    print("[measure]", dict(Counter(r["status"] for r in out)),
          dict(Counter(r.get("confidence") for r in out if r["status"] == "measured")))
    return out


# ------------------------------------------------------------------ export
def export(cfg, results, rows):
    from pyproj import Transformer
    C, H = cfg["cache"], cfg["hag"]
    meta = json.load(open(f"{C}/meta.json")); n = meta["n"]
    xyz = np.memmap(f"{C}/xyz.f4", np.float32, "r", shape=(n, 3))
    inst = np.asarray(np.memmap(f"{C}/inst.u4", np.uint32, "r", shape=(n,)))
    sem = np.asarray(np.memmap(f"{C}/sem.u1", np.uint8, "r", shape=(n,)))
    inten = np.asarray(np.memmap(f"{C}/intensity.u2", np.uint16, "r", shape=(n,)))
    hag = np.asarray(np.load(f"{H}/hag.npy", mmap_mode="r"))
    res = {r["tree"]: r for r in results}
    keys = [r[0] for r in rows]
    order = np.argsort(inst, kind="stable"); si = inst[order]
    bounds = np.searchsorted(si, np.arange(int(inst.max()) + 2))
    tf = Transformer.from_crs("EPSG:32618", "EPSG:4326", always_xy=True)
    rng = np.random.default_rng(0)
    table = []
    for k in keys:
        a, b = bounds[k], bounds[k + 1]
        if b <= a:
            continue
        idx = order[a:b]
        P = np.column_stack([np.asarray(xyz[idx, 0], np.float64) + meta["ox"],
                             np.asarray(xyz[idx, 1], np.float64) + meta["oy"],
                             np.asarray(xyz[idx, 2], np.float64) + meta["oz"]])
        hk = np.asarray(hag[idx], np.float64)
        lo = (hk >= 0.5) & (hk <= 2.5)
        cx, cy = (P[lo, :2].mean(0) if lo.sum() >= 8 else P[:, :2].mean(0))
        lon, lat = tf.transform(float(cx), float(cy))
        sel = rng.choice(len(idx), cfg["export_point_cap"], replace=False) if len(idx) > cfg["export_point_cap"] else np.arange(len(idx))
        sel.sort()
        hdr = laspy.LasHeader(point_format=1, version="1.2")
        hdr.offsets = [np.floor(P[sel, 0].min()), np.floor(P[sel, 1].min()), np.floor(P[sel, 2].min())]
        hdr.scales = [1e-4] * 3
        hdr.add_extra_dims([ExtraBytesParams(name="sem", type=np.uint8),
                            ExtraBytesParams(name="instance", type=np.uint32),
                            ExtraBytesParams(name="hag", type=np.float32)])
        las = laspy.LasData(hdr)
        las.x, las.y, las.z = P[sel, 0], P[sel, 1], P[sel, 2]
        las.intensity = np.asarray(inten[idx][sel], np.uint16)
        las.sem = np.asarray(sem[idx][sel], np.uint8)
        las.instance = np.full(len(sel), k, np.uint32)
        las.hag = hk[sel].astype(np.float32)
        r = res.get(k, {})
        conf = r.get("confidence") or "none"
        if r.get("status") == "refused":
            conf = "none"
        sub = conf if conf in ("high", "medium", "low") else "refused"
        d = os.path.join(cfg["trees"], sub); os.makedirs(d, exist_ok=True)
        las.write(os.path.join(d, f"tree_{k}.laz"))
        table.append(dict(id_arbol=k, latitud=round(lat, 7), longitud=round(lon, 7),
                          DBH_cm=r.get("d_ring_cm", ""), altura_m=r.get("top_m", ""),
                          confianza=conf, estado=r.get("status"),
                          arco=r.get("qc", {}).get("arc", ""), seed=r.get("seed_kind", ""),
                          inclinacion_deg=r.get("axis_tilt_deg", ""),
                          contrafuerte=r.get("buttress_suspected", ""),
                          n_pts=r.get("n_pts"), pts_banda=r.get("n_band"),
                          motivo=r.get("reason", ""), n_exportados=int(len(sel))))
    cols = list(table[0].keys())
    with open(os.path.join(cfg["out"], "tabla_arboles.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(table)
    print(f"[export] {len(table)} trees -> {cfg['trees']}/<confidence>/ ; table -> tabla_arboles.csv")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config"); ap.add_argument("--input"); ap.add_argument("--out")
    ap.add_argument("--steps", default="cache,ground,hag,canopy,measure,export")
    a = ap.parse_args()
    cfg = load_config(a)
    steps = set(a.steps.split(","))
    t0 = time.time()
    if "cache" in steps: build_cache(cfg)
    if "ground" in steps: build_ground(cfg)
    if "hag" in steps: build_hag(cfg)
    rows = canopy_summary(cfg) if "canopy" in steps else json.load(open(os.path.join(cfg["out"], "canopy_summary.json")))
    results = measure_all(cfg, rows) if "measure" in steps else json.load(open(os.path.join(cfg["out"], "results.json")))
    if "export" in steps: export(cfg, results, rows)
    print(f"[done] {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
