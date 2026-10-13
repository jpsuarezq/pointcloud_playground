"""Export one .laz per canopy tree (individualised, viewable one by one) and the
final per-tree table with lat/lon, DBH, height and QC.

Coordinates written in EPSG:32618 (UTM 18N); table adds WGS84 lat/lon.
"""
import os, json, time, numpy as np, laspy
from laspy import ExtraBytesParams
from pyproj import Transformer

C = "/tmp/opencode/dbh/out/cache"
H = "/tmp/opencode/dbh/out/hag"
O = "/tmp/opencode/dbh/out"
TREES = "/home/usuario/workspaces/NubePuntos/pointclouds/trees"
os.makedirs(TREES, exist_ok=True)
CAP = 400_000
meta = json.load(open(f"{C}/meta.json")); n = meta["n"]

xyz = np.memmap(f"{C}/xyz.f4", dtype=np.float32, mode="r", shape=(n, 3))
inst = np.asarray(np.memmap(f"{C}/inst_fm.u2", dtype=np.uint16, mode="r", shape=(n,)))
fm = np.asarray(np.memmap(f"{C}/sem_fm.u1", dtype=np.uint8, mode="r", shape=(n,)))
inten = np.asarray(np.memmap(f"{C}/intensity.u2", dtype=np.uint16, mode="r", shape=(n,)))
hag = np.asarray(np.load(f"{H}/hag.npy", mmap_mode="r"))

res = {r["tree"]: r for r in json.load(open(f"{O}/results_multi.json"))}
rows = json.load(open(f"{O}/canopy_summary.json"))
keys = [r[0] for r in rows]

t0 = time.time()
order = np.argsort(inst, kind="stable")
si = inst[order]
bounds = np.searchsorted(si, np.arange(int(inst.max()) + 2))

tf = Transformer.from_crs("EPSG:32618", "EPSG:4326", always_xy=True)
rng = np.random.default_rng(0)

table = []
for k in keys:
    a, b = bounds[k], bounds[k + 1]
    if b <= a:
        continue
    idx = order[a:b]
    P = np.column_stack([np.asarray(xyz[idx, 0], dtype=np.float64) + meta["ox"],
                         np.asarray(xyz[idx, 1], dtype=np.float64) + meta["oy"],
                         np.asarray(xyz[idx, 2], dtype=np.float64) + meta["oz"]])
    # centroid of the lower stem (0.5-2.5 m) for geolocation
    hk = np.asarray(hag[idx], dtype=np.float64)
    lo = (hk >= 0.5) & (hk <= 2.5)
    cx, cy = (P[lo, :2].mean(0) if lo.sum() >= 8 else P[:, :2].mean(0))
    lon, lat = tf.transform(float(cx), float(cy))

    # decimate for viewing
    if len(idx) > CAP:
        sel = rng.choice(len(idx), CAP, replace=False); sel.sort()
    else:
        sel = np.arange(len(idx))
    Pw = P[sel]; hw = hk[sel]

    hdr = laspy.LasHeader(point_format=1, version="1.2")
    hdr.offsets = [np.floor(Pw[:, 0].min()), np.floor(Pw[:, 1].min()), np.floor(Pw[:, 2].min())]
    hdr.scales = [1e-4, 1e-4, 1e-4]
    hdr.add_extra_dims([
        ExtraBytesParams(name="sem_fm", type=np.uint8),
        ExtraBytesParams(name="inst_fm", type=np.uint16),
        ExtraBytesParams(name="hag", type=np.float32),
    ])
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = Pw[:, 0], Pw[:, 1], Pw[:, 2]
    las.intensity = np.asarray(inten[idx][sel], dtype=np.uint16)
    las.sem_fm = np.asarray(fm[idx][sel], dtype=np.uint8)
    las.inst_fm = np.full(len(sel), k, dtype=np.uint16)
    las.hag = hw.astype(np.float32)
    las.write(f"{TREES}/tree_{k}.laz")

    r = res.get(k, {})
    table.append(dict(
        id_arbol=k, latitud=round(lat, 7), longitud=round(lon, 7),
        DBH_cm=r.get("d_ring_cm", ""), altura_m=r.get("top_m", ""),
        confianza=r.get("confidence", "refused"), estado=r.get("status"),
        arco=r.get("qc", {}).get("arc", ""), seed=r.get("seed_kind", ""),
        inclinacion_deg=r.get("axis_tilt_deg", ""), contrafuerte=r.get("buttress_suspected", ""),
        n_pts=r.get("n_pts"), pts_banda=r.get("n_band"),
        motivo=r.get("reason", ""), n_exportados=int(len(sel))))
    if (time.time() - t0) % 30 < 1:
        print(f"  {k} ({len(table)}/{len(keys)}) {time.time()-t0:.0f}s", flush=True)

import csv
cols = list(table[0].keys())
with open(f"{O}/tabla_arboles.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(table)
print("trees written:", len(table), "->", TREES)
print("table ->", f"{O}/tabla_arboles.csv", "%.0fs" % (time.time() - t0))
