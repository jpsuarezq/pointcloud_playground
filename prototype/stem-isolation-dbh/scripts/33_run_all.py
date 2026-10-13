"""Run the extended method over the 43-tree sample (for comparison) and over the
whole canopy (674 instances)."""
import json, sys, time, numpy as np
sys.path.insert(0, "/home/usuario/workspaces/NubePuntos/prototype/stem-isolation-dbh/scripts")
from dbh_ext import measure_tree_ext

C = "/tmp/opencode/dbh/out/cache"
H = "/tmp/opencode/dbh/out/hag"
O = "/tmp/opencode/dbh/out"
meta = json.load(open(f"{C}/meta.json")); n = meta["n"]

xyz = np.memmap(f"{C}/xyz.f4", dtype=np.float32, mode="r", shape=(n, 3))
inst = np.asarray(np.memmap(f"{C}/inst_fm.u2", dtype=np.uint16, mode="r", shape=(n,)))
fm = np.asarray(np.memmap(f"{C}/sem_fm.u1", dtype=np.uint8, mode="r", shape=(n,)))
inten = np.asarray(np.memmap(f"{C}/intensity.u2", dtype=np.uint16, mode="r", shape=(n,)))
cloud = np.asarray(np.memmap(f"{C}/cloud_index.u1", dtype=np.uint8, mode="r", shape=(n,)))
hag = np.asarray(np.load(f"{H}/hag.npy", mmap_mode="r"))

# ---- per-scan intensity thresholds (proxy: midpoint of FM wood/foliage medians)
tau = {}
for ci in (0, 1):
    s = cloud == ci
    tw = np.median(inten[s & (fm == 1)]); tf = np.median(inten[s & (fm == 2)])
    tau[ci] = float((tw + tf) / 2)
print("tau per scan:", tau, flush=True)
TAU = (tau[0], tau[1])

rows = json.load(open(f"{O}/canopy_summary.json"))          # (id, top, npts, band, wood)
band = {r[0]: r[3] for r in rows}
top = {r[0]: r[1] for r in rows}

# sample = the prototype's 43 trees
elig = [r for r in rows if r[3] >= 10 and r[2] >= 3000]
elig.sort(key=lambda r: r[1])
sample = []
for b in range(6):
    lo = int(b * len(elig) / 6); hi = int((b + 1) * len(elig) / 6)
    chunk = sorted(elig[lo:hi], key=lambda r: -r[3])
    sample += [r[0] for r in chunk[:6]]
sample += [r[0] for r in rows if 5 <= r[3] < 10 and r[2] >= 3000][:6]
if 2957 not in sample:
    sample.append(2957)
sample = list(dict.fromkeys(sample))
allc = [r[0] for r in rows]

print("sample", len(sample), "canopy", len(allc), flush=True)
t0 = time.time()
order = np.argsort(inst, kind="stable")
si = inst[order]
bounds = np.searchsorted(si, np.arange(int(inst.max()) + 2))
print("sorted %.1fs" % (time.time() - t0), flush=True)


def run(keys, tag):
    out = []
    for k in keys:
        a, b = bounds[k], bounds[k + 1]
        idx = order[a:b]
        P = np.column_stack([np.asarray(xyz[idx, 0], dtype=np.float64) + meta["ox"],
                             np.asarray(xyz[idx, 1], dtype=np.float64) + meta["oy"],
                             np.asarray(xyz[idx, 2], dtype=np.float64) + meta["oz"]])
        res = measure_tree_ext(P, hag[idx], fm[idx], inten[idx], cloud[idx], TAU, k)
        res["top_m"] = round(top[k], 1); res["n_band"] = int(band[k])
        out.append(res)
    json.dump(out, open(f"{O}/results_{tag}.json", "w"), indent=1, default=str)
    from collections import Counter
    st = Counter(r["status"] for r in out)
    conf = Counter(r.get("confidence", "-") for r in out if r["status"] == "measured")
    print(f"[{tag}] n={len(out)} {dict(st)} conf={dict(conf)}")
    return out


if "sample" in sys.argv:
    run(sample, "sample")
if "all" in sys.argv:
    run(allc, "all")
print("done", round(time.time() - t0, 1), "s")
