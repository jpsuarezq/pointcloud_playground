"""Run the multi-candidate method over the whole canopy."""
import json, sys, time, numpy as np
sys.path.insert(0, "/home/usuario/workspaces/NubePuntos/prototype/stem-isolation-dbh/scripts")
from dbh_multi import measure_tree_multi

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

tau = {}
for ci in (0, 1):
    s = cloud == ci
    tau[ci] = float((np.median(inten[s & (fm == 1)]) + np.median(inten[s & (fm == 2)])) / 2)
TAU = (tau[0], tau[1])
print("tau:", tau, flush=True)

rows = json.load(open(f"{O}/canopy_summary.json"))
band = {r[0]: r[3] for r in rows}; top = {r[0]: r[1] for r in rows}
keys = [r[0] for r in rows]
print("canopy", len(keys), flush=True)

t0 = time.time()
order = np.argsort(inst, kind="stable")
si = inst[order]
bounds = np.searchsorted(si, np.arange(int(inst.max()) + 2))

out = []
for j, k in enumerate(keys):
    a, b = bounds[k], bounds[k + 1]
    idx = order[a:b]
    P = np.column_stack([np.asarray(xyz[idx, 0], dtype=np.float64) + meta["ox"],
                         np.asarray(xyz[idx, 1], dtype=np.float64) + meta["oy"],
                         np.asarray(xyz[idx, 2], dtype=np.float64) + meta["oz"]])
    res = measure_tree_multi(P, hag[idx], fm[idx], inten[idx], cloud[idx], TAU, k)
    res["top_m"] = round(top[k], 1); res["n_band"] = int(band[k])
    out.append(res)
    if j % 50 == 0:
        print(f"  {j}/{len(keys)}  {time.time()-t0:.0f}s", flush=True)

json.dump(out, open(f"{O}/results_multi.json", "w"), indent=1, default=str)
from collections import Counter
st = Counter(r["status"] for r in out)
conf = Counter(r.get("confidence", "-") for r in out if r["status"] == "measured")
print("multi:", dict(st), dict(conf), "%.0fs" % (time.time() - t0))
