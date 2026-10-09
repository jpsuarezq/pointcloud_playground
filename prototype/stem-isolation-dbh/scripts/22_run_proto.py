"""Run the DBH prototype over a stratified sample of canopy instances."""
import numpy as np, json, sys, time
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE, OX, OY, N
import dbh_proto as D

HC = "/tmp/opencode/cloud-research/out/hag"
O = "/tmp/opencode/cloud-research/out/inst"
F = "/tmp/opencode/cloud-research/out/fig"

inst = np.asarray(np.memmap(f"{O}/instance.u2.npy", dtype=np.uint16, mode="r", shape=(N,)))
hag = np.asarray(np.load(f"{HC}/hag_smrf025.npy", mmap_mode="r"), dtype=np.float64)
ok = np.asarray(np.load(f"{HC}/ok_smrf025.npy", mmap_mode="r"))
fm = np.asarray(load_attr()["fm"])
X = load_xyz()

rows = json.load(open(f"{O}/canopy_summary.json"))
band = {r[0]: r[3] for r in rows}
top = {r[0]: r[1] for r in rows}
npts = {r[0]: r[2] for r in rows}

# stratified sample: trees with >=10 band points, spread across height
elig = [r for r in rows if r[3] >= 10 and r[2] >= 3000]
elig.sort(key=lambda r: r[1])
sample = []
nbins = 6
for b in range(nbins):
    lo = int(b * len(elig) / nbins); hi = int((b + 1) * len(elig) / nbins)
    chunk = elig[lo:hi]
    chunk.sort(key=lambda r: -r[3])            # best covered within the height band
    sample += [r[0] for r in chunk[:6]]
sample += [r[0] for r in rows if 5 <= r[3] < 10 and r[2] >= 3000][:6]
if 2957 not in sample:
    sample.append(2957)
sample = list(dict.fromkeys(sample))
print("sample trees:", len(sample))

# index groups by instance
t0 = time.time()
order = np.argsort(inst, kind="stable")
sorted_inst = inst[order]
bounds = np.searchsorted(sorted_inst, np.arange(int(inst.max()) + 2))
print("sorted in %.1fs" % (time.time() - t0))

out = []
for k in sample:
    a, b = bounds[k], bounds[k + 1]
    idx = order[a:b]
    P = np.column_stack([np.asarray(X["x"][idx], dtype=np.float64) * SCALE + OX,
                         np.asarray(X["y"][idx], dtype=np.float64) * SCALE + OY,
                         np.asarray(X["z"][idx], dtype=np.float64) * SCALE])
    h = hag[idx]; f = fm[idx]
    res = D.measure_tree(P, h, f, k)
    res["top_m"] = round(top[k], 1); res["n_total"] = int(npts[k]); res["n_band"] = int(band[k])
    out.append(res)
    print(f"  {k:5d} top {top[k]:5.1f} band {band[k]:5d} -> {res['status']:8s}"
          f" d_ring={res.get('d_ring_cm','-')} taper={res.get('d_taper_cm','-')}"
          f" tilt={res.get('axis_tilt_deg', float('nan'))}"
          f" arc={res.get('qc',{}).get('arc','-')} spread={res.get('span_spread_cm','-')}", flush=True)
    if k == 2957:
        np.savez(f"{F}/tree_2957.npz", P=P.astype(np.float32), h=h.astype(np.float32),
                 f=f.astype(np.uint8))

json.dump(out, open(f"{O}/proto_results.json", "w"), indent=1, default=str)
print("wrote", f"{O}/proto_results.json")
