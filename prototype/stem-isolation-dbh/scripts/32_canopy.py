"""Per-instance canopy summary + intensity calibration for the wood proxy."""
import json, numpy as np

C = "/tmp/opencode/dbh/out/cache"
H = "/tmp/opencode/dbh/out/hag"
n = json.load(open(f"{C}/meta.json"))["n"]

inst = np.asarray(np.memmap(f"{C}/inst_fm.u2", dtype=np.uint16, mode="r", shape=(n,)))
fm = np.asarray(np.memmap(f"{C}/sem_fm.u1", dtype=np.uint8, mode="r", shape=(n,)))
inten = np.asarray(np.memmap(f"{C}/intensity.u2", dtype=np.uint16, mode="r", shape=(n,)))
rn = np.asarray(np.memmap(f"{C}/return_number.u1", dtype=np.uint8, mode="r", shape=(n,)))
nr = np.asarray(np.memmap(f"{C}/number_of_returns.u1", dtype=np.uint8, mode="r", shape=(n,)))
cloud = np.asarray(np.memmap(f"{C}/cloud_index.u1", dtype=np.uint8, mode="r", shape=(n,)))
hag = np.asarray(np.load(f"{H}/hag.npy", mmap_mode="r"))

NI = int(inst.max()) + 1
npts = np.bincount(inst, minlength=NI)
top = np.full(NI, -1e9)
np.maximum.at(top, inst, hag)
band = np.bincount(inst[(hag >= 1.15) & (hag <= 1.45) & (inst > 0)], minlength=NI)
wood = np.bincount(inst[(hag >= 1.15) & (hag <= 1.45) & (inst > 0) & (fm == 1)], minlength=NI)

canopy = top > 12
canopy[0] = False
print(f"instances {NI-1}, canopy(top>12m) {int(canopy.sum())}")
print(f"canopy >=1 band pt  {int((canopy & (band>=1)).sum())}")
print(f"canopy >=20 band pt {int((canopy & (band>=20)).sum())}")
print(f"canopy >=50 band pt {int((canopy & (band>=50)).sum())}")
print(f"canopy with 0 FM wood in band: {int((canopy & (band>=1) & (wood==0)).sum())}")

rows = [(int(k), float(top[k]), int(npts[k]), int(band[k]), int(wood[k])) for k in np.where(canopy)[0]]
rows.sort(key=lambda r: -r[3])
json.dump(rows, open("/tmp/opencode/dbh/out/canopy_summary.json", "w"))

# ---- intensity calibration: FM wood vs foliage, per scan ----
print("\nintensity medians (all points, by FM class and scan):")
for ci in sorted(set(np.unique(cloud).tolist())):
    s = cloud == ci
    print(f"  scan {ci}: n={int(s.sum())}")
    for cls, nm in [(0, "ground"), (1, "wood"), (2, "foliage")]:
        m = s & (fm == cls)
        if m.sum():
            q = np.percentile(inten[m], [10, 50, 90]).round(0)
            print(f"    {nm:8s} n={int(m.sum()):9d} p10/50/90 = {q}")

# solid-return proxy stats
solid = (nr == 1) | (rn == nr)
print("\nsolid-return fraction by FM class:",
      {nm: round(float(solid[fm == c].mean()), 3) for c, nm in [(0, "ground"), (1, "wood"), (2, "foliage")]})
