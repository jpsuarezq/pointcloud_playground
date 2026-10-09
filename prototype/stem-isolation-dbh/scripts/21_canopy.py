"""Per-instance summary using PredInstance_FM + HAG (SMRF 0.25 m)."""
import numpy as np, json, sys
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_attr, N

HC = "/tmp/opencode/cloud-research/out/hag"
O = "/tmp/opencode/cloud-research/out/inst"
inst = np.asarray(np.memmap(f"{O}/instance.u2.npy", dtype=np.uint16, mode="r", shape=(N,)))
hag = np.asarray(np.load(f"{HC}/hag_smrf025.npy", mmap_mode="r"), dtype=np.float64)
ok = np.asarray(np.load(f"{HC}/ok_smrf025.npy", mmap_mode="r"))
fm = np.asarray(load_attr()["fm"])

NI = int(inst.max()) + 1
n_pts = np.bincount(inst, minlength=NI)
top = np.full(NI, -np.inf)
np.maximum.at(top, inst, hag)
band = np.bincount(inst[(hag >= 1.15) & (hag <= 1.45) & (inst > 0)], minlength=NI)
wood = np.bincount(inst[(hag >= 1.15) & (hag <= 1.45) & (inst > 0) & (fm == 1)], minlength=NI)

canopy = (top > 12) & (top < np.inf)
print(f"instances {NI-1}, canopy(top>12m) {int(canopy.sum())}")
print(f"canopy with >=1 band pt  {int((canopy & (band>=1)).sum())}")
print(f"canopy with >=20 band pt {int((canopy & (band>=20)).sum())}")
print(f"canopy with >=50 band pt {int((canopy & (band>=50)).sum())}")
print("top height pct of canopy:", np.percentile(top[canopy], [5, 25, 50, 75, 95]).round(1))

rows = [(int(k), float(top[k]), int(n_pts[k]), int(band[k]), int(wood[k])) for k in np.where(canopy)[0]]
rows.sort(key=lambda r: -r[3])
json.dump(rows, open(f"{O}/canopy_summary.json", "w"))
print("\ntop 25 by band points:\n inst   top_m   npts    band  wood")
for r in rows[:25]:
    print(f"{r[0]:5d} {r[1]:7.1f} {r[2]:8d} {r[3]:6d} {r[4]:5d}")
ks = [r[0] for r in rows]
if 2957 in ks:
    print("\nexemplar 2957:", rows[ks.index(2957)])
