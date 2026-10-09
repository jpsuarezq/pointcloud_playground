"""Stream the LAZ once, saving PredInstance_FM per point in cache order.

Output: out/inst/instance.u2.npy (uint16, 72.16M) and a per-instance summary.
"""
import os, time, numpy as np, laspy

SRC = "/home/usuario/workspaces/NubePuntos/pointclouds/segmentado_crop.laz"
OUT = "/tmp/opencode/cloud-research/out/inst"
os.makedirs(OUT, exist_ok=True)

f = laspy.open(SRC)
n = f.header.point_count
inst = np.memmap(f"{OUT}/instance.u2.npy", dtype=np.uint16, mode="w+", shape=(n,))
i = 0
t0 = time.time()
for pts in f.chunk_iterator(8_000_000):
    m = len(pts)
    inst[i:i + m] = np.asarray(pts["PredInstance_FM"], dtype=np.uint16)
    i += m
    print(f"  {i}/{n}  {time.time()-t0:.0f}s", flush=True)
f.close()
inst.flush()
assert i == n, (i, n)
print("done", time.time() - t0)
u, c = np.unique(np.asarray(inst), return_counts=True)
print("instances:", len(u), "max", u.max(), "points with inst 0:", int(c[u == 0].sum()))
print("points per instance pct:", np.percentile(c[c > 0], [10, 50, 90, 99]).round(0))
