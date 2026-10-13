"""Stream segmentado_crop.laz once into flat memmaps of every scalar field we
need for the extended DBH/height pipeline.

Local float32 XYZ (metres, relative to a local origin) keeps 0.1 mm resolution
over a 200 m tile; raw integer XYZ would need 64-bit scaling per the header.
"""
import os, json, time, numpy as np, laspy

SRC = "/home/usuario/workspaces/NubePuntos/pointclouds/segmentado_crop.laz"
OUT = "/tmp/opencode/dbh/out/cache"
os.makedirs(OUT, exist_ok=True)

f = laspy.open(SRC)
h = f.header
n = h.point_count
ox, oy, oz = h.offsets
print("points", n, "offsets", ox, oy, oz, "scales", h.scales, flush=True)

xyz = np.memmap(f"{OUT}/xyz.f4", dtype=np.float32, mode="w+", shape=(n, 3))
inst = np.memmap(f"{OUT}/inst_fm.u2", dtype=np.uint16, mode="w+", shape=(n,))
sat_inst = np.memmap(f"{OUT}/inst_sat.u2", dtype=np.uint16, mode="w+", shape=(n,))
fm = np.memmap(f"{OUT}/sem_fm.u1", dtype=np.uint8, mode="w+", shape=(n,))
sat_sem = np.memmap(f"{OUT}/sem_sat.u1", dtype=np.uint8, mode="w+", shape=(n,))
inten = np.memmap(f"{OUT}/intensity.u2", dtype=np.uint16, mode="w+", shape=(n,))
rn = np.memmap(f"{OUT}/return_number.u1", dtype=np.uint8, mode="w+", shape=(n,))
nr = np.memmap(f"{OUT}/number_of_returns.u1", dtype=np.uint8, mode="w+", shape=(n,))
cloud = np.memmap(f"{OUT}/cloud_index.u1", dtype=np.uint8, mode="w+", shape=(n,))
score = np.memmap(f"{OUT}/score_fm.f4", dtype=np.float32, mode="w+", shape=(n,))

i = 0
t0 = time.time()
for pts in f.chunk_iterator(8_000_000):
    m = len(pts)
    sl = slice(i, i + m)
    xyz[sl, 0] = np.asarray(pts.x, dtype=np.float64) - ox
    xyz[sl, 1] = np.asarray(pts.y, dtype=np.float64) - oy
    xyz[sl, 2] = np.asarray(pts.z, dtype=np.float64) - oz
    inst[sl] = np.asarray(pts["PredInstance_FM"], dtype=np.uint16)
    sat_inst[sl] = np.asarray(pts["PredInstance_SAT"], dtype=np.uint16)
    fm[sl] = np.asarray(pts["PredSemantic_FM"], dtype=np.uint8)
    sat_sem[sl] = np.asarray(pts["PredSemantic_SAT"], dtype=np.uint8)
    inten[sl] = np.asarray(pts.intensity, dtype=np.uint16)
    rn[sl] = np.asarray(pts.return_number, dtype=np.uint8)
    nr[sl] = np.asarray(pts.number_of_returns, dtype=np.uint8)
    cloud[sl] = np.asarray(pts["Original cloud index"], dtype=np.uint8)
    score[sl] = np.asarray(pts["PredScore_FM"], dtype=np.float32)
    i += m
    print(f"  {i}/{n}  {time.time()-t0:.0f}s", flush=True)
f.close()
for a in (xyz, inst, sat_inst, fm, sat_sem, inten, rn, nr, cloud, score):
    a.flush()
assert i == n, (i, n)
json.dump(dict(n=int(n), ox=float(ox), oy=float(oy), oz=float(oz),
               scales=[float(s) for s in h.scales]), open(f"{OUT}/meta.json", "w"), indent=1)
print("done", time.time() - t0)
