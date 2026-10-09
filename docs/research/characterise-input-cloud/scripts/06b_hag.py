"""Height above ground for every point, per candidate ground surface.

Usage: 06b_hag.py <surface> [<surface> ...]
Each surface writes out/hag/hag_<name>.npy (float32, metres) and out/hag/ok_<name>.npy.
"""
import numpy as np, os, sys, json, time
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, N, SCALE, OX, OY

G = "/tmp/opencode/cloud-research/out/grids"
O = "/tmp/opencode/cloud-research/out/dem"
HC = "/tmp/opencode/cloud-research/out/hag"
os.makedirs(HC, exist_ok=True)
meta = dict(l.split("=") for l in open(f"{G}/grid_2500.txt"))
RES = int(meta["res"]); X0 = int(meta["x0"]); Y0 = int(meta["y0"])
NX, NY = int(meta["nx"]), int(meta["ny"])

X = load_xyz()
x = np.asarray(X["x"], dtype=np.int64)
y = np.asarray(X["y"], dtype=np.int64)
z = np.asarray(X["z"], dtype=np.float64) * SCALE
fm = np.asarray(load_attr()["fm"])
del X
xr = x * SCALE + OX; yr = y * SCALE + OY
ix = np.clip(((x - X0) // RES).astype(np.int32), 0, NX - 1)
iy = np.clip(((y - Y0) // RES).astype(np.int32), 0, NY - 1)


def bilin(dem, rowf, colf):
    nx = dem.shape[1]
    i0 = np.floor(rowf).astype(np.int64); j0 = np.floor(colf).astype(np.int64)
    fx = rowf - i0; fy = colf - j0
    i0c = np.clip(i0, 0, dem.shape[0] - 1); i1c = np.clip(i0 + 1, 0, dem.shape[0] - 1)
    j0c = np.clip(j0, 0, nx - 1); j1c = np.clip(j0 + 1, 0, nx - 1)
    return (dem[i0c, j0c] * (1 - fx) * (1 - fy) + dem[i0c, j1c] * fx * (1 - fy) +
            dem[i1c, j0c] * (1 - fx) * fy + dem[i1c, j1c] * fx * fy)


GRIDS = {"gmin025": "dem_gmin025.npy", "pmf025": "dem_pmf025.npy", "fm0025": "dem_fm0025.npy",
         "gmn05": "dem_gmn05.npy", "open50": "dem_open50.npy", "open90": "dem_open90.npy",
         "cloth05": "dem_cloth05.npy"}
GRIDRES = {"gmin025": 2500, "pmf025": 2500, "fm0025": 2500, "gmn05": 5000,
           "open50": 5000, "open90": 5000, "cloth05": 5000}
TIFS = {"smrf025": "smrf025", "csf025": "csf025", "pdalpmf025": "pmf025", "csf050": "csf050"}

for name in sys.argv[1:]:
    t0 = time.time()
    if name in GRIDS:
        res = GRIDRES[name]
        gm = dict(l.split("=") for l in open(f"{G}/grid_{res}.txt"))
        gx0, gy0, gnx, gny = int(gm["x0"]), int(gm["y0"]), int(gm["nx"]), int(gm["ny"])
        dem = np.where(np.isfinite(np.load(f"{O}/{GRIDS[name]}")), np.load(f"{O}/{GRIDS[name]}"), np.nan)
        okgrid = np.isfinite(bilin(dem, (y - gy0) / res, (x - gx0) / res))
        hag = z - bilin(dem, (y - gy0) / res, (x - gx0) / res)
    elif name in TIFS:
        base = TIFS[name]
        j = json.load(open(f"{O}/{base}.json")); arr = np.load(f"{O}/{base}.npy")
        a = np.where(arr < -9000, np.nan, arr).astype(np.float64)
        # GeoTIFF row 0 is the NORTH edge -> gy measured downward from origin[1]
        gx = (xr - j["origin"][0]) / j["res"][0]; gy = (j["origin"][1] - yr) / j["res"][1]
        inb = (gx >= 0) & (gx < a.shape[1] - 1.001) & (gy >= 0) & (gy < a.shape[0] - 1.001)
        hap = bilin(a, gy, gx)
        rc = np.clip(gy, 0, a.shape[0] - 1).astype(np.int64)
        cc = np.clip(gx, 0, a.shape[1] - 1).astype(np.int64)
        okgrid = inb & np.isfinite(a[rc, cc])
        hag = z - hap
    elif name.startswith("coarse"):
        res = {"coarse1m": 10000, "coarse2m": 20000}[name]
        cm = dict(l.split("=") for l in open(f"{G}/grid_{res}.txt"))
        cx0, cy0 = int(cm["x0"]), int(cm["y0"])
        cnx, cny = int(cm["nx"]), int(cm["ny"])
        dem = np.load(f"{O}/dem_{'gmin1000' if res==10000 else 'gmin2000'}.npy")
        dem = np.where(np.isfinite(dem), dem, np.nan)
        gcx = (x - cx0) / res; gcy = (y - cy0) / res
        okgrid = np.isfinite(bilin(dem, gcy, gcx))
        hag = z - bilin(dem, gcy, gcx)
    else:
        raise SystemExit("unknown surface " + name)

    ok = okgrid & np.isfinite(hag) & (hag > -50) & (hag < 200)
    np.save(f"{HC}/hag_{name}.npy", hag.astype(np.float32))
    np.save(f"{HC}/ok_{name}.npy", ok)
    d = hag[ok]
    if d.size == 0:
        print(name, "EMPTY"); continue
    p = np.percentile(d, [0.5, 1, 5, 25, 50, 75, 95, 99, 99.5])
    print(f"{name:12s} cov {ok.mean()*100:6.2f}%  hag pct(0.5,1,5,25,50,75,95,99,99.5) "
          f"{np.round(p,2)}  frac<0 {np.mean(d<0)*100:5.2f}%  {time.time()-t0:.0f}s", flush=True)
print("done")