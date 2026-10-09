"""Convert PDAL writers.gdal GeoTIFFs to .npy + write a JSON sidecar with the
geotransform, so the nksr env (no rasterio) can use them. Run with DeepForest python.
"""
import sys, json, os
import numpy as np
import rasterio

OUT = "/tmp/opencode/cloud-research/out/dem"
for name in sys.argv[1:]:
    p = f"{OUT}/{name}.tif"
    with rasterio.open(p) as ds:
        a = ds.read(1)
        tr = ds.transform
        meta = dict(name=name, shape=list(a.shape), crs=str(ds.crs),
                    res=[tr.a, -tr.e], origin=[tr.c, tr.f],
                    nodata=ds.nodata, dtype=str(a.dtype))
    np.save(f"{OUT}/{name}.npy", a.astype(np.float64))
    json.dump(meta, open(f"{OUT}/{name}.json", "w"), indent=1)
    print(name, a.shape, a.dtype, "finite", int(np.isfinite(a).sum()), "min/max",
          float(np.nanmin(a)), float(np.nanmax(a)), meta["origin"], meta["res"])