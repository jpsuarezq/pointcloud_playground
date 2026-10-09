"""Shared helpers over the flat cache in out/cache."""
import numpy as np

CACHE = "/tmp/opencode/cloud-research/out/cache"
SCALE = 1e-4
OX, OY, OZ = 675000.0, 758000.0, 0.0
N = 72162549


def load_xyz():
    """int32 raw LAS coords (x,y,z)."""
    dt = np.dtype([("x", "<i4"), ("y", "<i4"), ("z", "<i4")])
    return np.memmap(f"{CACHE}/xyz.i4", dtype=dt, mode="r", shape=(N,))


def load_attr():
    dt = np.dtype([("fm", "u1"), ("cls", "u1"), ("intensity", "<u2")])
    return np.memmap(f"{CACHE}/attr.u1", dtype=dt, mode="r", shape=(N,))


def load_rgb():
    dt = np.dtype([("r", "<u2"), ("g", "<u2"), ("b", "<u2")])
    return np.memmap(f"{CACHE}/rgb.u2", dtype=dt, mode="r", shape=(N,))


def grid_index(x, y, res, x0=None, y0=None, nx=None, ny=None):
    """Return (ix, iy, x0, y0, nx, ny). Local coords in LAS raw units."""
    if x0 is None:
        x0 = int(np.floor(x.min() / res) * res)
        y0 = int(np.floor(y.min() / res) * res)
    if nx is None:
        nx = int(np.floor(x.max() / res)) - int(round(x0 / res)) + 1
    if ny is None:
        ny = int(np.floor(y.max() / res)) - int(round(y0 / res)) + 1
    ix = ((x.astype(np.int64) - x0) // res).astype(np.int32)
    iy = ((y.astype(np.int64) - y0) // res).astype(np.int32)
    np.clip(ix, 0, nx - 1, out=ix)
    np.clip(iy, 0, ny - 1, out=iy)
    return ix, iy, x0, y0, nx, ny


def cell_min_max(ix, iy, nx, ny):
    """min/max z per cell via np.minimum.at / maximum.at."""
    flat = iy.astype(np.int64) * nx + ix
    nc = nx * ny
    zmin = np.full(nc, np.inf)
    zmax = np.full(nc, -np.inf)
    np.minimum.at(zmin, flat, ix.astype(np.float64))  # placeholder, replaced below
    return zmin, zmax


def local_utm(x_raw, y_raw):
    return x_raw * SCALE + OX, y_raw * SCALE + OY


def z_real(z_raw):
    return z_raw * SCALE + OZ