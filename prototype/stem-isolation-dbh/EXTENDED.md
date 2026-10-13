# Extended run: every canopy tree, individualised, with a per-tree table

Built on top of the `README.md` prototype (ticket #11) to answer the driver's
follow-up: use the **other scalar fields** to recover the trees the prototype
refused, individualise **every** canopy tree into its own file, and emit a final
per-tree table.

Throwaway prototype. No accuracy claim, no ground truth, no tuning to truth.

## What changed

The v2 prototype measured the ring from **all** points near the axis — FM wood
was only used to **seed** the axis. So a tree FM left unlabelled is recoverable
if the axis can be seeded another way. Two additions:

**Multi-candidate seeding** (`scripts/dbh_multi.py`). The dominant failure was
`no_points_near_axis`: one seed landed in a gap, so the ring saw nothing. We now
generate several candidate axes, refine each, and keep the best by QC:

1. **FM-wood column** (v2 behaviour).
2. **Intensity-proxy column** — low intensity as a wood proxy. Threshold is the
   per-scan midpoint of the FM wood/foliage medians (scan 0 ≈ 10 880, scan 1 ≈
   11 008). Wood median intensity ~8.5–9.2 k vs foliage ~12.8–13.3 k, but the
   distributions overlap, so this is a *proxy*.
3. **Base anchor** — median XY of the lower-stem points (0.5–2.5 m), vertical.
4. **Band seed** — densest DBSCAN cluster of the 1.15–1.45 m points nearest the
   base, vertical.
5. **Refined axis** for each candidate.

Candidates are ranked by `(valid, ok, narrowest slab, arc)`, so a valid narrow
ring always wins. **Slab discipline:** the ±0.15 m slab is tried first and only
widened to ±0.30 m when the narrow slab yields no valid ring — 72% of
measurements stay on the narrow band.

Every measurement carries a **confidence** (`high` = FM-wood seed, arc ≥ 0.35,
d < 1 m; `medium` = proxy/base/band seed passing the same QC; `low` = emitted
anyway, flagged). Fields that did **not** help: `number_of_returns` /
`return_number` (wood "solid-return" fraction 0.74 ≈ foliage 0.77);
`PredSemantic_SAT` (binary, not a wood class). `Original cloud index` (two scans)
is cached and available.

## Inputs rebuilt from scratch

The old `/tmp/opencode/cloud-research/` cache no longer exists. Rebuilt here:

- `30_cache_fields.py` — one streaming pass of `segmentado_crop.laz` → memmaps of
  xyz, both instance partitions, both semantics, intensity, returns, scan index,
  score (1.9 GB, 16 s).
- `31_hag.py` — SMRF ground (`pdal filters.smrf`, 1 m 22 s, 6.4 GB RSS) →
  0.25 m grid-minimum DEM (hole-filled ≤ 2 m, 3×3 median) → bilinear HAG per
  point.
- `32_canopy.py` — per-instance canopy summary + intensity calibration.
- `33_run_all.py` — single-seed run (comparison).
- `35_run_multi.py` — multi-candidate run (the deliverable).
- `34_export.py` — one `.laz` per tree + the table.

## Results

**674 canopy instances** (FM top HAG > 12 m); 361 have ≥ 1 point in the 1.3 m
band; **142 have zero FM-labelled wood in the band**.

| Run | measured | refused | high | medium | low |
|---|---|---|---|---|---|
| 43-tree sample, v2 prototype | 14 | 29 | — | — | — |
| 43-tree sample, extended | 26 | 17 | 9 | 10 | 7 |
| whole canopy, single seed | 166 | 508 | 83 | 30 | 53 |
| **whole canopy, multi-candidate** | **312** | 362 | 78 | 193 | 41 |

Of the 362 refusals, **312 have zero points in the 1.3 m band** (genuinely
unmeasurable) and only **1 has ≥ 10**. New seeds that recovered trees: band 120,
wood 94, base 63, intensity 35.

Measured DBH cm p10/50/90: **high 16.9/26.0/46.1**, **medium 16.4/32.8/55.4**,
**low 13.3/33.5/79.0**. The `low` tail (up to 109 cm) is the visible failure.

## Deliverables

- **Table**: `results/tabla_arboles.csv` — `id_arbol, latitud, longitud, DBH_cm,
  altura_m, confianza, estado, arco, seed, inclinacion_deg, contrafuerte, n_pts,
  pts_banda, motivo, n_exportados`. Coordinates WGS84 (EPSG:4326).
- **Per-tree clouds**: `pointclouds/trees/tree_<id>.laz` — 674 files, EPSG:32618,
  extra dims `sem_fm`, `inst_fm`, `hag`; decimated to ≤ 400 k pts for viewing.
  (Git-ignored via `pointclouds`.)
- **Full results**: `results/results_multi.json`, `results/results_all.json`,
  `results/results_sample_ext.json`.

### Preview (tallest 12 measured)

| Id | lat | lon | DBH cm | altura m | conf | arco | seed |
|---|---|---|---|---|---|---|---|
| 2957 | 6.8628922 | -73.4127252 | 53.4 | 40.6 | medium | 0.38 | intensity |
| 2778 | 6.8623475 | -73.4128364 | 28.8 | 39.9 | medium | 1.0 | intensity |
| 2885 | 6.8630426 | -73.4127205 | 25.8 | 39.3 | high | 0.69 | wood |
| 2906 | 6.8637 | -73.4127826 | 53.4 | 36.2 | medium | 0.38 | intensity |
| 2265 | 6.8638718 | -73.4129301 | 19.0 | 36.1 | high | 0.56 | wood |
| 3103 | 6.8619142 | -73.4126482 | 95.0 | 36.1 | low | 0.25 | wood |
| 2666 | 6.862316 | -73.4128257 | 104.4 | 35.7 | low | 0.38 | base |
| 2238 | 6.8629651 | -73.412963 | 36.2 | 35.6 | medium | 0.38 | band |
| 3505 | 6.8641423 | -73.412437 | 22.8 | 35.4 | medium | 0.5 | intensity |
| 2367 | 6.8635548 | -73.4129685 | 67.3 | 35.2 | low | 0.25 | band |
| 2600 | 6.8640116 | -73.4128243 | 22.2 | 35.0 | medium | 0.56 | base |
| 3290 | 6.8635285 | -73.4125162 | 29.3 | 34.9 | medium | 0.38 | base |

## Caveats

- Confidence bounds *consistency and evidence*, never accuracy — there is no
  ground truth and none was used.
- The intensity proxy overlaps with foliage, and band-seeded axes can circle a
  fork or buttress: `medium` is *plausible*, not *verified*. `low` rows must not
  be read as measurements.
- The DEM is SMRF at 0.25 m rebuilt here; band membership may differ slightly
  from the earlier characterisation surface.
- 10 instances with thousands of band points carry zero FM wood — they are now
  measured, but "measured" is not "a trunk".
