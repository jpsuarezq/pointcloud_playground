# Extended run: every canopy tree, individualised, with a per-tree table

Built on top of the `README.md` prototype (ticket #11) to answer the driver's
follow-up: use the **other scalar fields** to recover the trees the prototype
refused, individualise **every** canopy tree into its own file, and emit a final
per-tree table.

Throwaway prototype. No accuracy claim, no ground truth, no tuning to truth.

## What changed

The v2 prototype measured the ring from **all** points near the axis — FM wood
was only used to **seed** the axis. So a tree FM left unlabelled is recoverable
if the axis can be seeded another way. The extended method
(`scripts/dbh_ext.py`) adds a fallback chain:

1. **FM-wood column search** (v2 behaviour).
2. **Intensity-proxy column search** — low intensity is a wood proxy. The
   threshold is the per-scan midpoint of the FM wood/foliage medians
   (scan 0 ≈ 10 880, scan 1 ≈ 11 008). Wood median intensity ~8.5–9.2 k vs
   foliage ~12.8–13.3 k, but the distributions overlap, so this is a *proxy*.
3. **Base anchor** — median XY of the lower-stem points (0.5–2.5 m), vertical
   axis. Uses no labels at all.

Every measurement carries a **confidence** (`high` = FM-wood seed, arc ≥ 0.35,
d < 1 m; `medium` = proxy/base seed passing the same QC; `low` = emitted anyway,
flagged). Fields that did **not** help: `number_of_returns`/`return_number`
(wood "solid-return" fraction 0.74 ≈ foliage 0.77), `PredSemantic_SAT` (binary,
not a wood class). `Original cloud index` (two scans) is cached and available.

## Inputs rebuilt from scratch

The old `/tmp/opencode/cloud-research/` cache no longer exists. Rebuilt here:

- `30_cache_fields.py` — one streaming pass of `segmentado_crop.laz` → memmaps of
  xyz, both instance partitions, both semantics, intensity, returns, scan index,
  score (1.9 GB, 16 s).
- `31_hag.py` — SMRF ground (`pdal filters.smrf`, 1 m 22 s, 6.4 GB RSS) →
  0.25 m grid-minimum DEM (hole-filled ≤ 2 m, 3×3 median) → bilinear HAG per
  point.
- `32_canopy.py` — per-instance canopy summary + intensity calibration.
- `33_run_all.py` — run over the 43-tree sample and the whole canopy.
- `34_export.py` — one `.laz` per tree + the table.

## Results

**674 canopy instances** (FM top HAG > 12 m); 361 have ≥ 1 point in the 1.3 m
band; **142 have zero FM-labelled wood in the band**.

| Set | measured | refused | high | medium | low |
|---|---|---|---|---|---|
| 43-tree sample (old) | 14 | 29 | — | — | — |
| 43-tree sample (extended) | **26** | 17 | 9 | 10 | 7 |
| whole canopy (674) | **166** | 508 | 83 | 30 | 53 |

Recovered by the new seeds: **43 intensity + 13 base = 56 trees** that the
FM-wood-only prototype refused. Of the 508 refusals, **313 have zero points in
the band** (genuinely unmeasurable) and 195 have ≥ 1 (103 with ≥ 10).

Measured DBH cm p10/25/50/75/90 = **15.8 / 22.7 / 29.8 / 46.1 / 65.0**; heights
16 / 25 / 33 m. `low`-confidence diameters are the implausible tail (median
46 cm, up to 116 cm) — visible, and flagged, as intended.

## Deliverables

- **Table**: `results/tabla_arboles.csv` — `id_arbol, latitud, longitud, DBH_cm,
  altura_m, confianza, estado, arco, seed, inclinacion_deg, contrafuerte, n_pts,
  pts_banda, motivo, n_exportados`. Coordinates WGS84 (EPSG:4326).
- **Per-tree clouds**: `pointclouds/trees/tree_<id>.laz` — 674 files, 314 MB,
  EPSG:32618, extra dims `sem_fm`, `inst_fm`, `hag`; decimated to ≤ 400 k pts
  for viewing. (Git-ignored via `pointclouds`.)
- **Full results**: `results/results_all.json`, `results/results_sample_ext.json`.

### Preview (tallest 15 measured)

| Id | lat | lon | DBH cm | altura m | conf | arco | seed |
|---|---|---|---|---|---|---|---|
| 2778 | 6.8623475 | -73.4128364 | 116.3 | 39.9 | low | 0.38 | wood |
| 2885 | 6.8630426 | -73.4127205 | 25.8 | 39.3 | high | 0.69 | wood |
| 2265 | 6.8638718 | -73.4129301 | 19.0 | 36.1 | high | 0.56 | wood |
| 2238 | 6.8629651 | -73.412963 | 100.0 | 35.6 | low | 0.25 | wood |
| 3505 | 6.8641423 | -73.412437 | 29.2 | 35.4 | high | 0.44 | wood |
| 3290 | 6.8635285 | -73.4125162 | 73.1 | 34.9 | low | 0.31 | wood |
| 2820 | 6.8643434 | -73.4127821 | 21.4 | 34.6 | high | 0.5 | wood |
| 2704 | 6.8634186 | -73.4126898 | 14.3 | 34.5 | high | 0.44 | wood |
| 1956 | 6.8638627 | -73.4130108 | 32.1 | 34.1 | high | 0.38 | wood |
| 2559 | 6.8622111 | -73.4128279 | 45.1 | 33.7 | high | 0.44 | wood |
| 3205 | 6.8629479 | -73.4126239 | 20.8 | 33.6 | high | 0.44 | wood |
| 2006 | 6.862235 | -73.4130184 | 33.8 | 33.6 | high | 0.88 | wood |
| 2963 | 6.8633067 | -73.412698 | 44.8 | 33.1 | high | 0.5 | wood |
| 2650 | 6.8614672 | -73.4128093 | 17.2 | 32.9 | high | 0.69 | wood |
| 3226 | 6.8640262 | -73.4125591 | 25.0 | 32.9 | medium | 0.44 | base |

## Caveats

- Confidence bounds *consistency and evidence*, never accuracy — there is no
  ground truth and none was used.
- The intensity proxy overlaps with foliage; `low` rows must not be read as
  measurements.
- The DEM is SMRF at 0.25 m rebuilt here; band membership may differ slightly
  from the earlier characterisation surface.
- The sample→canopy coverage jump is method, not a claim that all 166 are real
  trunks (10 instances with thousands of band points carry zero FM wood).
