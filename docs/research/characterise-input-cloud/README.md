# Characterising `segmentado_RGB_color.las`

Date: 2026-10-09. Branch: `research/characterise-input-cloud` (research only — never merged).
Input: `pointclouds/segmentado_RGB_color.las`, 72,162,549 points, LAS 1.2 point format 2, EPSG:32618.

Everything below is measured by the scripts in `scripts/`, in this order. Raw numbers are in
`results/`, figures in `figures/`. Nothing is estimated unless explicitly flagged as such.

---

## 0. What is actually in the file

| property | measured |
|---|---|
| points | 72,162,549 |
| header bbox | 199.88 m (X) x 407.69 m (Y), Z 860.43 - 975.55 m |
| **LAS `classification` byte** | **0 for 100% of points — the file has never been ground-classified** |
| `PredSemantic_FM` | int32 extra dim, values {0,1,2}: 12,332,906 (17.09%) / 2,082,633 (2.886%) / 57,747,010 (80.02%) |
| RGB | **B is 0 for 100% of points.** R > 0 for only 19.5% of points (FM0 97.9%, FM1 95.5%, FM2 0.0005%). Values are 8-bit channels shifted `<< 8` |
| occupied footprint | 47.6% of the bbox (the tile is a diamond, not a rectangle) |
| analysis area | 151,021 usable 0.5 m cells = **37,755 m² (3.78 ha)**; 70,810,287 valid points → **1,876 pts/m²** |
| relief (SMRF DEM) | 84.9 m |
| slope (SMRF DEM, 2.5 m smoothing) | median **18.8°**, p75 29.4°, p90 55.5°, p99 82.9° |

**Discrepancy to flag:** the ticket states median slope 27.7°, p90 39.8°, p99 54.2°. I measure
18.8° / 55.5° / 82.9° on the analysis area. The *shape* of the distribution differs (my tail is much
fatter). I have not reproduced the 27.7° figure with any slope estimator I tried (cell-to-cell
gradient at 0.5 m, 1.5 m, 2.5 m, 4.5 m, 8.5 m smoothing all give 17-19° median). Whichever is
right, the terrain is steep enough that both matter.

**Tooling corrections to the ticket's assumptions:**
* `pdal` 2.3.0 here **does** have `filters.csf` and `filters.pmf`, and `filters.smrf` *does* take
  options (`cell`, `slope`, `scalar`, `threshold`, `window`).
* `filters.csf` at `resolution: 0.25` on this tile **OOM-kills at 9.2 GB** after 5m47s. At
  `resolution: 0.5` it completes in 4m19s / 9.2 GB.
* **`writers.gdal` writes approximately the MINIMUM z per cell, not the mean.** Verified against
  per-cell min / mean / max on a 30x40 m crop (tif matches min to 0.1 m, mean by 17 m). Consequence:
  a DEM produced via `writers.gdal` is a *lower-envelope raster*, never an interpolated surface.
  This also explains the next point.
* `filters.smrf` and `filters.pmf` produce **bit-identical** ground rasters on this tile via
  `writers.gdal` (verified with an independent sequential re-run of PMF).

---

## 1. The vertical density profile and the 1.0-1.6 m blind spot — **CONFIRMED**

Density is normalised per 0.5 m ground cell, per metre of height above ground (HAG), over the
37,755 m² analysis area. The 0.25 m PDAL SMRF surface is used as the ground: it is the only surface
here with **0.11% of points below it** and 95.8% coverage, i.e. a genuine lower envelope.

| HAG (m) | pts / m² / m | | HAG (m) | pts / m² / m |
|---|---|---|---|---|
| 0.00-0.05 | 1096 | | 1.15-1.20 | 38 |
| 0.05-0.10 | 1780 | | 1.25-1.30 | 42 |
| 0.10-0.15 | 1389 | | 1.35-1.40 | 50 |
| 0.20-0.25 | 769 | | 1.45-1.50 | 59 |
| 0.30-0.35 | 452 | | 1.55-1.60 | 69 |
| 0.50-0.55 | 182 | | 1.75-1.80 | 95 |
| 0.75-0.80 | 65 | | 1.95-2.00 | 128 |
| 0.85-0.90 | 48 | | 2.15-2.20 | 163 |
| 0.95-1.00 | 38 | | 2.45-2.50 | 200 |
| **1.05-1.10 (min)** | **37** | | **2.65-2.70 (max)** | **208** |
| 2.85-2.90 | 200 | | 2.95-3.00 | 191 |

The trough is 4.9x below the immediate near-ground layer (0.50-0.55 m, 182 pts/m²/m).

**Minimum 37 pts/m²/m at 1.05-1.10 m; 5.73x below the 2.65 m canopy peak (208 pts/m²/m).** After slope
correction (mean cos θ = 0.916 over the analysis cells) the numbers all scale by 1.09 — the ratio is
unchanged. So the previously-reported 5-6x "blind spot" **reproduces against a trustworthy ground
surface**.

The trough is strongest on flat ground and weakens with slope:

| slope stratum | cells | trough HAG | d_min | peak HAG | d_max | ratio |
|---|---|---|---|---|---|---|
| 0-15° | 44,976 | 0.95 m | 22.2 | 2.60 m | 221.7 | **9.98** |
| 15-25° | 55,246 | 1.00 m | 28.3 | 2.60 m | 211.2 | **7.45** |
| 25-35° | 35,506 | 1.15 m | 48.6 | 2.95 m | 222.9 | 4.59 |
| 35-45° | 11,233 | 1.30 m | 73.6 | 2.65 m | 201.6 | 2.74 |
| 45°+ | 4,060 | — | — | — | — | too few points |

Stratifying by **canopy closure** (DSM - DTM per cell) is the sharpest result:

| closure stratum | area (m²) | trough HAG | d_min | peak HAG | d_max | ratio |
|---|---|---|---|---|---|---|
| < 2 m (ground directly visible) | 4,685 | — | 2.7 | 1.80 m | 60.0 | none, monotonic decay |
| 2-8 m | 9,886 | 1.00 m | 45.0 | 2.45 m | 402.8 | **8.95** |
| 8-15 m | 2,955 | 0.95 m | 24.2 | 2.45 m | 114.6 | 4.74 |
| > 15 m (closed canopy) | 21,057 | 0.95 m | 29.5 | 2.50 m | 160.5 | 5.45 |

The trough exists **only where the canopy is closed** and vanishes in open cells. That is the
signature of an opacity-limited capture, not of a height-dependent detector.

### Point budget at the 1.3 m DBH band

Measured on the full tile:

* band 1.15-1.45 m contains **533,105 points = 14.12 pts/m² of plot** over the 0.30 m band
  (47 pts/m²/m);
* median 16-NN radius in the band = **0.157 m** → surface point density **207 pts/m²**;
* band geometry: 41.2% linear (λ>0.6), 7.15% vertical (v>0.8), **2.67% woody (both)**, 9.77% leafy.

*Geometric* budget `n(d) = ρ_s · π · d · L`, with ρ_s = 207 pts/m² and L = 0.30 m:

| DBH | 10 cm | 15 cm | 20 cm | 30 cm | 50 cm | 80 cm |
|---|---|---|---|---|---|---|
| points on a 0.30 m segment | 19.5 | 29.3 | 39.0 | **58.5** | 97.5 | 156 |

*Direct* measurement (328 woody clusters detected tile-wide, single-linkage at 10 cm inside the
band; cluster diameter from its own point spread; count = band points inside radius d/2 + 3 cm):

| measured diameter | 3-7 cm | 7-10 cm | 10-13 cm | 13-16 cm | 16-20 cm |
|---|---|---|---|---|---|
| median points | 14 | 16 | 17 | 17 | 16 |

**Honest verdict:** the direct count is essentially **flat in diameter** (14-17 points) and no stem
wider than ~20 cm was detectable in the band at all. Fit `n = 771.7·d − 49.8` (R = 0.64) implies an
effective on-stem surface density of 819 pts/m², i.e. the 0.157 m median neighbourhood radius is set
by neighbouring foliage, not by the stem surface. **I cannot confirm the previously-stated d-scaling
of the point budget.** The 30/50/80 cm numbers above are extrapolations from ρ_s = 207 pts/m² and
should be treated as upper bounds. What *is* solid: a stem-sized object at 1.3 m in this cloud
carries **~15-20 points on a 0.30 m segment**, which is marginal but not hopeless.

---

## 2. Is the blind spot intrinsic, or an artefact of the ground surface? — **INTRINSIC**

Nine independent ground surfaces were built and each was asked for its vertical profile:

| surface | coverage | **% points below surface** | trough HAG | d_min | peak HAG | d_max | ratio |
|---|---|---|---|---|---|---|---|
| PDAL SMRF, cell 0.25 m | 95.83% | **0.11** | 1.05 m | 36.6 | 2.60 m | 209.3 | **5.73** |
| PDAL CSF, resolution 0.5 m | 95.84% | **0.02** | 1.30 m | 45.6 | 2.80 m | 205.2 | **4.50** |
| numpy grid-min of `PredSemantic_FM==0`, 0.25 m | 92.07% | 13.97 | none | — | 2.45 m | 157.8 | — |
| numpy slope-corrected opening 4.5 m, 0.5 m | 95.45% | 14.76 | none | — | 2.25 m | 173.6 | — |
| numpy slope-corrected opening 2.5 m, 0.5 m | 95.45% | 18.17 | none | — | 2.10 m | 158.9 | — |
| numpy grid-CSF, 0.5 m | 95.45% | 28.96 | none | — | 1.80 m | 99.1 | — |
| numpy grid-min, 0.5 m | 95.45% | 31.45 | none | — | 1.80 m | 88.4 | — |
| numpy PMF (Zhang 2003), 0.25 m | 95.86% | 35.86 | none | — | 1.80 m | 84.4 | — |
| numpy grid-min, 0.25 m (no morphology) | 95.86% | 36.04 | none | — | 1.80 m | 84.8 | — |

### Why the grid-minimum surfaces fail (and it is not a bug)

At 0.25 m the cell minimum is a **lower envelope of everything**, and in closed-canopy cells the
lowest point in the cell is a point *up in the canopy*: measured 3x3 range of the cell-minimum grid
has median 1.18 m but p95 = 24.6 m and 268,522 of 618,082 cells (43%) exceed 2 m. Blowing that up
gives 36% of all points *below* the "ground" surface. A morphological opening at 0.25 m (my PMF,
r = 1..12 cells) fixes only 1,351 cells (0.2%) — the opening is swamped by the spikes.

### Five independent arguments that the trough is not a surface artefact

1. **Lower-envelope test.** Only the two surfaces with ≤ 0.2% of points below them (SMRF, CSF) show
   the trough. Every surface with > 13% below it destroys it. `writers.gdal` writing cell-minima
   is exactly what makes SMRF/PMF useful here despite their otherwise odd behaviour.
2. **Two published, unrelated algorithms agree.** SMRF (Pingel et al. 2013) and CSF (Zhang et al.
   2016) put the trough at 1.05 m and 1.30 m with depth ratios 5.73 and 4.50.
3. **Offset invariance.** Shifting the SMRF surface by −0.25 / 0 / +0.25 m leaves the trough depth
   **exactly the same depth (36.6 pts/m²/m; 39.1 in the offset run, which used a marginally
   different cell mask)** in all three cases; only its apparent height moves (0.80 → 1.05 →
   1.30 m). A calibration artefact would move the depth too.
4. **Canopy points alone.** Restricting to `PredSemantic_FM == 2` (canopy) alone: trough 1.10 m,
   25.9 pts/m²/m, ratio **6.87**. The trough is a property of the canopy points, not of the
   near-ground mass bleeding across a bad surface.
5. **Closure stratification.** Present in every canopy-closed stratum (ratio 5-9), absent where the
   ground is directly visible (section 1).

### Honest caveat

The *existence* of the trough is robust; its **apparent severity is surface-sensitive**. Against a
surface that is even 0.5 m too low (the FM0 grid-min), the dip disappears entirely. So: if you are
reporting DBH statistics you need a lower-envelope ground surface at 0.25 m, or the blind band will
be silently absorbed into the near-ground mass.

### Notes on the methods

* My grid-level Cloth Simulation Filter (`scripts/10_numpy_ground.py`) **degenerates to the
  grid-minimum** when fed a min-z grid, because the cloth's resting height is
  `max(terrain in 3x3) − dmin`, and in open cells that equals the terrain minimum. It is retained
  in the table as a documented negative result, not as a recommendation. PDAL's per-point CSF is
  the real thing and it works, just not at 0.25 m on 72M points.
* Because `writers.gdal` writes cell-minima, **`writers.gdal` is the wrong tool for building a DEM
  in this pipeline.** Build the DEM in numpy (`dem_*.npy` in `scripts/03`, `10`) and compute HAG
  yourself. Round-tripping through PDAL's LAS writer would also silently drop `PredSemantic_FM`.

---

## 3. What do the three `PredSemantic_FM` values mean?

### Height-above-ground distribution (SMRF ground)

| class | share of cloud | HAG p1 / p5 / p25 / p50 / p75 / p95 / p99 (m) | % below 0.5 m | below 1.3 m | below 5 m |
|---|---|---|---|---|---|
| FM0 | 17.09% | 0.01 / 0.03 / 0.07 / **0.15** / 0.43 / 2.76 / 3.40 | 76.4% | 81.3% | 99.98% |
| FM1 | 2.89% | 0.46 / 2.22 / 8.12 / **15.72** / 21.05 / 26.76 / 30.52 | 1.1% | 2.9% | 15.3% |
| FM2 | 80.02% | 0.06 / 0.20 / 2.89 / **16.28** / 24.27 / 30.64 / 34.36 | 10.4% | 13.2% | 35.7% |

* **FM0 = ground surface + low understory, up to ~3.4 m.** It is *not* purely ground: 18.7% of FM0
  sits above 1.3 m and its p99 is 3.40 m. But it is overwhelmingly a near-ground class.
* **FM2 = canopy / foliage.** Median 16.3 m, but 13.2% of it is below 1.3 m — it is not a clean
  canopy mask either.
* **FM1 = woody material: branches and dead wood.** Evidence below.

### Geometry (KNN PCA, k = 16, 260,000-point spatially-uniform subsample)

| class | linearity | planarity | verticality | woody (lin>0.6 & vert>0.8) | leafy (lin<0.25) |
|---|---|---|---|---|---|
| FM0 | 0.428 | 0.452 | 0.231 | 0.66% | 17.05% |
| **FM1** | **0.503** | 0.301 | **0.583** | **13.79%** | **11.75%** |
| FM2 | 0.468 | 0.357 | 0.376 | 2.55% | 13.29% |

FM1 is **2.7x more woody than FM2 and 21x more woody than FM0**, is by far the most vertical class,
and is the least leafy. It sits 15.7 m above ground. That is branches and dead wood in the canopy,
not "dead ground". (The DBH band 1.0-1.6 m is also the single most stem-like band in the whole
cloud: 8.76% woody vs 0.03% at 0-0.3 m — consistent with real stems being there.)

### Colour (weak evidence — the RGB is broken)

B is zero everywhere and R is zero for all FM2 points, so colour is only usable for FM0/FM1. Using
R and G only, excess green (G−R)/(G+R): **FM0 +0.49 (green), FM1 −0.18 (brown/red)** — consistent with
bark. But do not build anything on the RGB in this file.

### **Can any threshold or combination of `PredSemantic_FM` isolate canopy stems at 1.3 m?**

**No. Directly and unambiguously no.**

Inside the 1.15-1.45 m band, 9.85% of points are "woody" by local geometry. Using `PredSemantic_FM`
as a predictor of woodiness:

| rule | share of woody points kept | precision | vs. keeping everything |
|---|---|---|---|
| keep all points | 100% | 9.85% | baseline |
| `FM == 2` (canopy) | 58.79% | 9.26% | **worse** |
| `FM == 0` (ground-ish) | 39.70% | 11.08% | within noise, < 40% recall |
| `FM == 1` | 1.51% | 6.82% | worse |

**AUC of the FM value as a continuous score for woodiness = 0.4763.** On a 0.5 scale, 0.5 is chance.
0.476 is *below* chance. The label carries **no information at all** about whether a point at 1.3 m
is on a stem or on a leaf.

This is not surprising in hindsight: FM is effectively a *height* label (FM0 ≈ below ~3 m, FM2 ≈
canopy, FM1 ≈ woody canopy). Stems and foliage at the same height receive the same label. **Stem
detection at the DBH band has to be geometric (PCA linearity + verticality + radius clustering),
not semantic.** The good news is that the geometric signal is there: 2.67% of band points are woody
and they cluster into ~530 stems per 3.78 ha.

---

## 4. Ground classification recommendation for this terrain

### Comparison

| method | coverage | % below surface | slope-robust? | verdict |
|---|---|---|---|---|
| PDAL SMRF, `cell: 0.25` | 95.83% | **0.11** | ratio 10.0 / 7.5 / 4.6 / 2.7 across 0-15° … 35-45° | **recommended** |
| PDAL CSF, `resolution: 0.5` | 95.84% | **0.02** | ratio 8.9 / 6.4 / 4.1 / 2.4 | **recommended cross-check** |
| coarse 1 m DEM (grid-min) | 95.45% | 31.45 | destroys the trough | reject |
| coarse 2 m DEM | — | — | — | reject |
| numpy slope-corrected opening | 95.45% | 14.76 | — | usable but 14% biased low |
| numpy PMF / grid-min / grid-CSF | 92-96% | 29-36% | — | reject |
| learned / pretrained | n/a | n/a | n/a | **none available, see below** |

### Recommendation

1. **Use `pdal` `filters.smrf` with `cell: 0.25`** as the primary ground classifier. It is the only
   method here that produces a genuine lower envelope, is fully slope-stratified-correct across the
   whole 0-45° range, needs no fitting, and runs in 2m04s / 4.7 GB.
2. **Cross-check with `filters.csf` at `resolution: 0.5`** (4m19s / 9.2 GB). If the two agree
   (they do: trough 1.05 vs 1.30 m, ratio 5.7 vs 4.5), you have confidence. **Do not run CSF at
   0.25 m on a 72M-point tile — it OOMs.**
3. **Never use a bare grid-minimum DEM at any resolution** on this terrain. It is 31-36% below the
   point cloud.
4. **Do not round-trip through PDAL's LAS writer.** It drops `PredSemantic_FM`. And
   `writers.gdal` writes cell-minima, so it cannot express an interpolated DEM.
5. **Coarse 1-2 m DEMs are fine as a *slope prior* and as the input to a morphological opening**,
   but not as the HAG reference. At 1 m the 0.25 m-scale trough is completely smeared out.
6. Optional improvement using data already in the file: seed with `PredSemantic_FM == 0 ∧ HAG_SMRF < 0.3 m`
   (a high-precision ground set) and re-fit. FM0 alone is *not* a ground surface (13.97% of all points
   lie below the FM0 grid-min), but `FM0` intersected with a low HAG is a good prior.
7. If you want a learned classifier later, MCC-LiDAR (Evans & Hudak 2007) is open-source classical
   code, scores 96.29% OA on OpenGF (vs PMF 90.63, CSF 93.07, PTD 94.82), and is not in PDAL.

### Pretrained / inference-only ground segmentation — what exists

* **Point-SCT** (Li et al., *IEEE TGRS* 63:5703118, 2025) is exactly this problem — ground
  filtering of ALS in complex mountainous terrain, using boundary-detector/curvature/average-elevation
  priors instead of RGB. **I found no public code or pretrained weights**; it does not appear to be
  released.
* **OpenGF** (Qin et al., CVPRW 2021) — 47 km² of finely labelled ground/non-ground ALS from 9 terrain
  scenes, plus training code for PointNet++ / DGCNN / KPConv / RandLA-Net / SCF-Net. Published test
  accuracies: **KPConv 97.79, PointNet++ 97.58, DGCNN 96.34, MCC 96.29, RandLA-Net 96.29,
  SCF-Net 95.75, PTD 94.82, CSF 93.07, PMF 90.63** (OA). This is the resource to use *if* fitting
  is allowed. It is not inference-only.
* PointNet / PointNet++ / RandLA-Net / StraightPCF release pretrained weights, but for indoor
  scenes, LiDAR semantic classes, or point *denoising* — none is a ground/non-ground binary
  classifier for forestry, and none can be applied zero-shot here.
* The 2025 survey *Towards intelligent ground filtering of large-scale topographic point clouds*
  states the blocker explicitly: DL ground filters "inherently [have] limitations in training … in
  forest regions" **"owing to the absence of ground points within areas covered by dense
  vegetation."** This tile is 21,057 m² (56% of the analysis area) of closed canopy with > 15 m
  closure. Fitting a ground model here would be fitting on labels that are absent by construction.

**Bottom line: there is no pretrained, inference-only ground-segmentation model for forestry /
photogrammetric point clouds that I could find and verify. Use SMRF at 0.25 m.**

---

## 5. How many trees are in this tile?

Three independent estimators over the 37,755 m² (3.78 ha) analysis area:

**(a) CHM local maxima** (max HAG per 0.5 m cell; CHM max 41.0 m, median 18.7 m, p90 29.8 m):

| window | 3 m | 4 m | 5 m | 6 m | 8 m |
|---|---|---|---|---|---|
| CHM > 3 m | 919 (243/ha) | 565 (150/ha) | 374 (99/ha) | 282 (75/ha) | 174 (46/ha) |
| CHM > 5 m | 733 (194/ha) | 482 (128/ha) | 330 (87/ha) | 249 (66/ha) | 163 (43/ha) |
| CHM > 8 m | 648 (172/ha) | 429 (114/ha) | 294 (78/ha) | 225 (60/ha) | 148 (39/ha) |

This varies by **5x** over a 3→8 m window and by 2.3x over the CHM threshold. Useless on its own:
cacao crowns are 2.5-3 m tall and fully overlapping while the CHM median is 18.7 m (tall shade trees
dominate). I report the 4-5 m window as the least-arbitrary reading: **374-565 stems, 99-150/ha.**

**(b) Crown connected components** above 1.6 m and 2.0 m: **5-7 components** for the whole tile. The
canopy is closed; crowns merge. Completely uninformative.

**(c) Direct stem clustering in the DBH band** — single-linkage at 10 cm over the 14,215 woody points
(2.67% of 533,105 band points) inside 1.15-1.45 m:

* 4,750 woody clusters total;
* **528 clusters with ≥ 6 points = 140 stems/ha = ~530 stems in 3.78 ha**;
* measured stem diameters: p10 6.2 cm, median 9.3 cm, p90 14.3 cm.

**(d) Consistency check on geometry.** 37,755 m² / 530 stems = 71 m² per stem = **8.4 m mean
spacing**. For a cacao agroforestry plot at 3.5-5 m nominal spacing you would expect 400-800
cacao stems/ha; 530 stems in 3.78 ha (140/ha) is therefore a **substantial undercount**, exactly as
expected given that the detector needs ≥ 6 woody points in a band where the median stem is 9 cm and
cacao at 10 cm DBH falls below the detection limit.

**Verdict: order 300-570 trees in the usable area, best point estimate ~530 detectable stems
(140/ha), and treat that as a lower bound on the true stem count.** Per-tree measurement is
tractable at this scale in aggregate (72M points / 530 trees ≈ 136,000 points per stem) but
**marginal at the DBH band**, where the budget is ~15-20 points per 0.30 m stem segment. Height,
crown position and canopy-layer metrics per tree are comfortable; DBH per tree will need either
multi-band stacking (several 0.3 m segments averaged, which cuts the noise by √n) or a vertical
integration over 0.5-1.3 m where there are 3-5x more points.

---

## Files

```
scripts/00_cache.py             stream LAS -> flat memmaps (5 s for 2.16 GB)
scripts/01_grids.py             0.25 m / 0.5 m per-cell min-z and counts
scripts/02_pdal_dems.sh         PDAL SMRF / CSF / PMF ground rasters
scripts/03_numpy_dems.py        numpy grid-min, PMF, FM0 surfaces at 0.25 m
scripts/04_tif2npy.py           GeoTIFF -> .npy bridge (DeepForest env)
scripts/05_coarse_grids.py      1 m and 2 m grids + DEMs
scripts/06b_hag.py              height above ground per surface
scripts/07b_profile.py          vertical density profile, slope-stratified
scripts/08_strata.py            profile stratified by slope and canopy closure
scripts/09_trough_sensitivity.py  trough vs ground-surface vertical offset
scripts/10_numpy_ground.py      numpy opening + grid Cloth Simulation Filter
scripts/11_surface_compare.py   all 9 surfaces side by side (main Q2 table)
scripts/12_fm_analysis.py       PredSemantic_FM semantics (Q3)
scripts/13_stems_trees.py       tree-count estimators + first stem budget attempt
scripts/14/15/16_stem_budget*.py stem point budget, three attempts
scripts/17_figures.py           all figures
figures/fig1_density_profile.png    profile per ground surface
figures/fig2_fm_semantics.png       FM composition vs height
figures/fig3_trough_offset.png      trough vs surface offset
figures/fig4_trees_and_profile.png  tree count vs CHM window
figures/fig5_dem_hillshade.png      hillshade of the candidate ground surfaces
results/*.json, results/*.txt       every number quoted above
```

Reproduce with: `/home/usuario/miniconda/envs/nksr/bin/python scripts/NN_*.py`
(only `04_tif2npy.py` needs `/home/usuario/miniconda/envs/DeepForest/bin/python` for rasterio).
Intermediate arrays live in `/tmp/opencode/cloud-research/out/` and are not committed.

## Open threads

* Reproduce or refute the ticket's slope figures (27.7° median) — I get 18.8°.
* `filters.smrf` and `filters.pmf` giving bit-identical rasters is worth a PDAL bug report; I could
  not reproduce a case where they differ at 0.25 m on this tile.
* The RGB channel is unusable (B ≡ 0, R ≡ 0 on all FM2 points). Worth raising with whoever produced
  the file, since the same pipeline probably feeds the other tiles.