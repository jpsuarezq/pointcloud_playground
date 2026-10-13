# Methodology: canopy-tree DBH and height from photogrammetric point clouds

Reference method for quantifying the structural inputs to above-ground biomass in
cacao agroforestry systems. This document describes the algorithm applied to the
point cloud, its parameters, its quality controls and its failure modes. It is the
paper-facing companion to `EXTENDED.md` and the code in `scripts/`.

> **Scope.** The pipeline measures, per canopy tree, **diameter at breast height
> (DBH)** and **total height**. These are the two predictors of standard
> pantropical allometric biomass equations (below). The pipeline does **not**
> measure biomass directly, and it makes **no accuracy claim**: no field truth was
> used at any point, for fitting or for tuning.

## 1. Design constraints

1. **No fitting, no tuning.** Field truth is deliberately withheld. Every threshold
   is set from principle or from a proxy (e.g. an intensity split derived from the
   given semantic labels), then locked in the config. Nothing is optimised against
   held-out truth, because tuning to error is fitting.
2. **Failure must be visible.** With no truth, a plausible number with a small
   residual that is actually circling a buttress is worse than an explicit refusal.
   Every measurement therefore carries QC metrics and a confidence; the method
   refuses rather than guesses.
3. **Consistency, not accuracy.** Internal cross-method agreement and per-tree QC
   bound *consistency*; they never bound accuracy.
4. **Research-first, single tile.** Scale, throughput and production hardening are
   out of scope.

## 2. Input and coordinate conventions

- Input: a LAS/LAZ point cloud carrying a **per-tree instance** dimension
  (`PredInstance_FM`) and a **semantic** dimension (`PredSemantic_FM`: 0 ground,
  1 wood, 2 foliage). Additional fields used: `intensity`, return number/count,
  and a scan index (`Original cloud index`, two merged scans).
- Projected CRS: **EPSG:32618** (UTM 18N). Output table also carries **WGS84
  (EPSG:4326)** latitude/longitude.
- **Height above ground (HAG)** is the **vertical** distance to the terrain
  surface (ADR-0001), not the distance perpendicular to the slope. DBH is defined
  at 1.3 m of HAG.

## 3. Pipeline

Automated end-to-end with `scripts/run_pipeline.py` (one cloud in, all outputs
out). Stages:

### 3.1 Field cache
One streaming pass (`laspy` + `lazrs`, 8 M-point chunks) writes flat memmaps of
x,y,z (local float32, mm) and every scalar field. 72.2 M points → 1.9 GB, ≈16 s.

### 3.2 Ground surface and HAG
Ground points are classified with **PDAL `filters.smrf`** at a 0.25 m cell
(≈1 m 22 s, 6.4 GB peak RSS). The ground DEM is a **grid-minimum** of those points
per 0.25 m cell, holes nearest-filled up to 2 m, then a 3×3 median (the PDAL GDAL
writer is not used: it emits cell minima, not means). HAG is the bilinear-sampled
DEM subtracted from each point's Z.

### 3.3 Instance partition and canopy separation
The per-tree instance dimension defines the trees. Canopy trees are those with
**max HAG > 12 m** (cacao is 3–5 m, so a height threshold separates the two
populations). For each canopy instance we tabulate point count, top HAG, and the
number of points in the **DBH band** (`1.15 ≤ HAG ≤ 1.45 m`).

### 3.4 Stem isolation — multi-candidate seeding
The stem axis is found by a **vertical-line Hough search**: candidate axes are
vertical lines through candidate points; an axis scores by the number of 0.5 m
height bins in 1–5 m that contain a candidate point within 0.30 m. The axis that
*persists upward* wins (not the densest cluster — that locks onto branches).
Candidate point sets, tried in order and refined, are:

1. **FM-wood** points (the semantic wood class).
2. **Intensity proxy** — points with intensity below a per-scan threshold τ, where
   τ is the midpoint of the FM wood/foliage intensity medians (a proxy, not a fit).
3. **Base anchor** — the median XY of the lower-stem points (0.5–2.5 m), vertical.
4. **Band seed** — the densest DBSCAN cluster of the 1.15–1.45 m points nearest
   the base, vertical.

Each candidate axis is refined by iterated trimming to points within 0.35 m,
refusing to roll past **35° of tilt**. The candidate whose 1.3 m ring best satisfies
the QC is kept.

### 3.5 Diameter measurement
The ring is measured **perpendicular to the fitted axis** (a horizontal cut of a
leaning stem yields an ellipse, inflating diameter ~1/cos θ). In a slab of half
width ±0.15 m about 1.3 m, points within `rmax` of the axis are projected onto the
plane normal to the axis. The radius is the **per-angular-bin 25th-percentile axis
distance** (the innermost coherent surface — robust to far clusters and to twigs
crossing the axis), cross-checked by a **radius-bounded RANSAC circle** (unbounded
RANSAC fits near-collinear arcs with hundreds-of-metres radii). Diameter = 2× the
median bin radius. The slab is widened to ±0.30 m only when the narrow slab yields
no valid ring.

### 3.6 Height
Total height = maximum HAG of the instance.

### 3.7 Quality control and confidence
Each measurement emits: **arc coverage** (fraction of angular bins filled), per-bin
radius spread, **axis tilt**, number of height bins the stem persists through, and
the spread of the diameter across four axis-fitting spans (stability). A
measurement is **valid** when arc ≥ 0.35 and diameter < 1.0 m. Confidence:

- **high** — seeded from FM wood and valid;
- **medium** — seeded from a proxy/base/band and valid;
- **low** — emitted but not valid (low arc, wide slab, diameter ≥ 1 m, or unstable);
- **none** — refused (no seed, or no points near the axis).

### 3.8 Outputs
A per-tree table (`tabla_arboles.csv`), one `.laz` per tree grouped by confidence
(`trees/{high,medium,low,refused}/tree_<id>.laz`), and a `results.json` with the
full per-tree record. Height is in metres; DBH in centimetres.

## 4. From structure to biomass

Biomass is obtained by applying a **pantropical allometric equation** to the
measured DBH and height, e.g. Chave et al. (2014)
`AGB = 0.0673 · (ρ · DBH² · H)^0.976`, with wood specific gravity ρ per species
(or a plot-level default). Because the pipeline produces a **per-tree** record with
a confidence, biomass can be:

- aggregated per tree, then per plot, and
- **stratified by confidence** (report a high/medium/low breakdown rather than a
  single number), which is the honest form given the absence of truth.

Species-specific ρ is **out of scope**: no species product exists on this dataset.

## 5. Algorithm behaviour (whole canopy, 72.2 M-point tile)

| Quantity | Value |
|---|---|
| Canopy instances (top HAG > 12 m) | 674 |
| Measured | 312 (46%) |
| Refused | 362 |
| Confidence: high / medium / low | 78 / 193 / 41 |
| Seed that isolated the axis: band / wood / base / intensity | 120 / 94 / 63 / 35 |
| Slab used: ±0.15 m / ±0.30 m | 224 / 88 |
| Refused with **zero** points in the band | 312 |
| Refused with ≥10 band points | **1** |
| DBH cm p10/25/50/75/90 (all) | 15.6 / 22.7 / 30.7 / 41.3 / 58.2 |
| DBH cm p10/50/90 (high) | 16.9 / 26.0 / 46.1 |
| Total height m p10/50/90 | 14.5 / 23.5 / 32.3 |
| Arc coverage p10/50/90 (all) | 0.31 / 0.50 / 0.88 |
| Buttress suspected (1.3 m ring > 1.35× the 3 m ring) | 36 |
| Median diameter spread across fitting spans | 5.6 cm |
| Runtime (6-core, whole tile) | ≈323 s |

Interpretation: measurability is governed by the **band point budget** — essentially
every tree with a usable number of points in the 1.3 m band is measured, and the
refusals are dominated by trees with **no points at 1.3 m** (a photogrammetric
mid-storey blind spot coincident with the DBH band). The remaining uncertainty is
captured by confidence, not hidden.

## 6. Figures

- `proto_fig1_exemplar.png` — stem isolation on the exemplar tree: side view with
  the fitted axis over the wood points.
- `proto_fig2_plan.png` — perpendicular plan views at 1.3/2/3/4 m showing partial
  arcs and off-axis branches.
- `proto_fig3_agreement.png` — naive horizontal slice vs axis-perpendicular (the
  horizontal slice inflates diameter several-fold), and taper extrapolation vs
  direct (no agreement).
- `proto_fig4_qc.png` — accepted-ring arc coverage, stability across fitting spans,
  and the outcome breakdown.
- `ext_fig1_map.png` — all canopy trees by confidence (marker size ∝ DBH).
- `ext_fig2_dbh.png` — DBH distribution by confidence; height vs DBH.
- `ext_fig3_arc_dbh.png` — QC: arc coverage vs DBH with the accept thresholds.
- `ext_fig4_outcomes.png` — outcome, confidence, seed and slab breakdowns.
- `ext_fig5_band_budget.png` — band point budget vs measurability.

## 7. Parameters (locked)

| Parameter | Value | Source |
|---|---|---|
| DBH height | 1.3 m HAG | ADR-0001 / forestry convention |
| Band half-width | 0.15 m (wide 0.30 m) | narrow-first, widen on failure |
| Ring radius cap `rmax` | 0.45 m (wide 0.65 m) | canopy-stem scale |
| RANSAC radius bounds | 0.015–0.80 m | excludes runaway circles |
| Valid arc | 0.35 | ≥ ~6 of 16 angular bins |
| Valid DBH | < 1.0 m | canopy stems are ≤ ~1 m |
| Axis tilt cap | 35° | beyond this the axis rolls onto a branch |
| Canopy threshold | top HAG > 12 m | separates canopy from cacao (3–5 m) |
| SMRF cell / DEM res | 0.25 m | input density |
| DEM hole fill | ≤ 2 m | SMRF fills ~44% of cells |
| Intensity proxy τ | per-scan wood/foliage midpoint | given semantic labels |

## 8. Limitations

- **No ground truth.** All values are consistent magnitudes, not verified ones.
  Confidence bounds evidence, never accuracy.
- **Partial arcs.** Most accepted rings rest on 40–60% angular coverage; the radius
  is an *effective* radius and inherits buttress bias on asymmetric stems.
- **Proxy seeds.** Band- and intensity-seeded axes can circle a fork or buttress;
  these are `medium`, i.e. plausible but unverified.
- **DEM dependence.** On steep ground the DEM error propagates into DBH; a
  slope-stratified confidence is a natural extension (ADR-0001).
- **"Measured" ≠ "a trunk".** Some instances carry thousands of band points and no
  labelled wood — they are measured, but may not be trees.

## 9. Reproducing / running on new clouds

```sh
python scripts/run_pipeline.py --config pipeline_config.example.json
# or
python scripts/run_pipeline.py --input /path/cloud.laz --out /path/run
```

For a new cloud, set `input` and, if its dimension names or class ids differ, the
`fields` / `*_class` entries. Everything else is locked. Outputs land under `out/`:
`tabla_arboles.csv`, `trees/<confidence>/tree_<id>.laz`, `results.json`,
`metrics`-style figures. To compare several clouds, run each into its own `out/`
and aggregate the tables.
