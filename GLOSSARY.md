# Glossary

Vocabulary for the NubePuntos effort: measuring canopy-tree diameter and height
from photogrammetric point clouds of cacao agroforestry systems.

## Trees and crops

**Canopy tree** — A tree of the overstory, 15–35 m tall, DBH roughly 30–100+ cm. The
subject of measurement in this effort. Not to be confused with cacao.

**Cacao** — The understory crop (*Theobroma cacao*), 3–5 m tall, DBH roughly 5–15 cm.
Out of scope as a measurement target: at this point density a cacao stem carries
~7 points at the DBH band, which cannot support a diameter estimate.

**Cacao agroforestry** — The land system: cacao grown beneath a canopy of taller
trees, on mountainous terrain. The canopy trees and the cacao are distinct
populations in the point cloud, not one class of "tree".

## Measurement

**DBH (diameter at breast height)** — Stem diameter measured at 1.3 m above ground.
In this project defined as *vertical* distance above the terrain surface, not
perpendicular to the slope. See ADR-0001.

**Height above ground (HAG)** — Vertical distance from a point to the terrain
surface beneath it. Requires a DEM; computed here from the point cloud itself.

**Total height** — Vertical distance from the terrain surface to the crown apex.

**Blind spot** — The collapse in point density at ~1.05–1.10 m above ground: 37
pts/m²/m, 5.7× below the canopy peak. Intrinsic to dense-matching reconstruction under
**closed canopy**, not a ground-surface artifact — reproduced against nine ground
surfaces and two published algorithms. Its apparent height shifts with the DEM; its depth
does not.

**Point budget** — Points on a stem-sized object in the DBH band: **~15–20 per 0.30 m
segment, and flat in diameter.** Wider stems do not collect more points, so canopy trees
are no better resolved than cacao. No stem wider than ~20 cm was detectable at all. The
governing constraint on method choice.

**Woody point** — A point of woody material, identified *geometrically* (kNN-PCA
linearity and verticality) rather than by any stored label. ~2.7% of points in the DBH
band. This is the only reliable stem signal available.

**Stem** — A woody vertical structure detected geometrically. ~530 estimated in the
3.78 ha analysis footprint (140/ha), a lower bound.

## Data and labels

**`PredSemantic_FM`** — Extra byte dimension in the input LAS: the *semantic* output
of **ForestMamba** as run by the 3dtrees.earth platform (FM = ForestMamba), 3 classes
— ground, wood, leaf. Verified by kNN-PCA geometry: wood is the most linear (0.555),
least planar (0.310) and most vertically-aligned class (14.4% of its neighbourhoods
elongated and vertical, vs 1.6% for ground). Wood is therefore the **2.9%** class, not
the 80% one.

**3dtrees.earth** — Open platform (Fraunhofer IPM, Univ. Freiburg, GFZ et al.) that
processes close-range LiDAR via a Galaxy pipeline. **The source of our LAS.** Its
documented default workflow already computes per-tree DBH, tree height, crown volume,
species and a CHM. See `Adopt or build` ticket.

**ForestFormer3D** — Joint semantic and instance segmentation model (Nguyen et al.,
ICCV 2025). Installed locally at `/home/usuario/3DFormer/ForestFormer3D/` with a
trained checkpoint, but **not run on our data**. Licence CC BY-NC 4.0.

**ForestMamba** — Successor to ForestFormer3D (Nguyen et al., BMVC 2026), same
codebase. Integrated into 3dtrees.earth since Aug 2026. CC BY-NC 4.0.

**SegmentAnyTree** — Tree instance segmentation model (Wielgosz et al., 2024) used by
3dtrees.earth. Its output is what makes per-tree structural metrics possible.

**CspStandSegmentation** — Cylindrical stand segmentation (Frey & Schindler, 2024),
the tool 3dtrees.earth uses to extract DBH, tree height and crown volume per tree.

**`PredInstance`** — Per-tree instance identifier written by the platform. We hold
only `PredSemantic_FM`; whether `PredInstance` can be re-fetched is open.

**ITSMe** — R package (`lmterryn/ITSMe`) of individual-tree structural metrics.
`dbh_pc()` fits a least-squares circle through a 6 mm slice at 1.3 m, with QC
metrics. The leading concrete method if we build rather than adopt.

**treeID** — ForestFormer3D's per-tree instance identifier (the local equivalent of
`PredInstance`).

**Wood point** — A point in the `FM=1` class. It is genuinely woody material — 21× more
woody than ground by geometry, verticality 0.583 vs 0.231, median HAG 15.7 m; branches
and dead wood. **But it carries no usable signal in the DBH band**: only 1.5% of points
there are `FM=1`, at 6.8% precision. Scored as a continuous woodiness predictor across
the band, `FM` gives **AUC 0.476 — below the 0.5 chance line**. It is effectively a
*height* label: stems and foliage at the same height get the same label. Stem detection
there must be geometric.

**Canopy stem** — Distinguished from cacao and from other wood by size and
position, not by class: both are labelled wood.

**Tile** — One processed point-cloud file. The repo LAS is a **stripped 72.2M-point tile
of 3dtrees.earth dataset 3479 "Agroforestal"** (207,364,864 pts, EPSG:32618). It kept
only `PredSemantic_FM`; the full file also carries `PredInstance_FM`, `PredScore_FM`,
`PredInstance_SAT`, `PredSemantic_SAT`, plus `ReturnNumber`, `NumberOfReturns`,
`Intensity` and `UserData`. `/home/usuario/3DFormer/cloud.laz` (168.9M pts, Z 248–740)
is a different dataset entirely.

**Field truth** — Tape-measured DBH and height from a field campaign. Deliberately
withheld from this effort and never used for fitting *or* tuning.

**Manual reference** — Expert digital measurement (e.g. in CloudCompare) used to
compare automatic methods against each other. Bounds consistency, not accuracy.

## Artifacts

**The map** — The wayfinding issue that tracks decisions for this effort.

**The spec** — The destination artifact: a full design for the measurement
pipeline, written after wayfinding ends. Nothing is built during wayfinding.