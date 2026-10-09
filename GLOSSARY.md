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

**Blind spot** — The measured collapse in point density at roughly 1.0–1.6 m above
ground (to ~72 pts/m²/m from ~390 below and ~130 above). An artifact of
photogrammetric dense matching in cluttered mid-storey, and coincident with the
DBH band.

**Point budget** — The number of points available on a stem within the DBH band.
Canopy stems carry ~20 (30 cm DBH) to ~34 (50 cm DBH) points on a 30 cm segment.
The governing constraint on method choice.

## Data and labels

**`PredSemantic_FM`** — Extra byte dimension in the input LAS holding the 3-class
*semantic* output of **ForestFormer3D** (Nguyen et al., ICCV 2025): ground, wood,
leaf. Verified by kNN-PCA geometry — wood is the most linear (0.555), least planar
(0.310) and most vertically-aligned class (14.4% of its neighbourhoods elongated and
vertical, vs 1.6% for ground). Wood is therefore the **2.9%** class, not the 80% one.

**ForestFormer3D** — Joint semantic and instance segmentation model for forest
point clouds. Outputs a per-point class *and* a per-tree `treeID`. Runs locally at
`/home/usuario/3DFormer/ForestFormer3D/` with a trained checkpoint. Licence
CC BY-NC 4.0 (non-commercial).

**ForestMamba** — Successor to ForestFormer3D (BMVC 2026), same codebase, with a
released checkpoint reported to beat the paper numbers. CC BY-NC 4.0.

**treeID** — ForestFormer3D's per-tree instance identifier. Present in the source
cloud but **discarded** by the script that produced the repo LAS. Recovering it
would give individual-tree segmentation essentially for free.

**Wood point** — A point in the wood class. Wood is a *sparse selective subset* of
the cloud: in the DBH band it is 0.10 pts/m² against 2.6 pts/m² for all classes,
so a 30 cm canopy stem may carry only 1–3 wood points. The label therefore serves
as a **stem seed generator**, not the surface a diameter is fitted to.

**Canopy stem** — Distinguished from cacao and from other wood by size and
position, not by class: both are labelled wood.

**Tile** — One processed point-cloud file. Note the source cloud and the repo LAS
are different tiles at different elevations; provenance per tile is not yet
established.

**Field truth** — Tape-measured DBH and height from a field campaign. Deliberately
withheld from this effort and never used for fitting *or* tuning.

**Manual reference** — Expert digital measurement (e.g. in CloudCompare) used to
compare automatic methods against each other. Bounds consistency, not accuracy.

## Artifacts

**The map** — The wayfinding issue that tracks decisions for this effort.

**The spec** — The destination artifact: a full design for the measurement
pipeline, written after wayfinding ends. Nothing is built during wayfinding.