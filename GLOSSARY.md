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

**`PredSemantic_FM`** — Extra byte dimension in the input LAS holding a 3-class
semantic label: **0 = ground, 1 = wood, 2 = foliage** (likely *follaje/madera*).
Verified by kNN-PCA geometry: wood has the highest linearity (0.555), lowest
planarity (0.310) and by far the most vertical principal axis (14.4% of its
neighbourhoods elongated and vertical, vs 1.6% for ground).

**Wood point** — A point labelled `FM=1`. Wood is a *sparse selective subset* of
the cloud: in the DBH band it is 0.10 pts/m² against 2.6 pts/m² for all classes,
so a 30 cm canopy stem may carry only 1–3 wood points. The label therefore serves
as a **stem seed generator**, not the surface a diameter is fitted to.

**Canopy stem** — Distinguished from cacao and from other wood by size and
position, not by class: both are labelled wood.

**Field truth** — Tape-measured DBH and height from a field campaign. Deliberately
withheld from this effort and never used for fitting *or* tuning.

**Manual reference** — Expert digital measurement (e.g. in CloudCompare) used to
compare automatic methods against each other. Bounds consistency, not accuracy.

## Artifacts

**The map** — The wayfinding issue that tracks decisions for this effort.

**The spec** — The destination artifact: a full design for the measurement
pipeline, written after wayfinding ends. Nothing is built during wayfinding.