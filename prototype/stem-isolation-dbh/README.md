# Prototype: stem isolation and DBH at breast height (ticket #11)

**Throwaway prototype.** Not production, not validated, no tests. Its only job was
to answer one question with a rough artifact:

> Can a canopy tree's main stem be **isolated and measured at 1.3 m** on this
> capture — and what method actually does it?

Verdict: **partly, for a minority of trees, and only with the failure made
visible.** A near-vertical woody-column search isolates a stem axis on ~a third
of well-sampled canopy trees; the naive horizontal slice is grossly wrong; and the
taper extrapolation the ticket hoped for does **not** hold up.

Branch: throwaway primary source, never merged. Companion issue: #11.

## Data and provenance

- Input: `pointclouds/segmentado_crop.laz` (72,162,549 pts), fields
  `PredInstance_FM` (per-tree id) and `PredSemantic_FM` (0 ground / 1 wood / 2 leaf).
- Height above ground (HAG): SMRF at 0.25 m, from the characterisation work on
  branch `research/characterise-input-cloud`; HAG is the *vertical* height above
  terrain (ADR-0001), so the DBH band is `1.15 <= HAG <= 1.45`.
- Canopy = FM instance with top HAG > 12 m: **598 instances**; 328 have >=1 band
  point, 184 >=20, 125 >=50.

> Tracker figures (40.7% coverage) and these (54.8%) differ because of the ground
> surface and instance-set used; both agree the band is sparse and the median
> count for covered trees is ~15.

## Method (what the prototype actually does)

`scripts/dbh_proto.py`:

1. **Isolate the stem axis — vertical-column search, not point clustering.**
   Candidate columns are vertical lines through wood points; a column scores by the
   number of 0.5 m bins in [1, 5] m that hold wood within 0.30 m. The winner is the
   structure that *continues upward*, not the densest cluster. (v1 used DBSCAN on
   wood and locked onto branches with axis tilt 50–87°.)
2. **Refine the axis** by iterated trimming to points within 0.35 m, refusing to
   roll past 35° of tilt.
3. **Measure perpendicular to the axis.** In a +/-0.15 m slab, project points onto
   the plane normal to the axis and estimate radius as the **per-angular-bin low
   quantile** of axis distance (innermost coherent surface), then cross-check with a
   **radius-bounded** RANSAC circle. (v1's unbounded RANSAC fitted near-collinear
   arcs with diameters of *hundreds of metres*.)
4. **QC every measurement**: arc coverage, per-bin radius spread, tilt, number of
   height bins the stem persists through, and stability of the diameter across four
   axis-fitting spans.
5. **Refuse explicitly** — never emit a bare number — when wood is insufficient, no
   persistent column exists, arc < 0.35, or diameter >= 1 m.

## Results (43 sampled canopy trees with >=10 band points)

| Outcome | n |
|---|---|
| measured at 1.3 m | 14 |
| refused | 29 |

Refusals: insufficient wood in the lower stem 10; no points near the axis 9; no
persistent wood column 4; insufficient arc 2; QC reject 4.

- Accepted ring diameter: **22 / 26 / 35 / 40 / 45 cm** (p10/25/50/75/90).
- Arc coverage of accepted rings: median **0.59** (p10 0.40, p90 0.79) — most are
  partial arcs, consistent with the capture's ~0.30–0.50 coverage in the band.
- **Naive horizontal slice vs axis-perpendicular: median ratio 3.25x (max 5.2x).**
  The failure is not the ~6–15% ellipse from lean; it is the slice catching
  buttress, fork and neighbours. See `figures/proto_fig3_agreement.png`.
- **Taper extrapolation from 2–5 m does not hold:** only 4 of 17 trees agree with
  the direct measurement within 25%; ratio spans 0.18x–13x. The "clean 2–5 m zone"
  is *not* clean on tall trees — the ring profile is contaminated by forks and
  branches, and is frequently non-monotonic with height.
- Buttress suspected (1.3 m ring > 1.35x the 3 m ring) on 3 of 14 measured trees.

### Figures

- `figures/proto_fig1_exemplar.png` — instance 2957 side view, fitted axis over the
  wood (red) trunk; the axis tracks the pole from the base into the crown.
- `figures/proto_fig2_plan.png` — 2957 perpendicular plan views at 1.3/2/3/4 m,
  showing partial arcs and a persistent secondary branch ~0.5–0.7 m off-axis.
- `figures/proto_fig3_agreement.png` — method agreement.
- `figures/proto_fig4_qc.png` — QC distributions and outcome counts.

## Failure modes (the point of the exercise)

1. **Confident wrong numbers are prevented, not avoided.** The pipeline refuses
   ~two thirds of sampled trees; the earlier "plausible 15 cm with 3 cm residual"
   is exactly what this refuses.
2. **The stem signal is often absent.** 10 instances with thousands of band points
   carry *zero* labelled wood in 1–5 m — those instances are not measuring a trunk
   (shrub/crown/understory, or a wood-labelling gap). Worth its own look.
3. **Partial arcs.** Most accepted rings rest on 40–60% angular coverage; the
   per-bin low-quantile radius is an *effective* radius, not a true circumference
   fit, and on asymmetric stems it inherits buttress bias.
4. **Forks.** A persistent off-axis wood column at 3–5 m (see 2957) means the
   "clean zone" used for taper is a fork, not a tapered pole.

## Reproduce

Environment: conda env `nksr` (laspy 2.5.4 + lazrs, numpy, scipy, open3d,
scikit-learn, matplotlib). `conda` itself is broken on this machine; call the env
python directly.

```sh
# one streaming pass keyed to the existing /tmp cache from research/characterise-input-cloud
python scripts/20_instances.py     # PredInstance_FM per point (cache order)
python scripts/21_canopy.py        # per-instance canopy summary
python scripts/22_run_proto.py     # run the method over the sample
python scripts/26_figures.py       # figures + summary
```

Intermediate state lived in `/tmp/opencode/cloud-research/out` (not committed).

## What this does NOT establish

- No accuracy claim. There is no ground truth and none was used.
- The diameters are **consistent magnitudes**, not verified values.
- Sample is 43 trees, biased to the best-covered instances.
