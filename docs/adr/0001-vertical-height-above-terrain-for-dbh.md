# DBH is measured as vertical height above terrain, not perpendicular to slope

The terrain is mountainous — median slope 27.7 degrees, p90 39.8, with 78.6 m of
relief across the tile — so the definition of "1.3 m above ground" is load-bearing
rather than a detail. We define height above ground (and therefore the DBH
measurement height) as **vertical distance above the terrain surface**, matching the
standard forestry convention and the fact that trees grow upward.

## Considered options

- **Perpendicular to the slope surface.** Rejected: on the median 27.7-degree slope
  this would place the measurement ~15% further up the trunk than intended, and the
  error would be *slope-correlated* — inflating diameters on one aspect and shrinking
  them on the opposite one. A systematic bias that looks like real ecological
  variation.
- **Vertical above terrain (chosen).** Truncates the stem perpendicular to gravity at
  1.3 m, which is how DBH is defined and measured in the field with a tape.

## Consequences

- Every height and diameter in the output is a vertical quantity. On slopes this
  differs materially from the along-surface distance a person would walk.
- The quality of the DEM now propagates directly into DBH. At 0.5 m cell size on a
  ~28-degree slope, horizontal DEM uncertainty converts to roughly 0.24 m of vertical
  error — about 18% of the entire 1.3 m measurement height.
- Therefore trees on steep ground must carry lower confidence than trees on flat
  ground, and the output schema needs slope-stratified QC rather than one uniform
  confidence field. See the output-schema decision.
- SMRF filled only ~35% of 0.5 m cells on this terrain, so the DEM method is a
  first-order risk to DBH quality, not a detail. It is settled by the
  characterisation ticket.