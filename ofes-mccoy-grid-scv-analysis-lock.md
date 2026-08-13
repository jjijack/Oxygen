# OFES McCoy-grid thermohaline SCV analysis lock

Frozen: 2026-08-14, before implementing the full-grid detector and before
inspecting its 56-event association result.

## Purpose and status

This document fixes a model-grid extension of the McCoy et al. (2020)
profile-based SCV method. It is a method lock written after exploratory OFES
work, not a claim of preregistration before any OFES result existed. The
existing OFES virtual-profile counts (9/56 center profiles, 19/56 any-event
profiles, and their velocity-confirmed subsets) are frozen comparison targets
only. They may not be used to tune this detector.

The detector is deliberately split into three labels:

1. `mccoy_seed_lens`: a three-dimensional thermohaline lens grown from a
   McCoy-compatible profile seed;
2. `grid_supported_scv`: a seed lens that also has an OFES-native
   anticyclonic velocity structure. This is the primary SCV estimand;
3. `persistent_grid_scv`: a grid-supported object that can be linked across
   daily outputs. Duration is descriptive and is not part of SCV identity.

DO and all event labels are excluded from every detector gate. DO is read only
after the instantaneous and tracked catalogs are complete.

The detector is complementary to Zhu et al. (2024). Zhu uses an absolute
Okubo-Weiss threshold and a vorticity-dominant 3-D object; this lock uses
relative thermohaline/PV structure and native velocity confirmation. The Zhu
zero association is retained as a resolution-limited negative control and is
not re-tuned.

## Frozen source and coordinate conventions

- The McCoy profile chain must call the existing OFES implementation and the
  archived public MATLAB source identified by `source_archive_sha256 =
  ee0abca349c4f090b3c99a695a88ee81df418fd3f756110d165f83b693b26406`. The
  background sampling geometry, 10-dbar interpolation, IQR construction,
  Gaussian search, weak-`N2` gate, and dynamic-height gate are copied without
  changing their values or no-op behavior.
- OFES temperature and salinity are converted with TEOS-10. The detector uses
  `sigma0`, spiciness, `N2`, and rescaled isopycnal PV; it does not use DO,
  SSH, META contours, or event masks.
- The fixed density nodes are
  `25.60, 25.65, ..., 27.50 kg m^-3` (step `0.05 kg m^-3`).
- Each water column is interpolated to these nodes only where the density
  crossing is unique, monotonic, finite, and bracketed by delivered z-levels.
  No extrapolation is permitted. A density node with more than one crossing
  or no valid bracket is invalid for that column.
- The valid physical depth window is `100--1000 m`. A valid isopycnal voxel
  must have its interpolated depth inside this window. The delivered bottom
  is retained as a censor boundary; it is never extrapolated.
- For a continuous event `target_sigma0`, the association node is determined
  before results are inspected by fixed Voronoi density cells around the
  nodes (midpoints between adjacent nodes; half-open intervals with the
  lower node owning a tie). The manifest records target density, selected
  node, cell bounds, and mismatch. There is no post-result nearest-layer
  choice or density tolerance.
- Horizontal areas use the delivered tracer-cell metric. Vertical voxel
  thickness uses interfaces at density-node midpoints only for catalog
  bookkeeping; the physical depth validity and object thickness use the
  interpolated z-level geometry.

## Frozen relative PV diagnostic

Rescaled isopycnal PV follows Ernst et al. (2023), Eq. (6), DOI
`10.1175/JTECH-D-22-0121.1`,

\[
  PV_r = (\nabla\times\mathbf U + f)\cdot\nabla Z(\rho),
\]

with `Z(rho)` a fixed reference-density profile constructed from the same
daily OFES background domain. The implementation must record the reference
profile construction and source/code hash. Ernst's Arabian-Sea smoothing,
threshold, and SSH-optimization values are not imported.

For each valid density node, the detector forms a local PV anomaly relative to
the same-day 120--240 km background annulus. Candidate PV cores are local
minima contained by a closed contour of the relative anomaly. There is no
absolute PV threshold and no absolute Rossby-number threshold. The contour
level is the largest closed contour that still encloses the local minimum,
with the deterministic grid contour implementation recorded in the manifest.
PV is a structural diagnostic and cannot by itself create a seed without a
thermohaline lens.

## Frozen McCoy-seeded thermohaline lens

### Candidate prefilter and recall

The vectorized prefilter is only a computational superset of the existing
McCoy gates. It must retain all 141 known McCoy-compatible event profiles
(`141/141` prefilter recall) before the final grid catalog is evaluated. Recall
is an engineering invariant, not a target for changing the detector.

Every retained seed is classified with the unchanged profile chain:

- profile QC, equatorial exclusion, and at least 61 valid background profiles;
- local along-isopycnal spice and `N2` IQR anomalies;
- the existing `1.5 * IQR` and minimum-level gates;
- pycnocline-below restriction and surface-connected anomaly rejection;
- Gaussian spice structure and its fixed height, amplitude, center, `R2`, and
  normalized-error gates;
- weak-`N2` core gate; and
- the existing dynamic-height internal-maximum gate after BC1 removal.

The detector records the seed's spicy/minty sign, peak density, source
position, and all failure stages. No DO value is available to this step.

### Three-dimensional connection

Seeds are connected across adjacent density nodes only when the horizontal
thermohaline core masks overlap in physical space, have the same spice sign,
and both nodes are valid. Missing nodes do not bridge an object. A connected
component must satisfy all of the following fixed identity gates:

- equivalent horizontal radius: `5--50 km`;
- at least three adjacent valid density nodes;
- physical vertical thickness at least `100 m`;
- volume-weighted centroid depth at least `200 m`;
- shallowest occupied interpolated depth strictly greater than `100 m`;
- no required connection to the surface; and
- object size and depth censor flags retained when a delivered boundary is
  touched.

The horizontal boundary at each occupied node is the outermost closed contour
of the `N2` anomaly enclosing the thermohaline core, following Barabinot et
al. (2025). If no such closed contour exists, that node is not part of a lens.
The boundary method is a thermohaline boundary definition, not a claim that
relative-vorticity contours are universally invalid.

The PV minimum, spice sign, and `N2` anomaly are saved as continuous fields.
No fixed absolute PV, `N2`, spice, or Ro amplitude is introduced beyond the
existing McCoy gates and the identity gates above.

## Frozen native-velocity confirmation

The second-layer `grid_supported_scv` label requires all of the following:

- OFES native `u/v` collocated to tracer centers with the validated four-corner
  interpolation and SI conversion;
- background translation removed using the fixed local background definition;
- the four Nencioli et al. (2010) vector-geometry conditions with fixed
  `a = 3` and `b = 2` native grid cells. Along the cardinal east-west and
  north-south axes, component signs must reverse across the candidate center
  and absolute component magnitudes must be non-decreasing over offsets
  `1..a`; the center speed must be the finite minimum in its `(2b+1)^2`
  neighborhood; and the vectors on the square perimeter at offset `a-1` must
  traverse all four quadrants with one constant sense of rotation and no pair
  separated by more than one quadrant. The candidate center is the valid
  velocity minimum nearest the lens volume centroid and must lie inside the
  thermohaline mask;
- a signed circulation integral sampled at 64 equally spaced points on the
  outermost closed `N2` boundary, with at least 80 percent finite samples and
  anticyclonic sign in the Northern Hemisphere;
- consistent anticyclonic sign on at least two adjacent occupied density
  nodes; and
- no absolute minimum Ro threshold.

The existing single-point Ro/strain classification is retained as a secondary
diagnostic only. A lens that passes the thermohaline identity but fails native
velocity confirmation remains in the seed-lens catalog and is not silently
discarded.

The technical catalog retains cyclonic candidates under the same gates. The
primary SCV table is the Northern Hemisphere anticyclonic subset.

## Engineering and analytic validation

Validation is performed after this lock and before any population association.
It can diagnose implementation errors but cannot amend a frozen rule:

- a stationary field yields zero objects;
- a synthetic solid-body anticyclone passes the velocity geometry and sign;
- a solid-body cyclone enters only the cyclone technical catalog;
- pure shear produces no closed grid-supported object;
- adding a uniform translation does not change the object or confirmation
  result after background removal;
- a heave-only density displacement without a thermohaline anomaly does not
  create a seed lens;
- a surface-connected lens is rejected by the subsurface identity gate; and
- density interpolation, PV units/sign, contour closure, and object-volume
  calculations satisfy independent analytic checks.

Validation failures require implementation repair or an explicitly separate
versioned sensitivity. They never authorize changing this lock after viewing
the 56-event table.

## Frozen event association and nulls

The primary association uses the 56 strict DO50 events, their formal peak day,
peak horizontal cell, and `target_sigma0`. The event target density is mapped
to the fixed density cell rule above. The primary binary estimand is whether
that event core voxel is contained by a same-day `grid_supported_scv` object.
There is no search radius, overlap cutoff, or post-result density tolerance.

Secondary fixed outputs are:

1. thermohaline seed-lens containment;
2. grid-supported SCV containment;
3. center-within-object-radius;
4. continuous peak-object overlap fraction;
5. spicy/minty and density-stratum summaries; and
6. cross-tabs against the frozen McCoy center and any-of-17 labels, with and
   without direct velocity confirmation.

The two nulls are both mandatory:

- primary: same-day, same density node, event-paired 120--240 km ring,
  area-weighted and event-equal;
- secondary: full-year wet-cell occupancy stratified by month, 1-degree
  latitude, and exact delivered density node.

The two nulls are not interchangeable. If their signs disagree, both results
and the conflict are reported; neither may replace the other after inspection.
Event-level bootstrap and paired testing are used, not particle- or voxel-level
pseudoreplication.

## Frozen daily tracking

Seed lenses, grid-supported SCVs, and cyclones are tracked separately on
consecutive daily outputs. A continuation must satisfy all of:

- center displacement at most `55.66 km`;
- relative `Cz`, radius, and thickness changes each below `0.60`;
- mutual nearest admissibility; and
- no missing-day bridging.

Splits and merges start new tracks. Tracks carry first/last model-day,
boundary-entry, and boundary-exit censor flags. Report observed durations and
the `>=3`, `>=5`, and `>=30` day classes as lower-bound descriptions. Do not
estimate complete lifetimes or formation rates from this one-year cropped
delivery.

## McCoy bridge and provenance

After the grid catalog is frozen, each object is sampled at its maximum-volume
object-day without DO. The existing center plus 0.35R/0.70R, eight-azimuth
sampling and McCoy chain are reused unchanged. Object-day and track-level
conversion rates are descriptive and cannot tune the detector.

Every run must serialize code, protocol, configuration, input-data and source
hashes, the density-node mapping, validity/censor counts, recall diagnostics,
and all three detector layers. `worker_count` and I/O concurrency are runtime
settings and are excluded from the scientific signature. A complete run is
success regardless of whether association is positive, null, or absent.

This lock is not amended by favorable or unfavorable results. Any changed
threshold, grid, density mapping, background definition, or association rule
must be a separately named sensitivity analysis.
