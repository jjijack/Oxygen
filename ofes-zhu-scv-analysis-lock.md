# OFES Zhu-SCV Analysis Lock

Frozen: 2026-08-13, before computing Zhu-SCV associations for the 56 strict
OFES events.

## Purpose and status

This document fixes the model-grid SCV analysis before its population result
is inspected. It is a pre-result method lock, not a claim that the analysis
was preregistered before any OFES exploration. The existing McCoy virtual-Argo
hits, five-day proxy, and persistent-carrier labels may be used only for
post-lock cross-method auditing; they must not be used to tune the detector.

The primary detector reproduces Zhu et al. (2024), which identified
three-dimensional SCVs in a one-year, approximately 1.5-km KOE simulation.
OFES uses a 1/30-degree z-level grid over 25.016-44.983 degrees N and
140.016-169.983 degrees E. Its meridional spacing is approximately 3.71 km;
zonal spacing varies from approximately 3.36 km at 25 degrees N to 2.62 km at
45 degrees N (approximately 3.04 km at 35 degrees N). Resolution differences
and all adaptations below must be reported explicitly.

One model year is sufficient for a daily object catalog, event association,
background occupancy, and the observed-duration thresholds used here. It is
not sufficient for an unbiased lifetime or generation-rate estimate. Zhu et
al. reported 260 tracks observed for more than 30 days in one year, while also
noting that many long-lived SCVs remained on the final output day and others
left the nested domain. OFES track durations must therefore be treated as
observed lower bounds when censored, rather than complete lifetimes.

The two simulations are process analogues, not date-matched realizations.
Zhu's nested ROMS run uses climatological surface and lateral forcing. The
NP30 output used by Hosoda et al. is a free-running, JRA25-reanalysis-forced
hindcast with no demonstrated ocean-state assimilation or event phase lock.
Neither its 2003 dates nor Zhu track dates may be matched to observed SCVs;
only methods and population properties may be compared.

## Frozen instantaneous detector

1. Collocate the native B-grid `u` and `v` to tracer centers using the already
   validated four-corner interpolation. Convert velocities to SI units.
2. Compute metric-aware horizontal derivatives at tracer centers and form
   `W = sn^2 + ss^2 - zeta^2` without using DO, temperature, salinity, SSH, or
   any existing event label.
3. Smooth W itself with Zhu's 0.1-degree by 0.1-degree centered boxcar. On the
   delivered 1/30-degree OFES grid this is a 3 by 3 kernel. A smoothed value is
   finite only when all nine source cells are finite; no coastal or delivery-
   edge padding is permitted.
4. The primary threshold is the Zhu constant
   `W0 = 5e-9 s^-2`; candidate cells satisfy `W <= -W0`. Mandatory threshold
   sensitivities are fixed at `2.5e-9` and `1.0e-8 s^-2` (0.5 and 2 times the
   primary value). The primary label and conclusion always use 1 times W0.
5. At each depth, identify eight-neighbour connected candidate regions. A
   region is a closed candidate only when it touches neither the delivered
   horizontal edge nor an invalid-cell boundary.
6. For each region, place an equal-area circle at its area-weighted centroid.
   Define shape error as the symmetric-difference area between the region and
   that circle divided by circle area. Retain shape error below 60 percent.
7. Reject a horizontal region spanning fewer than four distinct longitude
   indices or fewer than four distinct latitude indices.
8. Join regions on adjacent delivered z levels when their horizontal masks
   share at least one tracer cell. Connected components of this adjacency
   graph define three-dimensional objects. No missing-level bridge is allowed.
9. Remove every component containing the shallowest delivered level (2.5 m),
   reproducing Zhu's removal of surface-connected eddies.
10. Determine polarity from zeta at the occupied voxel nearest the object's
    volume-weighted centroid. The primary McCoy-facing catalog retains only
    Northern Hemisphere anticyclones (`zeta < 0`); cyclones remain in the
    technical catalog and are reported separately.

The 140.016-degree E delivered western edge lies immediately east of Zhu's
Izu-Ogasawara Ridge generation site near 140 degrees E. The delivered OFES
box therefore cannot support an unbiased analysis of that source region.
Horizontal-edge-touching components remain rejected as incomplete closed
objects. Tracks first seen within 0.5 degree of any delivered edge are marked
`boundary_entry_censored` and must not be counted as locally generated;
tracks last seen within 0.5 degree are marked `boundary_exit_censored`.
Birth-site and formation-rate inference, especially at the western edge, is
outside the primary analysis.

The code must pass analytic solid-body/isolated-vortex tests and reject a pure
shear field before any population association table is generated. The fixed
2003-04-05 Hosoda date may be used as a positive integration test, but its
object count or geometry must not be used to alter the frozen rules.

## Frozen z-coordinate adaptation and object properties

OFES has 75 fixed z-level centers rather than Zhu's terrain-following levels.
Layer interfaces are the midpoints of adjacent delivered depths; the uppermost
interface is fixed at 0 m and the bottom interface is extrapolated by half of
the final center spacing.

- Object volume and center are calculated with tracer-cell horizontal area
  multiplied by layer thickness.
- `Cz` is the volume-weighted depth center.
- At each occupied level, cross-sectional area uses latitude-dependent cell
  area. `Rz = sqrt(Az / pi)`.
- Zhu radius `R` is the maximum Rz over occupied levels.
- Thickness `H` is the deepest occupied lower interface minus the shallowest
  occupied upper interface, not a count of levels or a center-depth span.
- Objects containing the deepest delivered level (1076.94 m) are retained and
  marked `bottom_truncated`. They may enter the primary event-containment
  result, but are excluded from complete-H and complete-Cz summaries; a
  sensitivity result excludes them from event containment as well.

Spatial completeness is resolution limited. Zhu's cold-core mean radius of
approximately 9 km corresponds to a diameter of only about 5-7 zonal and 4-6
meridional OFES cells in this domain before the 0.1-degree boxcar is applied.
The four-by-four resolved-shape gate will therefore preferentially miss the
small end of the Zhu population. This is an a priori interpretation limit,
not permission to relax the size or W gates after seeing the detection rate.

## Frozen daily tracking

Track cyclones and anticyclones separately and compare consecutive days only.
An admissible continuation must meet all Zhu gates:

- great-circle center displacement at most 55.66 km (the physical equivalent
  of 0.5 degree latitude);
- `abs(Cz[k+1] - Cz[k]) / Cz[k] < 0.60`;
- `abs(R[k+1] - R[k]) / R[k] < 0.60`;
- `abs(H[k+1] - H[k]) / H[k] < 0.60`.

When exactly one predecessor and one successor are mutually admissible, use
the nearest-center match. If an object has multiple admissible predecessors or
successors, all members of that split/merge relation start new tracks. Missing
days are not bridged. Identity and duration remain separate outputs:
instantaneous object, tracked object, at least five consecutive days, and more
than 30 consecutive days. Neither five nor 30 days is required for the primary
instantaneous SCV identity.

Each track also carries four independent censor flags: first model day, last
model day, boundary entry, and boundary exit. Report observed duration and
the proportions censored. Do not estimate a mean complete lifetime from this
one-year, horizontally cropped delivery.

## Frozen event-to-SCV association

All association metrics use the event peak day and the formal DO50 peak daily
object. They match only anticyclonic, surface-disconnected Zhu objects from the
same day.

The primary binary estimand is `core_contained`: the event's maximum-DeltaDO
horizontal cell and nearest delivered peak-depth level belong to a Zhu
three-dimensional mask. It has no tunable search radius or overlap threshold.

Two fixed secondary descriptions must be reported rather than selected after
the result:

1. `center_within_Rz`: great-circle distance from the DeltaDO core to the Zhu
   object's area-weighted center at the event peak level is at most that
   cross-section's Rz. If several objects qualify, choose the smallest d/Rz.
2. `peak_voxel_overlap_fraction`: area-weighted fraction of all member pixels
   in the formal DO50 peak daily object whose horizontal cell and nearest
   diagnosed peak-depth level fall inside any same-day Zhu mask. Report this
   continuously; do not introduce a post-result overlap cutoff.

Report all three metrics, their agreement, and bottom-truncation sensitivity.
The headline prevalence and enrichment use `core_contained` only.

Surface disconnection is a detector property, not proof of altimetric
invisibility. Removal of objects whose OW mask reaches the shallowest level
shows subsurface intensification, but it does not establish weak SSH, weak
surface velocity, or META non-detection. Those require the separately planned
surface-eddy comparison and cannot be inferred from the Zhu gate alone.

Hydrographic sign is also separate from polarity. Both Zhu warm-core/spicy
and cold-core/minty classes may be anticyclonic. For event and Zhu-to-McCoy
tables, retain continuous core sigma0 and report the literature-anchored
density strata 26.2-26.8, 26.8-27.0, and 27.0-27.3 kg m^-3, plus an outside-
range category. Also stratify the McCoy bridge by its existing spicy/minty
label. These are descriptive strata only and never SCV detection gates.

## Frozen background estimands

The primary null is a same-day, same-nearest-depth, event-paired spatial null.
For each event, calculate Zhu-mask occupancy across every valid tracer cell
120-240 km from the event core, weighted by cell area. Exclude pixels belonging
to the focal DO50 daily object only; retain other naturally occurring DO or SCV
features. The event is the statistical unit. Report the event-equal mean of
`core_contained - ring_occupancy_fraction`, an event bootstrap confidence
interval, and a paired test. This is the paper's primary background number.

The secondary sanity check is the full-year wet-cell Zhu occupancy stratified
by calendar month, one-degree latitude bin, and exact delivered z level. Map
each event to its corresponding stratum and compare observed event prevalence
with expected occupancy. It is not interchangeable with the paired ring null
and is not used to replace it.

If the two nulls point in opposite directions, report the conflict and inspect
spatial/seasonal support. Do not select the more favorable number. A robust
general enrichment claim requires the signs to agree; otherwise retain the
paired-ring estimand as the predefined local contrast and state that the
broader occupancy comparison does not corroborate it.

## Frozen Zhu-to-McCoy bridge

McCoy profiles cannot tune Zhu detection. After the complete primary Zhu
catalog is frozen:

- For each Zhu track, choose the maximum-volume object-day without reference
  to DO. Sample the center plus eight azimuths at 0.35R and eight at 0.70R,
  then run the existing McCoy-method-compatible profile chain unchanged.
- Report trajectory-equal center and any-of-17 conversion rates. All
  object-days may be reported only as a descriptive secondary result.
- At the 56 formal event peaks, cross-tabulate Zhu core containment against
  the already frozen center and any-of-17 McCoy results, with and without the
  existing direct-velocity confirmation.
- Always retain all four Zhu/McCoy cells. Do not redefine grid SCV to force
  recovery of the six center or eleven footprint McCoy-plus-velocity events.

The interpretation is complementary rather than mathematically equivalent:
McCoy detects a vertically confined hydrographic/dynamic-height signature from
profiles, whereas Zhu detects a resolved three-dimensional vorticity-dominant
subsurface object. `Zhu+ / McCoy+` is joint confirmation; either single-positive
cell remains a scientifically interpretable method discrepancy.

## Execution and provenance gate

This lock does not alter `processing.yml` and therefore cannot contaminate the
current OFES run signatures before code exists. During implementation, every
value and definition above must be copied to a versioned config section and
serialized verbatim into the detector manifest. The implementation commit,
analytic validations, source-code hash, processing-config hash, and the locked
Zhu citation must exist before the 365-day catalog or 56-event association is
run. Population results do not authorize amendments to these gates; any later
change creates a separately tagged sensitivity analysis.
