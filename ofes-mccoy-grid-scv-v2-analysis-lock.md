# OFES McCoy-grid SCV detector v2 analysis lock

Frozen: 2026-08-14, before implementing the v2 detector and before inspecting
any v2 56-event association result.

## Status: prospective repair, not preregistration

The v1 detector (`ofes-mccoy-grid-scv-analysis-lock.md`, protocol hash
`fec8df0b20d110af05f212b60a5a1316576a89524c3dee6d9fcd96e7d51daaad`) was
audited on 2026-08-14 and failed on six counts, recorded in the project memory
note `grid-scv-surface-eddy-audit-state`:

1. the grid seed masks only approximated the McCoy chain (spice IQR + N2 IQR +
   PV local minimum) and `prefilter_recall` was only a finite-value mask, so
   the 141/141 recall number did not mean McCoy-chain recall;
2. the multiplicative identity gate (5--50 km and >=3 nodes and >=100 m) was
   implemented as a single binary `eligible_identity`, mixing resolution and
   morphology;
3. the cyclonic technical catalog in the v1 lock was never implemented (the
   velocity gate hard-codes anticyclonic signs);
4. the Nencioli implementation deviated from the frozen wording in four
   places (no constant rotation sense, an extra tangential-alignment gate, an
   area-weighted candidate center, and circulation on a circle instead of the
   closed N2 boundary);
5. the analytic validation implemented roughly 3.5 of the 8 frozen checks and
   wrote nothing to disk;
6. the annual catalog never completed (5/365 days, mixed fragment hashes).

This document therefore fixes a **prospective repair**: it claims to be a
repair lock written after the v1 audit failure, not the project's original
preregistration. Every rule below is frozen before any v2 56-event
association number is generated or inspected. The existing OFES virtual-profile
counts (9/56 center profiles, 19/56 any-event profiles, 6/56 and 11/56
velocity-confirmed) remain frozen comparison targets; they may not tune any
v2 gate.

The detector is deliberately layered:

- Tier 0 `mccoy_profile_seed`: full McCoy profile chain passes on a grid column;
- Tier 1 `closed_thermohaline_lens`: seeds grown into a three-dimensional
  closed-N2 thermohaline lens (the primary cross-dataset estimand);
- Tier 2 `native_anticyclonic_support`: OFES-native velocity confirmation,
  split into weak and strong; the strong subset is the only object that may be
  called a strict grid-supported SCV;
- Tier 3 `persistence`: tracking. Descriptive only, never identity.

Naming discipline (paper-safe):

- Tier 1 objects: "SCV-related thermohaline lens";
- Tier 2 weak: "partially velocity-supported anticyclonic lens";
- Tier 2 strong: "strict grid-supported SCV";
- cyclonic technical-catalog entries: "cyclonic thermohaline lens (technical)".

A Tier-1 object must never be called a strict SCV.

## Rejected v1 inputs

The following v1 outputs are **rejected** as formal science products and are
explicitly not read by any v2 step. They are retained on disk for provenance
and are marked with `SUPERSEDED.txt`:

- `event_catalog_final/` -- pre-fix code snapshot hash `684fd465` not present
  in git history; contains density-node counts of 143 (39 nodes exist) and
  radii to 103.4 km, i.e. a bug state.
- `annual_20030101_20031231_localrecert/` and
  `annual_20030101_20031231_localrecert_current/` -- incomplete (5/365 days)
  with mixed daily-fragment hashes and a `running` manifest.
- `/tmp/ofes_grid_scv_*` -- transient engineering snapshots and logs.

The authoritative v1 event-level result
`event_catalog_final_current/` (code `81868b60` = HEAD `3b29e28`, 56/56
complete, 0/56 identity and 0/56 grid-supported containment) is retained as
the historical v1 baseline for cross-version reporting only.

## Frozen source and coordinate conventions

- The McCoy profile chain must call the existing OFES implementation
  (`_ofes_mccoy_profile_properties`, `_ofes_mccoy_reference`,
  `_ofes_mccoy_anomalies`, `_ofes_mccoy_iqr_on_profile`,
  `_ofes_mccoy_gaussian_fit`, `_ofes_mccoy_dynamic_height_gate`, and
  `_ofes_mccoy_profile_classification`) unchanged, using the frozen OFES port
  parameters under `config/processing.yml` -> `mccoy_virtual_argo` and the
  archived public MATLAB source identified by
  `source_archive_sha256 = ee0abca349c4f090b3c99a695a88ee81df418fd3f756110d165f83b693b26406`
  (`data/mccoy2020_scv/ArgoSCVs-v1.3.zip`).
- OFES temperature and salinity are converted with TEOS-10. The detector uses
  `sigma0`, spiciness, `N2`, and native `u/v`; it never reads DO, event masks,
  event intensity, or event IDs in any detector gate.
- Fixed density nodes: `25.60, 25.65, ..., 27.50 kg m^-3` (step `0.05`,
  39 nodes). Each water column is interpolated to these nodes only where the
  crossing is unique, monotonic, finite, and bracketed by delivered z-levels;
  no extrapolation. Valid interpolated depth window: `100--1000 m`.
- Continuous event `target_sigma0` is mapped to the fixed node with the frozen
  Voronoi cell rule (midpoint edges, `searchsorted(..., side='right')-1`;
  the v1 lock's "lower node owning a tie" wording was not matched by the
  implementation, which resolves a midpoint tie upward; the mismatch was
  audited as zero impact on the 56 events and the implementation rule is the
  frozen one).
- Horizontal areas use the delivered tracer-cell metric.

## Tier 0: McCoy profile seed

A Tier-0 seed is a grid column at tracer-cell center `(lat_j, lon_i)` on a
given day where the full `_ofes_mccoy_profile_classification` chain passes
(`mccoy_profile_compatible == True`).

- Column extraction: the tracer column `temp[:, j, i]`, `salinity[:, j, i]`
  from the same-day snapshot, passed through the same
  `_ofes_mccoy_profile_properties` path as the virtual-Argo pipeline (same
  10-dbar grid, same 10-dbar kernel smoothing, same QC gates).
- Background pool: the same-day, same-position 120--240 km ring, sampled at
  the frozen virtual-Argo geometry -- radii `130, 155, 180, 205, 230 km` times
  16 azimuths = 80 control positions (`_ofes_mccoy_sampling_points` with the
  same settings). The reference and the per-density IQR pool are built from
  those controls exactly as in the virtual-Argo pipeline. A column whose ring
  contains fewer than `minimum_background_profiles` (61) QC-passing controls
  fails at the `background_iqr` stage and is recorded, never raised.
- The seed records: spicy/minty (`scv_type`), Gaussian center/height/
  amplitude/`R2`/`NRMSE`, weak-N2, dynamic-height fields, and the full
  `mccoy_failure_stage` at every failure. No failure point may collapse into a
  bare boolean.
- A column that fails the full chain cannot be promoted to a seed by PV,
  vorticity, or native velocity. PV and Ro never create a seed.

## Tier 1: closed thermohaline lens

Starting from Tier-0 seeds:

1. At the seed's mapped `sigma0` node, find the outermost closed contour of
   the `N2` anomaly field at level 0 (Barabinot et al. 2025 boundary) whose
   interior contains the seed cell. The contour machinery follows the
   deterministic grid implementation already used in v1
   (`_ofes_grid_closed_n2_masks`), with the seed as the anchor instead of the
   thermohaline component.
2. The contour interior is the layer mask. Because the seed itself carries the
   same-sign spice anomaly and the weak-N2 core, the interior contains them by
   construction; the implementation must verify and record the count of
   same-sign spice-anomaly voxels inside the mask and fail the layer if it is
   zero.
3. Extension to an adjacent `sigma0` node requires all of:
   - same spice sign as the seed (`scv_type` maps to the node-space anomaly
     sign: spicy = positive, minty = negative);
   - the adjacent layer has at least one same-sign spice-anomaly voxel
     enclosed by a closed `N2` anomaly=0 contour (its own layer mask);
   - the two layer masks intersect, or intersect after a single 8-neighbour
     one-cell dilation (grid discretisation at about 3.7 km and slight tilt);
     dilation is never applied twice;
   - no missing-density-node bridging: extension only between adjacent
     `sigma0` nodes, never across a gap.
4. Seeds of the same spice sign that fall inside the same connected
   three-dimensional lens merge into one object. Opposite-sign seeds never
   merge.

Identity is a morphology/resolution classification, not one multiplied gate:

- `profile_only`: has a Tier-0 seed but no closed lens on at least two
  adjacent density nodes;
- `grid_lens`: at least two adjacent nodes form a closed lens; volume-weighted
  centroid depth > `200 m`; shallowest occupied interpolated depth >
  `100 m`; equivalent radius (max over layers of `sqrt(area/pi)`) in
  `5--100 km`;
- `well_resolved_grid_lens`: a `grid_lens` with at least three nodes and
  physical vertical thickness at least `100 m`;
- `underresolved`: radius < `5 km` or a single occupied node; recorded, not a
  lens;
- `broad_structure`: radius > `100 km`; recorded, not an SCV lens;
- `boundary_censored`: touches the horizontal delivered boundary, the `100 m`
  upper bound, or the delivered bottom; retained with the censor flag and
  excluded from strict proportions.

Classical-size subset: radius `5--50 km` (reported separately as
classical-size). Large-lens subset: radius `50--100 km` (reported as
large-lens/ITE-like). Both are subsets of the same `grid_lens` identity.

PV, spice, and `N2` anomalies are saved as continuous fields at every node.

## Tier 2: native anticyclonic support

For each `grid_lens`, use OFES-native `u/v` (four-corner tracer interpolation,
SI conversion, as in v1) with the same-day 120--240 km ring translation
removed, then classify support in two levels. No fixed absolute Rossby-number
or Okubo-Weiss threshold is used anywhere; PV, Ro, strain, and circulation are
all saved as continuous variables.

`weak_native_support` requires:

- at least one occupied node with anticyclonic (Northern-Hemisphere negative)
  signed circulation on its outermost closed `N2` boundary, **or** at least
  one occupied node that passes the Nencioli vector geometry; **and**
- the object-core area-weighted `zeta` (relative vorticity) is anticyclonic
  (negative).

`strong_native_support` requires all of:

- the Nencioli vector geometry passes;
- the closed-`N2`-boundary circulation is anticyclonic;
- both agree on at least two adjacent occupied density nodes;
- finite boundary samples at least 80%.

Nencioli vector geometry (frozen wording, replacing the v1 deviations):

- offsets `a = 3`, `b = 2` native grid cells;
- along each cardinal axis, the absolute component magnitudes must be
  non-decreasing over offsets `1..a` and the component signs must reverse
  across the candidate center;
- the center speed must be the finite minimum in its `(2b+1)^2` neighborhood;
- the 8 vectors on the square perimeter at offset `b` must traverse all four
  quadrants with one constant sense of rotation, with no pair separated by
  more than one quadrant, and with at least 6 of 8 tangential alignments
  positive (the alignment gate is kept as an explicit recorded gate, separate
  from the four Nencioli conditions, and is reported per node);
- the candidate center is the valid velocity minimum nearest the lens volume
  centroid and must lie inside the thermohaline mask.

Circulation: 64 equally spaced samples on the outermost closed `N2` boundary
(at the contour vertices, not on a circle), signed anticyclonic (negative) in
the Northern Hemisphere, at least 80% finite samples. The v1 circle sampling
is replaced by contour sampling.

A lens that passes Tier 1 but fails Tier 2 remains in the Tier-1 catalog and
is never silently discarded.

Cyclonic technical catalog: objects whose lens identity passes but whose
`N2`-boundary circulation and/or core `zeta` are cyclonic are retained in a
separate technical catalog under the same gates (same Nencioli conditions,
cyclonic sign). The primary SCV table is the Northern-Hemisphere
anticyclonic subset.

## Tier 3: persistence

Persistence is a descriptive grade, never an identity gate. Track on
consecutive model days only, for each of the three detection families
(seed lens, weak-native, strong-native) separately:

- same spice sign;
- center displacement at most `55.66 km`;
- relative radius, thickness, and centroid-depth changes each at most `0.60`;
- adjacent days only, no missing-day bridging;
- splits and merges start new tracks;
- boundary censor flags retained (`0.5` degree margin at the delivered
  boundary, year start/end).

Report the `>=3`, `>=5`, and `>=30` day classes as lower-bound descriptions.
This one-year cropped delivery cannot estimate complete lifetimes or
formation rates.

## Prefilter and compute control

The full McCoy chain cannot run on every OFES column. A broad prefilter
produces the candidate pool, with two hard properties:

1. **Allowed prefilter inputs only**: valid unique `sigma0` crossing; the
   broad spice/N2 anomaly envelope; pycnocline-below restriction; finite
   vertical data coverage. PV local minimum, Ro, and DO are forbidden inputs.
2. **Recall definition**: "does every positive of the full McCoy classifier
   enter the candidate pool?" The finite-value mask of v1 is replaced.

The recall-gated prefilter gate (provable superset of the McCoy gates): the
column's profile properties pass QC; and there exists a density level below
the reference pycnocline where the profile-space spice anomaly is positive or
negative while the profile-space `N2` anomaly is negative within the frozen
`n2_density_tolerance_kg_m3` (0.2) of that level. Any McCoy-positive column
satisfies this by the `initial_spice_n2_candidate` and pycnocline gates, so
recall = 1 is a theorem; the audit below verifies the theorem on real data.
A cheap node-space shortlist may be used to skip profile construction for
quiet columns, but that shortlist is itself audited for positive retention
with the same 100% requirement, and both stages may only be widened for
recall.

Engineering gates (all before any DO association):

- the frozen 141 McCoy-compatible event profiles
  (`virtual_profile_diagnostics.parquet`, `sample_role == 'event'`,
  `mccoy_profile_compatible == True`) must be retained 141/141 by the
  prefilter at their mapped cells, and the full classifier at those positions
  must reproduce the stored classification (see controls);
- audit sample: 2003-01-15, 2003-04-15, 2003-07-15, 2003-10-15. On each date,
  draw 1000 QC-passing valid columns stratified by 2-degree latitude band and
  100-m band of the column's median valid-node depth, uniform per stratum,
  fixed random seed `20260814`, without replacement. Run the full classifier
  on these columns **without** the prefilter. Any McCoy positive among them
  must be in the candidate pool (100% recall); their Tier-0 results and
  failure stages are recorded.
- the prefilter may only be widened for recall; it may never be adjusted
  according to 56-event overlap.

## Analytic and engineering validation

All 14 checks below must run before any population association and every
result must be written to `detector_validation.json`, which is part of the
manifest gate. A failed check requires implementation repair or an explicit
separately named sensitivity; it never authorizes changing this lock after the
56-event table is viewed.

1. stationary ocean with no thermohaline anomaly: Tier 0--2 all zero;
2. tilted three-dimensional Gaussian spice + weak-N2 lens: Tier 0 and Tier 1
   pass;
3. the same lens shifted one tracer cell between layers: the one-cell
   dilation connects the layers;
4. a shift of more than one cell with no boundary overlap: layers do not
   connect;
5. pure heave without an along-isopycnal thermohaline anomaly: no Tier-0 seed;
6. an edge shear patch with high vorticity/OW but no thermohaline seed: no
   Tier-1 lens;
7. a thermohaline filament without a closed N2 contour: `profile_only`, not a
   `grid_lens`;
8. solid-body anticyclone: strong native gate passes;
9. solid-body cyclone: enters only the cyclonic technical catalog;
10. pure shear: no native anticyclonic support;
11. adding a uniform translation: object and support results identical after
    background removal;
12. surface-connected lens: technical record retained, strict subsurface
    identity rejected;
13. two adjacent lenses: no false merge from contour dilation;
14. PV units/sign, isopycnal interpolation, contour closure, volume, and
    circulation pass independent analytic checks.

## Positive and negative controls

Written before any DO association, under the v2 results root:

1. `seed_control_audit.parquet` -- for each of the 141 frozen McCoy positives:
   stored classification (reference), prefilter retention at the mapped cell,
   Tier-0 reproduction at the exact sample position (bilinear profile, frozen
   per-position ring) and at the nearest tracer-cell column, Tier-1
   extension, and the explicit reason when not extended. Tier-0 must reproduce
   141/141 at the exact positions and at the mapped cells; any failure is
   reported with its full failure stage and is a stop-and-report gate.
2. `matched_background_control.parquet` -- the frozen 80 per-event
   120--240 km background virtual profiles (4480 total). Grid Tier-0 is run
   at each control position (per-position ring) and compared with the stored
   control classification; report the Tier-0 false-positive rate (grid
   positive among stored negatives), the overall positive rate, and the
   failure-stage distribution. Background hits never delete results.
3. `tier_transition_summary.json` -- the attrition chain: 141 known positives
   -> prefilter retained -> exact Tier-0 -> profile_only/Tier-1 -> weak
   native -> strong native. The 9/56 center, 19/56 any-event, and 6/56 and
   11/56 velocity-confirmed counts must reproduce at the Tier-0 layer
   (exact-position classification); their Tier-1/Tier-2 upgrade rates are
   reported, not required.

## Run discipline and freeze ordering

- All v2 outputs live under the external directory
  `/mnt/w2/scratch/user3/Oxygen-cache/ofes_grid_scv_v2_results/`; v1 outputs
  are never overwritten.
- Required outputs: `manifest.json`, `detector_validation.json`,
  `seed_control_audit.parquet`, `matched_background_control.parquet`,
  `tier_transition_summary.json`, `object_inventory.parquet`,
  `object_voxels.parquet`, `event_association.parquet`,
  `gate_attrition.json`, `analysis_summary.json`, plus per-day
  day-level parquets.
- Daily runs write atomically. Resume accepts a day only when it is
  `complete` and carries the same code/protocol/config hashes. Any mixed-hash
  day is rejected, not reused.
- `worker_count` and I/O concurrency go only into execution provenance, never
  into the scientific signature.
- Data is read only through `./data/OFES_NP30` (the softlink root registered
  in `config/paths.yml`).
- `track.py` is not modified during a formal 56-event run or the annual run.
- Ordering is strict: (1) v2 lock written and committed; (2) detector
  implemented; (3) analytic validation and positive/background controls
  pass; (4) only bug fixes are allowed at this stage, never frozen-rule
  changes; (5) implementation committed with code/protocol/config/source
  hashes recorded; (6) then, and only then, the 56-event catalog is run; (7)
  the 365-day annual catalog and its month x 1-degree x sigma0 occupancy null
  start only if the 56-event catalog shows a nonzero and reasonable Tier-1
  background occupancy and no hash mixing.

## 56-event association and nulls

The association uses the frozen 56 strict DO50 events
(`population_peak_diagnostics.parquet`: threshold 50,
`population_diagnostic_passed`, `background_annulus_within_delivery_window`),
their formal peak day, peak cell, and `target_sigma0` mapped to the fixed
density node before results are inspected.

DO enters only at this stage. Per event and per Tier, save:

- whether the event core voxel is contained by an object of that Tier;
- event-core to object-center distance divided by object radius (`d/R`);
- the continuous overlap fraction of the event peak mask with the object
  mask;
- the `target_sigma0` mapped node and nearest-object distance;
- the four-level cross-tab: `profile_only`, Tier 1, weak-native,
  strong-native.

Primary estimand: Tier-1 `closed_thermohaline_lens` core containment.
Strict secondary: `strong_native_support` (strict grid-supported SCV)
containment.

Both nulls are mandatory and reported side by side; if their directions
disagree, both results and the conflict are reported:

- primary: same-day, same density node, event-paired 120--240 km ring,
  area-weighted and event-equal;
- secondary: full-year wet-cell occupancy stratified by month, 1-degree
  latitude, and exact delivered density node (computed by the annual catalog).

Event-level bootstrap (10,000 replicates, seed `20260729`) and paired
Wilcoxon are used; both the predeclared one-sided (`greater`) result and the
transparent two-sided result are saved and reported for every Tier. No Tier
is reported alone because its number looked best; all tiers are reported
together.

## Interpretation and stop conditions

- Result A: Tier-0 and Tier-1 controls are healthy and Tier-1/Tier-2 event
  overlap exceeds the background. Then OFES contains thermohaline/anticyclonic
  carriers continuous with the McCoy SCV spectrum; only the strong-native
  subset may be called strict grid-SCV.
- Result B: the positive controls form Tier-1 objects but the 56-event
  overlap is near zero and not above background. Accept the scientific null:
  these DO event peaks are usually not mature closed SCV cores; the OFES
  mechanism story becomes multi-mechanism (frontal/strain formation with
  anticyclonic organization).
- Result C: Tier-0 succeeds but most positive controls never form Tier-1
  objects and the annual background has almost no deep objects. Conclude the
  1/30-degree OFES delivery cannot test three-dimensional strict SCV;
  write "grid coherence unresolved", never "SCV absent".
- Result D: Tier-1 objects are numerous but native velocity barely supports
  them. Use "SCV-related thermohaline lens / weakly resolved carrier", never
  "strict SCV"; examine filaments, tilted water masses, and velocity
  resolution limits.

No gate may be changed on the basis of DO overlap. Any later rule change is a
named v2 sensitivity, with the main result retained.

## Commit plan

One compact local commit per node (no push):

1. `Stabilize grid SCV run state` -- done: verified-identical
   `_ofes_grid_connect_components` optimization (`d2c356e`), SUPERSEDED
   markers on rejected v1 outputs;
2. `Lock grid SCV detector repair` -- this document;
3. `Build McCoy-seeded grid lenses` -- Tier 0 + Tier 1 implementation;
4. `Validate grid SCV support tiers` -- analytic validation and controls;
5. `Evaluate grid SCV event overlap` -- 56-event association and nulls;
6. `Track grid SCV lifecycles` -- Tier 3 tracking, only after the 56-event
   gate passes.

Each commit runs `py_compile`, a diff check, and a small focused validation;
track.py/config/notebook edits follow the post-edit review checklist; config
changes bump the processing changelog; notebooks consume only
`status == complete` and hash-matched results, never `sorted(glob)[-1]`.

## Frozen input hash table

| Input | SHA-256 |
|---|---|
| McCoy MATLAB archive `data/mccoy2020_scv/ArgoSCVs-v1.3.zip` | `ee0abca349c4f090b3c99a695a88ee81df418fd3f756110d165f83b693b26406` |
| v1 protocol `ofes-mccoy-grid-scv-analysis-lock.md` | `fec8df0b20d110af05f212b60a5a1316576a89524c3dee6d9fcd96e7d51daaad` |
| `config/processing.yml` | `ed3e5bd164cdd66f5bd6453daa27f6c5f885ce6456d58401fe5adfeeb6552213` |
| `population_peak_diagnostics.parquet` (56-event authority) | `b91435f2fb3e6a7796e184058fb39a91f8a8fd1c7c3e706144fa12e974fe1fde` |
| `virtual_profile_diagnostics.parquet` (141 positives + controls) | `3a2b4e747ef7c5df044ab1f9ed9f0a76d4c727f1ab0845e71ca5d6209897cddf` |

This lock is not amended by favorable or unfavorable results.
