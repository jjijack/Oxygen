# OFES surface-eddy event-association lock

Status: pre-registered before reading any authoritative DO50 event result.
Lock date: 2026-08-14 (Asia/Shanghai)

This is the event-level companion to `ofes-surface-eddy-analysis-lock.md`.
It does not alter the PET detection or tracking signature. The PET catalog
must first have a complete 365-day manifest and complete tracking outputs.

## Inputs and identity gate

- PET manifest: explicitly passed `plot_outputs/do/ofes_np30_ke/surface_eddy/runs/<run_id>/manifest.json`.
- DO50 population root:
  `/mnt/w2/scratch/user3/Oxygen-ofes/plot_outputs/do/ofes_np30_ke/ofes_delta_do_catalog/20030101_20031231_cf957935d38a/`.
- Event-diagnostics root:
  `/mnt/w2/scratch/user3/Oxygen-ofes/plot_outputs/do/ofes_np30_ke/event_diagnostics/ofes_events_21efbe902ab7/`.
- McCoy-compatible event summary is selected only by an explicit path or one
  unambiguous `event_summary.parquet` below the diagnostics root.
- Every input file or explicit directory inventory is SHA-256 recorded in the
  association summary. No newest/lexicographically-last run is selected.
- The population authority is `population_peak_diagnostics.parquet`.
  Candidate selection is exactly `threshold == 50` (59 rows), followed by
  `population_diagnostic_passed == True` (strict 56 rows). Event ID, peak date,
  peak longitude, peak latitude, and daily object key are retained.
- The DO daily catalog must contain exactly one peak-day row per strict event,
  and its peak date/longitude/latitude must match the population authority.

## Primary PET association

- On the event peak day, use actual PET detections only (`is_virtual == False`).
- A PET row is association-eligible only when `filter_valid == True` and
  `boundary_censored == False`; no virtual observation enters the primary
  numerator.
- The primary indicator is whether the event peak-core grid-cell center is
  contained by any eligible actual PET effective contour, with `any polarity`
  reported first.
- Secondary fields include speed-contour containment, DO peak-footprint area
  overlap with effective/speed contours, nearest effective-contour/center
  distance, and distance divided by effective radius.
- The primary event denominator reports the strict 56, the filter/ocean-valid
  core count, the complete-background-ring eligible count, and the actual PET
  containment count separately.

## Local and annual nulls

- The primary local null is the same-day 120--240 km ring around each event
  core, using only filter-valid and ocean-valid cells. Event peak-footprint
  cells are excluded when available from the authoritative peak-pixel table.
- Ring occupancy is computed as a cell-area-weighted fraction and compared at
  the event level, never treating ring cells as independent events.
- Report the event-level mean of `core indicator - ring occupancy`, a seeded
  10,000-replicate event bootstrap 95% CI, and a paired Wilcoxon test.
- The secondary sanity null is the existing month × 1-degree latitude annual
  PET occupancy framework; it is reported separately and never substituted for
  the local ring estimate.

## Rotation and McCoy cross-tabs

- The rotation subset is exactly the existing 29 rows selected by
  `rotation_dominated == True` in the authoritative population diagnostics.
- Deep polarity is derived from the authoritative deep `rossby_number` using
  the existing convention: negative = anticyclonic, positive = cyclonic.
- Preserve the existing core-weighted surface-Ro expression and
  `surface_core_rotation_polarity_match`; do not redefine it using PET.
- PET rotation classes are `same`, `opposite`, and `no-PET`; mixed PET
  polarities or missing deep polarity are `ambiguous`.
- McCoy event-level fields must reproduce the locked counts:
  `center_profile_mccoy_compatible` 9/56,
  `center_profile_velocity_confirmed` 6/56,
  `any_event_profile_mccoy_compatible` 19/56, and
  `any_event_profile_velocity_confirmed` 11/56. A count mismatch aborts.
- The 161-row `threshold == 50` quality-eligible sensitivity population is
  reported separately and never merged into the strict 56 denominator.

## Fixed plots and outputs

The three representative cases are fixed to
`OFES_DO50_E000002`, `OFES_DO50_E000239`, and `OFES_DO50_E000176`; all three
must be drawn. Each case overlays the DO peak footprint/core, PET effective and
speed contours, PET center/polarity, and nearest-object normalized distance
when no effective containment exists.

The output names are fixed:

- `surface_eddy_event_association.parquet`
- `surface_eddy_event_nulls.parquet`
- `surface_eddy_surface_ro_crosstab.parquet`
- `surface_eddy_mccoy_crosstab.parquet`
- `surface_eddy_quality_eligible_161.parquet`
- `surface_eddy_summary.json`
- `surface_eddy_event_association_manifest.json`

## Frozen control hashes

- Detection lock is the current `ofes-surface-eddy-analysis-lock.md`; its
  SHA-256 is recorded in the PET detection manifest.
- Event postprocess config:
  `config/ofes_surface_eddy_postprocess.yml`
- Event postprocess config SHA-256:
  `3f0505cd292c25eb8964221d8e2ed81681732d689668a1eddc94074ba2474cde`
