# Oxygen Copilot Instructions

> **Canonical source: `CLAUDE.md`** — This file is kept for Copilot compatibility. In case of conflict, `CLAUDE.md` takes precedence.

## Repo Snapshot
- Majority of the logic lives in `track.py` covering geometry helpers, data ingestion, analytics, and plotting for eddy–Argo interactions.
- Configuration is entirely YAML-driven under `config/` (`paths.yml`, `processing.yml`, `regions.yml`); do not hardcode paths or thresholds.
- Data products are materialized into local folders referenced in `config/paths.yml` (e.g., `Argo_data/`, `META_tracks/`, `plot_outputs/`, `external/natural_earth/`).
- There are no tests; validate changes by running the specific workflow function you touched.

## Notebook Workflow & Comments
- Author new or updated functions in `track.py`, then import them into `GLORYS.ipynb` for experimentation; treat the notebook purely as a consumer of `track.py` logic.
- When updating docstrings or comments after making changes, describe the function for a first-time reader by focusing on its current behavior rather than a changelog, and ensure any functional modifications are reflected by synchronizing the function’s description accordingly.

## Config & Regions
- Always read paths via `_PATHS_CFG` accessors at the top of `track.py`; new features should expose knobs through the YAML instead of constants.
- Region bounds (`lonmin`, `lonmax`, `latmin`, `latmax`) are globals populated by `_load_region_config`; if you do not call `switch_region(...)`, code uses `default_region` from `config/regions.yml` (currently `global`). Call `switch_region('global')` (or another key) before workflows that depend on spatial filtering.
- In parallel code paths (e.g., Dask or `multiprocessing`), do not assume that region-related globals are preserved in worker processes—they may fall back to defaults. Make sure to reinitialize or propagate the required global state within each worker (e.g., by calling `switch_region(...)` or equivalent setup).
- Regions may cross the dateline; rely on `_region_lon_mask`, `_minimal_lon_diff_deg`, or `adaptive_distance_m` instead of naive `lon` comparisons.
- Basemap styling derives from `processing.yml:plot.basemap_colors`; reuse `_BASEMAP_COLORS` when adding new plots.

## META Pipelines
- Use `load_meta_data()` to open the official META 3.x NetCDF files defined in `config/paths.yml`.
- Convert META NetCDF to parquet/zarr via `export_meta_tracks(...)`; it streams by chunks, can run with `use_dask=True`, and writes under `META_tracks/<region>/<kind>_*`.
- `find_track(kind, track_id, include_contours=True)` is the canonical accessor for META tracks; it automatically fuses parquet daily records with contours from the matching zarr store.
- Prefer kind-string driven calls (`'acs'|'acl'|'cs'|'cl'`) across new workflows (`find_track`, `plot_track`, etc.); keep legacy list-based DS inputs only for backward compatibility.
- Keep chunk sizes and Dask worker counts configurable when extending `export_meta_tracks`—the code already autoscales chunking based on dataset length.

## Argo Ingestion
- Legacy `.mat` files go through `convert_mat_to_parquet(year, input_dir, output_dir)`; newer `.txt` drops are ingested by `process_argo_txt_to_yearly_parquet_dask(...)` which maps each file to parquet with Dask.
- `process_argo_txt_to_yearly_parquet_dask` is configuration-driven for origin/intermediate/output directories; in current code, txt input prefers `paths.argo_txt_input` and falls back to `./Argo_origin` when unset. Keep key naming consistent when evolving `paths.yml`.
- `load_argo_data(year, variable_selection=...)` is the normalized reader; it standardizes column names and picks adjusted variables. Extend variable handling by updating the `default_selection` dict.
- Massive joins should work on yearly parquet files instead of loading everything in memory—follow the pattern in `filtered_float_data` (lazy loading only the needed years).

## GLORYS & Vertical Diagnostics
- GLORYS fields are pulled with `get_track_area_glorys(...)` and interpolated via `get_vertical_glorys`; plotting entry points are split across track/Argo horizontal and vertical helpers (`plot_track_horizontal_glorys`, `plot_argo_horizontal_glorys`, `plot_track_vertical_glorys`, `plot_argo_vertical_glorys`, plus overview variants).
- Prefer passing kind strings (`'acs'|'acl'|'cs'|'cl'`) plus `track_id` and snapshot index; legacy ACS/ACL/CS/CL list-style inputs remain supported for compatibility.
- When adding new 2D/3D variables, update the `alias_map` + `var_dims` dictionaries near `get_vertical_glorys`; ensure interpolation uses `RegularGridInterpolator` with masked arrays, matching existing error handling.
- GLORYS diagnostic plots that depend on anomaly detection/projection config write under `plot_outputs/<method>/<region>/...`; method-independent replay plots from `plot_data_package(..., save_fig=True)` write under `plot_outputs/shared/<region>/plot_track_vertical_glorys/`. Keep filenames consistent with the current `{dataset}{id}_vertical_{var}_YYYYMMDD_k*b*.png` template so downstream scripts can glob predictably.

## Float/Eddy Analytics
- `filtered_float_data` is the shared matcher between META tracks and Argo profiles; it first date-joins, then applies both polygon containment (`is_point_in_contour`) and radius checks (`adaptive_distance_m`). Reuse these helpers whenever proximity logic is needed.
- Argo anomaly detection lives in `calculate_delta_do` and is configured through `DetectionConfig` / `make_detection_config(...)`; detection defaults and plotting defaults come from `processing.yml`. When changing defaults, expose them there and call `print_current_processing_defaults()` during debugging.
- `plot_track` orchestrates everything—loading a track, matching floats, computing Argo anomalies, and plotting; new visual embellishments should hook into the existing Matplotlib legend/grid patterns and honor `plot_radius`, `plot_unrelated_argo`, and colorbar settings.
- Regional anomaly maps are produced via `plot_argo_hotspots`, which spins up Dask tasks per year, writes plots to `plot_outputs/<method>/<region>/plot_argo_hotspots/`, and can persist anomalies to parquet; keep additions streaming-friendly.
- Method-independent outputs, such as baseline Argo/eddy interaction parquet files and generic profile/relative-position diagnostics, belong under `plot_outputs/shared/<region>/...`, not directly under `plot_outputs/<region>/...`.

## Geospatial Conventions
- Distances default to planar approximations adjusted by latitude (`approximate_degree_length`); switch to true great-circle calculations by toggling `force_great_circle` or the adaptive thresholds instead of reimplementing Haversine.
- Always normalize longitudes through `_normalize_lon_array` before polygon math, especially for datasets that straddle ±180°.
- Circular overlays on maps should be drawn via `ellipse_patch_for_eddy`, which already compensates for cos(lat) distortion.
- Use `local_xy_distance_m`/`great_circle_distance_m` helpers for any new metrics so that dateline handling stays consistent throughout the codebase.

## Running Things
- For quick tests, prefer `conda activate plot` and run in that environment instead of creating a new virtual environment.
- Workflows are invoked interactively; the typical pattern for exporting META tracks looks like:
  ```shell
  python - <<'PY'
  from track import load_meta_data, export_meta_tracks
  ACL = load_meta_data()[1]
  export_meta_tracks(ACL, kind='acl', use_dask=True, write_contours=True)
  PY
  ```
- Swap in `process_argo_txt_to_yearly_parquet_dask()` or `plot_argo_hotspots(...)` inside the same REPL snippet when testing ingestion or plotting; keep each invocation focused to avoid reloading multi-GB datasets unnecessarily.
- Prefer lightweight validation first (single kind, narrow date/year range, optional `show_fig=False`) before launching full-region or multi-year runs.
- Heavy Dask workflows will open dashboards—monitor saturation before tuning worker counts or `chunk_size`.
- Plots and parquet outputs are written relative to the configured paths; clean up `_tmp` folders only after downstream consumers finish reading them.
