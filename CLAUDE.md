# CLAUDE.md

## Project Overview

Physical oceanography research project analyzing mesoscale eddy–Argo float dissolved oxygen interactions, overlaid on GLORYS reanalysis fields. Detects subduction/anomalies of dissolved oxygen associated with eddies across ocean basins.

## Environment

```bash
conda activate plot
```

No `requirements.txt` exists. Key dependencies: NumPy, Pandas, Dask, xarray, Matplotlib, Cartopy, SciPy, PyArrow, netCDF4.

## Architecture

- **`track.py`** — monolithic core (single very large file): geometry helpers, data ingestion (Argo/META/GLORYS), analytics (DO/AOU/TRIM anomaly detection), and all plotting logic.
- **Notebooks** (`GLORYS.ipynb`, `GLORYS AOU.ipynb`, `GLORYS TRIM.ipynb`) — consumers of `track.py`; do not put new logic in notebooks.
- **`config/`** — YAML-driven configuration:
  - `paths.yml` — data directory layout
  - `regions.yml` — spatial region definitions (lon/lat bounds, dateline flags)
  - `processing.yml` — algorithm parameters, thresholds, plot colors (versioned; see its header changelog)

## Key Conventions

1. **Config over constants** — Always read paths via `_PATHS_CFG` at the top of `track.py`. New features expose knobs through YAML, not hardcoded values.
2. **Author in track.py, consume in notebooks** — New/updated functions go in `track.py`; notebooks import and call them.
3. **Region globals** — Call `switch_region(...)` before workflows that depend on spatial filtering. In Dask/multiprocessing, reinitialize region state in workers.
4. **Dateline safety** — Never use naive `lon` comparisons. Use `_region_lon_mask`, `_minimal_lon_diff_deg`, or `adaptive_distance_m`. Normalize longitudes with `_normalize_lon_array`. Basemap styling comes from `processing.yml:plot.basemap_colors`; reuse `_BASEMAP_COLORS` when adding new plots.
5. **Kind strings** — Use `'acs'|'acl'|'cs'|'cl'` strings (not legacy list-based DS inputs) for `find_track`, `plot_track`, etc.
6. **Geospatial helpers** — Distances via `approximate_degree_length` or `great_circle_distance_m`. Eddy circles via `ellipse_patch_for_eddy`. Reuse existing helpers; do not reimplement Haversine. For new distance metrics use `local_xy_distance_m`/`great_circle_distance_m`; toggle precise calculation via `force_great_circle` rather than reimplementing.
7. **Post-edit self-review** — After each batch of edits to `track.py`, invoke the post-edit-review skill to check for comment hygiene, stale docstrings, mismatched return types, and dead parameters.

## Pipelines

- **META**: `load_meta_data()` → `export_meta_tracks(...)` → `find_track(kind, track_id)` to access fused parquet + zarr contours. Keep chunk sizes and Dask worker counts configurable when extending `export_meta_tracks`.
- **Argo**: Legacy `.mat` files via `convert_mat_to_parquet(year, input_dir, output_dir)`; newer `.txt` drops via `process_argo_txt_to_yearly_parquet_dask(...)`. `load_argo_data(year, variable_selection=...)` to read; extend variable handling via `default_selection` dict. Prefer yearly parquet for joins, not full in-memory loads.
- **GLORYS**: `get_track_area_glorys(...)` + `get_vertical_glorys` for data; plotting via `plot_*_horizontal_glorys` / `plot_*_vertical_glorys`. New 2D/3D vars need `alias_map` + `var_dims` updates; interpolation must use `RegularGridInterpolator` with masked arrays. `plot_data_package(..., save_fig=True)` writes method-independent replay plots to `plot_outputs/shared/<region>/plot_track_vertical_glorys/`.
- **Float/Eddy matching**: `filtered_float_data` fuses META tracks and Argo profiles via date-join → polygon containment (`is_point_in_contour`) → radius check (`adaptive_distance_m`). Reuse these helpers for any new proximity logic. `plot_track` orchestrates the full pipeline: load track, match floats, compute anomalies, plot.
- **Anomaly detection**: `calculate_delta_do` via `DetectionConfig` / `make_detection_config(...)`. Defaults from `processing.yml`.
- **Hotspot maps**: `plot_argo_hotspots(...)` — Dask-parallel per year, writes to `plot_outputs/<method>/<region>/plot_argo_hotspots/`.
- **Argo 3D reconstruction**: `collect_argo_pool(...)` → `_build_argo_3d_field(...)` (depth-parallel Gaussian-kernel field; zarr cache via `_save_argo_3d_field` / `load_argo_3d_field`) → `slice_section_from_argo_field(...)`. Plot 2×2 vertical overviews via `plot_regional_vertical_argo_overview` (Eulerian box), `plot_track_vertical_argo_overview` (Lagrangian / eddy-following), `plot_argo_vertical_argo_overview` (Eulerian point snapshot). Defaults in `processing.yml:argo_reconstruction`; `h_bw` / `min_weight` trade coverage vs mesoscale resolution.

## Running & Validation

- There are no tests. Validate by running the specific workflow function you changed.
- Prefer lightweight validation first (single kind, narrow date/year range, `show_fig=False`) before full-region runs.
- Typical invocation pattern:
  ```python
  from track import load_meta_data, export_meta_tracks
  ACL = load_meta_data()[1]
  export_meta_tracks(ACL, kind='acl', use_dask=True, write_contours=True)
  ```
- Monitor Dask dashboards for saturation before tuning worker counts.
- Clean up `_tmp` folders only after downstream consumers finish reading them.

## Output Layout

- Method-specific plots → `plot_outputs/<method>/<region>/...`
- Method-independent/shared outputs → `plot_outputs/shared/<region>/...`
- Filenames follow `{dataset}{id}_vertical_{var}_YYYYMMDD_k*b*.png`
