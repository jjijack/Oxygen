# AGENTS.md

## Sync Rule

- `AGENTS.md` and `CLAUDE.md` are peer entrypoints for this repo. Keep them equivalent.
- If either file changes, sync the counterpart in the same edit.
- If repo instructions conflict with the latest user request in the current chat, the user request wins.

## Project Overview

Physical oceanography research project analyzing mesoscale eddy-Argo float dissolved oxygen interactions, overlaid on GLORYS reanalysis fields. Detects subduction and anomaly patterns of dissolved oxygen associated with eddies across ocean basins. Frame this as a physical oceanography project; dissolved oxygen is primarily a tracer here, not the main biogeochemical endpoint.

## AI Collaboration Surfaces

- Main repo instructions live in this file and `CLAUDE.md`.
- Oxygen hot memory is maintained in Claude project memory at `/home/user3/.claude/projects/-mnt-w2-scratch-user3-Oxygen/memory/`.
- When working on Oxygen in Codex, read and update the relevant files in that Claude memory directory directly instead of maintaining a separate Codex project memory set.
- Shared skill docs are mirrored in `.claude/skills/` and `.agents/skills/`. When one side changes, sync the peer file in the same edit.
- `claude-memory` is a Codex-only local bridge skill under `.agents/skills/`; it exists to route Oxygen durable memory work into Claude project memory and intentionally has no `.claude` peer.

## User And Workflow Preferences

- Reply to the user in Chinese unless they explicitly ask otherwise.
- The user commits manually. One commit should carry one headline feature. You may draft commit messages, but do not run `git commit`.
- New datasets should go under `data/<dataset>/` and be registered in `config/paths.yml`. Leave legacy top-level data directories in place.
- Treat `/home/user3/scratch/Oxygen` and `/mnt/w2/scratch/user3/Oxygen` as the same repo path.

## Environment

```bash
conda activate plot
```

No `requirements.txt` exists. Key dependencies: NumPy, Pandas, Dask, xarray, Matplotlib, Cartopy, SciPy, PyArrow, netCDF4.

GLORYS NetCDF is not truly local. It is mounted from SJTU HPC through `sshfs`. If GLORYS paths start failing with `Input/output error` or repeated SSH reset symptoms, consult the Claude memory note `/home/user3/.claude/projects/-mnt-w2-scratch-user3-Oxygen/memory/glorys-data-mount-sshfs.md` before retrying mounts.

## Architecture

- **`track.py`**: monolithic core for geometry helpers, data ingestion, analytics, and plotting.
- **Notebooks**: `GLORYS.ipynb`, `GLORYS AOU.ipynb`, `GLORYS TRIM.ipynb`, and `OFES.ipynb` are consumers and workflow entrypoints; do not put new scientific logic in notebooks.
- **OFES surface eddies**: `ofes_surface_eddy.py` contains the PET detector and `run_ofes_surface_eddy.py` is its standalone producer wrapper. Run the producer in the dedicated `ofes-pet` environment; ordinary `track.py` imports must not require PET.
- **`config/`**:
  - `paths.yml`: data directory layout
  - `regions.yml`: spatial region definitions
  - `processing.yml`: algorithm parameters, thresholds, plot colors, and versioned changelog

## Key Conventions

1. **Config over constants**: read paths via `_PATHS_CFG` at the top of `track.py`. Expose new knobs through YAML instead of hardcoding them.
2. **Author in `track.py`, consume in notebooks**: new or changed logic belongs in `track.py`; notebooks import and call it.
3. **Region globals**: call `switch_region(...)` before workflows that depend on spatial filtering. In Dask or multiprocessing code, reinitialize region state inside workers.
4. **Dateline safety**: never rely on naive longitude comparisons. Use `_region_lon_mask`, `_minimal_lon_diff_deg`, `adaptive_distance_m`, and `_normalize_lon_array`.
5. **Kind strings**: use `'acs'|'acl'|'cs'|'cl'` for workflows such as `find_track` and `plot_track`.
6. **Geospatial helpers**: reuse `approximate_degree_length`, `great_circle_distance_m`, `local_xy_distance_m`, and `ellipse_patch_for_eddy`; do not reimplement distance logic.
7. **Docstring house style**: public entry points need the full summary -> prose -> `参数:` -> `返回:` -> `输出:` -> `说明:` structure documented in the Claude memory note `/home/user3/.claude/projects/-mnt-w2-scratch-user3-Oxygen/memory/docstring-convention.md`.
8. **Post-edit self-review**: after each batch of edits to `track.py`, run the post-edit review checklist from the mirrored skill docs.

## Pipelines

- **META**: `load_meta_data()` -> `export_meta_tracks(...)` -> `find_track(kind, track_id)` for parquet plus zarr contours.
- **Argo**: legacy `.mat` via `convert_mat_to_parquet(...)`; newer `.txt` via `process_argo_txt_to_yearly_parquet_dask(...)`; read with `load_argo_data(...)`.
- **GLORYS**: `get_track_area_glorys(...)` plus `get_vertical_glorys`; plotting via `plot_*_horizontal_glorys` and `plot_*_vertical_glorys`.
- **Float/Eddy matching**: `filtered_float_data` handles date join, contour containment, and radius filtering.
- **Anomaly detection**: `calculate_delta_do` configured through `DetectionConfig` and `processing.yml`.
- **Hotspot maps**: `plot_argo_hotspots(...)` writes to `plot_outputs/<method>/<region>/plot_argo_hotspots/`.
- **Argo 3D reconstruction**: `collect_argo_pool(...)` -> `_build_argo_3d_field(...)` -> slice and overview plotting helpers.
- **OFES**: expensive public producers create fixed semantic outputs; lightweight loaders, reducers, and plotters consume them. `OFES.ipynb` keeps producer cells visible but unexecuted and retains executed lightweight summaries and figures.

## Running And Validation

- There are no formal tests. Validate by running the specific workflow you changed.
- Prefer lightweight validation first: one kind, narrow year range, `show_fig=False` where possible.
- Typical invocation pattern:

```python
from track import load_meta_data, export_meta_tracks
ACL = load_meta_data()[1]
export_meta_tracks(ACL, kind='acl', use_dask=True, write_contours=True)
```

- Watch Dask dashboard saturation before tuning worker counts.
- Clean `_tmp` folders only after downstream consumers finish reading them.

## Output Layout

- Method-specific plots: `plot_outputs/<method>/<region>/...`
- Method-independent and shared outputs: `plot_outputs/shared/<region>/...`
- Filenames follow the current `{dataset}{id}_vertical_{var}_YYYYMMDD_k*b*.png` pattern.
