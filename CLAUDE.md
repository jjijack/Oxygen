# CLAUDE.md

## Sync Rule

- `CLAUDE.md` and `AGENTS.md` are peer entrypoints for this repo. Keep them equivalent.
- If either file changes, sync the counterpart in the same edit.
- If repo instructions conflict with the latest user request in the current chat, the user request wins.

## Project Overview

Physical oceanography research project analyzing mesoscale eddy-Argo float dissolved oxygen interactions, overlaid on GLORYS reanalysis fields. Detects subduction and anomaly patterns of dissolved oxygen associated with eddies across ocean basins. Frame this as a physical oceanography project; dissolved oxygen is primarily a tracer here, not the main biogeochemical endpoint.

## AI Collaboration Surfaces

- Main repo instructions live in this file and `AGENTS.md`.
- Oxygen hot memory is maintained in Claude project memory at `/home/user3/.claude/projects/-mnt-w2-scratch-user3-Oxygen/memory/`.
- When working on Oxygen in Codex, read and update the relevant files in that Claude memory directory directly instead of maintaining a separate Codex project memory set.
- Shared skill docs (`draft-commit`, `post-edit-review`) are mirrored in `.claude/skills/` and `.agents/skills/`. When one side changes, sync the peer file in the same edit.
- `claude-memory` (routes Oxygen memory work into the Claude project memory) and `oxygen-background-runner` (detached long runs) are Codex-only skills under `.agents/skills/`; they intentionally have no `.claude` peer.

## User And Workflow Preferences

- Reply to the user in Chinese unless they explicitly ask otherwise.
- Git:
  - On any branch other than `main`, commit validated work yourself without asking, in the `draft-commit` style: one headline feature per commit. Unpushed commits may be amended or fixed up.
  - Stage only your own changes; leave another agent's in-flight edits alone.
  - Ask the user first before committing to or merging into `main`, pushing any branch (including force-pushes), deleting branches, or discarding uncommitted or untracked work.
  - Explore on a feature branch and commit each piece once it is confirmed. Commit `GLORYS.ipynb` once at the end of a direction; `OFES.ipynb` cells travel with their `track.py` commit.
  - Run Git mutations from `/mnt/w2/scratch/user3/Oxygen` (the home alias can fail on `.git/index.lock`).
  - For GitHub over Clash TUN, configure `IPQoS none` in the local SSH host settings and use `git push origin <branch>` normally. If port 22 is unavailable, use `GIT_SSH_COMMAND='ssh -4 -p 443 -o Hostname=ssh.github.com -o HostKeyAlias=github.com -o IPQoS=none' git push origin <branch>`.
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
- **Argo Letter**: producers (`build_argo_do_occurrence_table`, `audit_scv_profile_identity`, the root-level `run_scv_matched_effects.py`, `calculate_scv_core_alignment`, `calculate_scv_matched_level_sensitivity`, `build_argo_lens_case`) write fixed outputs under `plot_outputs/`; the journal-format `plot_*` figure functions and `export_scv_matched_effect_tables` consume them; `assemble_argo_paper_package.py` copies figures and tables into the untracked `argo_paper_package/`. The `GLORYS.ipynb` section keeps producers unexecuted and figures executed.

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
- When reworking a producer or plotter, write to a scratch directory first and compare with the existing formal output by content (rows, columns and values; pixels for figures) before replacing it. Do not add hash or checksum checks.

## Output Layout

- Method-specific outputs: `plot_outputs/<method>/<region>/<function name>/` via `DetectionConfig.output_dir(fn_name)`; name directories after what the function does, not after the paper that uses it.
- Method-independent and shared outputs: `plot_outputs/shared/<region>/...`
- Exploration and one-off checks: `plot_outputs/test/`.
- Vertical-section figures follow the `{dataset}{id}_vertical_{var}_YYYYMMDD_k*b*.png` pattern.
- Data, caches and results stay out of Git, and so do the paper packages (`argo_paper_package/`, `ofes_paper_package/`). Retired code goes to `old codes/`, retired data and outputs to `old data/`; do not create per-package `archive/` directories.
