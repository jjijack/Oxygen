"""Run the OFES NP30 surface-eddy producer in the dedicated PET environment.

The command-line runner keeps the expensive PET dependency separate from the
plotting environment.  It writes one fixed catalog directory, reuses complete
daily fragments, and leaves reading, association, and plotting to ``track.py``.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import multiprocessing
import os
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import yaml

import ofes_surface_eddy as pet


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (dt.date, dt.datetime)):
        return value.isoformat()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Cannot encode value of type {type(value).__name__}")


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, default=_json_default),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def _resolve_path(repo_root: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else repo_root / path


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        value = yaml.safe_load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a YAML mapping in {path}")
    return value


def _catalog_paths(root: Path) -> dict[str, Path]:
    return {
        "manifest": root / "manifest.json",
        "valid_domain": root / "surface_eddy_valid_domain.npz",
        "daily_objects": root / "surface_eddy_daily_objects.parquet",
        "track_observations": root / "surface_eddy_track_observations.parquet",
        "tracks": root / "surface_eddy_tracks.parquet",
        "summary": root / "surface_eddy_catalog_summary.json",
    }


def _build_request(
    config: Mapping[str, Any], start: dt.date, end: dt.date
) -> dict[str, Any]:
    return {
        "date_start": start.isoformat(),
        "date_end": end.isoformat(),
        "region": dict(config["input"].get("region", {})),
        "raw_eta_units": str(config["input"].get("raw_eta_units", "cm")),
        "filter": dict(config.get("filter", {})),
        "detection": dict(config.get("detection", {})),
        "tracking": dict(config.get("tracking", {})),
    }


def _prepare_configuration(
    repo_root: Path,
    paths_config: Path,
    processing_config: Path,
    start_date: str | dt.date | None,
    end_date: str | dt.date | None,
    eta_root: str | Path | None,
    output_dir: str | Path | None,
) -> tuple[dict[str, Any], dt.date, dt.date, Path, dict[str, Any]]:
    paths = _load_yaml(paths_config).get("paths", {})
    processing = _load_yaml(processing_config)
    surface = dict(processing.get("ofes", {}).get("surface_eddy", {}) or {})
    start = pet.parse_date(start_date or surface.get("date_start", "2003-01-01"))
    end = pet.parse_date(end_date or surface.get("date_end", "2003-12-31"))
    if end < start:
        raise ValueError("end_date precedes start_date")

    configured_eta_root = (
        eta_root
        if eta_root is not None
        else _resolve_path(
            repo_root,
            Path(paths.get("ofes_root", "./data/OFES_NP30"))
            / str(surface.get("eta_subdir", "eta")),
        )
    )
    resolved_eta_root = _resolve_path(repo_root, configured_eta_root).resolve()
    configured_output = (
        output_dir
        if output_dir is not None
        else paths.get(
            "ofes_surface_eddy_root",
            "./plot_outputs/do/ofes_np30_ke/surface_eddy",
        )
    )
    resolved_output = _resolve_path(repo_root, configured_output).resolve()
    config = {
        "input": {
            "eta_root": str(resolved_eta_root),
            "start_date": start.isoformat(),
            "end_date": end.isoformat(),
            "region": dict(surface.get("region", {})),
            "raw_eta_units": str(surface.get("raw_eta_units", "cm")),
        },
        "filter": dict(surface.get("filter", {})),
        "detection": dict(surface.get("detection", {})),
        "tracking": dict(surface.get("tracking", {})),
    }
    return (
        config,
        start,
        end,
        resolved_output,
        _build_request(config, start, end),
    )


def _write_valid_domain(
    config: Mapping[str, Any], start: dt.date, path: Path
) -> None:
    source = pet.eta_path_for_date(config["input"]["eta_root"], start)
    snapshot = pet.load_eta_snapshot(
        source,
        expected_date=start,
        raw_units=str(config["input"]["raw_eta_units"]),
    )
    prepared = pet.prepare_pet_grid(snapshot, config)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp.npz")
    np.savez_compressed(
        temporary,
        lon=np.asarray(snapshot.lon, dtype=float),
        lat=np.asarray(snapshot.lat, dtype=float),
        filter_valid=np.asarray(prepared.filter_valid_mask.T, dtype=bool),
        ocean_valid=np.asarray(prepared.ocean_valid_mask.T, dtype=bool),
    )
    os.replace(temporary, path)


def _read_fragment(path: Path, date_text: str) -> pd.DataFrame | None:
    if not path.is_file():
        return None
    try:
        frame = pd.read_parquet(path)
    except Exception:
        return None
    missing = sorted(set(pet.DAILY_COLUMNS).difference(frame.columns))
    if missing:
        return None
    if not frame.empty and set(frame["date"].astype(str)) != {date_text}:
        return None
    return frame.reindex(columns=list(pet.DAILY_COLUMNS))


def _detect_day(task: tuple[str, Mapping[str, Any]]) -> tuple[str, pd.DataFrame, dict[str, Any]]:
    date_text, config = task
    frame, _native, audit = pet.detect_ofes_eddies_day(
        date_text, config, config["input"]["eta_root"]
    )
    return date_text, frame.reindex(columns=list(pet.DAILY_COLUMNS)), audit


def _catalog_is_complete(
    root: Path,
    manifest: Mapping[str, Any],
    dates: list[dt.date],
) -> bool:
    if manifest.get("status") != "complete":
        return False
    paths = _catalog_paths(root)
    required = ("valid_domain", "daily_objects", "track_observations", "tracks")
    if not all(paths[key].is_file() for key in required):
        return False
    try:
        for key in ("daily_objects", "track_observations"):
            table = pd.read_parquet(paths[key])
            if not set(pet.DAILY_COLUMNS).issubset(table.columns):
                return False
        track_table = pd.read_parquet(paths["tracks"], columns=["track_id"])
        if "track_id" not in track_table.columns:
            return False
    except Exception:
        return False
    day_map = manifest.get("days", {})
    expected = [day.isoformat() for day in dates]
    if list(day_map) != expected:
        return False
    return all(
        day_map[date_text].get("status") == "complete"
        and _read_fragment(
            root / "daily" / f"{pd.Timestamp(date_text):%Y%m%d}.parquet",
            date_text,
        )
        is not None
        for date_text in expected
    )


def _result(root: Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    paths = _catalog_paths(root)
    return {
        "run_dir": root,
        "manifest": dict(manifest),
        "paths": {key: value for key, value in paths.items()},
    }


def run_catalog(
    repo_root: str | Path = ".",
    *,
    paths_config: str | Path = "config/paths.yml",
    processing_config: str | Path = "config/processing.yml",
    start_date: str | dt.date | None = None,
    end_date: str | dt.date | None = None,
    eta_root: str | Path | None = None,
    output_dir: str | Path | None = None,
    workers: int = 1,
) -> dict[str, Any]:
    """Produce the fixed OFES PET surface-eddy catalog from the command line.

    The producer reads daily OFES SSH snapshots, detects surface eddies, and
    links observations with the configured PET overlap rule.  Complete daily
    fragments are reused automatically so an interrupted production can be
    called again with the same scientific request.

    参数:
        - repo_root (str | Path): repository root used to resolve relative paths。
        - paths_config (str | Path): paths YAML containing OFES and output roots。
        - processing_config (str | Path): processing YAML containing PET settings。
        - start_date (str | date | None): inclusive first day; None uses the YAML value。
        - end_date (str | date | None): inclusive last day; None uses the YAML value。
        - eta_root (str | Path | None): optional explicit ETA directory。
        - output_dir (str | Path | None): optional explicit fixed catalog directory。
        - workers (int): number of daily PET worker processes。

    返回:
        - dict: fixed output directory, manifest, and stable product paths。

    输出:
        - `manifest.json`, daily fragments, valid-domain masks, and tracked
          parquet tables below the fixed output directory。

    说明:
        - 该 runner 不导入 `track.py`；Notebook 的 reducer、loader 和绘图仍在
          plotting 环境中运行。
        - 输出目录中的旧结果只按结构和请求日期复用；科学配置或核心算法改变时，
          由用户手动删除固定目录后重跑。
    """
    if int(workers) < 1:
        raise ValueError("workers must be >= 1")
    root = Path(repo_root).expanduser().resolve()
    paths_path = _resolve_path(root, paths_config).resolve()
    processing_path = _resolve_path(root, processing_config).resolve()
    config, start, end, output_root, request = _prepare_configuration(
        root,
        paths_path,
        processing_path,
        start_date,
        end_date,
        eta_root,
        output_dir,
    )
    dates = pet.date_range(start, end)
    paths = _catalog_paths(output_root)
    manifest = _read_json(paths["manifest"]) if paths["manifest"].is_file() else None
    if manifest is not None:
        existing_request = manifest.get("request")
        if existing_request is not None and existing_request != request:
            raise ValueError(
                "Existing surface-eddy output has a different scientific request; "
                "pass a new output_dir or remove the fixed directory explicitly."
            )
        if _catalog_is_complete(output_root, manifest, dates):
            return _result(output_root, manifest)

    output_root.mkdir(parents=True, exist_ok=True)
    day_root = output_root / "daily"
    day_root.mkdir(parents=True, exist_ok=True)
    if not paths["valid_domain"].is_file():
        _write_valid_domain(config, start, paths["valid_domain"])

    if manifest is None:
        manifest = {
            "schema_version": 1,
            "status": "running",
            "request": request,
            "worker_count": int(workers),
            "days": {},
        }
    manifest["status"] = "running"
    manifest["request"] = request
    manifest["worker_count"] = int(workers)
    old_days = manifest.get("days", {})
    manifest["days"] = {
        date_text: dict(old_days.get(date_text, {}))
        for date_text in (day.isoformat() for day in dates)
    }
    _write_json(paths["manifest"], manifest)

    frames_by_date: dict[str, pd.DataFrame] = {}
    pending: list[str] = []
    for day in dates:
        date_text = day.isoformat()
        fragment = day_root / f"{day:%Y%m%d}.parquet"
        candidate = _read_fragment(fragment, date_text)
        if manifest["days"][date_text].get("status") == "complete" and candidate is not None:
            frames_by_date[date_text] = candidate
        else:
            pending.append(date_text)
            manifest["days"][date_text] = {
                "status": "pending",
                "daily_path": str(fragment),
            }
    _write_json(paths["manifest"], manifest)

    def store_detected(
        date_text: str, frame: pd.DataFrame, audit: dict[str, Any]
    ) -> None:
        fragment = day_root / f"{pd.Timestamp(date_text):%Y%m%d}.parquet"
        frame.reindex(columns=list(pet.DAILY_COLUMNS)).to_parquet(
            fragment, index=False
        )
        frames_by_date[date_text] = frame
        manifest["days"][date_text] = {
            "status": "complete",
            "daily_path": str(fragment),
            "audit": audit,
        }
        _write_json(paths["manifest"], manifest)

    try:
        tasks = [(date_text, config) for date_text in pending]
        if int(workers) > 1 and len(tasks) > 1:
            context_name = (
                "fork"
                if "fork" in multiprocessing.get_all_start_methods()
                else "spawn"
            )
            with multiprocessing.get_context(context_name).Pool(
                min(int(workers), len(tasks))
            ) as pool:
                for date_text, frame, audit in pool.imap_unordered(_detect_day, tasks):
                    store_detected(date_text, frame, audit)
        else:
            for task in tasks:
                store_detected(*_detect_day(task))

        frames = [frames_by_date[day.isoformat()] for day in dates]
        actual, observations, tracks = pet.track_surface_eddies(frames, config)
        daily_columns = list(pet.DAILY_COLUMNS)
        actual[daily_columns].to_parquet(paths["daily_objects"], index=False)
        observations[daily_columns].to_parquet(
            paths["track_observations"], index=False
        )
        tracks.to_parquet(paths["tracks"], index=False)
        summary = {
            "request": request,
            "days": len(dates),
            "objects_actual": len(actual),
            "observations_actual_plus_virtual": len(observations),
            "tracks": len(tracks),
            "virtual_observations": int(observations["is_virtual"].astype(bool).sum())
            if "is_virtual" in observations
            else 0,
            "paths": {key: value for key, value in paths.items()},
        }
        _write_json(paths["summary"], summary)
        manifest.update(
            {
                "status": "complete",
                "catalog": {
                    "days_expected": len(dates),
                    "days_complete": len(dates),
                    "objects_total": len(actual),
                    "daily_output": str(paths["daily_objects"]),
                },
                "tracking": {
                    "status": "complete",
                    "paths": {
                        "daily": str(paths["daily_objects"]),
                        "track_observations": str(paths["track_observations"]),
                        "tracks": str(paths["tracks"]),
                    },
                },
            }
        )
        _write_json(paths["manifest"], manifest)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _write_json(paths["manifest"], manifest)
        raise
    return _result(output_root, manifest)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--paths-config", default="config/paths.yml")
    parser.add_argument("--processing-config", default="config/processing.yml")
    parser.add_argument("--start-date")
    parser.add_argument("--end-date")
    parser.add_argument("--eta-root")
    parser.add_argument("--output-dir")
    parser.add_argument("--workers", type=int, default=1)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    result = run_catalog(
        args.repo_root,
        paths_config=args.paths_config,
        processing_config=args.processing_config,
        start_date=args.start_date,
        end_date=args.end_date,
        eta_root=args.eta_root,
        output_dir=args.output_dir,
        workers=args.workers,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
