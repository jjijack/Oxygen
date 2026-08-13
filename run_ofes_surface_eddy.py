"""CLI orchestration for the standalone OFES PET surface-eddy workflow."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import os
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Mapping

import pandas as pd
import numpy as np
from netCDF4 import Dataset

from ofes_surface_eddy import (
    CONTOUR_COLUMNS,
    DAILY_COLUMNS,
    detect_ofes_eddies_day,
    date_range,
    empty_daily_frame,
    hash_file_inventory,
    load_yaml_config,
    parse_date,
    runtime_environment_manifest,
    sha256_bytes,
    sha256_file,
    track_surface_eddies,
    validate_basic_helpers,
    validate_tracking_helpers,
)

LOGGER = logging.getLogger("run_ofes_surface_eddy")
RUN_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    """Write JSON atomically so an interrupted run leaves the prior manifest."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def _relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _git_revision(repo_root: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _combined_code_hash(repo_root: Path) -> str:
    module_paths = [repo_root / "ofes_surface_eddy.py", repo_root / "run_ofes_surface_eddy.py"]
    payload = b"".join(
        path.name.encode() + b"\0" + path.read_bytes() + b"\0" for path in module_paths
    )
    return sha256_bytes(payload)


def _resolve_paths(
    repo_root: Path, config_path: Path
) -> tuple[dict[str, Any], Path, Path, Path]:
    config = load_yaml_config(config_path)
    eta_root = Path(config["input"]["eta_root"])
    if not eta_root.is_absolute():
        eta_root = repo_root / eta_root
    output_root = Path(config["output"]["root"])
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    protocol_path = repo_root / "ofes-surface-eddy-analysis-lock.md"
    return config, eta_root.resolve(), output_root.resolve(), protocol_path.resolve()


def _signature(
    repo_root: Path,
    config_path: Path,
    protocol_path: Path,
    eta_paths: list[Path],
    inventory_hash: str,
) -> dict[str, Any]:
    environment = runtime_environment_manifest()
    return {
        "code_sha256": _combined_code_hash(repo_root),
        "protocol_sha256": sha256_file(protocol_path),
        "config_sha256": sha256_file(config_path),
        "pet_module_sha256": environment["py_eddy_tracker_module_sha256"],
        "pet_version": environment["py_eddy_tracker_version"],
        "input_eta_inventory_sha256": inventory_hash,
        "date_window": {},
    }


def _make_scientific_signature(
    repo_root: Path,
    config_path: Path,
    protocol_path: Path,
    eta_paths: list[Path],
    inventory_hash: str,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    signature = _signature(
        repo_root, config_path, protocol_path, eta_paths, inventory_hash
    )
    signature["date_window"] = {
        "start": str(config["input"]["start_date"]),
        "end": str(config["input"]["end_date"]),
        "count": len(eta_paths),
    }
    signature["region"] = dict(config["input"]["region"])
    signature["filter"] = dict(config["filter"])
    signature["detection"] = dict(config["detection"])
    signature["tracking"] = dict(config["tracking"])
    return signature


def _manifest_path(run_dir: Path) -> Path:
    return run_dir / "manifest.json"


def _new_manifest(
    run_id: str,
    run_dir: Path,
    repo_root: Path,
    output_root: Path,
    signature: Mapping[str, Any],
    environment: Mapping[str, Any],
    inventory: list[dict[str, Any]],
    dates: list[dt.date],
    git_revision: str | None,
) -> dict[str, Any]:
    now = dt.datetime.now(dt.timezone.utc).isoformat()
    return {
        "schema_version": 1,
        "run_id": run_id,
        "status": "running",
        "created_at": now,
        "updated_at": now,
        "repo_root": str(repo_root),
        "run_dir": str(run_dir),
        "output_root": str(output_root),
        "git_revision_at_start": git_revision,
        "scientific_signature": dict(signature),
        "runtime_environment": dict(environment),
        "input_eta_inventory": inventory,
        "worker_count": None,
        "days": {
            day.isoformat(): {
                "status": "pending",
                "daily_path": str(run_dir / "daily" / f"{day:%Y%m%d}.parquet"),
                "native_paths": {
                    "anticyclonic": str(
                        run_dir / "native" / "anticyclonic" / f"{day:%Y%m%d}.nc"
                    ),
                    "cyclonic": str(
                        run_dir / "native" / "cyclonic" / f"{day:%Y%m%d}.nc"
                    ),
                },
            }
            for day in dates
        },
    }


def _validate_existing_manifest(
    manifest: Mapping[str, Any],
    signature: Mapping[str, Any],
    run_id: str,
    *,
    allow_code_change: bool = False,
) -> dict[str, str] | None:
    if manifest.get("run_id") != run_id:
        raise ValueError("run_id does not match existing manifest")
    existing_signature = dict(manifest.get("scientific_signature", {}))
    current_signature = dict(signature)
    if existing_signature == current_signature:
        return None
    changed_keys = {
        key
        for key in set(existing_signature) | set(current_signature)
        if existing_signature.get(key) != current_signature.get(key)
    }
    if not allow_code_change or changed_keys != {"code_sha256"}:
        raise ValueError(
            "Existing manifest scientific signature differs; only an explicit "
            "code_sha256-only resume is permitted"
        )
    return {
        "original_code_sha256": str(existing_signature["code_sha256"]),
        "resume_code_sha256": str(current_signature["code_sha256"]),
    }


def _detect_and_write_day(
    date_text: str, config_path_text: str, eta_root_text: str, run_dir_text: str
) -> dict[str, Any]:
    """Worker-safe daily detection and output writer."""

    config = load_yaml_config(config_path_text)
    frame, native, audit = detect_ofes_eddies_day(
        date_text, config, eta_root_text
    )
    run_dir = Path(run_dir_text)
    day = parse_date(date_text)
    daily_path = run_dir / "daily" / f"{day:%Y%m%d}.parquet"
    daily_path.parent.mkdir(parents=True, exist_ok=True)
    frame[list(DAILY_COLUMNS)].to_parquet(daily_path, index=False)
    native_paths: dict[str, str] = {}
    for polarity, observations in native.items():
        native_path = run_dir / "native" / polarity / f"{day:%Y%m%d}.nc"
        native_path.parent.mkdir(parents=True, exist_ok=True)
        observations.write_file(filename=str(native_path))
        native_paths[polarity] = str(native_path)
    return {
        "date": date_text,
        "audit": audit,
        "daily_path": str(daily_path),
        "native_paths": native_paths,
        "daily_sha256": sha256_file(daily_path),
        "daily_size": daily_path.stat().st_size,
    }


def _read_manifest(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"Manifest is not an object: {path}")
    return value


def _catalog_frames_from_manifest(manifest: Mapping[str, Any]) -> list[pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    for date_text, day_info in manifest["days"].items():
        if day_info.get("status") != "complete":
            raise RuntimeError(f"Day {date_text} is not complete")
        path = Path(day_info["daily_path"])
        if not path.is_file():
            raise FileNotFoundError(path)
        expected_hash = day_info.get("daily_sha256")
        if expected_hash and sha256_file(path) != expected_hash:
            raise RuntimeError(f"Daily catalog hash differs from manifest: {path}")
        frame = pd.read_parquet(path)
        if not frame.empty and set(frame["date"].astype(str)) != {date_text}:
            raise ValueError(f"Daily date mismatch in {path}")
        frames.append(frame)
    return frames


def _finite_sequence(value: Any) -> bool:
    """Return whether a parquet contour field contains only finite values."""

    try:
        array = np.asarray(value, dtype="f8")
    except (TypeError, ValueError):
        return False
    return bool(np.isfinite(array).all())


def _audit_complete_catalog(manifest: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    """Audit all daily parquet/native outputs before tracking is allowed."""

    expected_dates = [day.isoformat() for day in date_range(config["input"]["start_date"], config["input"]["end_date"])]
    day_map = manifest.get("days", {})
    if list(day_map) != expected_dates:
        raise RuntimeError("Catalog day inventory does not match the locked 365-day window")
    object_count = 0
    native_count = 0
    finite_contour_rows = 0
    for date_text in expected_dates:
        info = day_map[date_text]
        if info.get("status") != "complete":
            raise RuntimeError(f"Catalog day is not complete: {date_text}")
        daily_path = Path(info["daily_path"])
        frame = pd.read_parquet(daily_path)
        missing = sorted(set(DAILY_COLUMNS).difference(frame.columns))
        if missing:
            raise RuntimeError(f"Daily PET table lacks stable columns for {date_text}: {missing}")
        if not frame.empty:
            if set(frame["date"].astype(str)) != {date_text}:
                raise RuntimeError(f"Daily PET date mismatch: {daily_path}")
            if frame["object_id"].duplicated().any():
                raise RuntimeError(f"Duplicate daily object_id: {daily_path}")
            for column in CONTOUR_COLUMNS:
                if not frame[column].map(_finite_sequence).all():
                    raise RuntimeError(f"Non-finite parquet contour in {daily_path}: {column}")
            finite_contour_rows += len(frame)
            object_count += len(frame)
        for polarity, native_path_text in info.get("native_paths", {}).items():
            native_path = Path(native_path_text)
            if not native_path.is_file():
                raise FileNotFoundError(native_path)
            with Dataset(native_path, "r") as dataset:
                required = {
                    "effective_contour_latitude",
                    "effective_contour_longitude",
                    "speed_contour_latitude",
                    "speed_contour_longitude",
                }
                missing_native = sorted(required.difference(dataset.variables))
                if missing_native:
                    raise RuntimeError(f"Native PET file lacks contour variables: {native_path}: {missing_native}")
                for variable_name in required:
                    values = np.ma.asarray(dataset.variables[variable_name][:])
                    if not np.isfinite(values.compressed()).all():
                        raise RuntimeError(f"Non-finite native contour values: {native_path}: {variable_name}")
            native_count += 1
    return {
        "days_complete": len(expected_dates),
        "daily_objects": object_count,
        "native_files_read": native_count,
        "daily_rows_with_finite_contours": finite_contour_rows,
    }


def run_catalog(
    repo_root: str | Path,
    config_path: str | Path,
    run_id: str,
    workers: int = 1,
    resume: bool = True,
    allow_code_change_on_resume: bool = False,
) -> Path:
    """Run or resume the complete daily PET catalog for one explicit run id."""

    repo_root = Path(repo_root).resolve()
    config_path = Path(config_path).resolve()
    if not RUN_ID_RE.fullmatch(run_id):
        raise ValueError("run_id may contain only letters, numbers, '.', '_' and '-'")
    if not 1 <= workers <= 16:
        raise ValueError("workers must be in [1, 16]")
    config, eta_root, output_root, protocol_path = _resolve_paths(repo_root, config_path)
    dates = date_range(config["input"]["start_date"], config["input"]["end_date"])
    eta_paths = [eta_root / f"eta.{day:%m.%d.%Y}.nc" for day in dates]
    inventory_hash, inventory = hash_file_inventory(eta_paths)
    signature = _make_scientific_signature(
        repo_root, config_path, protocol_path, eta_paths, inventory_hash, config
    )
    run_dir = output_root / "runs" / run_id
    manifest_path = _manifest_path(run_dir)
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(
                f"Run already exists: {run_dir}; pass --resume or choose a new run_id"
            )
        manifest = _read_manifest(manifest_path)
        code_change = _validate_existing_manifest(
            manifest,
            signature,
            run_id,
            allow_code_change=allow_code_change_on_resume,
        )
        if code_change is not None:
            resume_audit = {
                "type": "code_hash_change_on_resume",
                "recorded_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "original_code_sha256": code_change["original_code_sha256"],
                "resume_code_sha256": code_change["resume_code_sha256"],
                "scientific_signature_non_code_fields_match": True,
                "explicit_cli_flag": "--allow-code-change-on-resume",
            }
            manifest.setdefault("resume_audits", []).append(resume_audit)
    else:
        run_dir.mkdir(parents=True, exist_ok=False)
        manifest = _new_manifest(
            run_id,
            run_dir,
            repo_root,
            output_root,
            signature,
            runtime_environment_manifest(),
            inventory,
            dates,
            _git_revision(repo_root),
        )
        _write_json(manifest_path, manifest)
    manifest["status"] = "running"
    manifest["worker_count"] = workers
    manifest["updated_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    _write_json(manifest_path, manifest)

    pending = [
        day
        for day in dates
        if not (
            manifest["days"][day.isoformat()].get("status") == "complete"
            and Path(manifest["days"][day.isoformat()]["daily_path"]).is_file()
        )
    ]
    try:
        if workers == 1:
            for day in pending:
                date_text = day.isoformat()
                manifest["days"][date_text]["status"] = "running"
                manifest["updated_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
                _write_json(manifest_path, manifest)
                result = _detect_and_write_day(
                    date_text, str(config_path), str(eta_root), str(run_dir)
                )
                manifest["days"][date_text].update(result)
                manifest["days"][date_text]["status"] = "complete"
                manifest["updated_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
                _write_json(manifest_path, manifest)
        elif pending:
            for day in pending:
                manifest["days"][day.isoformat()]["status"] = "running"
            _write_json(manifest_path, manifest)
            with ProcessPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(
                        _detect_and_write_day,
                        day.isoformat(),
                        str(config_path),
                        str(eta_root),
                        str(run_dir),
                    ): day
                    for day in pending
                }
                for future in as_completed(futures):
                    day = futures[future]
                    date_text = day.isoformat()
                    result = future.result()
                    manifest["days"][date_text].update(result)
                    manifest["days"][date_text]["status"] = "complete"
                    manifest["updated_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
                    _write_json(manifest_path, manifest)
    except Exception as exc:
        failed_date = next(
            (
                day_text
                for day_text, day_info in manifest["days"].items()
                if day_info.get("status") == "running"
            ),
            None,
        )
        if failed_date is not None:
            manifest["days"][failed_date]["status"] = "failed"
            manifest["days"][failed_date]["error"] = repr(exc)
        manifest["status"] = "failed"
        manifest["updated_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
        _write_json(manifest_path, manifest)
        raise

    catalog_audit = _audit_complete_catalog(manifest, config)
    frames = _catalog_frames_from_manifest(manifest)
    all_objects = (
        pd.concat(frames, ignore_index=True, sort=False)
        if frames
        else empty_daily_frame()
    )
    if not all_objects.empty:
        if all_objects["object_id"].duplicated().any():
            raise ValueError("Duplicate daily object_id in catalog")
        all_objects = all_objects[list(DAILY_COLUMNS)].sort_values(
            ["date", "polarity", "object_id"]
        )
    daily_output = output_root / "surface_eddy_daily_objects.parquet"
    daily_output.parent.mkdir(parents=True, exist_ok=True)
    all_objects.to_parquet(daily_output, index=False)
    manifest["status"] = "complete"
    manifest["catalog"] = {
        "days_expected": len(dates),
        "days_complete": sum(
            info.get("status") == "complete" for info in manifest["days"].values()
        ),
        "objects_total": int(len(all_objects)),
        "daily_output": str(daily_output),
        "audit": catalog_audit,
        "audit_code_sha256": _combined_code_hash(repo_root),
    }
    manifest["updated_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    _write_json(manifest_path, manifest)
    return manifest_path


def run_tracking(
    repo_root: str | Path, manifest_path: str | Path, config_path: str | Path
) -> dict[str, str]:
    """Track a complete catalog and write actual/virtual track tables."""

    repo_root = Path(repo_root).resolve()
    manifest_path = Path(manifest_path).resolve()
    config_path = Path(config_path).resolve()
    manifest = _read_manifest(manifest_path)
    if manifest.get("status") != "complete":
        raise RuntimeError("Tracking requires a catalog manifest with status=complete")
    config = load_yaml_config(config_path)
    catalog_audit = _audit_complete_catalog(manifest, config)
    frames = _catalog_frames_from_manifest(manifest)
    actual, observations, tracks = track_surface_eddies(frames, config)
    output_root = Path(manifest["output_root"]).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    paths = {
        "daily": output_root / "surface_eddy_daily_objects.parquet",
        "tracks": output_root / "surface_eddy_tracks.parquet",
        "track_observations": output_root / "surface_eddy_track_observations.parquet",
        "summary": output_root / "surface_eddy_summary.json",
    }
    actual[list(DAILY_COLUMNS)].to_parquet(paths["daily"], index=False)
    tracks.to_parquet(paths["tracks"], index=False)
    observations[list(DAILY_COLUMNS)].to_parquet(paths["track_observations"], index=False)
    summary = {
        "run_id": manifest["run_id"],
        "catalog_manifest": str(manifest_path),
        "objects_actual": int(len(actual)),
        "observations_actual_plus_virtual": int(len(observations)),
        "tracks": int(len(tracks)),
        "virtual_observations": int(observations["is_virtual"].astype(bool).sum()),
        "duration_class_counts": tracks["duration_class"].value_counts().to_dict()
        if not tracks.empty
        else {},
        "paths": {key: str(value) for key, value in paths.items()},
    }
    _write_json(paths["summary"], summary)
    manifest["tracking"] = {
        "status": "complete",
        "engine": "PET vertice_overlap effective-contour IoU",
        "minimum_overlap": config["tracking"]["minimum_overlap"],
        "maximum_search_days": config["tracking"]["maximum_search_days"],
        "maximum_consecutive_virtual": config["tracking"]["maximum_consecutive_virtual"],
        "paths": {key: str(value) for key, value in paths.items()},
        "objects_actual": int(len(actual)),
        "observations_actual_plus_virtual": int(len(observations)),
        "tracks": int(len(tracks)),
        "catalog_audit": catalog_audit,
        "code_sha256": _combined_code_hash(repo_root),
    }
    manifest["updated_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    _write_json(manifest_path, manifest)
    return {key: str(value) for key, value in paths.items()}


def run_validation(
    repo_root: str | Path,
    config_path: str | Path,
    output_dir: str | Path | None = None,
    workers: int = 1,
) -> dict[str, Any]:
    """Run fixed engineering dates and basic tests without event inputs."""

    repo_root = Path(repo_root).resolve()
    config_path = Path(config_path).resolve()
    config, eta_root, output_root, _ = _resolve_paths(repo_root, config_path)
    if not 1 <= workers <= 16:
        raise ValueError("workers must be in [1, 16]")
    validation_root = (
        Path(output_dir).resolve()
        if output_dir is not None
        else output_root / "validation"
    )
    validation_root.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "basic_helpers": validate_basic_helpers(),
        "tracking_helpers": validate_tracking_helpers(),
        "workers": workers,
        "dates": {},
    }
    validation_dates = [parse_date(value) for value in config["validation"]["fixed_dates"]]
    if workers == 1:
        day_results = [
            _detect_and_write_day(
                day.isoformat(), str(config_path), str(eta_root), str(validation_root)
            )
            for day in validation_dates
        ]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    _detect_and_write_day,
                    day.isoformat(),
                    str(config_path),
                    str(eta_root),
                    str(validation_root),
                ): day
                for day in validation_dates
            }
            day_results = [future.result() for future in as_completed(futures)]
    for day_result in sorted(day_results, key=lambda value: value["date"]):
        result["dates"][day_result["date"]] = day_result
    result_path = validation_root / "validation.json"
    _write_json(result_path, result)
    result["validation_path"] = str(result_path)
    return result


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-level", default="INFO")
    parser.add_argument("--log-file")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate")
    validate.add_argument("--repo-root", default=".")
    validate.add_argument("--config", default="config/ofes_surface_eddy.yml")
    validate.add_argument("--output-dir")
    validate.add_argument("--workers", type=int, default=1)

    catalog = subparsers.add_parser("catalog")
    catalog.add_argument("--repo-root", default=".")
    catalog.add_argument("--config", default="config/ofes_surface_eddy.yml")
    catalog.add_argument("--run-id", required=True)
    catalog.add_argument("--workers", type=int, default=1)
    catalog.add_argument("--no-resume", action="store_true")
    catalog.add_argument("--allow-code-change-on-resume", action="store_true")

    tracking = subparsers.add_parser("tracking")
    tracking.add_argument("--repo-root", default=".")
    tracking.add_argument("--config", default="config/ofes_surface_eddy.yml")
    tracking.add_argument("--manifest", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    log_handlers: list[logging.Handler] = [logging.StreamHandler()]
    if getattr(args, "log_file", None):
        log_path = Path(args.log_file).resolve()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_handlers.append(logging.FileHandler(log_path, encoding="utf-8"))
    logging.basicConfig(
        level=getattr(logging, str(args.log_level).upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        handlers=log_handlers,
    )
    if args.command == "validate":
        result = run_validation(args.repo_root, args.config, args.output_dir, args.workers)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    if args.command == "catalog":
        manifest = run_catalog(
            args.repo_root,
            args.config,
            args.run_id,
            workers=args.workers,
            resume=not args.no_resume,
            allow_code_change_on_resume=args.allow_code_change_on_resume,
        )
        print(manifest)
        return 0
    if args.command == "tracking":
        result = run_tracking(args.repo_root, args.manifest, args.config)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
