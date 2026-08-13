"""Post-process the locked OFES PET surface-eddy catalog.

The detector itself deliberately knows nothing about DO50 events.  This
module is the later, explicit-input gate: it consumes one complete catalog
manifest plus exact event tables, records their hashes, and writes the
event-level association/null/cross-tab products into the dedicated surface
eddy output tree.
"""

from __future__ import annotations

import argparse
import ast
import datetime as dt
import json
import logging
import math
import os
import re
import tempfile
from collections import Counter
from pathlib import Path as FilePath
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml

# Set plotting/runtime caches before importing matplotlib or PET.  The
# dedicated environment may be run with a read-only home directory.
_RUNTIME_CACHE = FilePath(tempfile.gettempdir()) / "ofes-pet-runtime-cache"
_RUNTIME_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_RUNTIME_CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_RUNTIME_CACHE / "xdg"))

from matplotlib.path import Path as PolygonPath
from scipy.stats import wilcoxon

from ofes_surface_eddy import (
    DAILY_COLUMNS,
    POLARITY_NAME,
    date_range,
    load_eta_snapshot,
    load_yaml_config,
    parse_date,
    prepare_pet_grid,
    sha256_bytes,
    sha256_file,
)

LOGGER = logging.getLogger("run_ofes_surface_eddy_postprocess")


EVENT_ALIASES: dict[str, tuple[str, ...]] = {
    "event_id": (
        "event_id",
        "eventid",
        "event",
        "event_name",
        "id",
        "name",
    ),
    "peak_date": (
        "peak_date",
        "event_peak_date",
        "date_peak",
        "peak_day",
        "date",
    ),
    "peak_lon": (
        "peak_lon",
        "peak_longitude",
        "event_lon",
        "center_lon",
        "lon",
        "longitude",
    ),
    "peak_lat": (
        "peak_lat",
        "peak_latitude",
        "event_lat",
        "center_lat",
        "lat",
        "latitude",
    ),
    "peak_core_lon": (
        "peak_core_lon",
        "peak_core_longitude",
        "core_lon",
        "core_longitude",
        "peak_lon",
        "event_lon",
        "center_lon",
        "lon",
    ),
    "peak_core_lat": (
        "peak_core_lat",
        "peak_core_latitude",
        "core_lat",
        "core_latitude",
        "peak_lat",
        "event_lat",
        "center_lat",
        "lat",
    ),
    "core_footprint_lon": (
        "core_footprint_lon",
        "peak_core_footprint_lon",
        "do_core_lon",
        "do_peak_lon",
        "core_points_lon",
    ),
    "core_footprint_lat": (
        "core_footprint_lat",
        "peak_core_footprint_lat",
        "do_core_lat",
        "do_peak_lat",
        "core_points_lat",
    ),
    "do_peak_footprint_lon": (
        "do_peak_footprint_lon",
        "peak_footprint_lon",
        "do_footprint_lon",
        "footprint_lon",
        "footprint_longitude",
    ),
    "do_peak_footprint_lat": (
        "do_peak_footprint_lat",
        "peak_footprint_lat",
        "do_footprint_lat",
        "footprint_lat",
        "footprint_latitude",
    ),
    "strict_eligible": (
        "strict_eligible",
        "strict_event",
        "strict_56",
        "main_56",
        "population_diagnostic_passed",
        "background_ring_complete",
        "ring_complete",
    ),
    "quality_eligible": (
        "quality_eligible",
        "quality_eligible_161",
        "quality_flag",
        "is_quality_eligible",
    ),
    "background_ring_complete": (
        "background_ring_complete",
        "ring_complete",
        "complete_ring",
    ),
    "rotation_dominated": (
        "rotation_dominated",
        "rotation_event",
        "is_rotation_dominated",
        "rotation_dominant",
        "deep_rotation_dominated",
    ),
    "deep_ro_sign": (
        "deep_ro_sign",
        "deep_core_ro_sign",
        "core_weighted_ro_sign",
        "rossby_number",
        "subsurface_rossby_number",
        "absolute_subsurface_rossby_number",
        "deep_ro_polarity",
        "deep_core_polarity",
        "ro_polarity",
    ),
    "surface_ro_expression": (
        "surface_ro_expression",
        "surface_ro_class",
        "surface_ro_status",
        "surface_rotation_expression",
        "surface_ro",
    ),
    "surface_ro_same_polarity": (
        "surface_ro_same_polarity",
        "same_polarity_surface_ro",
        "surface_ro_same_sign",
        "core_weighted_surface_ro_same_polarity",
        "surface_core_rotation_polarity_match",
    ),
    "surface_ro_polarity": (
        "surface_ro_polarity",
        "core_weighted_surface_ro_polarity",
        "surface_ro_sign",
        "surface_rotation_polarity",
        "surface_core_weighted_rossby_number",
    ),
    "mccoy_center": (
        "mccoy_center",
        "mccoy_center_match",
        "mccoy_center_compatible",
        "center_mccoy",
        "mccoy_method_center",
        "center_profile_mccoy_compatible",
    ),
    "mccoy_center_velocity": (
        "mccoy_center_velocity",
        "mccoy_center_plus_velocity",
        "center_plus_velocity",
        "mccoy_center_and_velocity",
        "center_profile_velocity_confirmed",
    ),
    "mccoy_any17": (
        "mccoy_any17",
        "mccoy_any_17",
        "any17_mccoy",
        "mccoy_any",
        "any_mccoy",
        "any_event_profile_mccoy_compatible",
    ),
    "mccoy_any17_velocity": (
        "mccoy_any17_velocity",
        "mccoy_any17_plus_velocity",
        "any17_plus_velocity",
        "mccoy_any_and_velocity",
        "any_event_profile_velocity_confirmed",
    ),
}

TABLE_SUFFIXES = (".parquet", ".pq", ".csv", ".json")
PRIMARY_PET_CATEGORIES = ("same", "opposite", "no-PET", "ambiguous")


def _load_postprocess_config(path: str | FilePath) -> tuple[dict[str, Any], FilePath]:
    """Load detection config plus the separate event-analysis controls."""

    config_path = FilePath(path).resolve()
    config = load_yaml_config(config_path)
    if "association" in config:
        return config, config_path
    postprocess_path = config_path.with_name("ofes_surface_eddy_postprocess.yml")
    if not postprocess_path.is_file():
        raise FileNotFoundError(
            f"Association config is absent from {config_path} and sibling "
            f"postprocess config is missing: {postprocess_path}"
        )
    with postprocess_path.open(encoding="utf-8") as stream:
        postprocess = yaml.safe_load(stream)
    if not isinstance(postprocess, dict) or not isinstance(postprocess.get("association"), dict):
        raise ValueError(f"Invalid association config: {postprocess_path}")
    config["association"] = postprocess["association"]
    return config, postprocess_path


def _require_association_lock(config_path: FilePath) -> tuple[FilePath, str]:
    """Require the pre-registered event-analysis lock before event reads."""

    lock_path = config_path.parent.parent / "ofes-surface-eddy-association-lock.md"
    if not lock_path.is_file():
        raise FileNotFoundError(f"Association lock is missing: {lock_path}")
    lock_text = lock_path.read_text(encoding="utf-8")
    if "Status: pre-registered" not in lock_text:
        raise RuntimeError(f"Association lock does not declare pre-registration: {lock_path}")
    return lock_path, sha256_file(lock_path)


def _write_json(path: FilePath, value: Mapping[str, Any]) -> None:
    """Write JSON atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    temporary.replace(path)


def _norm_column(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, (list, tuple, np.ndarray, dict)):
        return False
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _first_non_missing(values: Iterable[Any]) -> Any:
    for value in values:
        if not _is_missing(value):
            return value
    return np.nan


def _parse_literal(value: Any) -> Any:
    """Parse list-like parquet/csv fields without executing code."""

    if _is_missing(value):
        return None
    if isinstance(value, (list, tuple, np.ndarray)):
        return value
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        for loader in (json.loads, ast.literal_eval):
            try:
                return loader(text)
            except (ValueError, SyntaxError, json.JSONDecodeError, TypeError):
                pass
    return value


def _as_bool(value: Any) -> bool | float:
    if _is_missing(value):
        return np.nan
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        if np.isfinite(float(value)):
            return bool(float(value))
        return np.nan
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "t", "present", "same", "positive"}:
        return True
    if text in {"0", "false", "no", "n", "f", "absent", "none", "negative"}:
        return False
    return np.nan


def _sign(value: Any) -> int | float:
    if _is_missing(value):
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        number = float(value)
        if not np.isfinite(number) or number == 0:
            return np.nan
        return 1 if number > 0 else -1
    text = str(value).strip().lower()
    if text in {"1", "+1", "positive", "pos", "anticyclonic", "anticyclone", "a", "anti"}:
        return 1
    if text in {"-1", "negative", "neg", "cyclonic", "cyclone", "c", "cycl"}:
        return -1
    if "anti" in text or "positive" in text:
        return 1
    if "cycl" in text or "negative" in text:
        return -1
    return np.nan


def _read_table(path: FilePath) -> pd.DataFrame:
    """Read one explicitly selected input table."""

    suffix = path.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        frame = pd.read_parquet(path)
    elif suffix == ".csv":
        frame = pd.read_csv(path)
    elif suffix == ".json":
        frame = pd.read_json(path)
    else:
        raise ValueError(f"Unsupported event table format: {path}")
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise ValueError(f"Event table is empty or not tabular: {path}")
    return frame


def _resolve_input_file(value: str | FilePath, label: str) -> FilePath:
    """Resolve a file or an unambiguous input directory.

    Directory resolution is deliberately conservative: a directory with more
    than one candidate requires the caller to provide the exact table path.
    This prevents an arbitrary newest/lexicographically-last event table from
    silently changing the scientific denominator.
    """

    path = FilePath(value).expanduser().resolve()
    if path.is_file():
        return path
    if not path.is_dir():
        raise FileNotFoundError(f"{label} input does not exist: {path}")
    candidates = sorted(
        candidate
        for candidate in path.rglob("*")
        if candidate.is_file() and candidate.suffix.lower() in TABLE_SUFFIXES
    )
    if len(candidates) == 1:
        return candidates[0]
    preferred_names = {
        "population": (
            "population_peak_diagnostics.parquet",
            "population.parquet",
            "event_population.parquet",
            "events.parquet",
        ),
        "diagnostics": (
            "event_diagnostic_summary.parquet",
            "event_summary.parquet",
            "deep_sensitivity_ranking.parquet",
            "events.parquet",
            "event_diagnostics.parquet",
            "diagnostics.parquet",
        ),
        "catalog": (
            "quality_event_catalog.parquet",
            "event_catalog.parquet",
            "events.parquet",
            "catalog.parquet",
        ),
    }.get(label, ())
    # Directory inputs are allowed, but resolution is priority-based rather
    # than lexicographic: a diagnostics directory contains several valid
    # tables, while only one is the requested primary table.
    for preferred_name in preferred_names:
        preferred = [candidate for candidate in candidates if candidate.name == preferred_name]
        if len(preferred) == 1:
            return preferred[0]
        if len(preferred) > 1:
            raise ValueError(
                f"{label} input has multiple {preferred_name} files; pass the exact path"
            )
    listing = "\n".join(f"  - {candidate}" for candidate in candidates[:40])
    raise ValueError(
        f"{label} input directory is ambiguous; pass the exact table path.\n{listing}"
    )


def _alias_value(frame: pd.DataFrame, canonical: str) -> pd.Series | None:
    """Find one canonical field using normalized aliases."""

    normalized = {_norm_column(column): column for column in frame.columns}
    aliases = EVENT_ALIASES.get(canonical, (canonical,))
    for alias in aliases:
        original = normalized.get(_norm_column(alias))
        if original is not None:
            return frame[original]
    return None


def _canonicalize_table(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    """Add stable event fields while retaining source columns for audit."""

    result = frame.copy()
    for canonical in EVENT_ALIASES:
        if canonical in result.columns:
            continue
        series = _alias_value(result, canonical)
        if series is not None:
            result[canonical] = series
    if "event_id" not in result:
        raise ValueError(f"{source} table has no recognizable event-id column")
    result["event_id"] = result["event_id"].astype(str).str.strip()
    if result["event_id"].duplicated().any():
        duplicates = result.loc[result["event_id"].duplicated(), "event_id"].head(10).tolist()
        raise ValueError(f"{source} table has duplicate event ids: {duplicates}")
    if "peak_date" in result:
        result["peak_date"] = pd.to_datetime(result["peak_date"], errors="raise").dt.date.astype(str)
    for column in (
        "peak_lon",
        "peak_lat",
        "peak_core_lon",
        "peak_core_lat",
    ):
        if column in result:
            result[column] = pd.to_numeric(result[column], errors="coerce")
    for column in (
        "strict_eligible",
        "quality_eligible",
        "background_ring_complete",
        "rotation_dominated",
        "mccoy_center",
        "mccoy_center_velocity",
        "mccoy_any17",
        "mccoy_any17_velocity",
        "surface_ro_same_polarity",
    ):
        if column in result:
            result[column] = result[column].map(_as_bool)
    for column in ("deep_ro_sign", "surface_ro_polarity"):
        if column in result:
            result[column] = result[column].map(_sign)
    if "rotation_dominated" not in result:
        for column in result.columns:
            normalized = _norm_column(column)
            if "rotation" in normalized and "dominat" in normalized:
                result["rotation_dominated"] = result[column].map(_as_bool)
                break
    if "peak_core_lon" not in result and "peak_lon" in result:
        result["peak_core_lon"] = result["peak_lon"]
    if "peak_core_lat" not in result and "peak_lat" in result:
        result["peak_core_lat"] = result["peak_lat"]
    result["_source_table"] = source
    return result


def _coalesce_sources(
    population: pd.DataFrame,
    diagnostics: pd.DataFrame | None,
    catalog: pd.DataFrame | None,
) -> pd.DataFrame:
    """Left-join explicit diagnostics/catalog fields onto population rows."""

    base = _canonicalize_table(population, "population")
    for source_name, source in (("diagnostics", diagnostics), ("catalog", catalog)):
        if source is None:
            continue
        right = _canonicalize_table(source, source_name)
        merged = base.merge(
            right,
            on="event_id",
            how="left",
            suffixes=("", f"_{source_name}"),
            validate="one_to_one",
        )
        for column in right.columns:
            if column == "event_id":
                continue
            suffixed = f"{column}_{source_name}"
            if suffixed in merged.columns:
                if column in merged.columns:
                    merged[column] = merged[column].combine_first(merged[suffixed])
                    merged.drop(columns=[suffixed], inplace=True)
                else:
                    merged.rename(columns={suffixed: column}, inplace=True)
        base = merged
    return base


def _select_strict_events(
    frame: pd.DataFrame, expected_count: int
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Select only an explicitly flagged strict set, never by row order."""

    candidate_count = len(frame)
    selected = frame
    selection_method = "population_rows"
    if "strict_eligible" in frame:
        flags = frame["strict_eligible"]
        if flags.notna().any():
            selected = frame.loc[flags.fillna(False).astype(bool)].copy()
            selection_method = "strict_eligible_flag"
    if len(selected) != expected_count and "background_ring_complete" in frame:
        flags = frame["background_ring_complete"]
        if flags.notna().any():
            candidate = frame.loc[flags.fillna(False).astype(bool)].copy()
            if len(candidate) == expected_count:
                selected = candidate
                selection_method = "background_ring_complete_flag"
    if len(selected) != expected_count:
        raise ValueError(
            f"Strict event set has {len(selected)} rows, expected {expected_count}; "
            "an explicit strict/ring-complete flag is required when the source "
            "contains a different candidate count."
        )
    if selected["event_id"].duplicated().any():
        raise ValueError("Strict event set contains duplicate event ids")
    required = {"peak_date", "peak_core_lon", "peak_core_lat"}
    missing = sorted(required.difference(selected.columns))
    if missing:
        raise ValueError(f"Strict event table lacks required fields: {missing}")
    if selected[list(required)].isna().any().any():
        raise ValueError("Strict event table contains missing peak date/core coordinates")
    return selected.reset_index(drop=True), {
        "candidate_count": candidate_count,
        "strict_count": len(selected),
        "selection_method": selection_method,
    }


def _validate_identity_consistency(
    population: pd.DataFrame,
    *sources: pd.DataFrame | None,
) -> None:
    """Check event identity fields against the authoritative population table."""

    authority = _canonicalize_table(population, "population-authority")
    fields = ("peak_date", "peak_lon", "peak_lat", "peak_core_lon", "peak_core_lat")
    for source_index, source in enumerate(sources, start=1):
        if source is None:
            continue
        candidate = _canonicalize_table(source, f"source-{source_index}")
        missing_ids = sorted(
            set(authority["event_id"].astype(str))
            - set(candidate["event_id"].astype(str))
        )
        if missing_ids:
            raise ValueError(
                "Source is missing authoritative event IDs: "
                f"{missing_ids[:10]}"
            )
        common = authority.merge(
            candidate,
            on="event_id",
            how="inner",
            suffixes=("_authority", "_source"),
            validate="one_to_one",
        )
        for field in fields:
            left, right = f"{field}_authority", f"{field}_source"
            if left not in common or right not in common:
                continue
            if field == "peak_date":
                mismatch = common[left].notna() & common[right].notna() & common[left].ne(common[right])
            else:
                left_values = pd.to_numeric(common[left], errors="coerce")
                right_values = pd.to_numeric(common[right], errors="coerce")
                mismatch = left_values.notna() & right_values.notna() & ~np.isclose(
                    left_values, right_values, rtol=0.0, atol=1e-8
                )
            if bool(mismatch.any()):
                bad_ids = common.loc[mismatch, "event_id"].head(5).tolist()
                raise ValueError(
                    f"Authoritative population identity mismatch in {field}: {bad_ids}"
                )


def _hash_explicit_path(path: str | FilePath) -> tuple[str, list[dict[str, Any]]]:
    """Hash an explicit file or all parquet files below an explicit directory."""

    root = FilePath(path).expanduser().resolve()
    if root.is_file():
        files = [root]
    elif root.is_dir():
        files = sorted(
            candidate.resolve()
            for candidate in root.rglob("*.parquet")
            if candidate.is_file()
        )
        if not files:
            raise FileNotFoundError(f"No parquet file below explicit directory: {root}")
    else:
        raise FileNotFoundError(root)
    inventory = [
        {"path": str(file), "sha256": sha256_file(file), "size": file.stat().st_size}
        for file in files
    ]
    return sha256_bytes(json.dumps(inventory, sort_keys=True).encode()), inventory


def _find_unique_named(root: str | FilePath, filename: str) -> FilePath:
    """Find exactly one named parquet below an explicit path."""

    path = FilePath(root).expanduser().resolve()
    if path.is_file():
        if path.name != filename:
            raise ValueError(f"Expected {filename}, received {path}")
        return path
    if not path.is_dir():
        raise FileNotFoundError(path)
    candidates = sorted(candidate.resolve() for candidate in path.rglob(filename) if candidate.is_file())
    if not candidates:
        raise FileNotFoundError(f"{filename} not found below {path}")
    if len(candidates) != 1:
        raise ValueError(f"Ambiguous {filename} below {path}; pass the exact file path: {candidates}")
    return candidates[0]


def _require_complete_pet_manifest(
    manifest_path: str | FilePath,
    repo_root: str | FilePath | None = None,
    config_path: str | FilePath | None = None,
) -> tuple[FilePath, dict[str, Any]]:
    """Refuse event reads until the complete, current 365-day PET run exists."""

    path = FilePath(manifest_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete":
        raise RuntimeError("PET event association requires manifest status=complete")
    date_window = manifest.get("scientific_signature", {}).get("date_window", {})
    if int(date_window.get("count", 0)) != 365:
        raise RuntimeError("PET event association requires the locked 365-day catalog")
    if int(manifest.get("catalog", {}).get("days_complete", 0)) != 365:
        raise RuntimeError("PET manifest does not report all 365 days complete")
    days = manifest.get("days", {})
    if len(days) != 365 or any(info.get("status") != "complete" for info in days.values()):
        raise RuntimeError("PET manifest has a missing or incomplete daily status")
    tracking = manifest.get("tracking", {})
    if tracking.get("status") != "complete":
        raise RuntimeError("PET event association requires complete tracking output")
    for key in ("daily", "track_observations", "tracks"):
        output = tracking.get("paths", {}).get(key)
        if not output or not FilePath(output).is_file():
            raise FileNotFoundError(f"Missing tracked PET output {key}: {output}")
    if repo_root is not None:
        protocol_path = FilePath(repo_root).expanduser().resolve() / "ofes-surface-eddy-analysis-lock.md"
        expected_protocol = sha256_file(protocol_path)
        observed_protocol = manifest.get("scientific_signature", {}).get("protocol_sha256")
        if observed_protocol != expected_protocol:
            raise RuntimeError("PET manifest protocol hash differs from the current analysis lock")
    if config_path is not None:
        expected_config = sha256_file(FilePath(config_path).expanduser().resolve())
        observed_config = manifest.get("scientific_signature", {}).get("config_sha256")
        if observed_config != expected_config:
            raise RuntimeError("PET manifest config hash differs from the current detector config")
    return path, manifest


def _read_population_authoritative(
    population_path: str | FilePath, config: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, FilePath]:
    """Read the exact population diagnostics and derive only its flagged 56."""

    path = _find_unique_named(population_path, "population_peak_diagnostics.parquet")
    population = pd.read_parquet(path).copy()
    required = {
        "event_id",
        "threshold",
        "population_diagnostic_passed",
        "peak_date",
        "daily_object_key",
        "peak_depth_m",
        "peak_lon",
        "peak_lat",
    }
    missing = sorted(required - set(population.columns))
    if missing:
        raise KeyError(f"Authoritative population is missing columns: {missing}")
    population["event_id"] = population["event_id"].astype(str)
    population["threshold"] = pd.to_numeric(population["threshold"], errors="raise")
    population["population_diagnostic_passed"] = population[
        "population_diagnostic_passed"
    ].map(_as_bool)
    population["peak_date"] = pd.to_datetime(population["peak_date"], errors="raise").dt.normalize()
    association_cfg = config["association"]
    candidates = population.loc[np.isclose(population["threshold"], 50.0)].copy()
    expected_candidates = int(association_cfg["expected_candidate_events"])
    if len(candidates) != expected_candidates:
        raise ValueError(f"DO50 candidate count changed: expected {expected_candidates}, found {len(candidates)}")
    strict = candidates.loc[candidates["population_diagnostic_passed"].eq(True)].copy()
    strict = strict.sort_values(["peak_date", "event_id"], kind="mergesort").reset_index(drop=True)
    expected_strict = int(association_cfg["expected_strict_events"])
    if len(strict) != expected_strict:
        raise ValueError(f"Strict DO50 count changed: expected {expected_strict}, found {len(strict)}")
    if strict["event_id"].duplicated().any():
        raise ValueError("Strict population contains duplicate event_id values")
    return population, strict, path


def _read_do_catalog(
    manifest: Mapping[str, Any], delta_catalog_path: str | FilePath
) -> tuple[pd.DataFrame, FilePath]:
    """Read the explicit DO catalog daily objects and return its root."""

    root = FilePath(delta_catalog_path).expanduser().resolve()
    daily_path = _find_unique_named(root, "daily_objects.parquet")
    daily = pd.read_parquet(daily_path).copy()
    required = {"event_id", "date", "peak_lon", "peak_lat"}
    missing = sorted(required - set(daily.columns))
    if missing:
        raise KeyError(f"DO daily catalog is missing columns: {missing}")
    daily["event_id"] = daily["event_id"].astype(str)
    daily["date"] = pd.to_datetime(daily["date"], errors="raise").dt.normalize()
    for column in ("peak_lon", "peak_lat"):
        daily[column] = pd.to_numeric(daily[column], errors="raise")
    start = pd.Timestamp(manifest["scientific_signature"]["date_window"]["start"])
    end = pd.Timestamp(manifest["scientific_signature"]["date_window"]["end"])
    if daily["date"].min() < start or daily["date"].max() > end:
        raise ValueError("DO daily catalog extends outside the locked 2003 window")
    return daily, root


def _validate_peak_identity(
    strict: pd.DataFrame, do_daily: pd.DataFrame
) -> None:
    """Verify event id/date/lon/lat at the authoritative DO peak day."""

    expected = strict[["event_id", "peak_date", "peak_lon", "peak_lat"]].copy()
    expected["event_id"] = expected["event_id"].astype(str)
    expected["peak_date"] = pd.to_datetime(expected["peak_date"], errors="raise").dt.normalize()
    expected["peak_lon"] = pd.to_numeric(expected["peak_lon"], errors="raise")
    expected["peak_lat"] = pd.to_numeric(expected["peak_lat"], errors="raise")
    peak_rows = do_daily.merge(
        expected[["event_id", "peak_date"]],
        left_on=["event_id", "date"],
        right_on=["event_id", "peak_date"],
        how="inner",
    ).drop(columns=["peak_date"])
    if peak_rows["event_id"].duplicated().any():
        duplicates = peak_rows.loc[peak_rows["event_id"].duplicated(), "event_id"].head(10).tolist()
        raise ValueError(
            "DO catalog has multiple peak-day rows for strict events; "
            f"cannot choose one silently: {duplicates}"
        )
    missing = sorted(set(expected["event_id"]) - set(peak_rows["event_id"]))
    if missing:
        raise ValueError(f"DO catalog lacks peak-day rows for strict events: {missing[:10]}")
    observed = expected.merge(peak_rows[["event_id", "peak_lon", "peak_lat"]], on="event_id", suffixes=("_population", "_catalog"), validate="one_to_one")
    date_mismatch = observed["peak_date"].isna()
    coordinate_mismatch = (
        ~np.isclose(observed["peak_lon_population"], observed["peak_lon_catalog"], rtol=0.0, atol=1e-8)
        | ~np.isclose(observed["peak_lat_population"], observed["peak_lat_catalog"], rtol=0.0, atol=1e-8)
    )
    if bool(date_mismatch.any() or coordinate_mismatch.any()):
        bad = observed.loc[date_mismatch | coordinate_mismatch, "event_id"].head(10).tolist()
        raise ValueError(f"DO catalog peak identity differs from population authority: {bad}")


def _read_surface_outputs(manifest: Mapping[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    paths = manifest["tracking"]["paths"]
    daily = pd.read_parquet(paths["daily"]).copy()
    observations = pd.read_parquet(paths["track_observations"]).copy()
    tracks = pd.read_parquet(paths["tracks"]).copy()
    for frame in (daily, observations):
        frame["date"] = pd.to_datetime(frame["date"], errors="raise").dt.normalize()
        for column in ("is_virtual", "filter_valid", "boundary_censored", "coastal"):
            frame[column] = frame[column].map(_as_bool)
    return daily, observations, tracks


def _cell_area_m2(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Calculate spherical cell areas from ascending center coordinates."""

    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat_edges = np.empty(lat.size + 1, dtype=float)
    lon_edges = np.empty(lon.size + 1, dtype=float)
    lat_edges[1:-1] = (lat[:-1] + lat[1:]) / 2.0
    lon_edges[1:-1] = (lon[:-1] + lon[1:]) / 2.0
    lat_edges[0] = lat[0] - (lat[1] - lat[0]) / 2.0
    lat_edges[-1] = lat[-1] + (lat[-1] - lat[-2]) / 2.0
    lon_edges[0] = lon[0] - (lon[1] - lon[0]) / 2.0
    lon_edges[-1] = lon[-1] + (lon[-1] - lon[-2]) / 2.0
    lat_edges = np.clip(lat_edges, -90.0, 90.0)
    return (6_371_000.0**2) * np.diff(np.sin(np.deg2rad(lat_edges)))[:, None] * np.deg2rad(np.diff(lon_edges))[None, :]


def _load_surface_domain(manifest: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    """Reconstruct filter-valid/ocean-valid masks from the first ETA inventory file."""

    inventory = manifest.get("input_eta_inventory", [])
    if not inventory:
        raise ValueError("PET manifest has an empty ETA inventory")
    snapshot = load_eta_snapshot(
        inventory[0]["path"], raw_units=str(config["input"]["raw_eta_units"])
    )
    prepared = prepare_pet_grid(snapshot, config)
    lon_grid, lat_grid = np.meshgrid(snapshot.lon, snapshot.lat)
    return {
        "lon": snapshot.lon,
        "lat": snapshot.lat,
        "lon_grid": lon_grid,
        "lat_grid": lat_grid,
        "area_m2": _cell_area_m2(snapshot.lat, snapshot.lon),
        "filter_valid": np.asarray(prepared.filter_valid_mask.T, dtype=bool),
        "ocean_valid": np.asarray(prepared.ocean_valid_mask.T, dtype=bool),
    }


def _haversine_km(lon_a: Any, lat_a: Any, lon_b: Any, lat_b: Any) -> np.ndarray:
    """Great-circle distance in km."""

    lon_a, lat_a, lon_b, lat_b = [
        np.deg2rad(np.asarray(value, dtype=float))
        for value in (lon_a, lat_a, lon_b, lat_b)
    ]
    dlon = lon_b - lon_a
    dlat = lat_b - lat_a
    hav = np.sin(dlat / 2.0) ** 2 + np.cos(lat_a) * np.cos(lat_b) * np.sin(dlon / 2.0) ** 2
    return 2.0 * 6_371_000.0 * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0))) / 1000.0


def _as_float_array(value: Any) -> np.ndarray:
    parsed = _parse_literal(value)
    if parsed is None:
        return np.array([], dtype=float)
    array = np.asarray(parsed, dtype=float).ravel()
    return array[np.isfinite(array)]


def _polygon_from_row(row: Mapping[str, Any], contour: str = "effective") -> PolygonPath | None:
    lon = _as_float_array(row.get(f"{contour}_contour_lon"))
    lat = _as_float_array(row.get(f"{contour}_contour_lat"))
    if lon.size < 3 or lon.size != lat.size:
        return None
    vertices = np.column_stack([lon, lat])
    if not np.array_equal(vertices[0], vertices[-1]):
        vertices = np.vstack([vertices, vertices[0]])
    return PolygonPath(vertices, closed=True)


def _covered_by_rows(
    lon: np.ndarray, lat: np.ndarray, rows: pd.DataFrame, contour: str = "effective"
) -> np.ndarray:
    points = np.column_stack([np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)])
    covered = np.zeros(points.shape[0], dtype=bool)
    for row in rows.to_dict("records"):
        polygon = _polygon_from_row(row, contour)
        if polygon is not None:
            covered |= polygon.contains_points(points, radius=1e-10)
    return covered


def _valid_surface_rows(objects: pd.DataFrame, date: pd.Timestamp, virtual: bool | None = False) -> pd.DataFrame:
    rows = objects.loc[objects["date"].eq(date)].copy()
    if virtual is not None:
        rows = rows.loc[rows["is_virtual"].eq(virtual)]
    rows = rows.loc[rows["filter_valid"].eq(True) & ~rows["boundary_censored"].eq(True)]
    return rows.reset_index(drop=True)


def _peak_pixels(catalog_root: FilePath, population_row: Mapping[str, Any]) -> pd.DataFrame:
    date = pd.Timestamp(population_row["peak_date"]).normalize()
    path = catalog_root / "days" / f"peak_pixels_{date:%Y%m%d}.parquet"
    if not path.is_file():
        raise FileNotFoundError(path)
    pixels = pd.read_parquet(path)
    required = {"object_id_do50", "lat_index", "lon_index"}
    missing = sorted(required - set(pixels.columns))
    if missing:
        raise KeyError(f"Peak-pixel table missing columns: {missing}")
    key = str(population_row["daily_object_key"])
    try:
        object_id = int(key.rsplit("_", 1)[-1])
    except ValueError as error:
        raise ValueError(f"Cannot parse daily_object_key {key}") from error
    selected = pixels.loc[pixels["object_id_do50"].astype(int).eq(object_id)].copy()
    if selected.empty:
        raise ValueError(f"No peak pixels for {population_row['event_id']} / {key}")
    return selected


def _nearest_index(values: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(value))))


def _point_matches(rows: pd.DataFrame, lon: float, lat: float, contour: str = "effective") -> pd.DataFrame:
    point = np.array([float(lon), float(lat)], dtype=float)[None, :]
    matches = []
    for row in rows.to_dict("records"):
        polygon = _polygon_from_row(row, contour)
        if polygon is not None and bool(polygon.contains_points(point, radius=1e-10)[0]):
            matches.append(row)
    return pd.DataFrame(matches, columns=rows.columns)


def _event_row(
    population_row: Mapping[str, Any],
    objects: pd.DataFrame,
    domain: Mapping[str, Any],
    catalog_root: FilePath,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    event_id = str(population_row["event_id"])
    date = pd.Timestamp(population_row["peak_date"]).normalize()
    peak_lon, peak_lat = float(population_row["peak_lon"]), float(population_row["peak_lat"])
    lon_index, lat_index = _nearest_index(domain["lon"], peak_lon), _nearest_index(domain["lat"], peak_lat)
    core_filter_valid = bool(domain["filter_valid"][lat_index, lon_index])
    core_ocean_valid = bool(domain["ocean_valid"][lat_index, lon_index])
    eligible = core_filter_valid and core_ocean_valid
    rows = _valid_surface_rows(objects, date, virtual=False)
    effective = _point_matches(rows, peak_lon, peak_lat, "effective") if eligible else pd.DataFrame(columns=rows.columns)
    speed = _point_matches(rows, peak_lon, peak_lat, "speed") if eligible else pd.DataFrame(columns=rows.columns)
    pixels = _peak_pixels(catalog_root, population_row)
    pixel_lat = pixels["lat_index"].to_numpy(dtype=int)
    pixel_lon = pixels["lon_index"].to_numpy(dtype=int)
    inside = (
        (pixel_lat >= 0) & (pixel_lat < domain["lat"].size)
        & (pixel_lon >= 0) & (pixel_lon < domain["lon"].size)
    )
    pixel_lat, pixel_lon = pixel_lat[inside], pixel_lon[inside]
    pixel_filter_valid = domain["filter_valid"][pixel_lat, pixel_lon]
    pixel_ocean_valid = domain["ocean_valid"][pixel_lat, pixel_lon]
    pixel_valid = pixel_filter_valid & pixel_ocean_valid
    pixel_area = domain["area_m2"][pixel_lat, pixel_lon]
    pixel_covered = _covered_by_rows(domain["lon"][pixel_lon], domain["lat"][pixel_lat], rows, "effective")
    valid_peak_area = float(np.sum(pixel_area[pixel_valid]))
    peak_overlap = float(np.sum(pixel_area[pixel_valid & pixel_covered]) / valid_peak_area) if valid_peak_area else np.nan
    distance = _haversine_km(domain["lon_grid"], domain["lat_grid"], peak_lon, peak_lat)
    ring = (
        (distance >= float(config["association"]["ring_inner_km"]))
        & (distance <= float(config["association"]["ring_outer_km"]))
        & domain["filter_valid"] & domain["ocean_valid"]
    )
    for row_index, column_index in set(zip(pixel_lat.tolist(), pixel_lon.tolist())):
        ring[row_index, column_index] = False
    ring_points_covered = _covered_by_rows(domain["lon_grid"][ring], domain["lat_grid"][ring], rows, "effective")
    ring_area = float(np.sum(domain["area_m2"][ring]))
    ring_occupied_area = float(np.sum(domain["area_m2"][ring][ring_points_covered])) if ring_area else np.nan
    ring_fraction = ring_occupied_area / ring_area if ring_area else np.nan
    if rows.empty:
        nearest_center = nearest_vertex = distance_over_radius = np.nan
    else:
        center_distance = _haversine_km(rows["center_lon"].to_numpy(float), rows["center_lat"].to_numpy(float), peak_lon, peak_lat)
        nearest = int(np.nanargmin(center_distance))
        vertices = []
        for row in rows.to_dict("records"):
            vertices.append(np.nanmin(_haversine_km(_as_float_array(row["effective_contour_lon"]), _as_float_array(row["effective_contour_lat"]), peak_lon, peak_lat)))
        nearest_center = float(center_distance[nearest])
        nearest_vertex = float(np.nanmin(vertices))
        radius = float(rows.iloc[nearest]["effective_radius_km"])
        distance_over_radius = nearest_center / radius if radius > 0 else np.nan
    polarity_codes = sorted(set(effective["polarity_code"].astype(int).tolist())) if not effective.empty else []
    return {
        "event_id": event_id,
        "peak_date": date,
        "peak_lon": peak_lon,
        "peak_lat": peak_lat,
        "peak_depth_m": float(population_row["peak_depth_m"]),
        "peak_core_filter_valid": core_filter_valid,
        "peak_core_ocean_valid": core_ocean_valid,
        "peak_core_analysis_eligible": eligible,
        "peak_core_contained_by_actual_pet_effective_contour": bool(not effective.empty) if eligible else np.nan,
        "peak_core_contained_by_actual_pet_speed_contour": bool(not speed.empty) if eligible else np.nan,
        "effective_match_count": int(len(effective)),
        "effective_match_polarity_codes": json.dumps(polarity_codes),
        "effective_match_object_ids": json.dumps(sorted(effective["object_id"].astype(str).tolist())),
        "speed_match_count": int(len(speed)),
        "effective_peak_footprint_overlap_fraction": peak_overlap,
        "ring_occupancy_fraction": ring_fraction,
        "ring_valid_cell_count": int(np.count_nonzero(ring)),
        "ring_valid_area_m2": ring_area,
        "ring_occupied_area_m2": ring_occupied_area,
        "nearest_pet_center_distance_km": nearest_center,
        "nearest_pet_effective_contour_vertex_distance_km": nearest_vertex,
        "nearest_pet_center_distance_over_effective_radius": distance_over_radius,
    }


def _annual_occupancy(
    objects: pd.DataFrame,
    domain: Mapping[str, Any],
    dates: Sequence[pd.Timestamp],
) -> pd.DataFrame:
    """Compute the declared month x 1-degree latitude PET occupancy null."""

    valid = domain["filter_valid"] & domain["ocean_valid"]
    lat_bins = range(int(np.floor(domain["lat"].min())), int(np.floor(domain["lat"].max())) + 1)
    rows: list[dict[str, Any]] = []
    for value in dates:
        date = pd.Timestamp(value).normalize()
        actual = _valid_surface_rows(objects, date, virtual=False)
        covered = _covered_by_rows(
            domain["lon_grid"].ravel(), domain["lat_grid"].ravel(), actual, "effective"
        ).reshape(valid.shape)
        for lat_bin in lat_bins:
            latitude_mask = (domain["lat"] >= lat_bin) & (domain["lat"] < lat_bin + 1.0)
            mask = valid & latitude_mask[:, None]
            rows.append(
                {
                    "date": date,
                    "month": int(date.month),
                    "lat_bin_deg": int(lat_bin),
                    "denominator_cell_count": int(np.count_nonzero(mask)),
                    "occupied_cell_count": int(np.count_nonzero(mask & covered)),
                    "denominator_area_m2": float(np.sum(domain["area_m2"][mask])),
                    "occupied_area_m2": float(np.sum(domain["area_m2"][mask & covered])),
                }
            )
    daily = pd.DataFrame(rows)
    if daily.empty:
        return daily
    result = daily.groupby(["month", "lat_bin_deg"], as_index=False)[
        ["denominator_cell_count", "occupied_cell_count", "denominator_area_m2", "occupied_area_m2"]
    ].sum()
    result["occupancy_fraction"] = result["occupied_area_m2"] / result["denominator_area_m2"].replace(0.0, np.nan)
    return result


def _lifecycle_association(
    strict: pd.DataFrame,
    do_daily: pd.DataFrame,
    observations: pd.DataFrame,
) -> pd.DataFrame:
    """Count actual and virtual PET containment over each DO observed lifecycle."""

    records = []
    for population_row in strict.to_dict("records"):
        event_id = str(population_row["event_id"])
        days = do_daily.loc[do_daily["event_id"].eq(event_id)].sort_values("date")
        actual_effective = actual_speed = virtual_effective = virtual_speed = 0
        for day_row in days.to_dict("records"):
            date = pd.Timestamp(day_row["date"]).normalize()
            lon, lat = float(day_row["peak_lon"]), float(day_row["peak_lat"])
            actual = _valid_surface_rows(observations, date, virtual=False)
            virtual = _valid_surface_rows(observations, date, virtual=True)
            actual_effective += int(not _point_matches(actual, lon, lat, "effective").empty)
            actual_speed += int(not _point_matches(actual, lon, lat, "speed").empty)
            virtual_effective += int(not _point_matches(virtual, lon, lat, "effective").empty)
            virtual_speed += int(not _point_matches(virtual, lon, lat, "speed").empty)
        records.append(
            {
                "event_id": event_id,
                "observed_lifecycle_days": int(len(days)),
                "actual_pet_effective_contained_days": actual_effective,
                "actual_pet_speed_contained_days": actual_speed,
                "virtual_pet_effective_contained_days": virtual_effective,
                "virtual_pet_speed_contained_days": virtual_speed,
            }
        )
    return pd.DataFrame(records)


def _pet_peak_category(row: Mapping[str, Any], deep_code: int | float) -> str:
    """Classify a rotation-event PET peak as same/opposite/no/ambiguous."""

    if not bool(row.get("peak_core_analysis_eligible", False)):
        return "no-PET"
    contained = row.get("peak_core_contained_by_actual_pet_effective_contour")
    if _is_missing(contained) or not bool(contained):
        return "no-PET"
    if _is_missing(deep_code) or not np.isfinite(float(deep_code)):
        return "ambiguous"
    try:
        codes = json.loads(row.get("effective_match_polarity_codes", "[]"))
    except (TypeError, ValueError, json.JSONDecodeError):
        codes = []
    same = any(int(code) == int(deep_code) for code in codes)
    opposite = any(int(code) == -int(deep_code) for code in codes)
    if same and opposite:
        return "ambiguous"
    if same:
        return "same"
    if opposite:
        return "opposite"
    return "ambiguous"


def _rotation_crosstab(
    strict: pd.DataFrame, association: pd.DataFrame, config: Mapping[str, Any]
) -> pd.DataFrame:
    required = {
        "rotation_dominated",
        "rossby_number",
        "surface_core_rotation_polarity_match",
        "surface_core_weighted_rossby_number",
    }
    missing = sorted(required - set(strict.columns))
    if missing:
        raise KeyError(f"Population lacks locked rotation fields: {missing}")
    rotation = strict.loc[strict["rotation_dominated"].map(_as_bool).eq(True)].copy()
    expected = int(config["association"]["expected_rotation_events"])
    if len(rotation) != expected:
        raise ValueError(f"Expected {expected} rotation events, found {len(rotation)}")
    joined = rotation.merge(association, on="event_id", how="left", validate="one_to_one")
    deep_ro = pd.to_numeric(joined["rossby_number"], errors="coerce")
    # The source convention is Ro < 0 = anticyclonic, while PET uses +1 for
    # anticyclonic and -1 for cyclonic.
    joined["deep_polarity_code"] = np.where(deep_ro < 0, 1, np.where(deep_ro > 0, -1, np.nan))
    joined["surface_core_rotation_polarity_match"] = joined[
        "surface_core_rotation_polarity_match"
    ].map(_as_bool)
    joined["pet_peak_category"] = [
        _pet_peak_category(row, deep_code)
        for row, deep_code in zip(joined.to_dict("records"), joined["deep_polarity_code"])
    ]
    return joined[
        [
            "event_id",
            "peak_date",
            "peak_lon",
            "peak_lat",
            "deep_polarity_code",
            "surface_core_rotation_polarity_match",
            "surface_core_weighted_rossby_number",
            "pet_peak_category",
            "peak_core_contained_by_actual_pet_effective_contour",
            "effective_match_polarity_codes",
        ]
    ].sort_values("event_id", kind="mergesort").reset_index(drop=True)


def _mccoy_summary_path(
    diagnostics_path: str | FilePath, exact_path: str | FilePath | None = None
) -> FilePath:
    if exact_path is not None:
        path = FilePath(exact_path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    root = FilePath(diagnostics_path).expanduser().resolve()
    candidates = sorted(path.resolve() for path in root.rglob("event_summary.parquet") if path.is_file())
    mccoy = [path for path in candidates if any("mccoy" in part.lower() for part in path.parts)]
    if len(mccoy) == 1:
        return mccoy[0]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"No event_summary.parquet below {root}")
    raise ValueError(f"McCoy event summary is ambiguous; pass --mccoy-summary: {candidates}")


def _mccoy_crosstab(
    strict: pd.DataFrame,
    association: pd.DataFrame,
    diagnostics_path: str | FilePath,
    exact_path: str | FilePath | None = None,
) -> tuple[pd.DataFrame, FilePath]:
    path = _mccoy_summary_path(diagnostics_path, exact_path)
    mccoy = pd.read_parquet(path).copy()
    required = {
        "event_id",
        "center_profile_mccoy_compatible",
        "any_event_profile_mccoy_compatible",
        "center_profile_velocity_confirmed",
        "any_event_profile_velocity_confirmed",
    }
    missing = sorted(required - set(mccoy.columns))
    if missing:
        raise KeyError(f"McCoy event summary lacks columns: {missing}")
    mccoy["event_id"] = mccoy["event_id"].astype(str)
    if mccoy["event_id"].duplicated().any():
        raise ValueError("McCoy event summary contains duplicate event IDs")
    _validate_identity_consistency(strict, mccoy)
    joined = strict.merge(
        mccoy[list(required)], on="event_id", how="left", validate="one_to_one"
    ).merge(
        association[
            ["event_id", "pet_peak_category", "peak_core_contained_by_actual_pet_effective_contour"]
        ],
        on="event_id",
        how="left",
        validate="one_to_one",
    )
    flags = [
        "center_profile_mccoy_compatible",
        "any_event_profile_mccoy_compatible",
        "center_profile_velocity_confirmed",
        "any_event_profile_velocity_confirmed",
    ]
    for column in flags:
        joined[column] = joined[column].map(_as_bool)
        if joined[column].map(_is_missing).any():
            raise ValueError(f"McCoy event summary did not match all strict IDs for {column}")
    expected = {
        "center_profile_mccoy_compatible": 9,
        "center_profile_velocity_confirmed": 6,
        "any_event_profile_mccoy_compatible": 19,
        "any_event_profile_velocity_confirmed": 11,
    }
    for column, count in expected.items():
        observed = int(joined[column].sum())
        if observed != count:
            raise ValueError(f"McCoy count changed for {column}: expected {count}, found {observed}")
    return joined[
        [
            "event_id",
            "pet_peak_category",
            "peak_core_contained_by_actual_pet_effective_contour",
            *flags,
        ]
    ].sort_values("event_id", kind="mergesort").reset_index(drop=True), path


def _read_sensitivity_population(
    diagnostics_path: str | FilePath, config: Mapping[str, Any]
) -> pd.DataFrame:
    path = _find_unique_named(diagnostics_path, "deep_sensitivity_ranking.parquet")
    ranking = pd.read_parquet(path).copy()
    required = {"event_id", "threshold", "deep_sensitivity_eligible", "peak_date", "peak_lon", "peak_lat"}
    missing = sorted(required - set(ranking.columns))
    if missing:
        raise KeyError(f"Sensitivity ranking lacks columns: {missing}")
    selected = ranking.loc[
        np.isclose(pd.to_numeric(ranking["threshold"], errors="raise"), 50.0)
        & ranking["deep_sensitivity_eligible"].map(_as_bool).eq(True)
    ].copy()
    expected = int(config["association"]["expected_quality_eligible_events"])
    if len(selected) != expected:
        raise ValueError(f"Expected {expected} sensitivity events, found {len(selected)}")
    selected["event_id"] = selected["event_id"].astype(str)
    selected["peak_date"] = pd.to_datetime(selected["peak_date"], errors="raise").dt.normalize()
    return selected


def _sensitivity_association(ranking: pd.DataFrame, objects: pd.DataFrame) -> pd.DataFrame:
    records = []
    for row in ranking.to_dict("records"):
        date = pd.Timestamp(row["peak_date"]).normalize()
        actual = _valid_surface_rows(objects, date, virtual=False)
        effective = _point_matches(actual, float(row["peak_lon"]), float(row["peak_lat"]), "effective")
        speed = _point_matches(actual, float(row["peak_lon"]), float(row["peak_lat"]), "speed")
        records.append(
            {
                "event_id": str(row["event_id"]),
                "peak_date": date,
                "peak_lon": float(row["peak_lon"]),
                "peak_lat": float(row["peak_lat"]),
                "pet_effective_contained": bool(not effective.empty),
                "pet_speed_contained": bool(not speed.empty),
                "effective_match_count": int(len(effective)),
            }
        )
    return pd.DataFrame(records).sort_values("event_id", kind="mergesort").reset_index(drop=True)


def _bootstrap_mean(values: Sequence[float], replicates: int, seed: int) -> tuple[float, float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(int(seed))
    sample = rng.integers(0, array.size, size=(int(replicates), array.size))
    means = array[sample].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(array.mean()), float(low), float(high)


def _wilcoxon_greater(values: Sequence[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return np.nan, np.nan
    if np.allclose(array, 0.0, rtol=0.0, atol=0.0):
        return 0.0, 1.0
    result = wilcoxon(array, alternative="greater", zero_method="wilcox")
    return float(result.statistic), float(result.pvalue)


def _as_float(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return np.nan
    return number if np.isfinite(number) else np.nan


def _points_from_values(lon_value: Any, lat_value: Any) -> list[tuple[float, float]]:
    """Return finite ``(lon, lat)`` points from scalar/list-like fields."""

    lon_value, lat_value = _parse_literal(lon_value), _parse_literal(lat_value)
    if lon_value is None or lat_value is None:
        return []
    if isinstance(lon_value, dict) and lat_value is None:
        lon_value, lat_value = lon_value.get("lon"), lon_value.get("lat")
    try:
        lon_array = np.asarray(lon_value, dtype="f8")
        lat_array = np.asarray(lat_value, dtype="f8")
    except (TypeError, ValueError):
        return []
    if lon_array.ndim == 0 and lat_array.ndim == 0:
        pair = (_as_float(lon_array), _as_float(lat_array))
        return [pair] if np.isfinite(pair).all() else []
    lon_flat, lat_flat = lon_array.ravel(), lat_array.ravel()
    if lon_flat.size != lat_flat.size:
        return []
    return [
        (float(lon), float(lat))
        for lon, lat in zip(lon_flat, lat_flat)
        if np.isfinite(lon) and np.isfinite(lat)
    ]


def _points_from_row(row: Mapping[str, Any], prefix: str) -> list[tuple[float, float]]:
    """Read a footprint from separate lon/lat columns or a pair-list column."""

    lon_names = (f"{prefix}_lon", f"{prefix}_longitude")
    lat_names = (f"{prefix}_lat", f"{prefix}_latitude")
    for lon_name in lon_names:
        for lat_name in lat_names:
            if lon_name in row and lat_name in row:
                points = _points_from_values(row[lon_name], row[lat_name])
                if points:
                    return points
    for name in (
        prefix,
        f"{prefix}_points",
        f"{prefix}_polygon",
        f"{prefix}_footprint",
    ):
        if name not in row or _is_missing(row[name]):
            continue
        parsed = _parse_literal(row[name])
        if isinstance(parsed, dict):
            return _points_from_values(parsed.get("lon"), parsed.get("lat"))
        try:
            array = np.asarray(parsed, dtype="f8")
        except (TypeError, ValueError):
            continue
        if array.ndim == 2 and array.shape[1] >= 2:
            return [
                (float(lon), float(lat))
                for lon, lat in array[:, :2]
                if np.isfinite(lon) and np.isfinite(lat)
            ]
    return []


def _wrap_longitudes(longitudes: np.ndarray, reference: float) -> np.ndarray:
    return reference + ((longitudes - reference + 180.0) % 360.0) - 180.0


def _polygon_contains(points: Sequence[tuple[float, float]], lon: float, lat: float) -> bool:
    if len(points) < 3 or not np.isfinite([lon, lat]).all():
        return False
    vertices = np.asarray(points, dtype="f8")
    vertices[:, 0] = _wrap_longitudes(vertices[:, 0], lon)
    polygon = PolygonPath(vertices, closed=True)
    return bool(polygon.contains_point((lon, lat), radius=1.0e-10))


def _polygon_contains_many(
    points: Sequence[tuple[float, float]], locations: np.ndarray
) -> np.ndarray:
    if len(points) < 3 or locations.size == 0:
        return np.zeros(locations.shape[0], dtype=bool)
    vertices = np.asarray(points, dtype="f8")
    reference = float(np.nanmean(vertices[:, 0]))
    shifted = vertices.copy()
    shifted[:, 0] = _wrap_longitudes(shifted[:, 0], reference)
    locations = np.asarray(locations, dtype="f8").copy()
    locations[:, 0] = _wrap_longitudes(locations[:, 0], reference)
    return PolygonPath(shifted, closed=True).contains_points(locations, radius=1.0e-10)


def _haversine_km(lon_a: Any, lat_a: Any, lon_b: Any, lat_b: Any) -> np.ndarray:
    """Great-circle distance in km with NumPy broadcasting."""

    lon_a, lat_a, lon_b, lat_b = [
        np.deg2rad(np.asarray(value, dtype="f8"))
        for value in (lon_a, lat_a, lon_b, lat_b)
    ]
    dlon = (lon_b - lon_a + np.pi) % (2.0 * np.pi) - np.pi
    dlat = lat_b - lat_a
    hav = np.sin(dlat / 2.0) ** 2 + np.cos(lat_a) * np.cos(lat_b) * np.sin(dlon / 2.0) ** 2
    return 2.0 * 6371.0088 * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0)))


def _contour_rows_for_day(
    daily: pd.DataFrame, date_text: str, valid_only: bool = True
) -> list[dict[str, Any]]:
    if daily.empty:
        return []
    subset = daily.loc[daily["date"].astype(str).eq(date_text)].copy()
    if "is_virtual" in subset:
        subset = subset.loc[~subset["is_virtual"].astype(bool)]
    if valid_only and "filter_valid" in subset:
        subset = subset.loc[subset["filter_valid"].astype(bool)]
    if valid_only and "boundary_censored" in subset:
        subset = subset.loc[~subset["boundary_censored"].astype(bool)]
    return subset.to_dict("records")


def _covered_rows(
    rows: Sequence[Mapping[str, Any]], lon: float, lat: float, contour_kind: str
) -> list[Mapping[str, Any]]:
    covered = []
    for row in rows:
        points = _points_from_row(row, contour_kind)
        if _polygon_contains(points, lon, lat):
            covered.append(row)
    return covered


def _nearest_contour_distance_km(
    rows: Sequence[Mapping[str, Any]], lon: float, lat: float, contour_kind: str
) -> float:
    distances: list[float] = []
    for row in rows:
        points = _points_from_row(row, contour_kind)
        if not points:
            continue
        values = _haversine_km(
            lon,
            lat,
            np.asarray([point[0] for point in points]),
            np.asarray([point[1] for point in points]),
        )
        if values.size:
            distances.append(float(np.nanmin(values)))
    return float(np.nanmin(distances)) if distances else np.nan


def _nearest_row(rows: Sequence[Mapping[str, Any]], lon: float, lat: float) -> Mapping[str, Any] | None:
    if not rows:
        return None
    distances = []
    for row in rows:
        center_lon, center_lat = _as_float(row.get("center_lon")), _as_float(row.get("center_lat"))
        if np.isfinite([center_lon, center_lat]).all():
            distances.append((float(_haversine_km(lon, lat, np.asarray([center_lon]), np.asarray([center_lat]))[0]), row))
    return min(distances, key=lambda item: item[0])[1] if distances else None


def _footprint_coverage(
    rows: Sequence[Mapping[str, Any]], points: Sequence[tuple[float, float]], contour_kind: str
) -> tuple[float, int]:
    if not points:
        return np.nan, 0
    locations = np.asarray(points, dtype="f8")
    covered = np.zeros(locations.shape[0], dtype=bool)
    for row in rows:
        covered |= _polygon_contains_many(_points_from_row(row, contour_kind), locations)
    return float(covered.mean()), int(locations.shape[0])


def _grid_point_index(
    lon: float, lat: float, grid_lon: np.ndarray, grid_lat: np.ndarray
) -> tuple[int, int]:
    x = int(np.abs(grid_lon - lon).argmin())
    y = int(np.abs(grid_lat - lat).argmin())
    return x, y


def _grid_ring_metrics(
    lon: float,
    lat: float,
    prepared: Any,
    inner_km: float,
    outer_km: float,
) -> dict[str, Any]:
    """Calculate a point-centered, cell-area-weighted local ring null."""

    grid_lon = np.asarray(prepared.lon, dtype="f8")
    grid_lat = np.asarray(prepared.lat, dtype="f8")
    lon_grid, lat_grid = np.meshgrid(grid_lon, grid_lat, indexing="ij")
    distance = _haversine_km(lon, lat, lon_grid, lat_grid)
    ring = (distance >= inner_km) & (distance <= outer_km)
    filter_valid = np.asarray(prepared.filter_valid_mask, dtype=bool)
    ocean_valid = np.asarray(prepared.ocean_valid_mask, dtype=bool)
    valid = ring & filter_valid & ocean_valid
    # Regular longitude/latitude grid: the cos(latitude) factor is enough for
    # relative area weights and keeps this calculation independent of a
    # separate geodesic-cell package.
    weights = np.cos(np.deg2rad(lat_grid))
    total_weight = float(weights[ring].sum())
    ocean_weight = float(weights[ring & ocean_valid].sum())
    valid_weight = float(weights[valid].sum())
    x, y = _grid_point_index(lon, lat, grid_lon, grid_lat)
    return {
        "ring_total_cells": int(ring.sum()),
        "ring_ocean_cells": int((ring & ocean_valid).sum()),
        "ring_filter_valid_ocean_cells": int(valid.sum()),
        "ring_ocean_fraction": ocean_weight / total_weight if total_weight else np.nan,
        "ring_filter_valid_fraction": valid_weight / ocean_weight if ocean_weight else np.nan,
        "ring_valid_fraction": valid_weight / total_weight if total_weight else np.nan,
        "ring_valid_weight": valid_weight,
        "ring_mask": ring,
        "valid_mask": valid,
        "grid_lon": grid_lon,
        "grid_lat": grid_lat,
        # The event core is not part of the 120--240 km ring.  Its primary
        # eligibility therefore uses the point masks directly, not the ring
        # mask (using ``valid[x, y]`` would censor every event by construction).
        "grid_point_valid": bool(filter_valid[x, y] and ocean_valid[x, y]),
        "grid_point_ocean_valid": bool(ocean_valid[x, y]),
        "grid_point_filter_valid": bool(filter_valid[x, y]),
    }


def _exclude_footprint_cells(
    metrics: dict[str, Any], footprint_points: Sequence[tuple[float, float]]
) -> dict[str, Any]:
    """Exclude DO peak-footprint cells from the local background ring."""

    if not footprint_points:
        return metrics
    valid = np.asarray(metrics["valid_mask"], dtype=bool).copy()
    for point_lon, point_lat in footprint_points:
        x, y = _grid_point_index(
            float(point_lon), float(point_lat), metrics["grid_lon"], metrics["grid_lat"]
        )
        valid[x, y] = False
    metrics = dict(metrics)
    metrics["valid_mask"] = valid
    metrics["ring_excluded_footprint_cells"] = int(
        np.count_nonzero(np.asarray(metrics["ring_mask"], dtype=bool))
        - np.count_nonzero(valid)
    )
    metrics["ring_filter_valid_ocean_cells"] = int(np.count_nonzero(valid))
    return metrics


def _ring_occupancy(
    rows: Sequence[Mapping[str, Any]], metrics: Mapping[str, Any]
) -> float:
    valid = np.asarray(metrics["valid_mask"], dtype=bool)
    if not valid.any():
        return np.nan
    lon_grid, lat_grid = np.meshgrid(metrics["grid_lon"], metrics["grid_lat"], indexing="ij")
    locations = np.column_stack([lon_grid[valid], lat_grid[valid]])
    covered = np.zeros(locations.shape[0], dtype=bool)
    for row in rows:
        covered |= _polygon_contains_many(_points_from_row(row, "effective_contour"), locations)
    weights = np.cos(np.deg2rad(locations[:, 1]))
    return float(weights[covered].sum() / weights.sum()) if weights.sum() else np.nan


def _load_valid_grid(
    date_text: str,
    config: Mapping[str, Any],
    eta_root: FilePath,
    cache: dict[str, Any],
) -> Any:
    if date_text not in cache:
        snapshot = load_eta_snapshot(
            FilePath(eta_root) / f"eta.{parse_date(date_text):%m.%d.%Y}.nc",
            expected_date=date_text,
            raw_units=str(config["input"]["raw_eta_units"]),
        )
        cache[date_text] = prepare_pet_grid(snapshot, config)
    return cache[date_text]


def _extract_primary_contour_properties(
    rows: Sequence[Mapping[str, Any]],
    core_lon: float,
    core_lat: float,
) -> dict[str, Any]:
    effective = _covered_rows(rows, core_lon, core_lat, "effective_contour")
    speed = _covered_rows(rows, core_lon, core_lat, "speed_contour")
    nearest = _nearest_row(rows, core_lon, core_lat)
    effective_ids = [str(row.get("object_id")) for row in effective]
    speed_ids = [str(row.get("object_id")) for row in speed]
    all_polarities = sorted({str(row.get("polarity")) for row in effective})
    primary = nearest
    if effective:
        primary = min(
            effective,
            key=lambda row: float(
                _haversine_km(
                    core_lon,
                    core_lat,
                    np.asarray([_as_float(row.get("center_lon"))]),
                    np.asarray([_as_float(row.get("center_lat"))]),
                )[0]
            ),
        )
    return {
        "effective_core_contained": bool(effective),
        "speed_core_contained": bool(speed),
        "effective_core_contained_anticyclonic": any(row.get("polarity") == "anticyclonic" for row in effective),
        "effective_core_contained_cyclonic": any(row.get("polarity") == "cyclonic" for row in effective),
        "effective_object_ids": json.dumps(effective_ids, ensure_ascii=False),
        "speed_object_ids": json.dumps(speed_ids, ensure_ascii=False),
        "effective_polarities": json.dumps(all_polarities, ensure_ascii=False),
        "primary_object_id": str(primary.get("object_id")) if primary else None,
        "primary_track_id": primary.get("track_id") if primary else None,
        "primary_polarity": primary.get("polarity") if primary else None,
        "nearest_effective_contour_distance_km": _nearest_contour_distance_km(
            rows, core_lon, core_lat, "effective_contour"
        ),
        "nearest_effective_center_distance_km": (
            float(
                _haversine_km(
                    core_lon,
                    core_lat,
                    np.asarray([_as_float(nearest.get("center_lon"))]),
                    np.asarray([_as_float(nearest.get("center_lat"))]),
                )[0]
            )
            if nearest
            else np.nan
        ),
        "nearest_effective_radius_km": _as_float(nearest.get("effective_radius_km")) if nearest else np.nan,
    }


def _event_surface_ro_category(row: Mapping[str, Any]) -> str:
    value = row.get("surface_ro_expression")
    if not _is_missing(value):
        text = str(value).strip().lower()
        if text in {"none", "no", "absent", "no-expression", "no_expression"}:
            return "none"
        if "same" in text or "present" in text or "express" in text:
            return "expressed"
    same = _as_bool(row.get("surface_ro_same_polarity"))
    if same is True:
        return "same-polarity"
    polarity = _sign(row.get("surface_ro_polarity"))
    return "expressed" if np.isfinite(polarity) else "unknown"


def _read_manifest(path: FilePath) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"Manifest is not an object: {path}")
    return value


def _validate_complete_manifest(
    manifest_path: FilePath, config: Mapping[str, Any]
) -> dict[str, Any]:
    manifest = _read_manifest(manifest_path)
    if manifest.get("status") != "complete":
        raise RuntimeError("Event association requires manifest status=complete")
    expected_dates = [day.isoformat() for day in date_range(config["input"]["start_date"], config["input"]["end_date"])]
    days = manifest.get("days", {})
    if list(days) != expected_dates:
        raise RuntimeError("Manifest date inventory does not exactly match the locked 365-day window")
    incomplete = [date_text for date_text in expected_dates if days[date_text].get("status") != "complete"]
    if incomplete:
        raise RuntimeError(f"Manifest has incomplete dates: {incomplete[:10]}")
    for date_text in expected_dates:
        info = days[date_text]
        daily_path = FilePath(info["daily_path"])
        if not daily_path.is_file():
            raise FileNotFoundError(daily_path)
        expected_hash = info.get("daily_sha256")
        if expected_hash and sha256_file(daily_path) != expected_hash:
            raise RuntimeError(f"Daily input hash mismatch for tracking: {daily_path}")
    if manifest.get("tracking", {}).get("status") != "complete":
        raise RuntimeError("Run track_surface_eddies first; association requires complete tracking outputs")
    tracking_paths = manifest["tracking"].get("paths", {})
    for key in ("daily", "tracks", "track_observations"):
        output = tracking_paths.get(key)
        if not output or not FilePath(output).is_file():
            raise FileNotFoundError(f"Missing tracked PET output {key}: {output}")
    return manifest


def _load_daily_objects(manifest: Mapping[str, Any]) -> pd.DataFrame:
    output_root = FilePath(manifest["output_root"])
    aggregate = output_root / "surface_eddy_daily_objects.parquet"
    if aggregate.is_file():
        frame = pd.read_parquet(aggregate)
    else:
        frames = [pd.read_parquet(info["daily_path"]) for info in manifest["days"].values()]
        frame = pd.concat(frames, ignore_index=True, sort=False)
    missing = sorted(set(DAILY_COLUMNS).difference(frame.columns))
    if missing:
        raise ValueError(f"Daily PET table lacks stable fields: {missing}")
    frame = frame.loc[~frame["is_virtual"].astype(bool)].copy()
    frame["date"] = frame["date"].astype(str)
    return frame[list(DAILY_COLUMNS)].sort_values(["date", "polarity", "object_id"]).reset_index(drop=True)


def _annual_point_occupancy(
    daily: pd.DataFrame,
    event_date: str,
    lon: float,
    lat: float,
) -> tuple[float, int, str]:
    """Eulerian month/latitude-stratified null at the event core location."""

    date_value = parse_date(event_date)
    # Match the existing Zhu annual occupancy convention: month x 1-degree
    # latitude strata, used only as the secondary Eulerian sanity null.
    lat_band = f"{math.floor(lat):g}_{math.floor(lat) + 1:g}N"
    month = date_value.month
    dates = sorted(
        {
            str(value)
            for value in daily.loc[
                pd.to_datetime(daily["date"]).dt.month.eq(month), "date"
            ]
        }
    )
    covered_days = 0
    for day in dates:
        rows = _contour_rows_for_day(daily, day, valid_only=True)
        if _covered_rows(rows, lon, lat, "effective_contour"):
            covered_days += 1
    return (
        covered_days / len(dates) if dates else np.nan,
        len(dates),
        f"month={month:02d};latitude_band={lat_band}",
    )


def _event_record(
    event: Mapping[str, Any],
    daily: pd.DataFrame,
    config: Mapping[str, Any],
    eta_root: FilePath,
    grid_cache: dict[str, Any],
) -> dict[str, Any]:
    event_id = str(event["event_id"])
    date_text = str(event["peak_date"])
    core_lon, core_lat = _as_float(event["peak_core_lon"]), _as_float(event["peak_core_lat"])
    rows = _contour_rows_for_day(daily, date_text, valid_only=True)
    pet = _extract_primary_contour_properties(rows, core_lon, core_lat)
    prepared = _load_valid_grid(date_text, config, eta_root, grid_cache)
    association = config["association"]
    footprint_points = _points_from_row(event, "do_peak_footprint")
    if not footprint_points:
        footprint_points = _points_from_row(event, "core_footprint")
    ring = _grid_ring_metrics(
        core_lon,
        core_lat,
        prepared,
        float(association["ring_inner_km"]),
        float(association["ring_outer_km"]),
    )
    ring = _exclude_footprint_cells(ring, footprint_points)
    ring_occupancy = _ring_occupancy(rows, ring)
    annual_occupancy, annual_days, annual_stratum = _annual_point_occupancy(
        daily, date_text, core_lon, core_lat
    )
    footprint_effective_fraction, footprint_n = _footprint_coverage(
        rows, footprint_points, "effective_contour"
    )
    footprint_speed_fraction, _ = _footprint_coverage(rows, footprint_points, "speed_contour")
    event_dict: dict[str, Any] = {
        "event_id": event_id,
        "peak_date": date_text,
        "peak_lon": _as_float(event.get("peak_lon")),
        "peak_lat": _as_float(event.get("peak_lat")),
        "peak_core_lon": core_lon,
        "peak_core_lat": core_lat,
        "strict_eligible": True,
        "population_background_ring_complete": event.get("background_ring_complete", np.nan),
        "quality_eligible": event.get("quality_eligible", np.nan),
        "rotation_dominated": event.get("rotation_dominated", np.nan),
        "deep_ro_sign": event.get("deep_ro_sign", np.nan),
        "surface_ro_expression": _event_surface_ro_category(event),
        "surface_ro_polarity": event.get("surface_ro_polarity", np.nan),
        "surface_ro_same_polarity": event.get("surface_ro_same_polarity", np.nan),
        "mccoy_center": event.get("mccoy_center", np.nan),
        "mccoy_center_velocity": event.get("mccoy_center_velocity", np.nan),
        "mccoy_any17": event.get("mccoy_any17", np.nan),
        "mccoy_any17_velocity": event.get("mccoy_any17_velocity", np.nan),
        **pet,
        **{key: value for key, value in ring.items() if not isinstance(value, np.ndarray)},
        "ring_occupancy_fraction": ring_occupancy,
        "annual_month_latitude_occupancy_fraction": annual_occupancy,
        "annual_null_days": annual_days,
        "annual_null_stratum": annual_stratum,
        "do_peak_footprint_n": footprint_n,
        "do_peak_footprint_effective_overlap_fraction": footprint_effective_fraction,
        "do_peak_footprint_speed_overlap_fraction": footprint_speed_fraction,
        "peak_core_filter_valid": ring["grid_point_filter_valid"],
        "peak_core_ocean_valid": ring["grid_point_ocean_valid"],
        "peak_core_primary_eligible": bool(
            ring["grid_point_valid"] and ring["ring_filter_valid_ocean_cells"] > 0
        ),
        "ring_excluded_footprint_cells": int(
            ring.get("ring_excluded_footprint_cells", 0)
        ),
    }
    event_dict["effective_core_contained_raw"] = event_dict[
        "effective_core_contained"
    ]
    event_dict["speed_core_contained_raw"] = event_dict["speed_core_contained"]
    if not event_dict["peak_core_primary_eligible"]:
        event_dict["effective_core_contained"] = np.nan
        event_dict["speed_core_contained"] = np.nan
    event_dict["core_minus_ring_occupancy"] = (
        float(event_dict["effective_core_contained"]) - ring_occupancy
        if np.isfinite(ring_occupancy)
        and _as_bool(event_dict["peak_core_primary_eligible"]) is True
        else np.nan
    )
    event_dict["core_minus_annual_occupancy"] = (
        float(event_dict["effective_core_contained"]) - annual_occupancy
        if np.isfinite(annual_occupancy)
        and _as_bool(event_dict["peak_core_primary_eligible"]) is True
        else np.nan
    )
    primary_polarity = event_dict.get("primary_polarity")
    deep_sign = _sign(event_dict.get("deep_ro_sign"))
    primary_code = {"anticyclonic": 1, "cyclonic": -1}.get(primary_polarity, np.nan)
    if not np.isfinite(deep_sign) or not np.isfinite(primary_code):
        event_dict["rotation_pet_category"] = "no-PET" if not event_dict["effective_core_contained"] else "ambiguous"
    elif primary_code == deep_sign:
        event_dict["rotation_pet_category"] = "same"
    else:
        event_dict["rotation_pet_category"] = "opposite"
    if np.isfinite(event_dict["nearest_effective_radius_km"]) and event_dict["nearest_effective_radius_km"] > 0:
        event_dict["nearest_effective_center_distance_over_radius"] = (
            event_dict["nearest_effective_center_distance_km"]
            / event_dict["nearest_effective_radius_km"]
        )
    else:
        event_dict["nearest_effective_center_distance_over_radius"] = np.nan
    return event_dict


def _bootstrap_mean_ci(values: pd.Series, replicates: int, seed: int) -> tuple[float, float, float]:
    clean = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype="f8")
    if clean.size == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    samples = rng.choice(clean, size=(replicates, clean.size), replace=True).mean(axis=1)
    return float(clean.mean()), float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def _paired_wilcoxon(values: pd.Series) -> dict[str, Any]:
    clean = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype="f8")
    if clean.size < 2 or np.allclose(clean, 0.0):
        return {"n": int(clean.size), "statistic": np.nan, "pvalue": np.nan}
    result = wilcoxon(clean, zero_method="wilcox", alternative="two-sided", method="auto")
    return {"n": int(clean.size), "statistic": float(result.statistic), "pvalue": float(result.pvalue)}


def _validate_mccoy_counts(frame: pd.DataFrame) -> dict[str, int]:
    """Require the four pre-existing McCoy event-level counts."""

    expected = {
        "mccoy_center": 9,
        "mccoy_center_velocity": 6,
        "mccoy_any17": 19,
        "mccoy_any17_velocity": 11,
    }
    observed: dict[str, int] = {}
    for column, target in expected.items():
        if column not in frame:
            raise KeyError(
                f"McCoy event-level field {column!r} is absent; pass the explicit "
                "McCoy event_summary source or a population table containing it"
            )
        values = frame[column].map(_as_bool)
        count = int(values.eq(True).sum())
        if count != target:
            raise ValueError(
                f"McCoy count changed for {column}: expected {target}, found {count}"
            )
        observed[column] = count
    return observed


def _distribution(values: pd.Series) -> dict[str, Any]:
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if clean.empty:
        return {"n": 0, "p05": np.nan, "median": np.nan, "p95": np.nan}
    return {
        "n": int(clean.size),
        "p05": float(clean.quantile(0.05)),
        "median": float(clean.median()),
        "p95": float(clean.quantile(0.95)),
    }


def _crosstab_with_categories(frame: pd.DataFrame, definition: str) -> pd.DataFrame:
    values = frame[definition].map(_as_bool)
    pet = frame["rotation_pet_category"].where(frame["effective_core_contained"], "no-PET")
    output = pd.crosstab(values, pet, dropna=False).reindex(index=[False, True], fill_value=0)
    output = output.reindex(columns=PRIMARY_PET_CATEGORIES, fill_value=0)
    output.index = ["false", "true"]
    output.index.name = "mccoy_positive"
    output.insert(0, "definition", definition)
    return output.reset_index()


def run_association(
    repo_root: str | FilePath,
    config_path: str | FilePath,
    manifest_path: str | FilePath,
    population_input: str | FilePath,
    diagnostics_input: str | FilePath,
    catalog_input: str | FilePath,
) -> dict[str, str]:
    """Run the locked 56-event PET association and write all event tables."""

    repo_root = FilePath(repo_root).resolve()
    config_path = FilePath(config_path).resolve()
    manifest_path = FilePath(manifest_path).resolve()
    config, postprocess_config_path = _load_postprocess_config(config_path)
    association_lock_path, association_lock_sha256 = _require_association_lock(config_path)
    manifest = _validate_complete_manifest(manifest_path, config)
    population_path = _resolve_input_file(population_input, "population")
    diagnostics_path = _resolve_input_file(diagnostics_input, "diagnostics")
    catalog_path = _resolve_input_file(catalog_input, "catalog")
    population = _read_table(population_path)
    diagnostics = _read_table(diagnostics_path)
    catalog = _read_table(catalog_path)
    expected_candidates = int(config["association"]["expected_candidate_events"])
    if len(population) != expected_candidates:
        raise ValueError(
            f"Authoritative population candidate count is {len(population)}, "
            f"expected {expected_candidates}"
        )
    _validate_identity_consistency(population, diagnostics, catalog)
    merged = _coalesce_sources(population, diagnostics, catalog)
    mccoy_path: FilePath | None = None
    try:
        mccoy_path = _mccoy_summary_path(population_input)
    except FileNotFoundError:
        # Some already materialized population tables contain the McCoy
        # flags directly.  In that case no separate nested summary is needed.
        mccoy_path = None
    if mccoy_path is not None and mccoy_path not in {
        population_path,
        diagnostics_path,
        catalog_path,
    }:
        mccoy = _read_table(mccoy_path)
        _validate_identity_consistency(population, mccoy)
        merged = _coalesce_sources(merged, mccoy, None)
    strict, selection_audit = _select_strict_events(
        merged, int(config["association"]["expected_strict_events"])
    )
    mccoy_counts = _validate_mccoy_counts(strict)
    daily = _load_daily_objects(manifest)
    eta_root = FilePath(config["input"]["eta_root"])
    if not eta_root.is_absolute():
        eta_root = repo_root / eta_root
    grid_cache: dict[str, Any] = {}
    records = [
        _event_record(row, daily, config, eta_root.resolve(), grid_cache)
        for row in strict.to_dict("records")
    ]
    association = pd.DataFrame.from_records(records)
    association = association.sort_values("event_id").reset_index(drop=True)

    diagnostics_root = FilePath(diagnostics_input).expanduser().resolve()
    if diagnostics_root.is_file():
        diagnostics_root = diagnostics_root.parent
    sensitivity = _read_sensitivity_population(diagnostics_root, config)
    sensitivity_association = _sensitivity_association(sensitivity, daily)
    output_root = FilePath(manifest["output_root"]).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    event_association_path = output_root / "surface_eddy_event_association.parquet"
    association.to_parquet(event_association_path, index=False)

    null_rows = []
    for row in association.to_dict("records"):
        for null_type, value, comparison in (
            (
                "local_ring_120_240km",
                row["ring_occupancy_fraction"],
                row["core_minus_ring_occupancy"],
            ),
            (
                "annual_month_latitude_point_occupancy",
                row["annual_month_latitude_occupancy_fraction"],
                row["core_minus_annual_occupancy"],
            ),
        ):
            null_rows.append(
                {
                    "event_id": row["event_id"],
                    "peak_date": row["peak_date"],
                    "peak_core_lon": row["peak_core_lon"],
                    "peak_core_lat": row["peak_core_lat"],
                    "null_type": null_type,
                    "occupancy_fraction": value,
                    "core_indicator": row["effective_core_contained"],
                    "core_minus_null": comparison,
                    "ring_valid_fraction": row["ring_valid_fraction"],
                    "null_stratum": row["annual_null_stratum"] if "annual" in null_type else "same-day local ring",
                }
            )
    nulls = pd.DataFrame(null_rows)
    null_path = output_root / "surface_eddy_event_nulls.parquet"
    nulls.to_parquet(null_path, index=False)

    rotation = association.loc[association["rotation_dominated"].map(_as_bool).eq(True)].copy()
    expected_rotation = int(config["association"]["expected_rotation_events"])
    if len(rotation) != expected_rotation:
        raise ValueError(
            f"Rotation-dominated cross-tab has {len(rotation)} rows, expected {expected_rotation}; "
            "the authoritative diagnostics must provide the locked 29-event flag."
        )
    surface_ro = rotation[
        [
            "event_id",
            "deep_ro_sign",
            "surface_ro_expression",
            "surface_ro_polarity",
            "surface_ro_same_polarity",
            "effective_core_contained",
            "primary_polarity",
            "rotation_pet_category",
        ]
    ].copy()
    surface_ro["surface_ro_same_polarity"] = surface_ro["surface_ro_same_polarity"].map(_as_bool)
    surface_ro.insert(
        1,
        "pet_surface_expression",
        surface_ro["effective_core_contained"].map({True: "PET-effective-core", False: "no-PET-effective-core"}),
    )
    surface_ro_path = output_root / "surface_eddy_surface_ro_crosstab.parquet"
    surface_ro.to_parquet(surface_ro_path, index=False)

    mccoy_definitions = (
        "mccoy_center",
        "mccoy_center_velocity",
        "mccoy_any17",
        "mccoy_any17_velocity",
    )
    mccoy_frames = [_crosstab_with_categories(association, definition) for definition in mccoy_definitions]
    mccoy = pd.concat(mccoy_frames, ignore_index=True)
    mccoy_path = output_root / "surface_eddy_mccoy_crosstab.parquet"
    mccoy.to_parquet(mccoy_path, index=False)
    sensitivity_path = output_root / "surface_eddy_161_event_sensitivity.parquet"
    sensitivity_association.to_parquet(sensitivity_path, index=False)

    bootstrap_replicates = int(config["association"]["bootstrap_replicates"])
    seed = int(config["association"]["random_seed"])
    local_mean, local_low, local_high = _bootstrap_mean_ci(
        association["core_minus_ring_occupancy"], bootstrap_replicates, seed
    )
    annual_mean, annual_low, annual_high = _bootstrap_mean_ci(
        association["core_minus_annual_occupancy"], bootstrap_replicates, seed + 1
    )
    tracks_path = output_root / "surface_eddy_tracks.parquet"
    tracks = pd.read_parquet(tracks_path) if tracks_path.is_file() else pd.DataFrame()
    observations_path = output_root / "surface_eddy_track_observations.parquet"
    observations = pd.read_parquet(observations_path) if observations_path.is_file() else pd.DataFrame()
    daily_object_counts = [
        {"date": str(row["date"]), "polarity": str(row["polarity"]), "objects": int(row["objects"])}
        for row in daily.groupby(["date", "polarity"]).size().reset_index(name="objects").to_dict("records")
    ]
    duration_counts = (
        {str(key): int(value) for key, value in tracks["duration_class"].value_counts().items()}
        if "duration_class" in tracks
        else {}
    )
    duration_total = sum(duration_counts.values())
    summary: dict[str, Any] = {
        "run_id": manifest["run_id"],
        "manifest": str(manifest_path),
        "event_input_files": {
            "population": {"path": str(population_path), "sha256": sha256_file(population_path), "size": population_path.stat().st_size},
            "diagnostics": {"path": str(diagnostics_path), "sha256": sha256_file(diagnostics_path), "size": diagnostics_path.stat().st_size},
            "catalog": {"path": str(catalog_path), "sha256": sha256_file(catalog_path), "size": catalog_path.stat().st_size},
        },
        "selection": selection_audit,
        "primary_denominator": {
            "strict_input": int(len(association)),
            "core_filter_valid_and_ocean_valid": int(association["peak_core_filter_valid"].astype(bool).mul(association["peak_core_ocean_valid"].astype(bool)).sum()),
            "primary_eligible": int(association["peak_core_primary_eligible"].astype(bool).sum()),
            "effective_core_contained": int(association["effective_core_contained"].astype(bool).sum()),
            "speed_core_contained": int(association["speed_core_contained"].astype(bool).sum()),
        },
        "primary_estimate": {
            "core_minus_local_ring_mean": local_mean,
            "core_minus_local_ring_bootstrap_95ci": [local_low, local_high],
            "paired_wilcoxon_core_minus_local_ring": _paired_wilcoxon(association["core_minus_ring_occupancy"]),
            "core_minus_annual_stratified_mean": annual_mean,
            "core_minus_annual_stratified_bootstrap_95ci": [annual_low, annual_high],
            "paired_wilcoxon_core_minus_annual_stratified": _paired_wilcoxon(association["core_minus_annual_occupancy"]),
        },
        "rotation_events": {
            "n": int(len(rotation)),
            "pet_category_counts": rotation["rotation_pet_category"].value_counts(dropna=False).to_dict(),
            "surface_ro_expression_counts": rotation["surface_ro_expression"].value_counts(dropna=False).to_dict(),
        },
        "mccoy_crosstab": {
            definition: {
                "available": bool(association[definition].notna().any()),
                "positive_count": int(association[definition].map(_as_bool).eq(True).sum()),
            }
            for definition in mccoy_definitions
        },
        "mccoy_locked_counts": mccoy_counts,
        "sensitivity_161": {
            "n": int(len(sensitivity_association)),
            "effective_core_contained": int(
                sensitivity_association["pet_effective_contained"].sum()
            ),
            "speed_core_contained": int(
                sensitivity_association["pet_speed_contained"].sum()
            ),
        },
        "daily_object_counts": daily_object_counts,
        "object_distributions": {
            "effective_radius_km": _distribution(daily["effective_radius_km"]),
            "speed_radius_km": _distribution(daily["speed_radius_km"]),
            "amplitude_m": _distribution(daily["amplitude_m"]),
        },
        "tracks": {
            "count": int(len(tracks)),
            "duration_class_counts": duration_counts,
            "duration_class_proportions": {
                key: (value / duration_total if duration_total else np.nan)
                for key, value in duration_counts.items()
            },
            "virtual_observation_fraction": float(observations["is_virtual"].astype(bool).mean()) if "is_virtual" in observations and len(observations) else np.nan,
            "boundary_censored_fraction": float(daily["boundary_censored"].astype(bool).mean()) if len(daily) else np.nan,
        },
        "paths": {
            "association": str(event_association_path),
            "nulls": str(null_path),
            "surface_ro_crosstab": str(surface_ro_path),
            "mccoy_crosstab": str(mccoy_path),
            "sensitivity_161": str(sensitivity_path),
        },
    }
    summary_path = output_root / "surface_eddy_summary.json"
    _write_json(summary_path, summary)
    association_manifest = {
        "schema_version": 1,
        "status": "complete",
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "catalog_manifest": str(manifest_path),
        "catalog_run_id": manifest["run_id"],
        "protocol_sha256": sha256_file(repo_root / "ofes-surface-eddy-analysis-lock.md"),
        "config_sha256": sha256_file(config_path),
        "postprocessor_sha256": sha256_file(FilePath(__file__)),
        "event_input_files": summary["event_input_files"],
        "selection": selection_audit,
        "outputs": summary["paths"],
    }
    association_manifest_path = output_root / "surface_eddy_event_association_manifest.json"
    _write_json(association_manifest_path, association_manifest)
    summary["paths"]["manifest"] = str(association_manifest_path)
    _write_json(summary_path, summary)
    return {key: str(value) for key, value in {"association": event_association_path, "nulls": null_path, "surface_ro_crosstab": surface_ro_path, "mccoy_crosstab": mccoy_path, "summary": summary_path, "manifest": association_manifest_path}.items()}


def run_authoritative_association(
    repo_root: str | FilePath,
    config_path: str | FilePath,
    manifest_path: str | FilePath,
    population_input: str | FilePath,
    diagnostics_input: str | FilePath,
    catalog_input: str | FilePath,
    mccoy_summary: str | FilePath | None = None,
) -> dict[str, str]:
    """Run association using the known OFES DO50 artifact schemas.

    The population table supplies the locked 59 -> 56 selection and the
    rotation fields.  The DO catalog supplies peak footprints and lifecycle
    positions.  Event-diagnostics supplies the precomputed McCoy-compatible
    event flags and the 161-event sensitivity ranking.  None of these inputs
    is opened until the complete PET manifest gate passes.
    """

    repo_root = FilePath(repo_root).resolve()
    config_path = FilePath(config_path).resolve()
    config, postprocess_config_path = _load_postprocess_config(config_path)
    association_lock_path, association_lock_sha256 = _require_association_lock(config_path)
    manifest_file, manifest = _require_complete_pet_manifest(
        manifest_path,
        repo_root=repo_root,
        config_path=config_path,
    )
    population, strict, population_file = _read_population_authoritative(population_input, config)
    do_daily, catalog_root = _read_do_catalog(manifest, catalog_input)
    _validate_peak_identity(strict, do_daily)
    objects, observations, tracks = _read_surface_outputs(manifest)
    domain = _load_surface_domain(manifest, config)

    event_rows = [
        _event_row(row, objects, domain, catalog_root, config)
        for row in strict.to_dict("records")
    ]
    association = pd.DataFrame(event_rows).sort_values("event_id", kind="mergesort").reset_index(drop=True)
    association["peak_date"] = pd.to_datetime(association["peak_date"]).dt.normalize()
    association["peak_core_contained_by_actual_pet_effective_contour"] = association[
        "peak_core_contained_by_actual_pet_effective_contour"
    ]
    association["peak_core_contained_by_actual_pet_speed_contour"] = association[
        "peak_core_contained_by_actual_pet_speed_contour"
    ]
    association["peak_core_containment_any_polarity"] = association[
        "peak_core_contained_by_actual_pet_effective_contour"
    ]
    association["pet_peak_category"] = np.where(
        association["peak_core_contained_by_actual_pet_effective_contour"].map(_as_bool).eq(True),
        "PET-present",
        "no-PET",
    )

    lifecycle = _lifecycle_association(strict, do_daily, observations)
    association = association.merge(lifecycle, on="event_id", how="left", validate="one_to_one")

    dates = [pd.Timestamp(day) for day in date_range(config["input"]["start_date"], config["input"]["end_date"])]
    annual = _annual_occupancy(objects, domain, dates)
    annual_path = FilePath(manifest["output_root"]) / "surface_eddy_annual_occupancy.parquet"
    annual.to_parquet(annual_path, index=False)
    annual_lookup = {
        (int(row["month"]), int(row["lat_bin_deg"])): float(row["occupancy_fraction"])
        for row in annual.to_dict("records")
    }
    association["annual_month_latitude_occupancy_fraction"] = [
        annual_lookup.get((date.month, int(np.floor(lat))), np.nan)
        for date, lat in zip(association["peak_date"], association["peak_lat"])
    ]
    association["core_minus_local_ring_occupancy"] = (
        association["peak_core_contained_by_actual_pet_effective_contour"].astype(float)
        - association["ring_occupancy_fraction"]
    )
    association["core_minus_annual_stratified_occupancy"] = (
        association["peak_core_contained_by_actual_pet_effective_contour"].astype(float)
        - association["annual_month_latitude_occupancy_fraction"]
    )
    association["annual_null_stratum"] = [
        f"month={date.month:02d};latitude_band={int(np.floor(lat))}--{int(np.floor(lat)) + 1}N"
        for date, lat in zip(association["peak_date"], association["peak_lat"])
    ]

    output_root = FilePath(manifest["output_root"]).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    association_path = output_root / "surface_eddy_event_association.parquet"
    association.to_parquet(association_path, index=False)

    null_rows: list[dict[str, Any]] = []
    for row in association.to_dict("records"):
        for null_type, value, difference, stratum in (
            (
                "local_ring_120_240km",
                row["ring_occupancy_fraction"],
                row["core_minus_local_ring_occupancy"],
                "same-day 120--240 km filter-valid/ocean-valid ring",
            ),
            (
                "annual_month_latitude_occupancy",
                row["annual_month_latitude_occupancy_fraction"],
                row["core_minus_annual_stratified_occupancy"],
                row["annual_null_stratum"],
            ),
        ):
            null_rows.append(
                {
                    "event_id": row["event_id"],
                    "peak_date": row["peak_date"],
                    "peak_lon": row["peak_lon"],
                    "peak_lat": row["peak_lat"],
                    "null_type": null_type,
                    "occupancy_fraction": value,
                    "core_indicator": row["peak_core_contained_by_actual_pet_effective_contour"],
                    "core_minus_null": difference,
                    "ring_valid_cell_count": row["ring_valid_cell_count"],
                    "ring_valid_area_m2": row["ring_valid_area_m2"],
                    "stratum": stratum,
                }
            )
    null_path = output_root / "surface_eddy_event_nulls.parquet"
    pd.DataFrame(null_rows).to_parquet(null_path, index=False)

    rotation = _rotation_crosstab(strict, association, config)
    association = association.drop(columns=["pet_peak_category"]).merge(
        rotation[["event_id", "pet_peak_category"]], on="event_id", how="left", validate="one_to_one"
    )
    association.to_parquet(association_path, index=False)
    surface_ro_path = output_root / "surface_eddy_surface_ro_crosstab.parquet"
    rotation.to_parquet(surface_ro_path, index=False)

    resolved_mccoy_summary = mccoy_summary
    if resolved_mccoy_summary is None:
        population_root = FilePath(population_input).expanduser().resolve()
        if population_root.is_file():
            population_root = population_root.parent
        resolved_mccoy_summary = _mccoy_summary_path(population_root)
    mccoy_frame, mccoy_path_input = _mccoy_crosstab(
        strict, association, diagnostics_input, exact_path=resolved_mccoy_summary
    )
    mccoy_path = output_root / "surface_eddy_mccoy_crosstab.parquet"
    mccoy_frame.to_parquet(mccoy_path, index=False)

    sensitivity = _read_sensitivity_population(diagnostics_input, config)
    _validate_identity_consistency(strict, sensitivity)
    sensitivity_result = _sensitivity_association(sensitivity, objects)
    sensitivity_path = output_root / "surface_eddy_quality_eligible_161.parquet"
    sensitivity_result.to_parquet(sensitivity_path, index=False)
    plot_paths = _write_plot_outputs(
        output_root,
        config,
        strict,
        association,
        objects,
        observations,
        tracks,
        domain,
        catalog_root,
    )

    local_values = association["core_minus_local_ring_occupancy"].to_numpy(dtype=float)
    annual_values = association["core_minus_annual_stratified_occupancy"].to_numpy(dtype=float)
    local_mean, local_low, local_high = _bootstrap_mean(
        local_values, int(config["association"]["bootstrap_replicates"]), int(config["association"]["random_seed"])
    )
    annual_mean, annual_low, annual_high = _bootstrap_mean(
        annual_values, int(config["association"]["bootstrap_replicates"]), int(config["association"]["random_seed"]) + 1
    )
    eligible = association["peak_core_analysis_eligible"].map(_as_bool).eq(True)
    core_effective = association["peak_core_contained_by_actual_pet_effective_contour"].map(_as_bool)
    duration_counts = (
        {str(key): int(value) for key, value in tracks["duration_class"].value_counts().items()}
        if "duration_class" in tracks
        else {}
    )
    input_hashes: dict[str, Any] = {}
    for label, path in (
        ("population", population_input),
        ("event_diagnostics", diagnostics_input),
        ("ofes_delta_do_catalog", catalog_input),
    ):
        digest, inventory = _hash_explicit_path(path)
        input_hashes[label] = {"sha256": digest, "files": inventory}
    summary = {
        "run_id": manifest["run_id"],
        "manifest": str(manifest_file),
        "postprocess_config": {
            "path": str(postprocess_config_path),
            "sha256": sha256_file(postprocess_config_path),
        },
        "association_lock": {
            "path": str(association_lock_path),
            "sha256": association_lock_sha256,
        },
        "input_hashes": input_hashes,
        "population_input_file": str(population_file),
        "mccoy_input_file": str(mccoy_path_input),
        "candidate_events": int(len(population.loc[np.isclose(pd.to_numeric(population["threshold"], errors="coerce"), 50.0)])),
        "strict_events": int(len(association)),
        "strict_primary_analysis_eligible": int(eligible.sum()),
        "peak_core_effective_contained_any_polarity": int(core_effective.fillna(False).sum()),
        "peak_core_speed_contained_any_polarity": int(association["peak_core_contained_by_actual_pet_speed_contour"].map(_as_bool).eq(True).sum()),
        "effective_peak_footprint_overlap": {
            "n": int(association["effective_peak_footprint_overlap_fraction"].notna().sum()),
            "median": float(association["effective_peak_footprint_overlap_fraction"].median()),
        },
        "local_ring_null": {
            "core_minus_ring_mean": local_mean,
            "bootstrap_95ci": [local_low, local_high],
            "wilcoxon_greater": _wilcoxon_greater(local_values),
        },
        "annual_month_latitude_null": {
            "core_minus_annual_mean": annual_mean,
            "bootstrap_95ci": [annual_low, annual_high],
            "wilcoxon_greater": _wilcoxon_greater(annual_values),
            "occupancy_table": str(annual_path),
        },
        "rotation_29": {
            "n": int(len(rotation)),
            "pet_peak_category_counts": {str(key): int(value) for key, value in rotation["pet_peak_category"].value_counts().items()},
            "surface_ro_same_polarity_counts": {str(key): int(value) for key, value in rotation["surface_core_rotation_polarity_match"].value_counts(dropna=False).items()},
        },
        "mccoy": {
            "n": int(len(mccoy_frame)),
            "positive_counts": {
                column: int(mccoy_frame[column].map(_as_bool).eq(True).sum())
                for column in (
                    "center_profile_mccoy_compatible",
                    "center_profile_velocity_confirmed",
                    "any_event_profile_mccoy_compatible",
                    "any_event_profile_velocity_confirmed",
                )
            },
        },
        "quality_eligible_161": {
            "n": int(len(sensitivity_result)),
            "effective_contained": int(sensitivity_result["pet_effective_contained"].sum()),
            "speed_contained": int(sensitivity_result["pet_speed_contained"].sum()),
            "path": str(sensitivity_path),
        },
        "daily_objects": {
            "total": int(len(objects)),
            "by_polarity": {str(key): int(value) for key, value in objects["polarity"].value_counts().items()},
            "effective_radius_km": _distribution(objects["effective_radius_km"]),
            "speed_radius_km": _distribution(objects["speed_radius_km"]),
            "amplitude_m": _distribution(objects["amplitude_m"]),
            "boundary_censored_fraction": float(objects["boundary_censored"].mean()),
            "filter_valid_fraction": float(objects["filter_valid"].mean()),
        },
        "tracks": {
            "count": int(len(tracks)),
            "duration_class_counts": duration_counts,
            "virtual_observation_fraction": float(observations["is_virtual"].mean()),
        },
        "filter_valid_domain_by_day": [
            {
                "date": date_text,
                "filter_valid_cells": int(info.get("audit", {}).get("filter_valid_cells", 0)),
                "filter_valid_fraction": float(info.get("audit", {}).get("filter_valid_fraction", np.nan)),
            }
            for date_text, info in manifest["days"].items()
        ],
        "paths": {
            "association": str(association_path),
            "nulls": str(null_path),
            "surface_ro_crosstab": str(surface_ro_path),
            "mccoy_crosstab": str(mccoy_path),
            "annual_occupancy": str(annual_path),
            "quality_eligible_161": str(sensitivity_path),
            "plots": plot_paths,
        },
    }
    summary_path = output_root / "surface_eddy_summary.json"
    _write_json(summary_path, summary)
    manifest["event_association"] = {
        "status": "complete",
        "strict_events": int(len(association)),
        "primary_analysis_eligible": int(eligible.sum()),
        "postprocess_config_path": str(postprocess_config_path),
        "postprocess_config_sha256": sha256_file(postprocess_config_path),
        "association_lock_path": str(association_lock_path),
        "association_lock_sha256": association_lock_sha256,
        "input_hashes": input_hashes,
        "paths": {key: str(value) for key, value in summary["paths"].items()},
    }
    manifest["updated_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    _write_json(manifest_file, manifest)
    return {key: str(value) for key, value in {"association": association_path, "nulls": null_path, "surface_ro_crosstab": surface_ro_path, "mccoy_crosstab": mccoy_path, "summary": summary_path}.items()}


def _plot_contour(ax: Any, row: Mapping[str, Any], name: str, **kwargs: Any) -> None:
    lon = _as_float_array(row.get(f"{name}_contour_lon"))
    lat = _as_float_array(row.get(f"{name}_contour_lat"))
    if lon.size >= 3 and lon.size == lat.size:
        if lon[0] != lon[-1] or lat[0] != lat[-1]:
            lon, lat = np.r_[lon, lon[0]], np.r_[lat, lat[0]]
        ax.plot(lon, lat, **kwargs)


def _write_plot_outputs(
    output_root: FilePath,
    config: Mapping[str, Any],
    strict: pd.DataFrame,
    association: pd.DataFrame,
    objects: pd.DataFrame,
    observations: pd.DataFrame,
    tracks: pd.DataFrame,
    domain: Mapping[str, Any],
    catalog_root: FilePath,
) -> dict[str, str]:
    """Write the predeclared diagnostic plots, including all three cases."""

    import matplotlib.pyplot as plt

    plot_root = output_root / "plots"
    plot_root.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    region = config["input"]["region"]
    xlim = (float(region["lon_min"]), float(region["lon_max"]))
    ylim = (float(region["lat_min"]), float(region["lat_max"]))

    for date_text in config["validation"]["fixed_dates"]:
        date = pd.Timestamp(date_text).normalize()
        rows = objects.loc[objects["date"].eq(date)].copy()
        fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=True, sharey=True)
        for axis, polarity in zip(axes, ("anticyclonic", "cyclonic")):
            subset = rows.loc[rows["polarity"].eq(polarity)]
            for row in subset.to_dict("records"):
                _plot_contour(
                    axis,
                    row,
                    "effective",
                    color="#c2410c" if polarity == "anticyclonic" else "#2563eb",
                    lw=1.0,
                    alpha=0.9 if bool(row.get("filter_valid")) else 0.35,
                )
                _plot_contour(
                    axis,
                    row,
                    "speed",
                    color="black",
                    lw=0.7,
                    ls="--",
                    alpha=0.8 if bool(row.get("filter_valid")) else 0.3,
                )
                axis.plot(row["center_lon"], row["center_lat"], "o", ms=2, color="black")
            axis.set_title(f"{date_text} {polarity}; n={len(subset)}")
            axis.set_xlim(*xlim)
            axis.set_ylim(*ylim)
            axis.set_xlabel("longitude (deg E)")
            axis.grid(alpha=0.2)
        axes[0].set_ylabel("latitude (deg N)")
        fig.suptitle("PET model-SSH effective (solid) / speed (dashed) contours")
        fig.tight_layout()
        path = plot_root / f"pet_contours_{date:%Y%m%d}.png"
        fig.savefig(path, dpi=160, bbox_inches="tight")
        plt.close(fig)
        paths[f"pet_contours_{date:%Y%m%d}"] = str(path)

    fig, axis = plt.subplots(figsize=(9, 6))
    if not tracks.empty:
        for track_id, group in observations.groupby("track_id", sort=True):
            group = group.sort_values("date")
            axis.plot(group["center_lon"], group["center_lat"], lw=0.45, alpha=0.35)
    censored = observations.loc[observations["boundary_censored"].astype(bool)]
    if not censored.empty:
        axis.scatter(censored["center_lon"], censored["center_lat"], s=4, c="#dc2626", label="boundary censored")
    axis.set_xlim(*xlim)
    axis.set_ylim(*ylim)
    axis.set_xlabel("longitude (deg E)")
    axis.set_ylabel("latitude (deg N)")
    axis.set_title(f"PET tracks; n={len(tracks)}")
    axis.grid(alpha=0.2)
    if not censored.empty:
        axis.legend()
    path = plot_root / "pet_tracks_boundary_censor.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    paths["pet_tracks_boundary_censor"] = str(path)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(
        association["ring_occupancy_fraction"],
        association["peak_core_contained_by_actual_pet_effective_contour"].astype(float),
        c=association["peak_lat"],
        cmap="viridis",
        s=26,
    )
    axes[0].set_xlabel("120--240 km ring occupancy")
    axes[0].set_ylabel("actual PET effective core indicator")
    axes[0].set_title("56-event peak-core containment vs local ring null")
    axes[0].grid(alpha=0.2)
    comparison = association[
        [
            "peak_core_contained_by_actual_pet_effective_contour",
            "peak_core_contained_by_actual_pet_speed_contour",
        ]
    ].astype(float).mean()
    axes[1].bar(["effective", "speed"], comparison.to_numpy(), color=["#0f766e", "#7c3aed"])
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("fraction of strict events")
    axes[1].set_title("Effective vs speed contour")
    fig.tight_layout()
    path = plot_root / "event_containment_and_contour_comparison.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    paths["event_containment_and_contour_comparison"] = str(path)

    rotation = association.loc[
        association["pet_peak_category"].isin(PRIMARY_PET_CATEGORIES),
        "pet_peak_category",
    ].value_counts()
    fig, axis = plt.subplots(figsize=(6, 4))
    axis.bar(rotation.index.astype(str), rotation.to_numpy(), color="#0f766e")
    axis.set_ylabel("events")
    axis.set_title("Rotation-event PET effective expression")
    axis.grid(axis="y", alpha=0.2)
    path = plot_root / "rotation_pet_expression.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    paths["rotation_pet_expression"] = str(path)

    case_ids = tuple(
        config["association"].get(
            "representative_event_ids",
            ("OFES_DO50_E000002", "OFES_DO50_E000239", "OFES_DO50_E000176"),
        )
    )
    strict_by_id = strict.set_index("event_id")
    for event_id in case_ids:
        if event_id not in strict_by_id.index:
            raise ValueError(f"Required representative event is absent: {event_id}")
        population_row = strict_by_id.loc[event_id]
        date = pd.Timestamp(population_row["peak_date"]).normalize()
        rows = _valid_surface_rows(objects, date, virtual=False)
        peak_lon, peak_lat = float(population_row["peak_lon"]), float(population_row["peak_lat"])
        event_result = association.loc[association["event_id"].astype(str).eq(event_id)]
        if event_result.empty:
            raise ValueError(f"Association result is missing required case: {event_id}")
        event_result = event_result.iloc[0]
        fig, axis = plt.subplots(figsize=(8, 6))
        pixels = _peak_pixels(catalog_root, population_row)
        pixel_lat = pixels["lat_index"].to_numpy(dtype=int)
        pixel_lon = pixels["lon_index"].to_numpy(dtype=int)
        valid_pixels = (
            (pixel_lat >= 0)
            & (pixel_lat < len(domain["lat"]))
            & (pixel_lon >= 0)
            & (pixel_lon < len(domain["lon"]))
        )
        axis.scatter(
            domain["lon"][pixel_lon[valid_pixels]],
            domain["lat"][pixel_lat[valid_pixels]],
            s=6,
            color="#a16207",
            alpha=0.45,
            label="DO peak footprint",
        )
        matches = _point_matches(rows, peak_lon, peak_lat, "effective")
        for row in rows.to_dict("records"):
            _plot_contour(axis, row, "effective", color="#0f766e" if str(row["object_id"]) in set(matches["object_id"].astype(str)) else "#94a3b8", lw=1.4)
            _plot_contour(axis, row, "speed", color="#7c3aed", lw=0.9, ls="--")
            axis.plot(row["center_lon"], row["center_lat"], "o", ms=3, color="#111827")
        axis.plot(peak_lon, peak_lat, "*", ms=11, color="#dc2626", label="DO peak core")
        if matches.empty and not rows.empty:
            nearest = _nearest_row(rows.to_dict("records"), peak_lon, peak_lat)
            if nearest is not None:
                axis.plot(nearest["center_lon"], nearest["center_lat"], "o", ms=8, mfc="none", mec="#dc2626", label="nearest PET center")
        pet_polarity = event_result.get("primary_polarity")
        normalized_distance = event_result.get("nearest_pet_center_distance_over_effective_radius")
        effective_flag = _as_bool(
            event_result.get("peak_core_contained_by_actual_pet_effective_contour")
        )
        speed_flag = _as_bool(
            event_result.get("peak_core_contained_by_actual_pet_speed_contour")
        )
        effective_text = str(effective_flag) if isinstance(effective_flag, bool) else "NA"
        speed_text = str(speed_flag) if isinstance(speed_flag, bool) else "NA"
        distance_text = (
            f"{float(normalized_distance):.2f}"
            if np.isfinite(_as_float(normalized_distance))
            else "unavailable"
        )
        axis.text(
            0.02,
            0.98,
            f"PET effective: {effective_text}\n"
            f"PET speed: {speed_text}\n"
            f"polarity: {pet_polarity or 'none'}\n"
            f"nearest d/R: {distance_text}",
            transform=axis.transAxes,
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
        axis.set_xlim(max(xlim[0], peak_lon - 4), min(xlim[1], peak_lon + 4))
        axis.set_ylim(max(ylim[0], peak_lat - 4), min(ylim[1], peak_lat + 4))
        axis.set_xlabel("longitude (deg E)")
        axis.set_ylabel("latitude (deg N)")
        axis.set_title(f"{event_id} | {date:%Y-%m-%d} | PET effective/speed")
        axis.grid(alpha=0.2)
        axis.legend(loc="best", fontsize=8)
        fig.tight_layout()
        path = plot_root / f"case_{event_id}.png"
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        paths[f"case_{event_id}"] = str(path)
    return paths


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--config", default="config/ofes_surface_eddy.yml")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--population", required=True, help="Exact population table or an unambiguous population directory")
    parser.add_argument("--diagnostics", required=True, help="Exact event-diagnostics table or an unambiguous directory")
    parser.add_argument("--catalog", required=True, help="Exact DO50 catalog table or an unambiguous directory")
    parser.add_argument("--mccoy-summary", help="Exact McCoy event_summary.parquet when diagnostics contains multiple summaries")
    parser.add_argument("--log-level", default="INFO")
    parser.add_argument("--log-file")
    parser.add_argument("command", choices=("associate",))
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if args.log_file:
        log_path = FilePath(args.log_file).resolve()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path, encoding="utf-8"))
    logging.basicConfig(
        level=getattr(logging, str(args.log_level).upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        handlers=handlers,
    )
    result = run_authoritative_association(
        args.repo_root,
        args.config,
        args.manifest,
        args.population,
        args.diagnostics,
        args.catalog,
        args.mccoy_summary,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
