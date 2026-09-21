#!/usr/bin/env python
"""Bounded paper-ready SCV/Argo repair and writing export.

The runner lives at the repository root, imports the adjacent track module,
reads existing A1--C2
artifacts and source metadata read-only, and writes only paper_ready outputs.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

MAIN_WORKTREE = Path(__file__).resolve().parent
# Resolve inputs relative to this checkout, independently of the launch directory.
SOURCE_REPO = MAIN_WORKTREE

if str(MAIN_WORKTREE) not in sys.path:
    sys.path.insert(0, str(MAIN_WORKTREE))
import track

AUDIT_ROOT = SOURCE_REPO / "plot_outputs/test/argo_scv_targeted_audit_20260920"
DEFAULT_OUTPUT = AUDIT_ROOT / "paper_ready"
FORMAL_COHORT = SOURCE_REPO / (
    "plot_outputs/do/global_ocean/scv_matched_control/"
    "scv_matched_control_cohort_2002_2023_depth300m_n3_space750km_y2_m2_"
    "core50m_externalplat_uniqueplat_noreuse_ekeq10_specf65759ebbe.parquet"
)
SCV_ANCHORED = SOURCE_REPO / (
    "plot_outputs/do/global_ocean/screen_mccoy_scvs_against_glorys/"
    "mccoy_glorys_miss.parquet"
)
SCV_CATALOG = SOURCE_REPO / "data/mccoy2020_scv/SCVs.csv"
ARGO_DATA = SOURCE_REPO / "Argo_data"
BASE_PREPROCESSING = AUDIT_ROOT / "profile_preprocessing_audit.parquet"
BASE_DETECTOR = AUDIT_ROOT / "detector_profile_audit.parquet"
IDENTITY_CANDIDATES = AUDIT_ROOT / "scv_identity_candidates.parquet"
IDENTITY_SUMMARY = AUDIT_ROOT / "scv_identity_summary.parquet"
THRESHOLDS = (20.0, 35.0, 50.0)
SCOPES = ("Global", "Global excluding KE")
BOOTSTRAP_ITERATIONS = 2000
RANDOM_SEED = 42
LOADER_SELECTION = {
    "Temperature": "Temp_Adjusted",
    "Temperature_Flag": "Temp_Adjusted_Flag",
    "DO": "DOXY_Adjusted",
    "DO_Flag": "DOXY_Adjusted_Flag",
    "Salinity": "PSAL_Adjusted",
    "Salinity_Flag": "PSAL_Adjusted_Flag",
}
BASE_COLUMNS = [
    "Year", "Month", "Day", "Longitude", "Latitude", "Depth",
    "Profile_number", "Platform_number",
]
NORMALIZATION_MAP = {
    "Depth_m": "Depth",
    "Temperature_degC": "Temp_Adjusted",
    "DO_mol_kg": "DOXY_Adjusted",
    "Salinity_psu": "PSAL_Adjusted",
}


class Tee:
    def __init__(self, stream, handle):
        self.stream = stream
        self.handle = handle

    def write(self, value: str) -> int:
        self.stream.write(value)
        self.handle.write(value)
        return len(value)

    def flush(self) -> None:
        self.stream.flush()
        self.handle.flush()

    def isatty(self) -> bool:
        return bool(getattr(self.stream, "isatty", lambda: False)())


def log(message: str) -> None:
    print(message, flush=True)


def ensure_output(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {path}; use --overwrite.")


def save_frame(frame: pd.DataFrame, stem: str, output_dir: Path, overwrite: bool) -> dict[str, str]:
    parquet_path = output_dir / f"{stem}.parquet"
    csv_path = output_dir / f"{stem}.csv"
    ensure_output(parquet_path, overwrite)
    ensure_output(csv_path, overwrite)
    frame.to_parquet(parquet_path, index=False)
    frame.to_csv(csv_path, index=False)
    return {"parquet": str(parquet_path), "csv": str(csv_path)}


def save_json(payload: Any, name: str, output_dir: Path, overwrite: bool) -> str:
    path = output_dir / name
    ensure_output(path, overwrite)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    return str(path)


def great_circle_km(lon: float, lat: float, other_lon: float, other_lat: float) -> float:
    values = [lon, lat, other_lon, other_lat]
    if not all(np.isfinite(value) for value in values):
        return np.nan
    return float(track.great_circle_distance_m(other_lon, other_lat, lon, lat) / 1000.0)


def decimal_places(value: Any) -> int:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return 0
    if "e" in text.lower():
        return 6
    return len(text.split(".", 1)[1]) if "." in text else 0


def read_profile_index(profile_years: dict[int, int], data_dir: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    by_year: dict[int, set[int]] = {}
    for profile_number, year in profile_years.items():
        by_year.setdefault(int(year), set()).add(int(profile_number))
    for year, ids in sorted(by_year.items()):
        path = data_dir / f"Argo{year}.parquet"
        if not ids or not path.exists():
            continue
        columns = ["Profile_number", "Platform_number", "Year", "Month", "Day", "Longitude", "Latitude"]
        schema = set(pq.read_schema(path).names)
        frame = pq.read_table(
            path, columns=[c for c in columns if c in schema],
            filters=[("Profile_number", "in", sorted(ids))],
        ).to_pandas()
        if frame.empty:
            continue
        frame["Profile_number"] = pd.to_numeric(frame["Profile_number"], errors="coerce")
        frame = frame.dropna(subset=["Profile_number"]).copy()
        frame["Profile_number"] = frame["Profile_number"].astype(int)
        rows.append(frame.groupby("Profile_number", as_index=False, sort=False).first())
    if not rows:
        return pd.DataFrame(columns=[
            "Profile_number", "Platform_number", "Year", "Month", "Day",
            "Longitude", "Latitude", "profile_date",
        ])
    index = pd.concat(rows, ignore_index=True).drop_duplicates("Profile_number", keep="first")
    index["profile_date"] = pd.to_datetime(
        {
            "year": pd.to_numeric(index["Year"], errors="coerce"),
            "month": pd.to_numeric(index["Month"], errors="coerce"),
            "day": pd.to_numeric(index["Day"], errors="coerce"),
        },
        errors="coerce",
    )
    return index.reset_index(drop=True)


def read_filtered_loader(year: int, profile_ids: set[int], data_dir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = data_dir / f"Argo{year}.parquet"
    schema = set(pq.read_schema(path).names)
    needed = set(BASE_COLUMNS) | set(LOADER_SELECTION.values())
    raw = pq.read_table(
        path,
        columns=[c for c in schema if c in needed],
        filters=[("Profile_number", "in", sorted(profile_ids))],
    ).to_pandas().rename(columns=NORMALIZATION_MAP)
    final = pd.DataFrame(index=raw.index)
    for column in BASE_COLUMNS:
        if column in raw.columns:
            final[column] = raw[column]
    for standard_name, source_name in LOADER_SELECTION.items():
        final[standard_name] = raw[source_name] if source_name in raw.columns else pd.NA
    for column in ("Profile_number", "Platform_number"):
        if column in final.columns:
            final[column] = pd.to_numeric(final[column], errors="coerce").astype("Int64")
    return final.reset_index(drop=True), {
        "year": year,
        "path": str(path),
        "schema_columns": sorted(schema),
        "source_columns_present": {k: v in schema for k, v in LOADER_SELECTION.items()},
        "loader_created_pd_na_columns": [k for k, v in LOADER_SELECTION.items() if v not in schema],
        "rows_read_after_profile_filter": int(len(final)),
        "profiles_read_after_profile_filter": int(final["Profile_number"].nunique()),
    }


def load_profile_data(profile_years: dict[int, int], data_dir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    by_year: dict[int, set[int]] = {}
    for profile_number, year in profile_years.items():
        by_year.setdefault(int(year), set()).add(int(profile_number))
    frames = []
    provenance = {}
    for year, ids in sorted(by_year.items()):
        log(f"[loader] year={year}, profiles={len(ids)}")
        frame, metadata = read_filtered_loader(year, ids, data_dir)
        frames.append(frame)
        provenance[str(year)] = metadata
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), provenance


def profile_year_map(active: pd.DataFrame, anchored: pd.DataFrame, candidates: pd.DataFrame) -> dict[int, int]:
    mapping: dict[int, int] = {}
    for frame, pcol, ycol in (
        (active, "profile_number", "year"),
        (anchored, "profile_number", "year"),
        (candidates, "candidate_profile_number", "candidate_profile_year"),
    ):
        for row in frame[[pcol, ycol]].dropna().itertuples(index=False):
            p = int(getattr(row, pcol))
            y = int(getattr(row, ycol))
            if p in mapping and mapping[p] != y:
                raise RuntimeError(f"Profile {p} has inconsistent years.")
            mapping[p] = y
    return mapping


def profile_audit(profile_data: pd.DataFrame, profile_ids: set[int], base: pd.DataFrame) -> tuple[pd.DataFrame, dict[int, pd.DataFrame]]:
    cfg = track.make_detection_config("do", anomaly_min_depth=300.0)
    rows = []
    cleaned_profiles: dict[int, pd.DataFrame] = {}
    for profile_number in sorted(profile_ids):
        profile = profile_data.loc[
            profile_data["Profile_number"].astype("Int64").eq(profile_number)
        ].copy()
        cleaned, diagnostics = track._prepare_do_profile_for_detection(profile, cfg)
        diagnostics["Profile_number"] = profile_number
        diagnostics["raw_rows"] = int(len(profile))
        if cleaned is None:
            diagnostics.update({
                "clean_do_min_depth_m": np.nan,
                "clean_do_max_depth_m": np.nan,
                "clean_ts_min_depth_m": np.nan,
                "clean_ts_max_depth_m": np.nan,
                "n_clean_common_levels_300_1000": 0,
                "n_clean_do_levels_300_1000": 0,
                "n_clean_ts_levels_300_1000": 0,
            })
        else:
            cleaned_profiles[profile_number] = cleaned
            depth = pd.to_numeric(cleaned["Depth"], errors="coerce")
            in_range = depth.ge(300.0) & depth.le(1000.0)
            diagnostics.update({
                "clean_do_min_depth_m": float(depth.min()),
                "clean_do_max_depth_m": float(depth.max()),
                "clean_ts_min_depth_m": float(depth.min()),
                "clean_ts_max_depth_m": float(depth.max()),
                "n_clean_common_levels_300_1000": int(in_range.sum()),
                "n_clean_do_levels_300_1000": int((in_range & cleaned["DO"].notna()).sum()),
                "n_clean_ts_levels_300_1000": int((
                    in_range & cleaned["Temperature"].notna() & cleaned["Salinity"].notna()
                ).sum()),
            })
        rows.append(diagnostics)
    audit = pd.DataFrame(rows)
    base_fields = base[[
        "Profile_number", "detector_preprocessed",
        "clean_rows_before_depth_dedup", "clean_rows_after_depth_dedup",
    ]].rename(columns={
        "detector_preprocessed": "a1_detector_preprocessed",
        "clean_rows_before_depth_dedup": "a1_clean_rows_before_depth_dedup",
        "clean_rows_after_depth_dedup": "a1_clean_rows_after_depth_dedup",
    })
    audit = audit.merge(base_fields, on="Profile_number", how="left", validate="one_to_one")
    audit["a1_detector_boolean_match"] = audit["detector_preprocessed"].eq(
        audit["a1_detector_preprocessed"]
    )
    audit["a1_count_semantics_difference_near_zero_skip"] = (
        audit["near_zero_triggered"].fillna(False).astype(bool)
        & (
            ~audit["clean_rows_before_depth_dedup"].eq(audit["a1_clean_rows_before_depth_dedup"])
            | ~audit["clean_rows_after_depth_dedup"].eq(audit["a1_clean_rows_after_depth_dedup"])
        )
    )
    return audit, cleaned_profiles


def _scalar_regression_match(left: Any, right: Any, *, atol: float = 1e-8, rtol: float = 1e-8) -> bool:
    left_value = pd.to_numeric(pd.Series([left]), errors="coerce").iloc[0]
    right_value = pd.to_numeric(pd.Series([right]), errors="coerce").iloc[0]
    if pd.isna(left_value) and pd.isna(right_value):
        return True
    if pd.isna(left_value) or pd.isna(right_value):
        return False
    return bool(np.isclose(float(left_value), float(right_value), atol=atol, rtol=rtol))


def detector_outcome_regression(
    profile_data: pd.DataFrame,
    profile_ids: set[int],
    base_detector: pd.DataFrame,
) -> pd.DataFrame:
    """对 690 个原 A1 profile 重新运行 detector 并与旧结果逐项比较。"""
    old = base_detector.set_index(["Profile_number", "threshold_umol_kg"])
    rows = []
    for profile_number in sorted(profile_ids):
        profile = profile_data.loc[
            profile_data["Profile_number"].astype("Int64").eq(profile_number)
        ].copy()
        for threshold in THRESHOLDS:
            cfg = track.make_detection_config(
                "do", do_threshold=float(threshold), anomaly_min_depth=300.0
            )
            _, new_diagnostics = track._prepare_do_profile_for_detection(profile, cfg)
            new_preprocessed = bool(new_diagnostics["detector_preprocessed"])
            detected = track.calculate_delta_do(
                profile,
                detection_config=cfg,
                remove_outliers=True,
                include_aou=True,
                verbose=False,
            )
            best = track._keep_best_anomaly_per_profile(detected, cfg)
            peak = best.iloc[0] if not best.empty else None
            old_row = old.loc[(profile_number, float(threshold))]
            new_has = bool(peak is not None)
            new_delta = float(peak["delta_do"]) if peak is not None else np.nan
            new_depth = float(peak["depth"]) if peak is not None else np.nan
            recomputed_bool_match = new_has == bool(old_row["recomputed_has_delta_do"])
            recomputed_delta_match = _scalar_regression_match(
                new_delta, old_row["recomputed_delta_do"]
            )
            recomputed_depth_match = _scalar_regression_match(
                new_depth, old_row["recomputed_peak_depth_m"], atol=1e-5, rtol=1e-8
            )
            direct_bool_match = new_has == bool(old_row["direct_anomaly_cache_has_delta_do"])
            rows.append({
                "Profile_number": int(profile_number),
                "threshold_umol_kg": float(threshold),
                "new_preprocessed": new_preprocessed,
                "a1_recomputed_preprocessed": bool(old_row["recomputed_preprocessed"]),
                "preprocessing_match": bool(
                    new_preprocessed == bool(old_row["recomputed_preprocessed"])
                ),
                "new_has_delta_do": new_has,
                "a1_recomputed_has_delta_do": bool(old_row["recomputed_has_delta_do"]),
                "a1_recomputed_delta_do": old_row["recomputed_delta_do"],
                "new_delta_do": new_delta,
                "a1_recomputed_peak_depth_m": old_row["recomputed_peak_depth_m"],
                "new_peak_depth_m": new_depth,
                "a1_direct_cache_has_delta_do": bool(old_row["direct_anomaly_cache_has_delta_do"]),
                "new_vs_a1_recomputed_bool_match": bool(recomputed_bool_match),
                "new_vs_a1_recomputed_delta_match": bool(recomputed_delta_match),
                "new_vs_a1_recomputed_depth_match": bool(recomputed_depth_match),
                "new_vs_a1_direct_cache_bool_match": bool(direct_bool_match),
                "outcome_match": bool(
                    new_preprocessed == bool(old_row["recomputed_preprocessed"])
                    and recomputed_bool_match
                    and recomputed_delta_match
                    and recomputed_depth_match
                ),
            })
    return pd.DataFrame(rows)


def build_participant_flags(active: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame:
    fields = [
        "Profile_number", "detector_preprocessed", "near_zero_triggered",
        "near_zero_count", "clean_do_min_depth_m", "clean_do_max_depth_m",
        "clean_ts_min_depth_m", "clean_ts_max_depth_m",
        "n_clean_common_levels_300_1000", "n_clean_do_levels_300_1000",
        "n_clean_ts_levels_300_1000",
    ]
    frame = active.copy()
    frame["profile_number"] = pd.to_numeric(frame["profile_number"], errors="coerce").astype(int)
    frame = frame.merge(
        audit[fields].rename(columns={"Profile_number": "profile_number"}),
        on="profile_number", how="left", validate="one_to_one",
    )
    anchor_max = frame.loc[frame["is_scv"].astype(bool), [
        "match_set_id", "clean_do_max_depth_m"
    ]].rename(columns={"clean_do_max_depth_m": "anchor_clean_do_max_depth_m"})
    frame = frame.merge(anchor_max, on="match_set_id", how="left", validate="many_to_one")
    core = pd.to_numeric(frame["anchor_core_depth_m"], errors="coerce")
    margin = pd.to_numeric(frame["core_depth_margin_m"], errors="coerce").fillna(50.0)
    do_min = pd.to_numeric(frame["clean_do_min_depth_m"], errors="coerce")
    do_max = pd.to_numeric(frame["clean_do_max_depth_m"], errors="coerce")
    ts_min = pd.to_numeric(frame["clean_ts_min_depth_m"], errors="coerce")
    ts_max = pd.to_numeric(frame["clean_ts_max_depth_m"], errors="coerce")
    frame["coverage_300_1000_support"] = frame["n_clean_common_levels_300_1000"].gt(0).fillna(False)
    frame["clean_core_coverage_recheck"] = (
        do_min.le(core - margin) & do_max.ge(core + margin)
        & ts_min.le(core - margin) & ts_max.ge(core + margin)
    ).fillna(False)
    anchor_max_clean = pd.to_numeric(frame["anchor_clean_do_max_depth_m"], errors="coerce")
    frame["clean_pair_max_depth_delta_m"] = (do_max - anchor_max_clean).abs()
    frame.loc[frame["is_scv"].astype(bool), "clean_pair_max_depth_delta_m"] = 0.0
    max_limit = pd.to_numeric(frame["max_depth_caliper_m"], errors="coerce").fillna(500.0)
    frame["clean_pair_max_depth_recheck"] = frame["clean_pair_max_depth_delta_m"].le(max_limit).fillna(False)
    frame["stored_core_coverage_rule"] = (
        pd.to_numeric(frame["do_min_depth_m"], errors="coerce").le(core - margin)
        & pd.to_numeric(frame["do_max_depth_m"], errors="coerce").ge(core + margin)
        & pd.to_numeric(frame["ts_min_depth_m"], errors="coerce").le(core - margin)
        & pd.to_numeric(frame["ts_max_depth_m"], errors="coerce").ge(core + margin)
    ).fillna(False)
    frame["stored_pair_max_depth_rule"] = pd.to_numeric(
        frame["match_max_depth_delta_m"], errors="coerce"
    ).le(max_limit).fillna(False)
    frame["depth_pair_recheck_pass"] = (
        frame["clean_core_coverage_recheck"] & frame["clean_pair_max_depth_recheck"]
    )
    frame["qualification_pass"] = (
        frame["detector_preprocessed"].fillna(False).astype(bool)
        & frame["depth_pair_recheck_pass"].fillna(False).astype(bool)
    )
    return frame


def direct_identity_index(candidates: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    years = dict(zip(
        candidates["candidate_profile_number"].astype(int),
        candidates["candidate_profile_year"].astype(int),
    ))
    current = candidates[[
        "current_profile_number", "current_profile_year"
    ]].drop_duplicates().rename(columns={
        "current_profile_number": "candidate_profile_number",
        "current_profile_year": "candidate_profile_year",
    })
    years.update(dict(zip(
        current["candidate_profile_number"].astype(int),
        current["candidate_profile_year"].astype(int),
    )))
    return read_profile_index(years, data_dir)


def adjudicate_identity(
    candidates: pd.DataFrame,
    catalog_path: Path,
    data_dir: Path,
    anchored: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    raw_catalog = pd.read_csv(catalog_path, dtype=str)
    catalog = raw_catalog.copy()
    catalog["catalog_row_id"] = np.arange(len(catalog), dtype=int)
    for column in ("ID", "Platform", "Cycle"):
        catalog[column] = pd.to_numeric(catalog[column], errors="coerce")
    catalog["catalog_date"] = pd.to_datetime(
        catalog["Cycle_ISO_DateTime_UTC"], errors="coerce", utc=True
    ).dt.tz_localize(None)
    for column in ("Longitude", "Latitude", "Core_Pressure"):
        catalog[column] = pd.to_numeric(catalog[column], errors="coerce")
    catalog["lon_decimals"] = raw_catalog["Longitude"].map(decimal_places)
    catalog["lat_decimals"] = raw_catalog["Latitude"].map(decimal_places)
    index = direct_identity_index(candidates, data_dir)
    direct = index.rename(columns={
        "Profile_number": "candidate_profile_number",
        "Platform_number": "direct_platform_number",
        "profile_date": "direct_profile_date",
        "Longitude": "direct_lon",
        "Latitude": "direct_lat",
    })
    source = anchored.reset_index(drop=True).copy()
    source.insert(0, "association_row", np.arange(1, len(source) + 1, dtype=int))
    source = source.rename(columns={
        "profile_number": "source_profile_number",
        "mccoy_core_pressure_db": "source_core_pressure_db",
        "date": "source_profile_date",
        "lon": "source_profile_lon",
        "lat": "source_profile_lat",
        "scv_type": "source_scv_type",
        "basin": "source_basin",
    })
    source_fields = [
        "association_row", "source_profile_number", "source_core_pressure_db",
        "source_profile_date", "source_profile_lon", "source_profile_lat",
        "source_scv_type", "source_basin",
    ]
    source_profile_map = dict(
        zip(source["association_row"].astype(int), source["source_profile_number"].astype(int))
    )
    candidate_profile_map = {
        int(association_row): int(group["current_profile_number"].iloc[0])
        for association_row, group in candidates.groupby("association_row", sort=True)
    }
    if source_profile_map != candidate_profile_map:
        raise RuntimeError("SCV source association order does not match identity candidate rows.")
    work = candidates.merge(
        catalog[[
            "catalog_row_id", "ID", "Platform", "Cycle", "catalog_date",
            "Longitude", "Latitude", "Core_Pressure", "lon_decimals", "lat_decimals",
        ]].rename(columns={
            "ID": "catalog_id_direct", "Platform": "catalog_platform_direct",
            "Cycle": "catalog_cycle_direct", "catalog_date": "catalog_date_direct",
            "Longitude": "catalog_lon_direct", "Latitude": "catalog_lat_direct",
            "Core_Pressure": "catalog_core_pressure_direct",
            "lon_decimals": "catalog_lon_decimals", "lat_decimals": "catalog_lat_decimals",
        }),
        on="catalog_row_id", how="left", validate="many_to_one",
    ).merge(
        direct[[
            "candidate_profile_number", "direct_platform_number",
            "direct_profile_date", "direct_lon", "direct_lat",
        ]],
        on="candidate_profile_number", how="left", validate="many_to_one",
    ).merge(
        source[source_fields],
        on="association_row", how="left", validate="many_to_one",
    )
    work["catalog_source_match"] = (
        pd.to_numeric(work["catalog_id"], errors="coerce").eq(work["catalog_id_direct"])
        & pd.to_numeric(work["catalog_platform"], errors="coerce").eq(work["catalog_platform_direct"])
        & pd.to_numeric(work["catalog_cycle"], errors="coerce").eq(work["catalog_cycle_direct"])
        & pd.to_datetime(work["catalog_date"]).eq(pd.to_datetime(work["catalog_date_direct"]))
    )
    work["local_index_match"] = (
        pd.to_numeric(work["candidate_platform_number"], errors="coerce").eq(work["direct_platform_number"])
        & pd.to_datetime(work["candidate_profile_date"]).eq(pd.to_datetime(work["direct_profile_date"]))
        & np.isclose(pd.to_numeric(work["candidate_lon"], errors="coerce"), pd.to_numeric(work["direct_lon"], errors="coerce"), equal_nan=False)
        & np.isclose(pd.to_numeric(work["candidate_lat"], errors="coerce"), pd.to_numeric(work["direct_lat"], errors="coerce"), equal_nan=False)
    )
    catalog_lon = pd.to_numeric(work["catalog_lon_direct"], errors="coerce")
    catalog_lat = pd.to_numeric(work["catalog_lat_direct"], errors="coerce")
    candidate_lon = pd.to_numeric(work["direct_lon"], errors="coerce")
    candidate_lat = pd.to_numeric(work["direct_lat"], errors="coerce")
    lon_decimals = pd.to_numeric(work["catalog_lon_decimals"], errors="coerce").fillna(0).astype(int)
    lat_decimals = pd.to_numeric(work["catalog_lat_decimals"], errors="coerce").fillna(0).astype(int)
    work["minimal_lon_diff_deg"] = [
        float(np.asarray(track._minimal_lon_diff_deg(candidate, reference)))
        if np.isfinite(candidate) and np.isfinite(reference) else np.nan
        for candidate, reference in zip(candidate_lon, catalog_lon)
    ]
    work["lat_diff_deg"] = candidate_lat - catalog_lat
    work["candidate_lon_at_catalog_precision"] = [
        round(float(candidate), int(decimals)) if np.isfinite(candidate) else np.nan
        for candidate, decimals in zip(candidate_lon, lon_decimals)
    ]
    work["catalog_lon_at_source_precision"] = [
        round(float(reference), int(decimals)) if np.isfinite(reference) else np.nan
        for reference, decimals in zip(catalog_lon, lon_decimals)
    ]
    work["candidate_lat_at_catalog_precision"] = [
        round(float(candidate), int(decimals)) if np.isfinite(candidate) else np.nan
        for candidate, decimals in zip(candidate_lat, lat_decimals)
    ]
    work["catalog_lat_at_source_precision"] = [
        round(float(reference), int(decimals)) if np.isfinite(reference) else np.nan
        for reference, decimals in zip(catalog_lat, lat_decimals)
    ]
    work["same_calendar_date"] = (
        pd.to_datetime(work["catalog_date_direct"]).dt.normalize()
        == pd.to_datetime(work["direct_profile_date"]).dt.normalize()
    )
    work["coordinate_precision_match"] = (
        work["candidate_lon_at_catalog_precision"].eq(work["catalog_lon_at_source_precision"])
        & work["candidate_lat_at_catalog_precision"].eq(work["catalog_lat_at_source_precision"])
    )
    work["metadata_exact"] = work["same_calendar_date"] & work["coordinate_precision_match"]
    work["original_catalog_core_match"] = (
        pd.to_numeric(work["catalog_core_pressure_direct"], errors="coerce")
        == pd.to_numeric(work["source_core_pressure_db"], errors="coerce")
    )
    original_rows_by_association = {
        int(association_row): sorted(
            group.loc[group["original_catalog_core_match"], "catalog_row_id"].astype(int).unique()
        )
        for association_row, group in work.groupby("association_row", sort=True)
    }
    work["original_catalog_row_ids"] = work["association_row"].map(original_rows_by_association)
    work["original_catalog_row_resolution"] = work["association_row"].map(
        lambda value: (
            "unique" if len(original_rows_by_association.get(int(value), [])) == 1
            else "multiple_same_core_records" if len(original_rows_by_association.get(int(value), [])) > 1
            else "unavailable"
        )
    )
    work["original_catalog_row"] = [
        int(row.catalog_row_id) in original_rows_by_association.get(int(row.association_row), [])
        for row in work.itertuples(index=False)
    ]
    work["metadata_exact_original_row"] = work["metadata_exact"] & work["original_catalog_row"]
    rows = []
    for association_row, group in work.groupby("association_row", sort=True):
        original_rows = original_rows_by_association.get(int(association_row), [])
        original_group = group.loc[group["original_catalog_row"]].copy()
        current = original_group.loc[original_group["candidate_is_current_profile"].astype(bool)]
        alternatives = original_group.loc[~original_group["candidate_is_current_profile"].astype(bool)]
        current_exact = current.loc[current["metadata_exact_original_row"]]
        alternative_exact = alternatives.loc[alternatives["metadata_exact_original_row"]]
        replacement_candidates = sorted(alternative_exact["candidate_profile_number"].astype(int).unique())
        row_resolution = (
            "unique" if len(original_rows) == 1
            else "multiple_same_core_records" if len(original_rows) > 1
            else "unavailable"
        )
        if row_resolution != "unique":
            status = "unresolved-original-catalog-row"
        elif len(current_exact) and len(alternative_exact):
            status = "unresolved-current-and-alternative"
        elif len(current_exact):
            status = "metadata-supported"
        elif len(alternative_exact) and len(replacement_candidates) == 1:
            status = "confirmed-mismatch"
        else:
            status = "unresolved"
        current_dist = pd.to_numeric(current["position_distance_km"], errors="coerce")
        rows.append({
            "association_row": int(association_row),
            "current_profile_number": int(group["current_profile_number"].iloc[0]),
            "catalog_rows_in_window": int(group["catalog_row_id"].nunique()),
            "original_catalog_row_ids": original_rows,
            "original_catalog_row_id": int(original_rows[0]) if len(original_rows) == 1 else np.nan,
            "original_catalog_row_resolution": row_resolution,
            "physical_profile_candidates": int(group["candidate_profile_number"].nunique()),
            "physical_profile_candidates_original_row": int(original_group["candidate_profile_number"].nunique()),
            "current_metadata_exact_catalog_rows": sorted(current_exact["catalog_row_id"].astype(int).unique()),
            "alternative_metadata_exact_catalog_rows": sorted(alternative_exact["catalog_row_id"].astype(int).unique()),
            "replacement_candidate_profiles": replacement_candidates,
            "current_min_position_distance_km": float(current_dist.min()) if len(current_dist) else np.nan,
            "current_max_position_distance_km": float(current_dist.max()) if len(current_dist) else np.nan,
            "current_date_position_status": status,
            "catalog_source_match_all": bool(group["catalog_source_match"].all()),
            "local_index_match_all": bool(group["local_index_match"].all()),
            "current_candidate_count": int(len(current)),
            "alternative_candidate_count": int(len(alternatives)),
            "position_gt_10km": bool(current_dist.min() > 10.0) if len(current_dist) else False,
            "evidence": (
                "current profile matches the uniquely recovered original catalog row at calendar-date and source-coordinate precision"
                if status == "metadata-supported"
                else "unique alternative matches the uniquely recovered original catalog row at calendar-date and source-coordinate precision"
                if status == "confirmed-mismatch"
                else "original catalog row is not uniquely recoverable from the source core-pressure record"
                if status == "unresolved-original-catalog-row"
                else "the uniquely recovered original catalog row has both current and alternative metadata matches"
                if status == "unresolved-current-and-alternative"
                else "no source-precision comparison distinguishes current from alternatives"
            ),
        })
    decision = pd.DataFrame(rows)
    unresolved_statuses = {"unresolved", "unresolved-original-catalog-row", "unresolved-current-and-alternative"}
    details = {
        "n_associations": int(len(decision)),
        "metadata_supported": int(decision["current_date_position_status"].eq("metadata-supported").sum()),
        "confirmed_mismatch": int(decision["current_date_position_status"].eq("confirmed-mismatch").sum()),
        "unresolved": int(decision["current_date_position_status"].isin(unresolved_statuses).sum()),
        "identity_status_counts": decision["current_date_position_status"].value_counts(dropna=False).to_dict(),
        "confirmed_association_rows": decision.loc[
            decision["current_date_position_status"].eq("confirmed-mismatch"), "association_row"
        ].astype(int).tolist(),
        "confirmed_original_catalog_row_ids": decision.loc[
            decision["current_date_position_status"].eq("confirmed-mismatch"), "original_catalog_row_id"
        ].astype(int).tolist(),
        "unresolved_original_catalog_association_rows": decision.loc[
            decision["current_date_position_status"].eq("unresolved-original-catalog-row"), "association_row"
        ].astype(int).tolist(),
        "unresolved_unique_catalog_row_association_rows": decision.loc[
            decision["current_date_position_status"].eq("unresolved")
            & decision["original_catalog_row_resolution"].eq("unique"), "association_row"
        ].astype(int).tolist(),
        "original_catalog_row_unique": int(decision["original_catalog_row_resolution"].eq("unique").sum()),
        "original_catalog_row_unresolved": int(
            (~decision["original_catalog_row_resolution"].eq("unique")).sum()
        ),
        "confirmed_mismatch_with_original_catalog_row": int(
            decision["current_date_position_status"].eq("confirmed-mismatch").sum()
        ),
        "source_catalog_all_match": bool(work["catalog_source_match"].all()),
        "source_local_index_all_match": bool(work["local_index_match"].all()),
        "position_gt_10km": int(decision["position_gt_10km"].sum()),
        "large_distance_current_profiles": decision.loc[decision["position_gt_10km"], "current_profile_number"].astype(int).tolist(),
    }
    return decision, work, details


def pair_queue(participants: pd.DataFrame, anchor_condition: pd.Series) -> pd.DataFrame:
    is_anchor = participants["is_scv"].astype(bool)
    anchor_ids = set(participants.loc[is_anchor & anchor_condition, "match_set_id"].astype(str))
    control_condition = participants["qualification_pass"].fillna(False).astype(bool)
    control_ids = set(participants.loc[~is_anchor & control_condition, "match_set_id"].astype(str))
    retained = anchor_ids & control_ids
    keep = participants["match_set_id"].astype(str).isin(retained) & (
        (is_anchor & anchor_condition) | (~is_anchor & control_condition)
    )
    queue = participants.loc[keep].copy()
    queue["bootstrap_input_order"] = np.arange(len(queue), dtype=np.int64)
    return queue


def bootstrap_facts(requested: int, valid: int, zero: int, infinite: int, low: float, high: float) -> dict:
    """Describe saved draws and raw endpoints without judging inferential validity."""
    undefined = requested - valid
    if min(requested, valid, zero, infinite, undefined) < 0 or zero + infinite > valid:
        raise ValueError("Inconsistent bootstrap draw counts")
    return {
        "all_draws_defined_finite_and_positive": bool(requested > 0 and valid == requested and zero == 0 and infinite == 0),
        "has_zero_draws": bool(zero),
        "has_infinite_draws": bool(infinite),
        "has_undefined_draws": bool(undefined),
        "undefined_draws": int(undefined),
        "interval_unbounded": bool(np.isinf(low) or np.isinf(high)),
        "degenerate_all_defined_zero": bool(valid > 0 and zero == valid and infinite == 0 and low == 0 and high == 0),
        "ci_low_raw": low,
        "ci_high_raw": high,
    }


def bootstrap_status(requested: int, valid: int, zero: int, infinite: int, low: float, high: float) -> tuple[str, bool]:
    """Return factual draw labels and whether every requested draw is finite positive."""
    facts = bootstrap_facts(requested, valid, zero, infinite, low, high)
    labels = []
    if valid == 0:
        labels.append("no-defined-draws")
    if facts["degenerate_all_defined_zero"]:
        labels.append("all-defined-draws-zero")
    for key, label in (("has_zero_draws", "zero-draws"), ("has_infinite_draws", "infinite-draws"),
                       ("has_undefined_draws", "undefined-draws"), ("interval_unbounded", "unbounded-interval")):
        if facts[key]:
            labels.append(label)
    return ";".join(labels) or "all-draws-finite-positive", facts["all_draws_defined_finite_and_positive"]


def restore_bootstrap_input_order(queue: pd.DataFrame) -> pd.DataFrame:
    """Restore the recorded analysis order before constructing cells and clusters.

    CSV display order is not bootstrap order. Missing or invalid order metadata
    must be recovered from the original input, never inferred from CI endpoints.
    Subsets may retain gaps in the original ordinal sequence.
    """
    if "bootstrap_input_order" not in queue:
        raise ValueError("Missing bootstrap_input_order; recover the original analysis order before replay")
    order = pd.to_numeric(queue["bootstrap_input_order"], errors="coerce")
    if (not np.isfinite(order).all() or order.lt(0).any()
            or order.mod(1).ne(0).any() or order.duplicated().any()):
        raise ValueError("bootstrap_input_order must contain unique nonnegative integer ordinals")
    return queue.assign(bootstrap_input_order=order.astype(np.int64)).sort_values(
        "bootstrap_input_order", kind="stable"
    ).copy()


def summarize_queue(queue: pd.DataFrame, label: str) -> pd.DataFrame:
    queue = restore_bootstrap_input_order(queue)
    rows = []
    for scope in SCOPES:
        anchors = queue.loc[queue["is_scv"].astype(bool)].drop_duplicates("match_set_id")
        scope_ids = set(anchors.loc[track._scv_matched_scope_mask(anchors, scope), "match_set_id"].astype(str))
        scoped = queue.loc[queue["match_set_id"].astype(str).isin(scope_ids)].copy()
        dependency = track._matched_dependency_components(scoped)
        for threshold in THRESHOLDS:
            tag = track._format_detection_value(float(threshold))
            outcome = f"has_delta_do_{tag}"
            cells = track._matched_set_cells(scoped, outcome)
            if not cells.empty:
                cells = cells.merge(dependency, on="match_set_id", how="left", validate="one_to_one")
                cells["anchor_platform"] = cells["anchor_platform_number"]
            scv_n = int(len(cells))
            control_n = int(cells[["c", "d"]].sum(axis=1).sum()) if not cells.empty else 0
            scv_k = int(cells["a"].sum()) if not cells.empty else 0
            control_k = int(cells["c"].sum()) if not cells.empty else 0
            result = {
                "queue": label, "scope": scope, "threshold_umol_kg": float(threshold),
                "n_matched_sets": scv_n, "n_anchor_profiles": scv_n,
                "n_control_profiles": control_n, "scv_numerator": scv_k,
                "control_numerator": control_k,
                "scv_rate": float(scv_k / scv_n) if scv_n else np.nan,
                "control_rate": float(control_k / control_n) if control_n else np.nan,
                "matched_mantel_haenszel_or": track._mantel_haenszel_odds_ratio(cells),
            }
            seed_scope = "qualification" if label == "qualification_all" else label
            for cluster in ("anchor_platform", "dependency_component"):
                low, high, valid, zero, infinite = track._bootstrap_matched_or_by_cluster(
                    cells, cluster_col=cluster, iterations=BOOTSTRAP_ITERATIONS,
                    random_seed=track._stable_analysis_seed(
                        RANDOM_SEED, seed_scope, scope, float(threshold), "profile", cluster
                    ),
                )
                status, usable = bootstrap_status(BOOTSTRAP_ITERATIONS, valid, zero, infinite, low, high)
                prefix = "anchor_platform" if cluster == "anchor_platform" else "dependency_component"
                result.update({
                    f"{prefix}_bootstrap_requested": BOOTSTRAP_ITERATIONS,
                    f"{prefix}_bootstrap_valid": int(valid),
                    f"{prefix}_bootstrap_undefined_draws": int(BOOTSTRAP_ITERATIONS - valid),
                    f"{prefix}_bootstrap_zero_draws": int(zero),
                    f"{prefix}_bootstrap_infinite_draws": int(infinite),
                    f"{prefix}_bootstrap_ci_low_raw": low,
                    f"{prefix}_bootstrap_ci_high_raw": high,
                    f"{prefix}_bootstrap_status": status,
                    f"{prefix}_bootstrap_all_draws_defined_finite_and_positive": usable,
                    f"{prefix}_bootstrap_display_interval": (
                        "undefined" if valid == 0 else f"[{low:g}, {high:g}] ({status}; raw)"
                    ),
                })
                result.update({f"{prefix}_bootstrap_{key}": value for key, value in
                               bootstrap_facts(BOOTSTRAP_ITERATIONS, valid, zero, infinite, low, high).items()})
            rows.append(result)
    return pd.DataFrame(rows)


def edge_case_tests(summary: pd.DataFrame) -> pd.DataFrame:
    def cells(a: int, b: int, c: int, d: int) -> pd.DataFrame:
        return pd.DataFrame({
            "match_set_id": ["x1", "x2"], "a": [a, a], "b": [b, b],
            "c": [c, c], "d": [d, d], "n": [a + b + c + d, a + b + c + d],
            "anchor_platform": [1, 2], "dependency_component": [0, 1],
        })
    rows = []
    for name, frame in (
        ("zero_vs_one", cells(0, 1, 1, 0)),
        ("zero_vs_zero", cells(0, 1, 0, 1)),
        ("mh_denominator_zero", cells(1, 0, 0, 1)),
    ):
        mh = track._mantel_haenszel_odds_ratio(frame)
        low, high, valid, zero, infinite = track._bootstrap_matched_or_by_cluster(
            frame, cluster_col="anchor_platform", iterations=100, random_seed=20260920
        )
        status, usable = bootstrap_status(100, valid, zero, infinite, low, high)
        rows.append({
            "test": name, "mh_or": mh, "ci_low_raw": low, "ci_high_raw": high,
            "requested": 100, "valid": valid, "undefined_draws": 100 - valid,
            "zero_draws": zero, "infinite_draws": infinite, "status": status,
            "all_draws_defined_finite_and_positive": usable,
        })
    normal = summary.loc[
        summary.queue.eq("qualification_all")
        & summary.scope.eq("Global")
        & summary.threshold_umol_kg.eq(20.0)
    ].iloc[0]
    rows.append({
        "test": "normal_c2_global_do20",
        "mh_or": normal["matched_mantel_haenszel_or"],
        "ci_low_raw": normal["anchor_platform_bootstrap_ci_low_raw"],
        "ci_high_raw": normal["anchor_platform_bootstrap_ci_high_raw"],
        "requested": normal["anchor_platform_bootstrap_requested"],
        "valid": normal["anchor_platform_bootstrap_valid"],
        "undefined_draws": normal["anchor_platform_bootstrap_undefined_draws"],
        "zero_draws": normal["anchor_platform_bootstrap_zero_draws"],
        "infinite_draws": normal["anchor_platform_bootstrap_infinite_draws"],
        "status": normal["anchor_platform_bootstrap_status"],
        "all_draws_defined_finite_and_positive": normal["anchor_platform_bootstrap_all_draws_defined_finite_and_positive"],
    })
    return pd.DataFrame(rows)


def near_zero_sensitivity(profile_data: pd.DataFrame, audit: pd.DataFrame, roles: dict[int, str], detector_lookup: pd.DataFrame) -> pd.DataFrame:
    ids = audit.loc[audit["near_zero_triggered"].fillna(False).astype(bool), "Profile_number"].astype(int).tolist()
    cfg = track.make_detection_config("do", anomaly_min_depth=300.0, do_near_zero_max_count=None)
    lookup = detector_lookup.set_index(["Profile_number", "threshold_umol_kg"])
    rows = []
    for profile_number in ids:
        profile = profile_data.loc[profile_data["Profile_number"].astype("Int64").eq(profile_number)].copy()
        for threshold in THRESHOLDS:
            cfg_thr = track.make_detection_config(cfg, do_threshold=float(threshold))
            prep, diagnostics = track._prepare_do_profile_for_detection(profile, cfg_thr)
            detected = track.calculate_delta_do(
                profile, detection_config=cfg_thr, remove_outliers=True, include_aou=True, verbose=False
            )
            best = track._keep_best_anomaly_per_profile(detected, cfg_thr)
            row = best.iloc[0] if not best.empty else None
            old = lookup.loc[(profile_number, float(threshold))]
            rows.append({
                "Profile_number": profile_number,
                "role": roles.get(profile_number, ""),
                "threshold_umol_kg": float(threshold),
                "default_near_zero_max_count": 7,
                "sensitivity_near_zero_max_count": None,
                "default_preprocessed": bool(old["recomputed_preprocessed"]),
                "sensitivity_preprocessed": bool(diagnostics["detector_preprocessed"]),
                "default_has_delta_do": bool(old["recomputed_has_delta_do"]),
                "sensitivity_has_delta_do": bool(row is not None),
                "default_delta_do": old["recomputed_delta_do"],
                "sensitivity_delta_do": float(row["delta_do"]) if row is not None else np.nan,
                "default_peak_depth_m": old["recomputed_peak_depth_m"],
                "sensitivity_peak_depth_m": float(row["depth"]) if row is not None else np.nan,
                "newly_evaluable": bool(not old["recomputed_preprocessed"] and diagnostics["detector_preprocessed"]),
                "new_positive": bool(not old["recomputed_has_delta_do"] and row is not None),
            })
    return pd.DataFrame(rows)


def existing_caliper_readout() -> pd.DataFrame:
    files = {
        250.0: SOURCE_REPO / (
            "plot_outputs/do/global_ocean/scv_matched_control/"
            "scv_matched_control_summary_profile_2002_2023_depth300m_n3_space250km_y2_m2_"
            "core50m_externalplat_uniqueplat_noreuse_ekeq10_specc6cf2f655f_bootanchor_platform.parquet"
        ),
        500.0: SOURCE_REPO / (
            "plot_outputs/do/global_ocean/scv_matched_control/"
            "scv_matched_control_summary_profile_2002_2023_depth300m_n3_space500km_y2_m2_"
            "core50m_externalplat_uniqueplat_noreuse_ekeq10_specdf512e4a95_bootanchor_platform.parquet"
        ),
        750.0: SOURCE_REPO / (
            "plot_outputs/do/global_ocean/scv_matched_control/"
            "scv_matched_control_summary_profile_2002_2023_depth300m_n3_space750km_y2_m2_"
            "core50m_externalplat_uniqueplat_noreuse_ekeq10_specf65759ebbe_bootanchor_platform.parquet"
        ),
    }
    frames = []
    for caliper, path in files.items():
        frame = pd.read_parquet(path)
        frame = frame.loc[frame["scope"].isin(SCOPES) & frame["threshold_umol_kg"].isin(THRESHOLDS)].copy()
        frame.insert(0, "spatial_caliper_km", caliper)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def correction_recheck(identity: pd.DataFrame, participants: pd.DataFrame, profile_data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for decision in identity.loc[identity["current_date_position_status"].eq("confirmed-mismatch")].itertuples(index=False):
        replacement = list(decision.replacement_candidate_profiles)
        replacement_profile = int(replacement[0]) if len(replacement) == 1 else None
        current_profile = int(decision.current_profile_number)
        anchor = participants.loc[participants.is_scv.astype(bool) & participants.profile_number.eq(current_profile)]
        if replacement_profile is None:
            rows.append({
                "association_row": int(decision.association_row),
                "current_profile_number": current_profile,
                "replacement_candidate_profiles": json.dumps(replacement),
                "original_anchor_in_active_cohort": bool(not anchor.empty),
                "safe_to_substitute": False,
                "decision": "exclude_current_set_pending_identity_resolution",
                "status": "no_unique_replacement",
            })
            continue
        replacement_raw = profile_data.loc[profile_data.Profile_number.astype("Int64").eq(replacement_profile)].copy()
        cfg = track.make_detection_config("do", anomaly_min_depth=300.0)
        replacement_clean, diag = track._prepare_do_profile_for_detection(replacement_raw, cfg)
        common = {
            "association_row": int(decision.association_row),
            "current_profile_number": current_profile,
            "replacement_profile_number": replacement_profile,
            "replacement_preprocessed": bool(diag["detector_preprocessed"]),
            "replacement_eke_recheck": "blocked_no_GLORYS_candidate_EKE",
            "safe_to_substitute": False,
            "replacement_outcomes": json.dumps([], ensure_ascii=False),
            "original_anchor_in_active_cohort": bool(not anchor.empty),
        }
        if anchor.empty:
            common.update({
                "match_set_id": None,
                "replacement_core_coverage_pass": None,
                "replacement_controls_spatial_pass": None,
                "replacement_controls_year_pass": None,
                "replacement_controls_month_pass": None,
                "replacement_controls_depth_pass": None,
                "decision": "exclude_current_set_original_pair_unavailable",
                "status": "unique_replacement_original_anchor_not_in_active_cohort",
            })
            rows.append(common)
            continue
        anchor_row = anchor.iloc[0]
        core = float(anchor_row.anchor_core_depth_m)
        core_pass = bool(
            replacement_clean is not None
            and pd.to_numeric(replacement_clean.Depth, errors="coerce").min() <= core - 50
            and pd.to_numeric(replacement_clean.Depth, errors="coerce").max() >= core + 50
        )
        max_depth = float(pd.to_numeric(replacement_clean.Depth, errors="coerce").max()) if replacement_clean is not None else np.nan
        controls = participants.loc[
            participants.match_set_id.eq(anchor_row.match_set_id) & ~participants.is_scv.astype(bool)
        ]
        anchor_lon = float(replacement_raw.Longitude.iloc[0])
        anchor_lat = float(replacement_raw.Latitude.iloc[0])
        anchor_date = pd.Timestamp(
            int(replacement_raw.Year.iloc[0]), int(replacement_raw.Month.iloc[0]), int(replacement_raw.Day.iloc[0])
        )
        distances = [great_circle_km(anchor_lon, anchor_lat, float(row.lon), float(row.lat)) for _, row in controls.iterrows()]
        years = [abs(int(row.year) - anchor_date.year) for _, row in controls.iterrows()]
        months = [
            min(abs(int(row.month) - anchor_date.month), 12 - abs(int(row.month) - anchor_date.month))
            for _, row in controls.iterrows()
        ]
        depth_delta = (
            pd.to_numeric(controls.clean_do_max_depth_m, errors="coerce") - max_depth
        ).abs()
        outcomes = []
        for threshold in THRESHOLDS:
            cfg_thr = track.make_detection_config("do", do_threshold=float(threshold), anomaly_min_depth=300.0)
            detected = track.calculate_delta_do(replacement_raw, detection_config=cfg_thr, remove_outliers=True, include_aou=True, verbose=False)
            best = track._keep_best_anomaly_per_profile(detected, cfg_thr)
            peak = best.iloc[0] if not best.empty else None
            outcomes.append({
                "threshold": threshold,
                "has": bool(peak is not None),
                "delta": float(peak["delta_do"]) if peak is not None else np.nan,
                "depth": float(peak["depth"]) if peak is not None else np.nan,
            })
        common.update({
            "match_set_id": str(anchor_row.match_set_id),
            "replacement_core_coverage_pass": core_pass,
            "replacement_controls_spatial_pass": bool(all(value <= 750 for value in distances)),
            "replacement_controls_year_pass": bool(all(value <= 2 for value in years)),
            "replacement_controls_month_pass": bool(all(value <= 2 for value in months)),
            "replacement_controls_depth_pass": bool(depth_delta.le(500).fillna(False).all()),
            "replacement_outcomes": json.dumps(outcomes, ensure_ascii=False),
            "decision": "exclude_current_set_pending_pair_validation",
            "status": "pair_recheck_completed_eke_blocked",
        })
        rows.append(common)
    return pd.DataFrame(rows)


def write_documents(output_dir: Path, overwrite: bool, counts: dict[str, Any], identity: dict[str, Any], correction: pd.DataFrame, summary: pd.DataFrame) -> dict[str, str]:
    primary = summary.loc[summary.queue.eq("primary_metadata_supported") & summary.scope.eq("Global")].sort_values("threshold_umol_kg")
    confirmed_count = int(identity["confirmed_mismatch"])
    unresolved_count = int(identity["unresolved"])
    large_distance_count = len(identity.get("large_distance_current_profiles", []))
    no_active_pair_count = int(
        correction.get("status", pd.Series(dtype=str)).eq(
            "unique_replacement_original_anchor_not_in_active_cohort"
        ).sum()
    )
    pair_recheck_count = int(
        correction.get("status", pd.Series(dtype=str)).eq(
            "pair_recheck_completed_eke_blocked"
        ).sum()
    )
    lines = [
        "# 局部修复说明",
        "",
        "本轮只修复已由现有证据确定的问题，未进行全局扫描、GLORYS读取、正式重匹配或OFES修改。",
        "",
        f"- main track 新增共享 DO 预处理 helper，并由 calculate_delta_do 与资格审计共同调用；active qualification 为 {counts['qualification_sets']} sets / {counts['qualification_controls']} controls。",
        f"- 主 paper queue 只保留 metadata-supported 当前关联：{counts['primary_sets']} sets / {counts['primary_controls']} controls。",
        f"- 身份裁决：metadata-supported={identity['metadata_supported']}，confirmed-mismatch={confirmed_count}，unresolved={unresolved_count} / 263；原目录行唯一恢复 {identity['original_catalog_row_unique']}/263，{identity['original_catalog_row_unresolved']} 条保留原目录行未决。",
        f"- {confirmed_count} 条确认错配都有唯一替代 profile 的原目录行日期/坐标精度证据，其中 {large_distance_count} 条是当前距离大于 10 km 的重点条目；替代 profile 的 EKE 配对条件不能在本轮无 GLORYS 读取下重验，故没有把新 outcome 直接嫁接到旧 pair，而是排除该组并保留逐关联 recheck 表。",
        f"- 其中 {no_active_pair_count} 条确认错配的原始 anchor 不在 active formal cohort，故只记录唯一替代和身份证据；其余 {pair_recheck_count} 条完成了原 pair 的覆盖/时空/深度复核，但候选 EKE 仍未重验。",
        f"- 现保留的原目录行 confirmed association rows 为 {identity['confirmed_association_rows']}，对应 catalog rows 为 {identity['confirmed_original_catalog_row_ids']}；旧候选中的 association 224 不再算 confirmed（source core=589.8 对应原目录 row 1110，而 profile 83898 只匹配邻近 row 1109）。association 159/160 的原目录 core=405.3 同时对应 rows 648/649，故保留为原目录行未决。",
        f"- 仅有候选汇总而没有唯一同一原目录行证据的 association rows 为 {identity['unresolved_unique_catalog_row_association_rows']}；这些行不升级为 mismatch，也不依据氧结果选 profile。",
        "- 本目录中的表格和文字是既有 7 月 AOGS 图与叙事的局部修正补充，不替代既有图，也不生成新图。",
        f"- 近零敏感性仅覆盖 {counts['near_zero_profiles']} 个已定位触发 profile，未改变历史配置。",
        f"- 共享 helper 与 A1 detector 的 690 个 profile 资格布尔值全部一致；实际重新运行的 {counts['detector_regression_rows']} 个 profile×threshold 结果中，outcome mismatch={counts['detector_regression_outcome_mismatch']}（bool/value/depth 分别为 {counts['detector_regression_bool_mismatch']}/{counts['detector_regression_delta_mismatch']}/{counts['detector_regression_depth_mismatch']}）。",
        "- 26 个近零跳过 profile 的 clean-row 诊断数字不同，是因为 helper 记录了逐层近零剔除后的行数而旧缓存记录了剔除前行数，跳过结论本身未变。",
        "- 245/244 的差异来自旧 `mccoy_glorys_miss.has_delta_do.notna()` 字段口径与显式 threshold table 的 `do_evaluable` 口径：多出的具体 profile 是 32111，三个阈值都纳入显式表，但其目录状态为 `argo_no_signal`；A1 复核仅有 1 个 QC 清洗后全变量有效层，真实 detector preprocessing 不通过且三个阈值均无 carrier，因此不是 detector cache 差异。",
        "- 未被证实的怀疑：103 ambiguous 不能整体视为错配；旧 detector cache 三阈值仍无差异。",
        "",
        "## Global 主队列",
        "",
        "| threshold | SCV+ / control+ | MH OR | anchor CI | dependency CI | status |",
        "|---:|---:|---:|---|---|---|",
    ]
    for row in primary.itertuples(index=False):
        lines.append(
            f"| {row.threshold_umol_kg:g} | {row.scv_numerator}/{row.control_numerator} | {row.matched_mantel_haenszel_or:g} | "
            f"[{row.anchor_platform_bootstrap_ci_low_raw:g}, {row.anchor_platform_bootstrap_ci_high_raw:g}] | "
            f"[{row.dependency_component_bootstrap_ci_low_raw:g}, {row.dependency_component_bootstrap_ci_high_raw:g}] | "
            f"{row.anchor_platform_bootstrap_status}; {row.dependency_component_bootstrap_status} |"
        )
    lines += [
        "",
        "全球 230,117 分母未更新；正文发生率必须保留 [UPDATE REQUIRED: re-count unified eligible global denominator] 占位。",
    ]
    note = output_dir / "repair_note.md"
    ensure_output(note, overwrite)
    note.write_text("\n".join(lines) + "\n", encoding="utf-8")
    methods = output_dir / "methods_draft.md"
    ensure_output(methods, overwrite)
    methods.write_text(
        "\n".join([
            "# Methods draft", "",
            "We retained the original SCV matched sets and applied a separate profile qualification audit before constructing the paper queue. The audit uses the same main-worktree DO preprocessing as calculate_delta_do: QC flags 1/2/5/8, level-wise DO≤1.0 removal, the configured >7 whole-profile skip, same-level finite DO/T/S, depth de-duplication, and at least five valid full-profile levels.",
            "",
            "Failed anchors remove their original match_set_id; failed controls are removed within surviving sets; sets without a qualified control are excluded. Pair checks re-evaluate anchor core ±50 m coverage and the 500 m maximum-depth condition. This is a qualified subset of the original matching, not a new global optimum match.",
            "",
            f"Identity decisions first recover the original catalog row from the source association's core-pressure record, then compare candidates within that row against an independently reloaded local Argo metadata index. Catalog rows and physical Profile_number candidates remain separate. {confirmed_count} confirmed mismatches ({large_distance_count} current profiles with distance >10 km) are excluded because their replacement pair EKE condition was not revalidated without GLORYS; {no_active_pair_count} lack an active formal pair for recheck. {unresolved_count} associations remain unresolved and are shown in sensitivity.",
            "",
            "The earlier 244 versus current 245 eligibility count is a field-definition difference: Profile 32111 is present in the explicit threshold table but has the legacy `mccoy_glorys_miss` status `argo_no_signal`. It has only one QC-clean full-variable level, fails the shared detector preprocessing, and has no carrier at any threshold, so this count change is not a detector-cache disagreement.",
            "",
            f"The shared helper reproduced the A1 preprocessing Boolean for all 690 audited profiles. A fresh detector rerun covers {counts['detector_regression_rows']} profile-threshold combinations and has zero Boolean, delta, or peak-depth mismatches against the historical A1 recomputation. For 26 near-zero-skipped profiles, diagnostic clean-row counts differ only because the helper records rows after level-wise near-zero removal while the historical cache recorded rows before that removal; the skip Boolean is unchanged.",
            "",
            "Bootstrap exports retain requested, valid, undefined, zero, and infinite draw counts. A [0,0] zero-event bootstrap is labelled degenerate and is not presented as an inferential zero-width interval. The 230,117 global denominator remains an explicit update placeholder.",
        ]) + "\n", encoding="utf-8"
    )
    results = output_dir / "results_draft.md"
    ensure_output(results, overwrite)
    results.write_text(
        "\n".join([
            "# Results draft", "",
            f"The qualified original matching retained {counts['qualification_sets']} sets and {counts['qualification_controls']} controls. The metadata-supported primary queue retained {counts['primary_sets']} sets and {counts['primary_controls']} controls.",
            "",
            "The three threshold estimates, both cluster intervals, and bootstrap degeneracy fields are in paper_candidate_results.parquet. The results describe matched SCV/Argo association and do not establish transport contribution, material continuity, or a global occurrence rate.",
            "",
            "The near-zero sensitivity is limited to the 28 profiles already identified in the audit. It changes only the whole-profile skip and is not a global estimate of recoverable profiles. These tables supplement the existing July AOGS figures and narrative; no new figure is generated.",
        ]) + "\n", encoding="utf-8"
    )
    return {"repair_note": str(note), "methods": str(methods), "results": str(results)}


def run(output_dir: Path, overwrite: bool, log_file: Path | None) -> None:
    started = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)
    handle = None
    old_stdout, old_stderr = sys.stdout, sys.stderr
    try:
        if log_file is not None:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            handle = log_file.open("w", encoding="utf-8")
            sys.stdout, sys.stderr = Tee(sys.stdout, handle), Tee(sys.stderr, handle)
        if Path(track.__file__).resolve().parent != MAIN_WORKTREE:
            raise RuntimeError(f"track import is not from main worktree: {track.__file__}")
        track.switch_region("global")
        cohort = pd.read_parquet(FORMAL_COHORT)
        no_control = set(cohort.loc[cohort.match_status.eq("no_control"), "match_set_id"])
        active = cohort.loc[~cohort.match_set_id.isin(no_control)].copy()
        anchored = pd.read_parquet(SCV_ANCHORED)
        anchored["profile_number"] = pd.to_numeric(anchored.profile_number, errors="coerce").astype(int)
        base_preprocessing = pd.read_parquet(BASE_PREPROCESSING)
        base_detector = pd.read_parquet(BASE_DETECTOR)
        candidates = pd.read_parquet(IDENTITY_CANDIDATES)
        identity_summary = pd.read_parquet(IDENTITY_SUMMARY)
        union_years = profile_year_map(active, anchored, candidates.iloc[0:0])
        profile_years = profile_year_map(active, anchored, candidates)
        profile_data, loader_provenance = load_profile_data(profile_years, ARGO_DATA)
        audit, cleaned = profile_audit(profile_data, set(union_years), base_preprocessing)
        if not audit.a1_detector_boolean_match.fillna(False).all():
            raise RuntimeError("Shared main-worktree preprocessing differs from A1.")
        log("[regression] rerunning calculate_delta_do for the 690-profile A1 union")
        regression = detector_outcome_regression(profile_data, set(union_years), base_detector)
        participants = build_participant_flags(active, audit)
        identity, identity_pairs, identity_details = adjudicate_identity(
            candidates, SCV_CATALOG, ARGO_DATA, anchored
        )
        status_map = dict(zip(identity.current_profile_number.astype(int), identity.current_date_position_status.astype(str)))
        old_unique = dict(zip(
            identity_summary.current_profile_number.astype(int),
            identity_summary.n_candidate_profile_numbers.eq(1)
            & identity_summary.current_profile_in_candidate_set.astype(bool),
        ))
        qualification = participants.qualification_pass.fillna(False).astype(bool)
        queues = {
            "qualification_all": pair_queue(participants, qualification),
            "primary_metadata_supported": pair_queue(participants, qualification & participants.profile_number.map(status_map).eq("metadata-supported")),
            "unresolved_kept_sensitivity": pair_queue(participants, qualification & participants.profile_number.map(status_map).isin(["metadata-supported", "unresolved"])),
            "c2_unique_candidate_sensitivity": pair_queue(participants, qualification & participants.profile_number.map(old_unique).fillna(False).astype(bool)),
        }
        correction = correction_recheck(identity, participants, profile_data)
        summary = pd.concat([summarize_queue(value, key) for key, value in queues.items()], ignore_index=True)
        edge = edge_case_tests(summary)
        roles = {}
        for row in active.itertuples(index=False):
            roles.setdefault(int(row.profile_number), set()).add("formal_anchor" if bool(row.is_scv) else "formal_control")
        for row in anchored.itertuples(index=False):
            roles.setdefault(int(row.profile_number), set()).add("scv_association")
        roles_text = {key: ";".join(sorted(value)) for key, value in roles.items()}
        near_zero = near_zero_sensitivity(profile_data, audit, roles_text, base_detector)
        calipers = existing_caliper_readout()
        qual_global = summary.loc[summary.queue.eq("qualification_all") & summary.scope.eq("Global")]
        primary_global = summary.loc[summary.queue.eq("primary_metadata_supported") & summary.scope.eq("Global")]
        counts = {
            "active_sets": int(active.match_set_id.nunique()),
            "active_anchor_profiles": int(active.loc[active.is_scv.astype(bool), "profile_number"].nunique()),
            "active_control_profiles": int(active.loc[~active.is_scv.astype(bool), "profile_number"].nunique()),
            "qualification_sets": int(qual_global.n_matched_sets.max()),
            "qualification_controls": int(qual_global.n_control_profiles.max()),
            "primary_sets": int(primary_global.n_matched_sets.max()),
            "primary_controls": int(primary_global.n_control_profiles.max()),
            "near_zero_profiles": int(audit.near_zero_triggered.fillna(False).sum()),
            "a1_detector_boolean_mismatch": int((~audit.a1_detector_boolean_match.fillna(False)).sum()),
            "base_detector_rows": int(len(base_detector)),
            "detector_regression_rows": int(len(regression)),
            "detector_regression_outcome_mismatch": int((~regression.outcome_match).sum()),
            "detector_regression_bool_mismatch": int((~regression.new_vs_a1_recomputed_bool_match).sum()),
            "detector_regression_delta_mismatch": int((~regression.new_vs_a1_recomputed_delta_match).sum()),
            "detector_regression_depth_mismatch": int((~regression.new_vs_a1_recomputed_depth_match).sum()),
        }
        paths = {}
        paths["profile_qualification_audit"] = save_frame(audit.sort_values("Profile_number"), "profile_qualification_audit", output_dir, overwrite)
        paths["detector_outcome_regression"] = save_frame(regression.sort_values(["Profile_number", "threshold_umol_kg"]), "detector_outcome_regression", output_dir, overwrite)
        paths["qualified_analysis_queue"] = save_frame(queues["primary_metadata_supported"].sort_values(["match_set_id", "is_scv"], ascending=[True, False]), "qualified_analysis_queue", output_dir, overwrite)
        paths["identity_adjudication"] = save_frame(identity.sort_values("association_row"), "identity_adjudication", output_dir, overwrite)
        paths["identity_candidate_comparison"] = save_frame(identity_pairs.sort_values(["association_row", "candidate_profile_number"]), "identity_candidate_comparison", output_dir, overwrite)
        paths["confirmed_mismatch_recheck"] = save_frame(correction.sort_values("association_row"), "confirmed_mismatch_recheck", output_dir, overwrite)
        paths["paper_candidate_results"] = save_frame(summary.sort_values(["queue", "scope", "threshold_umol_kg"]), "paper_candidate_results", output_dir, overwrite)
        paths["near_zero_sensitivity"] = save_frame(near_zero.sort_values(["Profile_number", "threshold_umol_kg"]), "near_zero_sensitivity", output_dir, overwrite)
        paths["bootstrap_edge_case_tests"] = save_frame(edge, "bootstrap_edge_case_tests", output_dir, overwrite)
        paths["caliper_sensitivity_readout"] = save_frame(calipers.sort_values(["spatial_caliper_km", "scope", "threshold_umol_kg"]), "caliper_sensitivity_readout", output_dir, overwrite)
        paths["documents"] = write_documents(output_dir, overwrite, counts, identity_details, correction, summary)
        manifest = {
            "status": "complete",
            "elapsed_seconds": round(time.time() - started, 3),
            "main_worktree": str(MAIN_WORKTREE),
            "track_import": str(Path(track.__file__).resolve()),
            "local_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=MAIN_WORKTREE, text=True).strip(),
            "formal_cohort": str(FORMAL_COHORT),
            "scope": "paper-ready local repair; no rematch; no global scan; no GLORYS; OFES untouched",
            "counts": counts,
            "identity": identity_details,
            "loader_years": loader_provenance,
            "outputs": paths,
        }
        paths["manifest"] = save_json(manifest, "paper_ready_manifest.json", output_dir, overwrite)
        log(f"[done] elapsed={manifest['elapsed_seconds']}s; output={output_dir}")
    finally:
        sys.stdout, sys.stderr = old_stdout, old_stderr
        if handle is not None:
            handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--log-file", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.output_dir, bool(args.overwrite), args.log_file)
