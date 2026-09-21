#!/usr/bin/env python
"""Validate the bounded McCoy-to-Argo matcher repair without GLORYS reads."""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

MAIN_WORKTREE = Path(__file__).resolve().parent
# This runner now lives at the repository root, alongside track.py.
SOURCE_REPO = MAIN_WORKTREE
AUDIT_ROOT = SOURCE_REPO / "plot_outputs/test/argo_scv_targeted_audit_20260920"
IDENTITY_TABLE = AUDIT_ROOT / "paper_ready/identity_adjudication.parquet"
SCV_CATALOG = SOURCE_REPO / "data/mccoy2020_scv/SCVs.csv"
ARGO_DATA = SOURCE_REPO / "Argo_data"

if str(MAIN_WORKTREE) not in sys.path:
    sys.path.insert(0, str(MAIN_WORKTREE))
import track


def decimal_places(value: object) -> int:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return 0
    if "e" in text.lower():
        return 8
    return len(text.split(".", 1)[1].rstrip("0")) if "." in text else 0


def synthetic_validation() -> list[dict[str, object]]:
    rows = []

    def add_profile(profile: int, platform: int, date: str, lon: float, lat: float) -> None:
        day = pd.Timestamp(date)
        rows.append({
            "Profile_number": profile, "Platform_number": platform,
            "Year": day.year, "Month": day.month, "Day": day.day,
            "Longitude": lon, "Latitude": lat,
        })

    add_profile(210, 2, "2020-01-01", 1.0, 2.0)
    add_profile(211, 2, "2020-01-02", 9.0, 9.0)
    add_profile(310, 3, "2020-05-05", 10.0, 20.0)
    add_profile(311, 3, "2020-05-05", 11.0, 21.0)
    add_profile(410, 4, "2020-06-06", 30.0, 40.0)
    add_profile(411, 4, "2020-06-06", 31.0, 41.0)
    add_profile(510, 5, "2019-12-31", 50.0, 60.0)
    add_profile(610, 6, "2020-08-08", -180.001, 5.0)
    add_profile(611, 6, "2020-08-08", 170.0, 5.0)
    add_profile(710, 7, "2021-09-09", 70.0, 8.0)
    index = pd.DataFrame(rows)

    catalog = []

    def add_catalog(
        row_id: int, platform: int, date: str, lon: float, lat: float,
    ) -> None:
        catalog.append({
            "catalog_row_id": row_id, "ID": row_id + 1,
            "Platform": platform, "Cycle": row_id + 1,
            "date": pd.Timestamp(date), "Longitude": lon, "Latitude": lat,
            "Core_Pressure": 100.0, "Shallow_Pressure": 50.0,
            "Deep_Pressure": 150.0, "SCV_Type": "S", "basin": "test",
            "catalog_lon_decimals": 3, "catalog_lat_decimals": 3,
        })

    add_catalog(0, 2, "2020-01-01 23:30", 1.0, 2.0)
    add_catalog(1, 3, "2020-05-05 12:00", 11.0, 21.0)
    add_catalog(2, 4, "2020-06-06 12:00", 32.0, 42.0)
    add_catalog(3, 5, "2020-01-01 12:00", 50.0, 60.0)
    add_catalog(4, 6, "2020-08-08 14:00", 179.999, 5.0)
    add_catalog(5, 7, "2021-09-09 00:01", 70.0, 8.0)
    catalog_frame = pd.DataFrame(catalog)

    with tempfile.TemporaryDirectory(prefix="matcher_source_fix_", dir="/tmp") as temp_dir:
        argo_dir = Path(temp_dir)
        for year, frame in index.groupby("Year"):
            frame.to_parquet(argo_dir / f"Argo{year}.parquet", index=False)
        result = track._match_mccoy_scvs_to_profiles(catalog_frame, argo_dir, 1)

    expected = {0: 210, 1: 311, 2: None, 3: 510, 4: 610, 5: 710}
    expected_status = {
        0: "unique_same_calendar_day",
        1: "same_day_coordinate_supported",
        2: "unresolved_same_day_multiple_profiles",
        3: "window_coordinate_supported",
        4: "same_day_coordinate_supported",
        5: "unique_same_calendar_day",
    }
    checks = []
    for row in result.itertuples(index=False):
        actual = None if pd.isna(row.profile_number) else int(row.profile_number)
        checks.append({
            "test": f"synthetic_row_{row.catalog_row_id}",
            "passed": bool(actual == expected[int(row.catalog_row_id)]
                           and row.match_status == expected_status[int(row.catalog_row_id)]),
            "match_status": row.match_status,
            "selected_profile": actual,
            "expected_profile": expected[int(row.catalog_row_id)],
        })
    return checks


def real_263_validation() -> tuple[list[dict[str, object]], dict[str, object]]:
    identity = pd.read_parquet(IDENTITY_TABLE)

    def as_list(value: object) -> list[object]:
        return list(value) if isinstance(value, (list, tuple, np.ndarray, pd.Series)) else []

    catalog_ids = sorted({
        int(value)
        for values in identity["original_catalog_row_ids"]
        for value in as_list(values)
    })
    raw = pd.read_csv(SCV_CATALOG, dtype=str)
    catalog = pd.read_csv(SCV_CATALOG)
    catalog.insert(0, "catalog_row_id", np.arange(len(catalog), dtype=int))
    catalog["date"] = pd.to_datetime(
        catalog["Cycle_ISO_DateTime_UTC"], errors="coerce"
    ).dt.tz_localize(None)
    catalog["basin"] = [
        track._residual_basin_of(lat, lon)
        for lat, lon in zip(catalog["Latitude"], catalog["Longitude"])
    ]
    catalog["catalog_lon_decimals"] = raw["Longitude"].map(decimal_places)
    catalog["catalog_lat_decimals"] = raw["Latitude"].map(decimal_places)
    subset = catalog.loc[catalog["catalog_row_id"].isin(catalog_ids)].copy()
    matched = track._match_mccoy_scvs_to_profiles(subset, ARGO_DATA, 5)
    selected = dict(zip(matched["catalog_row_id"].astype(int), matched["profile_number"]))
    source_rows = matched.set_index("catalog_row_id", drop=False)
    checks = []
    for row in identity.loc[identity["current_date_position_status"].isin(
        ["confirmed-mismatch", "metadata-supported"]
    )].itertuples(index=False):
        catalog_row_id = int(row.original_catalog_row_id)
        actual = selected.get(catalog_row_id)
        source = source_rows.loc[catalog_row_id]
        expected = (
            int(row.replacement_candidate_profiles[0])
            if row.current_date_position_status == "confirmed-mismatch"
            else int(row.current_profile_number)
        )
        checks.append({
            "test": f"real_association_{row.association_row}",
            "passed": bool(pd.notna(actual) and int(actual) == expected),
            "match_status": source.match_status,
            "identity_status": row.current_date_position_status,
            "catalog_row_id": catalog_row_id,
            "selected_profile": None if pd.isna(actual) else int(actual),
            "expected_profile": expected,
            "source_coordinate_match": bool(source.coordinate_match),
            "source_calendar_day_delta": source.calendar_day_delta,
            "source_minimal_lon_diff_deg": source.minimal_lon_diff_deg,
            "source_lat_diff_deg": source.lat_diff_deg,
            "source_n_profile_candidates": int(source.n_profile_candidates),
        })
    special = matched.loc[matched["catalog_row_id"].isin([648, 649, 1110])]
    special_map = dict(zip(special["catalog_row_id"].astype(int), special["profile_number"]))
    details = {
        "catalog_rows_tested": int(len(subset)),
        "matcher_rows": int(len(matched)),
        "status_counts": matched["match_status"].value_counts().to_dict(),
        "confirmed_mismatch_checks": int(identity["current_date_position_status"].eq("confirmed-mismatch").sum()),
        "confirmed_mismatch_failures": int(sum(
            not item["passed"] for item in checks
            if item["identity_status"] == "confirmed-mismatch"
        )),
        "metadata_supported_checks": int(identity["current_date_position_status"].eq("metadata-supported").sum()),
        "metadata_supported_failures": int(sum(
            not item["passed"] for item in checks
            if item["identity_status"] == "metadata-supported"
        )),
        "confirmed_mismatch_source_status_counts": pd.Series([
            item["match_status"] for item in checks
            if item["identity_status"] == "confirmed-mismatch"
        ]).value_counts().to_dict(),
        "metadata_supported_source_status_counts": pd.Series([
            item["match_status"] for item in checks
            if item["identity_status"] == "metadata-supported"
        ]).value_counts().to_dict(),
        "special_rows_648_649_1110": special_map,
    }
    return checks, details


def screen_compatibility_validation() -> dict[str, object]:
    """Exercise the screen consumer with one resolved and one unresolved row."""
    with tempfile.TemporaryDirectory(prefix="matcher_screen_compatibility_", dir="/tmp") as temp_dir:
        root = Path(temp_dir)
        argo_dir = root / "argo"
        argo_dir.mkdir()
        pd.DataFrame([
            {"Profile_number": 810, "Platform_number": 8, "Year": 2020,
             "Month": 1, "Day": 1, "Longitude": 10.0, "Latitude": 20.0},
            {"Profile_number": 811, "Platform_number": 9, "Year": 2020,
             "Month": 2, "Day": 2, "Longitude": 30.0, "Latitude": 40.0},
            {"Profile_number": 812, "Platform_number": 9, "Year": 2020,
             "Month": 2, "Day": 2, "Longitude": 31.0, "Latitude": 41.0},
        ]).to_parquet(argo_dir / "Argo2020.parquet", index=False)
        catalog_path = root / "scv.csv"
        pd.DataFrame([
            {"ID": 1, "Platform": 8, "Cycle": 1,
             "Cycle_ISO_DateTime_UTC": "2020-01-01T15:00:00Z",
             "Longitude": 10.0, "Latitude": 20.0, "Core_Pressure": 100.0,
             "Shallow_Pressure": 50.0, "Deep_Pressure": 150.0, "SCV_Type": "S"},
            {"ID": 2, "Platform": 9, "Cycle": 2,
             "Cycle_ISO_DateTime_UTC": "2020-02-02T12:00:00Z",
             "Longitude": 35.0, "Latitude": 45.0, "Core_Pressure": 110.0,
             "Shallow_Pressure": 55.0, "Deep_Pressure": 160.0, "SCV_Type": "S"},
        ]).to_csv(catalog_path, index=False)
        output_dir = root / "screen_output"
        original = track.glorys_reproduces_mccoy_scv

        def fake_glorys(*args: object, **kwargs: object) -> dict[str, object]:
            return {
                "lon": 10.0, "lat": 20.0, "date": "2020-01-01",
                "core_depth_m": 100.0, "argo_dpi_core": 1.0,
                "argo_lens_sign": "warm-salty", "glorys_dpi_core": 0.2,
                "glorys_dpi_ratio": 0.2, "glorys_n2_core": 1.0,
                "glorys_misses": False, "status": "ok", "skip_reason": None,
            }

        track.glorys_reproduces_mccoy_scv = fake_glorys
        try:
            summary = track.screen_mccoy_scvs_against_glorys(
                mccoy_csv=catalog_path,
                argo_data_dir=argo_dir,
                output_dir=output_dir,
                match_days=1,
                return_details=True,
            )
        finally:
            track.glorys_reproduces_mccoy_scv = original
        rows = summary["rows"]
        unresolved = rows.loc[rows["status"].eq("unresolved_profile_match")]
        passed = bool(
            summary["n_catalog"] == 2
            and summary["n_matched"] == 1
            and summary["n_unresolved_profile_match"] == 1
            and summary["n_evaluated"] == 1
            and len(unresolved) == 1
            and unresolved.iloc[0]["match_status"] == "unresolved_same_day_multiple_profiles"
        )
        return {
            "passed": passed,
            "n_catalog": int(summary["n_catalog"]),
            "n_matched": int(summary["n_matched"]),
            "n_unresolved_profile_match": int(summary["n_unresolved_profile_match"]),
            "n_evaluated": int(summary["n_evaluated"]),
        }


def main(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    checks = synthetic_validation()
    real_checks, details = real_263_validation()
    checks.extend(real_checks)
    screen_details = screen_compatibility_validation()
    frame = pd.DataFrame(checks)
    frame.to_csv(output_dir / "matcher_source_fix_validation.csv", index=False)
    failed_checks = int((~frame["passed"].astype(bool)).sum())
    if not screen_details["passed"]:
        failed_checks += 1
    synthetic = [row for row in checks if row["test"].startswith("synthetic_")]
    payload = {
        "track_import": str(Path(track.__file__).resolve()),
        "synthetic_checks": len(synthetic),  # Legacy total-count field.
        "synthetic_total": len(synthetic),
        "synthetic_passed": sum(bool(row["passed"]) for row in synthetic),
        "real_checks": int(len(real_checks)),
        "failed_checks": failed_checks,
        "real_details": details,
        "screen_compatibility": screen_details,
        "scope": "source matcher only; no GLORYS and no historical cache refresh",
    }
    (output_dir / "matcher_source_fix_validation.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    note = "\n".join([
        "# McCoy–Argo source matcher repair",
        "",
        f"- Import checked: `{payload['track_import']}`.",
        "- Rule: compare catalog and Argo at calendar-day precision, then use source-precision independent coordinates; unresolved multiple candidates remain explicit.",
        f"- Synthetic cases: {payload['synthetic_passed']}/{payload['synthetic_total']} passed (afternoon timestamp, same-day coordinate tie-break, unresolved tie, cross-year, dateline, unique match).",
        f"- Existing catalog subset: {details['catalog_rows_tested']} original rows, {payload['failed_checks']} failed checks across matcher and screen checks.",
        f"- Status counts: `{json.dumps(details['status_counts'], ensure_ascii=False, sort_keys=True)}`.",
        f"- Confirmed old mismatches retained the independently supported replacement: {details['confirmed_mismatch_checks']} checked, {details['confirmed_mismatch_failures']} failures.",
        f"- Confirmed-mismatch source statuses: `{json.dumps(details['confirmed_mismatch_source_status_counts'], ensure_ascii=False, sort_keys=True)}`.",
        f"- Metadata-supported current associations retained: {details['metadata_supported_checks']} checked, {details['metadata_supported_failures']} failures.",
        f"- Diagnostic rows 648/649/1110 select profiles `{details['special_rows_648_649_1110']}`.",
        f"- Screen compatibility: `{json.dumps(screen_details, ensure_ascii=False, sort_keys=True)}` (GLORYS consumer was monkeypatched; no GLORYS read).",
        "- Scope: source matcher validation only; no oxygen selection, GLORYS read, historical cache refresh, or cohort statistic update.",
        "",
    ])
    (output_dir / "matcher_source_fix_note.md").write_text(note, encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2, default=str))
    if payload["failed_checks"]:
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=AUDIT_ROOT / "paper_ready")
    main(parser.parse_args().output_dir)
