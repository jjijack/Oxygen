#!/usr/bin/env python
"""Estimate SCV matched effects from the formal matched-control cohort.

The runner lives at the repository root and imports the adjacent track module;
inputs and outputs resolve through the configured data paths, relative to the
launch directory, as in track.py.
It re-qualifies the original matched sets with the shared DO preprocessing,
keeps anchors whose catalogue association is metadata-supported according to
``track.audit_scv_profile_identity``, and writes the qualified queue, matched
Mantel-Haenszel ORs with cluster-bootstrap intervals (whole-profile and
core-near maxima, globally, outside and inside the KE), and the near-zero and
caliper sensitivities to ``plot_outputs/do/global_ocean/scv_matched_effects``.
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

CODE_ROOT = Path(__file__).resolve().parent

if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))
import track

DEFAULT_OUTPUT = track.make_detection_config("do").output_dir("scv_matched_effects", "global_ocean")
MATCHED_CONTROL_DIR = track.make_detection_config("do").output_dir("scv_matched_control", "global_ocean")
FORMAL_COHORT = MATCHED_CONTROL_DIR / (
    "scv_matched_control_cohort_2002_2023_depth300m_n3_space750km_y2_m2_"
    "core50m_externalplat_uniqueplat_noreuse_ekeq10_specf65759ebbe.parquet"
)
SCV_ANCHORED = track._default_mccoy_anchored_path("global_ocean")
ARGO_DATA = track.argo_path
THRESHOLDS = (20.0, 35.0, 50.0)
SCOPES = ("Global", "Global excluding KE")
PRIMARY_QUEUE = "primary_metadata_supported"


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




def load_profile_data(profile_years: dict[int, int], data_dir: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    by_year: dict[int, set[int]] = {}
    for profile_number, year in profile_years.items():
        by_year.setdefault(int(year), set()).add(int(profile_number))
    frames = []
    for year, ids in sorted(by_year.items()):
        log(f"[loader] year={year}, profiles={len(ids)}")
        frames.append(track.load_argo_data(year, data_dir, profile_ids=ids))
    profile_data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return profile_data, {str(year): len(ids) for year, ids in sorted(by_year.items())}


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


def profile_audit(profile_data: pd.DataFrame, profile_ids: set[int]) -> tuple[pd.DataFrame, dict[int, pd.DataFrame]]:
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
    return pd.DataFrame(rows), cleaned_profiles


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
        status, usable = track._bootstrap_draw_status(100, valid, zero, infinite, low, high)
        rows.append({
            "test": name, "mh_or": mh, "ci_low_raw": low, "ci_high_raw": high,
            "requested": 100, "valid": valid, "undefined_draws": 100 - valid,
            "zero_draws": zero, "infinite_draws": infinite, "status": status,
            "all_draws_defined_finite_and_positive": usable,
        })
    normal = summary.loc[
        summary.queue.eq("qualification_all")
        & summary.outcome.eq("profile")
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


def near_zero_sensitivity(profile_data: pd.DataFrame, audit: pd.DataFrame, roles: dict[int, str]) -> pd.DataFrame:
    """Rerun the profiles skipped by the whole-profile near-zero rule with the rule disabled."""
    ids = audit.loc[audit["near_zero_triggered"].fillna(False).astype(bool), "Profile_number"].astype(int).tolist()
    default_cfg = track.make_detection_config("do", anomaly_min_depth=300.0)
    sensitivity_cfg = track.make_detection_config(default_cfg, do_near_zero_max_count=None)

    def detect(profile: pd.DataFrame, cfg) -> tuple[bool, pd.Series | None]:
        _, diagnostics = track._prepare_do_profile_for_detection(profile, cfg)
        detected = track.calculate_delta_do(
            profile, detection_config=cfg, remove_outliers=True, include_aou=True, verbose=False
        )
        best = track._keep_best_anomaly_per_profile(detected, cfg)
        return bool(diagnostics["detector_preprocessed"]), (best.iloc[0] if not best.empty else None)

    rows = []
    for profile_number in ids:
        profile = profile_data.loc[profile_data["Profile_number"].astype("Int64").eq(profile_number)].copy()
        for threshold in THRESHOLDS:
            default_prep, default_row = detect(profile, track.make_detection_config(default_cfg, do_threshold=float(threshold)))
            sens_prep, sens_row = detect(profile, track.make_detection_config(sensitivity_cfg, do_threshold=float(threshold)))
            rows.append({
                "Profile_number": profile_number,
                "role": roles.get(profile_number, ""),
                "threshold_umol_kg": float(threshold),
                "default_near_zero_max_count": int(default_cfg.do_near_zero_max_count),
                "sensitivity_near_zero_max_count": None,
                "default_preprocessed": default_prep,
                "sensitivity_preprocessed": sens_prep,
                "default_has_delta_do": default_row is not None,
                "sensitivity_has_delta_do": sens_row is not None,
                "default_delta_do": float(default_row["delta_do"]) if default_row is not None else np.nan,
                "sensitivity_delta_do": float(sens_row["delta_do"]) if sens_row is not None else np.nan,
                "default_peak_depth_m": float(default_row["depth"]) if default_row is not None else np.nan,
                "sensitivity_peak_depth_m": float(sens_row["depth"]) if sens_row is not None else np.nan,
                "newly_evaluable": bool(not default_prep and sens_prep),
                "new_positive": bool(default_row is None and sens_row is not None),
            })
    return pd.DataFrame(rows)


def existing_caliper_readout() -> pd.DataFrame:
    files = {
        250.0: MATCHED_CONTROL_DIR / (
            "scv_matched_control_summary_profile_2002_2023_depth300m_n3_space250km_y2_m2_"
            "core50m_externalplat_uniqueplat_noreuse_ekeq10_specc6cf2f655f_bootanchor_platform.parquet"
        ),
        500.0: MATCHED_CONTROL_DIR / (
            "scv_matched_control_summary_profile_2002_2023_depth300m_n3_space500km_y2_m2_"
            "core50m_externalplat_uniqueplat_noreuse_ekeq10_specdf512e4a95_bootanchor_platform.parquet"
        ),
        750.0: MATCHED_CONTROL_DIR / (
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
        replacement = [int(value) for value in decision.replacement_candidate_profiles]
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
        if Path(track.__file__).resolve().parent != CODE_ROOT:
            raise RuntimeError(f"track import is not from {CODE_ROOT}: {track.__file__}")
        track.switch_region("global")
        cohort = pd.read_parquet(FORMAL_COHORT)
        no_control = set(cohort.loc[cohort.match_status.eq("no_control"), "match_set_id"])
        active = cohort.loc[~cohort.match_set_id.isin(no_control)].copy()
        anchored = pd.read_parquet(SCV_ANCHORED)
        anchored["profile_number"] = pd.to_numeric(anchored.profile_number, errors="coerce").astype(int)
        identity_tables = track.load_scv_profile_identity()
        candidates = identity_tables["candidates"]
        identity_summary = identity_tables["summary"]
        identity = identity_tables["adjudication"]
        union_years = profile_year_map(active, anchored, candidates.iloc[0:0])
        profile_years = profile_year_map(active, anchored, candidates)
        profile_data, loader_provenance = load_profile_data(profile_years, ARGO_DATA)
        audit, cleaned = profile_audit(profile_data, set(union_years))
        participants = build_participant_flags(active, audit)
        status_map = dict(zip(identity.current_profile_number.astype(int), identity.current_date_position_status.astype(str)))
        old_unique = dict(zip(
            identity_summary.current_profile_number.astype(int),
            identity_summary.n_candidate_profile_numbers.eq(1)
            & identity_summary.current_profile_in_candidate_set.astype(bool),
        ))
        qualification = participants.qualification_pass.fillna(False).astype(bool)
        queues = {
            "qualification_all": pair_queue(participants, qualification),
            PRIMARY_QUEUE: pair_queue(participants, qualification & participants.profile_number.map(status_map).eq("metadata-supported")),
            "unresolved_kept_sensitivity": pair_queue(participants, qualification & participants.profile_number.map(status_map).isin(["metadata-supported", "unresolved"])),
            "c2_unique_candidate_sensitivity": pair_queue(participants, qualification & participants.profile_number.map(old_unique).fillna(False).astype(bool)),
        }
        correction = correction_recheck(identity, participants, profile_data)
        summaries = [
            track._summarize_scv_matched_queue(
                value, key, scopes=SCOPES,
                seed_label="qualification" if key == "qualification_all" else key,
            )
            for key, value in queues.items()
        ]
        summaries.append(track._summarize_scv_matched_queue(queues[PRIMARY_QUEUE], PRIMARY_QUEUE, scopes=("KE",)))
        summaries.append(track._summarize_scv_matched_queue(
            queues[PRIMARY_QUEUE], PRIMARY_QUEUE, scopes=SCOPES + ("KE",), outcome="core_aligned",
        ))
        summary = pd.concat(summaries, ignore_index=True)
        edge = edge_case_tests(summary)
        roles = {}
        for row in active.itertuples(index=False):
            roles.setdefault(int(row.profile_number), set()).add("formal_anchor" if bool(row.is_scv) else "formal_control")
        for row in anchored.itertuples(index=False):
            roles.setdefault(int(row.profile_number), set()).add("scv_association")
        roles_text = {key: ";".join(sorted(value)) for key, value in roles.items()}
        near_zero = near_zero_sensitivity(profile_data, audit, roles_text)
        calipers = existing_caliper_readout()
        qual_global = summary.loc[summary.queue.eq("qualification_all") & summary.scope.eq("Global")]
        primary_global = summary.loc[summary.queue.eq(PRIMARY_QUEUE) & summary.scope.eq("Global")]
        counts = {
            "active_sets": int(active.match_set_id.nunique()),
            "active_anchor_profiles": int(active.loc[active.is_scv.astype(bool), "profile_number"].nunique()),
            "active_control_profiles": int(active.loc[~active.is_scv.astype(bool), "profile_number"].nunique()),
            "qualification_sets": int(qual_global.n_matched_sets.max()),
            "qualification_controls": int(qual_global.n_control_profiles.max()),
            "primary_sets": int(primary_global.n_matched_sets.max()),
            "primary_controls": int(primary_global.n_control_profiles.max()),
            "near_zero_profiles": int(audit.near_zero_triggered.fillna(False).sum()),
        }
        paths = {}
        paths["profile_qualification_audit"] = save_frame(audit.sort_values("Profile_number"), "profile_qualification_audit", output_dir, overwrite)
        paths["qualified_analysis_queue"] = save_frame(queues[PRIMARY_QUEUE].sort_values(["match_set_id", "is_scv"], ascending=[True, False]), "qualified_analysis_queue", output_dir, overwrite)
        paths["confirmed_mismatch_recheck"] = save_frame(correction.sort_values("association_row"), "confirmed_mismatch_recheck", output_dir, overwrite)
        paths["scv_matched_effects"] = save_frame(summary.sort_values(["queue", "outcome", "scope", "threshold_umol_kg"]), "scv_matched_effects", output_dir, overwrite)
        paths["near_zero_sensitivity"] = save_frame(near_zero.sort_values(["Profile_number", "threshold_umol_kg"]), "near_zero_sensitivity", output_dir, overwrite)
        paths["bootstrap_edge_case_tests"] = save_frame(edge, "bootstrap_edge_case_tests", output_dir, overwrite)
        paths["caliper_sensitivity_readout"] = save_frame(calipers.sort_values(["spatial_caliper_km", "scope", "threshold_umol_kg"]), "caliper_sensitivity_readout", output_dir, overwrite)
        manifest = {
            "status": "complete",
            "elapsed_seconds": round(time.time() - started, 3),
            "code_root": str(CODE_ROOT),
            "data_root": str(Path.cwd()),
            "track_import": str(Path(track.__file__).resolve()),
            "local_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=CODE_ROOT, text=True).strip(),
            "formal_cohort": str(FORMAL_COHORT),
            "counts": counts,
            "identity_status_counts": identity.current_date_position_status.value_counts().to_dict(),
            "loader_years": loader_provenance,
            "outputs": paths,
        }
        paths["manifest"] = save_json(manifest, "scv_matched_effects_manifest.json", output_dir, overwrite)
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
