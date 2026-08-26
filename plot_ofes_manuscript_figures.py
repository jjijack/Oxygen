"""Render the six-panel OFES manuscript figure package from completed results.

This module only assembles existing parquet/JSON products. It does not read
the OFES NetCDF fields, rerun event detection, or define new scientific
statistics. Inputs are explicit so the same figures can be regenerated in a
different result directory without copying an old run path into the code.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

RED = "#d62728"
BLUE = "#4c78a8"
TEAL = "#72b7b2"
PURPLE = "#7b61a8"
ORANGE = "#d95f59"
INK = "#111827"
GRAY = "#9ca3af"


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise KeyError(f"{label} lacks required columns: {missing}")


def _true_mask(values: pd.Series) -> pd.Series:
    return values.eq(True).fillna(False)  # noqa: E712 - explicit nullable-bool match


def _read_table(tables_root: Path, name: str) -> pd.DataFrame:
    path = tables_root / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_parquet(path)


def _read_json(path: Path) -> Mapping[str, Any]:
    return json.loads(path.expanduser().resolve().read_text())


def _style_axes(axes: Any) -> None:
    for axis in np.ravel(axes):
        axis.grid(alpha=0.2, zorder=0)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)


def _save(fig: Any, output_dir: Path, name: str, dpi: int) -> Path:
    path = output_dir / name
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def _global_categories(thresholds: pd.DataFrame) -> pd.DataFrame:
    """Create the mutually exclusive DO-evaluable McCoy SCV categories."""

    _require_columns(
        thresholds,
        ("stable_scv_id", "threshold_umol_kg", "do_evaluable", "has_delta_do", "lat", "lon"),
        "global McCoy threshold table",
    )
    base = thresholds.loc[
        thresholds["threshold_umol_kg"].eq(20)
        & _true_mask(thresholds["do_evaluable"])
    ].copy()
    # The published OR denominator preserves the legacy exclusion of P32111.
    base = base.loc[base["stable_scv_id"].ne("2017:32111")].copy()
    if len(base) != 244 or base["stable_scv_id"].nunique() != 244:
        raise ValueError(f"Expected 244 unique DO-evaluable SCVs, found {len(base)}")
    for threshold in (20, 35, 50):
        rows = thresholds.loc[thresholds["threshold_umol_kg"].eq(threshold)]
        lookup = rows.set_index("stable_scv_id")["has_delta_do"].fillna(False).astype(bool)
        base[f"has_{threshold}"] = base["stable_scv_id"].map(lookup).fillna(False).astype(bool)
    base["category"] = "below DO20"
    base.loc[base["has_20"], "category"] = "DO20–35 carriers"
    base.loc[base["has_35"], "category"] = "DO35–50 carriers"
    base.loc[base["has_50"], "category"] = "DO50 carriers"
    counts = base["category"].value_counts().to_dict()
    expected = {
        "DO50 carriers": 17,
        "DO35–50 carriers": 22,
        "DO20–35 carriers": 28,
        "below DO20": 177,
    }
    if counts != expected or sum(counts.values()) != 244:
        raise ValueError(f"Mutually exclusive global categories changed: {counts}")
    if int(base["has_50"].sum()) != 17 or int((base["has_35"] | base["has_50"]).sum()) != 39:
        raise ValueError("DO50/DO35+ category assertions failed")
    if int((base["has_20"] | base["has_35"] | base["has_50"]).sum()) != 67:
        raise ValueError("DO20+ category assertion failed")
    return base


def _plot_figure1(
    thresholds: pd.DataFrame,
    sweep: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> Path:
    """Render the global DO-evaluable SCV, OR, and KE/OFES panels."""

    categories = _global_categories(thresholds)
    _require_columns(
        sweep,
        ("threshold_umol_kg", "group", "odds_ratio_vs_all_argo", "odds_ratio_ci_low", "odds_ratio_ci_high"),
        "SCV threshold sweep",
    )
    scv = sweep.loc[sweep["group"].eq("DO-evaluable McCoy SCVs")].sort_values("threshold_umol_kg")
    meta = sweep.loc[sweep["group"].eq("META-matched DO-evaluable Argo")].sort_values("threshold_umol_kg")
    if len(scv) != 3 or len(meta) != 3:
        raise ValueError("The formal OR table must contain three SCV and META rows")

    colors = {
        "DO50 carriers": RED,
        "DO35–50 carriers": ORANGE,
        "DO20–35 carriers": PURPLE,
        "below DO20": "#d9dee7",
    }
    fig, axes = plt.subplots(
        1, 3, figsize=(17.2, 5.3),
        gridspec_kw={"width_ratios": (1.25, 0.95, 0.9)},
        constrained_layout=True,
    )

    # Panel a: mutually exclusive global categories.
    for label, color in colors.items():
        rows = categories.loc[categories["category"].eq(label)]
        axes[0].scatter(
            rows["lon"], rows["lat"], s=15 if "below" in label else 24,
            color=color, alpha=0.8 if "below" not in label else 0.45,
            linewidths=0, label=label,
        )
    axes[0].set_xlim(-180, 180)
    axes[0].set_ylim(-70, 70)
    axes[0].set_xlabel("Longitude (°)")
    axes[0].set_ylabel("Latitude (°)")
    axes[0].set_title("244 DO-evaluable McCoy SCVs")
    axes[0].legend(frameon=False, fontsize=7.5, loc="lower left")
    axes[0].add_patch(Rectangle((140, 25), 40, 20, fill=False, edgecolor=INK, linewidth=1.4))
    axes[0].text(141, 46, "KE", fontsize=8, color=INK)

    # Panel b: ORs are read directly from the formal threshold sweep.
    y = np.arange(3)[::-1]
    for yi, (_, row) in zip(y, scv.iterrows()):
        axes[1].errorbar(
            float(row["odds_ratio_vs_all_argo"]), yi,
            xerr=[[float(row["odds_ratio_vs_all_argo"] - row["odds_ratio_ci_low"])],
                  [float(row["odds_ratio_ci_high"] - row["odds_ratio_vs_all_argo"])]],
            fmt="o", color=RED, capsize=3, label="McCoy SCV" if yi == y[0] else None,
        )
    for yi, (_, row) in zip(y, meta.iterrows()):
        axes[1].errorbar(
            float(row["odds_ratio_vs_all_argo"]), yi + 0.16,
            xerr=[[float(row["odds_ratio_vs_all_argo"] - row["odds_ratio_ci_low"])],
                  [float(row["odds_ratio_ci_high"] - row["odds_ratio_vs_all_argo"])]],
            fmt="s", color=TEAL, capsize=3, label="META" if yi == y[0] else None,
        )
    axes[1].axvline(1, color=INK, linestyle="--", linewidth=1)
    axes[1].set_xscale("log")
    axes[1].set_yticks(y + 0.08, ["DO50", "DO35", "DO20"])
    axes[1].set_xlabel("Odds ratio versus all DO-evaluable Argo")
    axes[1].set_title("Among the tested eddy classes")
    axes[1].legend(frameon=False, fontsize=8)

    # Panel c: the KE concentration is descriptive, not a regional probability test.
    ke = categories.loc[
        categories["lon"].between(140, 180) & categories["lat"].between(25, 45)
    ]
    if int((categories["category"] == "DO50 carriers").sum()) != 17 or int(
        (ke["category"] == "DO50 carriers").sum()
    ) != 16:
        raise ValueError("The frozen KE concentration count changed")
    for label, color in colors.items():
        rows = ke.loc[ke["category"].eq(label)]
        axes[2].scatter(rows["lon"], rows["lat"], s=26, color=color, alpha=0.85, label=label)
    axes[2].add_patch(Rectangle((140, 25), 40, 20, fill=False, edgecolor=INK, linewidth=1.4))
    axes[2].add_patch(Rectangle((140, 30), 35, 12, fill=False, edgecolor=BLUE, linewidth=1.4, linestyle="--"))
    axes[2].text(141, 40.8, "OFES analysis area", fontsize=8, color=BLUE)
    axes[2].set_xlim(138, 181)
    axes[2].set_ylim(23, 47)
    axes[2].set_xlabel("Longitude (°)")
    axes[2].set_ylabel("Latitude (°)")
    axes[2].set_title("KE concentration and OFES domain")
    axes[2].text(
        0.03, 0.04, "16 of 17 DO50 carriers\noccurred in the Kuroshio Extension",
        transform=axes[2].transAxes, fontsize=8.5, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none"},
    )
    _style_axes(axes)
    fig.suptitle("Global DO-evaluable SCV association and its KE concentration", fontsize=14)
    return _save(fig, output_dir, "Figure1_global_scv_ke.png", dpi)


def _water_mass_fraction(summary: pd.DataFrame) -> pd.Series:
    wm = pd.to_numeric(summary["peak_water_mass_do_contrast"], errors="coerce").abs()
    heave = pd.to_numeric(summary["peak_heave_do_contribution"], errors="coerce").abs()
    return wm / (wm + heave)


def _plot_figure2(
    quality_summary: pd.DataFrame,
    mccoy_event: pd.DataFrame,
    mccoy_audit: Mapping[str, Any],
    output_dir: Path,
    dpi: int,
) -> Path:
    """Render quality-event water-mass decomposition and McCoy enrichment."""

    _require_columns(
        quality_summary,
        (
            "event_id", "peak_water_mass_do_contrast", "peak_heave_do_contribution",
            "peak_same_sigma_theta_contrast", "peak_same_sigma_salinity_contrast",
        ),
        "quality-event diagnostic summary",
    )
    _require_columns(
        mccoy_event,
        (
            "event_id", "event_profile_mccoy_compatible_fraction",
            "background_control_mccoy_compatible_fraction",
            "any_event_profile_mccoy_compatible",
        ),
        "McCoy event summary",
    )
    strict_ids = set(mccoy_event["event_id"])
    if len(quality_summary) != 161 or len(strict_ids) != 56:
        raise ValueError("Expected 161 quality events and 56 strict events")
    strict = quality_summary["event_id"].isin(strict_ids)
    fraction = _water_mass_fraction(quality_summary)
    strict_fraction = fraction.loc[strict]
    quality_fraction = fraction
    strict_dominated = int((strict_fraction > 0.5).fillna(False).sum())
    quality_dominated = int((quality_fraction > 0.5).fillna(False).sum())
    if strict_dominated != 54 or quality_dominated != 138:
        raise ValueError("Frozen water-mass dominated counts changed")
    if not np.isclose(float(strict_fraction.median()), 0.8621853205, atol=2e-3):
        raise ValueError("Strict water-mass fraction median changed")
    if not np.isclose(float(quality_fraction.median()), 0.7819122081, atol=2e-3):
        raise ValueError("Quality water-mass fraction median changed")

    event_fraction = mccoy_event["event_profile_mccoy_compatible_fraction"].to_numpy(float)
    control_fraction = mccoy_event["background_control_mccoy_compatible_fraction"].to_numpy(float)
    if int(_true_mask(mccoy_event["any_event_profile_mccoy_compatible"]).sum()) != int(
        mccoy_audit["any_event_profile_mccoy_compatible_count"]
    ):
        raise ValueError("McCoy event-compatible count disagrees with the audit JSON")

    fig, axes = plt.subplots(2, 2, figsize=(13.4, 10.0), constrained_layout=True)
    ax = axes[0, 0]
    ax.scatter(
        quality_summary.loc[~strict, "peak_water_mass_do_contrast"],
        quality_summary.loc[~strict, "peak_heave_do_contribution"],
        s=17, color=GRAY, alpha=0.45, label="quality events",
    )
    ax.scatter(
        quality_summary.loc[strict, "peak_water_mass_do_contrast"],
        quality_summary.loc[strict, "peak_heave_do_contribution"],
        s=27, color=RED, alpha=0.8, label="strict DO50 events",
    )
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.axvline(0, color=INK, linewidth=0.8)
    ax.set_xlabel("Peak water-mass DO contrast (µmol kg⁻¹)")
    ax.set_ylabel("Peak heave DO contribution (µmol kg⁻¹)")
    ax.set_title("Quality-event decomposition")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[0, 1]
    finite_quality = quality_fraction.dropna()
    finite_strict = strict_fraction.dropna()
    ax.hist(finite_quality, bins=np.linspace(0, 1, 18), color=GRAY, alpha=0.55, label="quality (161)")
    ax.hist(finite_strict, bins=np.linspace(0, 1, 18), color=RED, alpha=0.62, label="strict (56)")
    ax.axvline(float(strict_fraction.median()), color=RED, linestyle="--", linewidth=1.4)
    ax.axvline(float(quality_fraction.median()), color=INK, linestyle=":", linewidth=1.4)
    ax.set_xlabel("Absolute water-mass contribution fraction")
    ax.set_ylabel("Events")
    ax.set_xlim(0, 1)
    ax.set_title("Water-mass contribution dominates")
    ax.text(
        0.03, 0.96,
        f"strict median {strict_fraction.median():.1%}; {strict_dominated}/56 dominated\n"
        f"quality median {quality_fraction.median():.1%}; {quality_dominated}/161 dominated",
        transform=ax.transAxes, va="top", fontsize=8.5,
        bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none"},
    )
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[1, 0]
    ax.scatter(
        quality_summary.loc[~strict, "peak_same_sigma_theta_contrast"],
        quality_summary.loc[~strict, "peak_same_sigma_salinity_contrast"],
        s=17, color=GRAY, alpha=0.42,
    )
    ax.scatter(
        quality_summary.loc[strict, "peak_same_sigma_theta_contrast"],
        quality_summary.loc[strict, "peak_same_sigma_salinity_contrast"],
        s=27, color=PURPLE, alpha=0.8,
    )
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.axvline(0, color=INK, linewidth=0.8)
    ax.set_xlabel("Same-σ₀ θ contrast (°C)")
    ax.set_ylabel("Same-σ₀ salinity contrast")
    ax.set_title("Joint thermohaline state (θ–S)")
    ax.text(
        0.03, 0.04, "θ and S jointly express one thermohaline state;\nnot three independent tests",
        transform=ax.transAxes, fontsize=8.5, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none"},
    )

    ax = axes[1, 1]
    order = np.argsort(event_fraction - control_fraction)
    for rank, idx in enumerate(order):
        ax.plot(
            [control_fraction[idx], event_fraction[idx]], [rank, rank],
            color=GRAY, alpha=0.34, linewidth=0.8,
        )
    ax.scatter(control_fraction[order], np.arange(len(order)), color=BLUE, s=16, label="same-event control")
    ax.scatter(event_fraction[order], np.arange(len(order)), color=RED, s=16, label="event core")
    mean = float(mccoy_audit["event_equal_mean_pass_fraction_difference"])
    low, high = (float(x) for x in mccoy_audit["event_equal_mean_difference_bootstrap95"])
    ax.axvline(0, color=INK, linestyle="--", linewidth=0.8)
    ax.axvline(mean, color=RED, linewidth=2, alpha=0.7)
    ax.set_yticks([])
    ax.set_xlabel("McCoy-compatible profile fraction")
    ax.set_title("Event-core enrichment in McCoy-compatible profiles")
    ax.text(
        0.03, 0.04,
        f"19/56 events have ≥1 compatible profile\nΔ event-equal = {mean:.3f} [{low:.3f}, {high:.3f}]\n"
        f"one-sided p={mccoy_audit['event_equal_paired_wilcoxon_p']:.2g}; "
        f"two-sided p={2 * float(mccoy_audit['event_equal_paired_wilcoxon_p']):.2g}",
        transform=ax.transAxes, fontsize=8.5, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
    )
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _style_axes(axes)
    fig.suptitle(
        "OFES oxygen-anomaly cores are enriched in McCoy-compatible profile signatures",
        fontsize=14,
    )
    return _save(fig, output_dir, "Figure2_ofes_water_mass_mccoy.png", dpi)


def _plot_figure3(
    ventilation_group: pd.DataFrame,
    trajectory: pd.DataFrame,
    transition: pd.DataFrame,
    walong: pd.DataFrame,
    strict_count: int,
    output_dir: Path,
    dpi: int,
) -> Path:
    """Render paired ventilation contrasts and distinct vertical diagnostics."""

    _require_columns(
        ventilation_group,
        (
            "analysis_subset", "horizon_days", "control_group", "metric", "event_count",
            "event_equal_mean_difference", "bootstrap95_low", "bootstrap95_high",
            "paired_wilcoxon_p",
        ),
        "ventilation group comparison",
    )
    _require_columns(
        trajectory,
        (
            "event_id", "direction_evaluable", "observed_vertical_motion",
            "resolved_vertical_pathway_supported", "resolved_downward_pathway_supported",
            "resolved_upward_pathway_supported",
        ),
        "trajectory population classification",
    )
    _require_columns(transition, ("event_id", "resolved_downward"), "transition phase table")
    _require_columns(
        walong,
        ("event_id", "aggregation", "region", "core_w_along_mean", "ring_w_along_mean"),
        "w_along formal audit",
    )
    rows = ventilation_group.loc[
        ventilation_group["analysis_subset"].eq("all_events")
        & ventilation_group["horizon_days"].eq(30)
    ].copy()
    metrics = {
        "ever_direct_mld_contact_fraction": "direct MLD contact",
        "ever_near_mld_contact_fraction": "within 25 m of MLD",
        "ever_outcrop_opportunity_fraction": "isopycnal outcrop opportunity",
    }
    rows = rows.loc[rows["metric"].isin(metrics)]
    if set(rows["control_group"]) != {"hydrographic_control", "kinematic_control"}:
        raise ValueError("30-day ventilation table lacks both formal controls")

    direction = trajectory.loc[_true_mask(trajectory["direction_evaluable"])]
    downward = int(direction["observed_vertical_motion"].eq("downward").sum())
    upward = int(direction["observed_vertical_motion"].eq("upward").sum())
    if len(direction) != 24 or (downward, upward) != (21, 3):
        raise ValueError("Frozen displacement-classifiable counts changed")
    resolved = trajectory.loc[_true_mask(trajectory["resolved_vertical_pathway_supported"])]
    resolved_down = int(_true_mask(resolved["resolved_downward_pathway_supported"]).sum())
    resolved_up = int(_true_mask(resolved["resolved_upward_pathway_supported"]).sum())
    if len(resolved) != 18 or (resolved_down, resolved_up) != (15, 3):
        raise ValueError("Frozen strict pathway counts changed")
    pathway_ids = set(transition.loc[_true_mask(transition["resolved_downward"]), "event_id"])
    if len(pathway_ids) != 19:
        raise ValueError("Frozen w_along population subset changed")
    daily = walong.loc[
        walong["aggregation"].eq("daily_mean") & walong["event_id"].isin(pathway_ids)
    ]
    core = daily.loc[daily["region"].eq("core"), "core_w_along_mean"].dropna()
    ring = daily.loc[daily["region"].eq("ring"), "ring_w_along_mean"].dropna()
    if len(core) != 19 or len(ring) != 19:
        raise ValueError("w_along subset is incomplete")
    if not np.isclose(float(core.mean()), 4.8095, atol=0.02):
        raise ValueError("Formal w_along daily-mean value changed")

    fig, axes = plt.subplots(
        1, 2, figsize=(15.2, 6.0),
        gridspec_kw={"width_ratios": (1.2, 0.9)},
        constrained_layout=True,
    )
    ax = axes[0]
    colors = {"hydrographic_control": BLUE, "kinematic_control": TEAL}
    offsets = {"hydrographic_control": 0.16, "kinematic_control": -0.16}
    for idx, metric in enumerate(metrics):
        for control in ("hydrographic_control", "kinematic_control"):
            row = rows.loc[(rows["metric"].eq(metric)) & rows["control_group"].eq(control)]
            if len(row) != 1:
                raise ValueError(f"Missing 30-day ventilation row for {metric}/{control}")
            row = row.iloc[0]
            yi = idx + offsets[control]
            ax.errorbar(
                float(row["event_equal_mean_difference"]), yi,
                xerr=[[float(row["event_equal_mean_difference"] - row["bootstrap95_low"])],
                      [float(row["bootstrap95_high"] - row["event_equal_mean_difference"])]],
                fmt="o", color=colors[control], capsize=3,
                label=(
                    "hydrographic control (n=28)"
                    if control == "hydrographic_control" and idx == 0
                    else "kinematic control (n=27)"
                    if control == "kinematic_control" and idx == 0
                    else None
                ),
            )
    ax.axvline(0, color=INK, linestyle="--", linewidth=1)
    ax.set_yticks(np.arange(3), list(metrics.values()))
    ax.set_xlabel("Event-equal anomaly − control fraction")
    ax.set_title("Event-equal ventilation contrasts in the\ntrajectory-complete paired subset")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.text(
        0.03, 0.03,
        "Each event is contrasted first, then weighted equally.\n"
        "Metrics are correlated; controls share most anomaly events.\n"
        "The six points are not six independent tests.",
        transform=ax.transAxes, fontsize=8.3, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
    )

    ax = axes[1]
    box_specs = [
        (
            "Displacement-classifiable",
            f"{len(direction)}/{strict_count}",
            f"{downward} downward · {upward} upward",
            BLUE,
        ),
        (
            "Strict resolved pathways",
            str(len(resolved)),
            f"{resolved_down} downward · {resolved_up} upward",
            RED,
        ),
        (
            "w_along subset",
            "n=19",
            f"daily mean +{core.mean():.1f} m d⁻¹ · nominal p=0.040",
            PURPLE,
        ),
    ]
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(box_specs) - 0.5)
    ax.axis("off")
    for idx, (title, count, detail, color) in enumerate(box_specs):
        y = len(box_specs) - idx - 1
        ax.add_patch(
            Rectangle(
                (0.04, y - 0.34), 0.92, 0.68,
                facecolor=color, alpha=0.12, edgecolor=color, linewidth=1.4,
            )
        )
        ax.text(0.08, y + 0.08, title, fontsize=10, weight="bold", color=color)
        ax.text(0.08, y - 0.12, f"{count}   {detail}", fontsize=9.4, color=INK)
    ax.text(
        0.04, -0.48,
        "Related but non-independent pathway diagnostics with distinct eligibility criteria.\n"
        "The w_along p value is nominal and exploratory; boxes are not additive or strictly nested.",
        fontsize=8.3, va="top",
    )
    _style_axes(axes[:1])
    fig.suptitle("OFES ventilation history and resolved downward pathways", fontsize=14)
    return _save(fig, output_dir, "Figure3_ventilation_downward.png", dpi)


def _contiguous_spans(dates: pd.Series, mask: pd.Series) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    dates = pd.to_datetime(dates).reset_index(drop=True)
    mask = mask.reset_index(drop=True).astype(bool)
    spans: list[tuple[pd.Timestamp, pd.Timestamp]] = []
    start: pd.Timestamp | None = None
    previous: pd.Timestamp | None = None
    for date, active in zip(dates, mask):
        if active and start is None:
            start = date
        if not active and start is not None:
            spans.append((start, previous or start))
            start = None
        previous = date
    if start is not None:
        spans.append((start, previous or start))
    return spans


def _plot_figure4(ventilation_event: pd.DataFrame, output_dir: Path, dpi: int) -> Path:
    """Render the E000073 backward-integrated pathway in chronological display order."""

    _require_columns(
        ventilation_event,
        (
            "date", "particle_id", "particle_group", "depth_m", "lat", "lon", "mld_m",
            "particle_sigma0", "direct_mld_contact",
        ),
        "E000073 ventilation daily diagnostics",
    )
    data = ventilation_event.copy()
    data["date"] = pd.to_datetime(data["date"])
    data = data.sort_values(["particle_group", "particle_id", "date"]).reset_index(drop=True)
    groups = {"anomaly": RED, "hydrographic_control": BLUE}
    if set(groups) - set(data["particle_group"]):
        raise ValueError("E000073 case is missing an anomaly or hydrographic-control ensemble")
    rep_ids = {
        "anomaly": "z+0_r1_a270",
        "hydrographic_control": "hydrographic_control_z+0_r1_a270",
    }
    reps: dict[str, pd.DataFrame] = {}
    for group, rep_id in rep_ids.items():
        rep = data.loc[
            data["particle_group"].eq(group) & data["particle_id"].eq(rep_id)
        ].copy()
        if len(rep) != 51:
            raise ValueError(f"Expected 51 representative trajectory days for {group}, found {len(rep)}")
        reps[group] = rep.sort_values("date").reset_index(drop=True)

    anomaly_rep = reps["anomaly"]
    control_rep = reps["hydrographic_control"]
    direct_days = int(_true_mask(control_rep["direct_mld_contact"]).sum())
    if direct_days != 0:
        raise ValueError("Hydrographic control direct-contact count changed")
    anomaly_direct = _true_mask(anomaly_rep["direct_mld_contact"])
    direct_count = int(anomaly_direct.sum())
    direct_spans = _contiguous_spans(anomaly_rep["date"], anomaly_direct)
    initial_span = direct_spans[0]
    initial_days = int((initial_span[1] - initial_span[0]).days + 1)
    if (initial_days, direct_count) != (13, 20):
        raise ValueError("E000073 MLD contact counts changed")

    fig, axes = plt.subplots(2, 2, figsize=(14.4, 10.0), constrained_layout=True)
    peak_date = pd.to_datetime(anomaly_rep["date"]).max()

    # a. Horizontal ensemble trajectories, with one representative particle emphasized.
    ax = axes[0, 0]
    for group, color in groups.items():
        group_data = data.loc[data["particle_group"].eq(group)]
        for _, particle in group_data.groupby("particle_id", sort=False):
            particle = particle.sort_values("date")
            ax.plot(particle["lon"], particle["lat"], color=color, alpha=0.08, linewidth=0.45)
        rep = reps[group]
        ax.plot(rep["lon"], rep["lat"], color=color, linewidth=2.2, label=f"{group} representative")
        ax.scatter(rep["lon"].iloc[0], rep["lat"].iloc[0], color=color, s=30, marker="o")
        ax.scatter(rep["lon"].iloc[-1], rep["lat"].iloc[-1], color=color, s=55, marker="*", zorder=4)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title("E000073 horizontal reconstructed histories")
    ax.legend(frameon=False, fontsize=8)
    ax.text(
        0.03, 0.03,
        "Trajectories were initialized at the event peak and integrated backward;\n"
        "reconstructed histories are displayed in chronological order.",
        transform=ax.transAxes, fontsize=8.2, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
    )

    # b. Depth and MLD timeline with the exact contact intervals.
    ax = axes[0, 1]
    for group, color in groups.items():
        group_data = data.loc[data["particle_group"].eq(group)]
        summary = group_data.groupby("date")["depth_m"].agg(
            ["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]
        )
        summary.columns = ["median", "q25", "q75"]
        dates = pd.to_datetime(summary.index)
        ax.plot(dates, summary["median"], color=color, linewidth=2, label=f"{group} depth median")
        ax.fill_between(dates, summary["q25"], summary["q75"], color=color, alpha=0.14)
    mld = data.loc[data["particle_group"].eq("anomaly")].groupby("date")["mld_m"].agg(
        ["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]
    )
    mld.columns = ["median", "q25", "q75"]
    mld_dates = pd.to_datetime(mld.index)
    ax.plot(mld_dates, mld["median"], color=INK, linestyle="--", linewidth=1.5, label="anomaly MLD median")
    ax.fill_between(mld_dates, mld["q25"], mld["q75"], color=INK, alpha=0.08)
    for start, end in direct_spans:
        ax.axvspan(start, end + pd.Timedelta(days=1), color=RED, alpha=0.10)
    ax.axvline(pd.Timestamp("2003-02-07"), color=ORANGE, linestyle=":", linewidth=1.3)
    ax.text(
        pd.Timestamp("2003-02-07"), 0.93, "final detachment", rotation=90,
        fontsize=8, va="top", transform=ax.get_xaxis_transform(),
    )
    ax.invert_yaxis()
    ax.set_ylabel("Depth (m; positive downward)")
    ax.set_title("Depth, mixed-layer depth, and direct contact")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.text(
        0.03, 0.03,
        f"Initial continuous contact: {initial_days} d (Jan 1–13)\n"
        f"Cumulative direct contact: {direct_count} d\n"
        "Re-encounters: Jan 26 and Feb 1–6\n"
        f"Hydrographic control: {direct_days} of {len(control_rep)} trajectory days with direct MLD contact",
        transform=ax.transAxes, fontsize=8.1, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
    )

    # c. Density histories.
    ax = axes[1, 0]
    for group, color in groups.items():
        group_data = data.loc[data["particle_group"].eq(group)]
        summary = group_data.groupby("date")["particle_sigma0"].agg(
            ["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]
        )
        summary.columns = ["median", "q25", "q75"]
        dates = pd.to_datetime(summary.index)
        ax.plot(dates, summary["median"], color=color, linewidth=2, label=f"{group} σ₀ median")
        ax.fill_between(dates, summary["q25"], summary["q75"], color=color, alpha=0.14)
        rep = reps[group]
        ax.plot(rep["date"], rep["particle_sigma0"], color=color, linewidth=1.0, alpha=0.75)
    ax.axvline(peak_date, color=INK, linestyle="--", linewidth=1, label="event peak")
    ax.set_ylabel("σ₀ (kg m⁻³)")
    ax.set_title("Density-matched trajectories under event-peak conditions")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.text(
        0.03, 0.03,
        "The controls match the event-peak condition; their earliest reconstructed\n"
        "positions are not a common natural release date or depth.",
        transform=ax.transAxes, fontsize=8.2, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
    )

    # d. Density-depth phase space.
    ax = axes[1, 1]
    for group, color in groups.items():
        group_data = data.loc[data["particle_group"].eq(group)]
        summary = group_data.groupby("date").agg(
            sigma=("particle_sigma0", "median"), depth=("depth_m", "median")
        )
        ax.plot(summary["sigma"], summary["depth"], color=color, linewidth=2, label=f"{group} ensemble median")
        rep = reps[group]
        ax.plot(rep["particle_sigma0"], rep["depth_m"], color=color, linewidth=1.1, alpha=0.75)
        ax.scatter(rep["particle_sigma0"].iloc[0], rep["depth_m"].iloc[0], color=color, s=30)
        ax.scatter(rep["particle_sigma0"].iloc[-1], rep["depth_m"].iloc[-1], color=color, s=55, marker="*")
    ax.invert_yaxis()
    ax.set_xlabel("σ₀ (kg m⁻³)")
    ax.set_ylabel("Depth (m; positive downward)")
    ax.set_title("Density–depth phase space")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.text(
        0.03, 0.03,
        "Anomaly earliest reconstructed depth: "
        f"{anomaly_rep['depth_m'].iloc[0]:.0f} m\n"
        "Control earliest reconstructed depth: "
        f"{control_rep['depth_m'].iloc[0]:.0f} m\n"
        "Event-peak matched condition is at the common endpoint.",
        transform=ax.transAxes, fontsize=8.2, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
    )
    for axis in (axes[0, 1], axes[1, 0]):
        axis.xaxis.set_major_locator(mdates.AutoDateLocator())
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        axis.tick_params(axis="x", rotation=35)
    _style_axes(axes)
    fig.suptitle(
        "E000073: contrasting ventilation histories under a matched event-peak condition",
        fontsize=14,
    )
    return _save(fig, output_dir, "Figure4_E000073_case.png", dpi)


def _plot_figure5(
    lifecycle: pd.DataFrame,
    stage: pd.DataFrame,
    retention: pd.DataFrame,
    transition_audit: Mapping[str, Any],
    fixed_water_audit: Mapping[str, Any],
    moving_water_audit: Mapping[str, Any],
    output_dir: Path,
    dpi: int,
) -> Path:
    """Render nested rotational organization, McCoy stages, and retention."""

    _require_columns(
        lifecycle,
        (
            "event_id", "persistent_anticyclonic_rotational_carrier",
            "scv_compatible", "surface_obscured_scv_compatible",
        ),
        "lifecycle event summary",
    )
    _require_columns(stage, ("event_id", "stage", "event_compatible_count"), "per-stage McCoy diagnostics")
    _require_columns(
        retention, ("event_id", "persistent_carrier", "post_peak_decay_slope"), "retention comparison"
    )
    if len(lifecycle) != 56:
        raise ValueError("Expected 56 lifecycle events")
    carrier = int(_true_mask(lifecycle["persistent_anticyclonic_rotational_carrier"]).sum())
    scv = int(_true_mask(lifecycle["scv_compatible"]).sum())
    obscured = int(_true_mask(lifecycle["surface_obscured_scv_compatible"]).sum())
    if (carrier, scv, obscured) != (27, 6, 1):
        raise ValueError(f"Frozen rotational hierarchy changed: {(carrier, scv, obscured)}")
    stage_counts = stage.groupby("stage").apply(
        lambda frame: int(frame["event_compatible_count"].gt(0).sum()),
        include_groups=False,
    ).to_dict()
    if stage_counts != {"start": 10, "peak": 19, "last": 9}:
        raise ValueError(f"Frozen per-stage McCoy counts changed: {stage_counts}")
    if len(retention) != 56:
        raise ValueError("Retention comparison must contain 56 events")

    primary = transition_audit["retention_comparison"]["primary_persistent_carrier"]["post_peak_decay_slope"]
    fixed = fixed_water_audit["carrier_retention"]["wm_normalized_decay_slope"]
    moving = moving_water_audit["carrier_retention"]["wm_normalized_decay_slope"]
    forest_rows = [
        ("ΔDO proxy / fixed site", primary),
        ("water mass / fixed site", fixed),
        ("water mass / moving core", moving),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.7), constrained_layout=True)
    ax = axes[0]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        (0.05, 0.08, 0.90, 0.84, "56 strict events", INK),
        (0.16, 0.20, 0.68, 0.60, f"{carrier} persistent anticyclonic\nrotational carriers", RED),
        (0.28, 0.32, 0.44, 0.36, f"{scv} SCV-compatible", PURPLE),
        (0.38, 0.40, 0.24, 0.20, f"{obscured} surface-obscured\nSCV-compatible", TEAL),
    ]
    for x, y, width, height, label, color in boxes:
        ax.add_patch(Rectangle((x, y), width, height, fill=False, edgecolor=color, linewidth=2.0))
        ax.text(x + width / 2, y + height - 0.06, label, ha="center", va="top", fontsize=9.2, color=color)
    ax.text(
        0.06, -0.08,
        "Nested definitions: 6/6 are within the 27 carriers;\n"
        "1/1 is within the six SCV-compatible events.\n"
        "These are containment relations, not independent tests.",
        fontsize=8.3, va="bottom",
    )
    ax.set_ylim(-0.12, 1)
    ax.set_title("Rotational organization hierarchy")

    ax = axes[1]
    stages = ["start", "peak", "last"]
    counts = [stage_counts[s] for s in stages]
    bars = ax.bar(np.arange(3), counts, color=[BLUE, RED, BLUE], width=0.55)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, count + 0.5, f"{count}/56", ha="center", fontsize=10)
    ax.set_xticks(np.arange(3), ["start", "peak", "last"])
    ax.set_ylim(0, 24)
    ax.set_ylabel("Events with ≥1 compatible profile")
    ax.set_title("McCoy-compatible profile expression by stage")
    ax.text(
        0.03, 0.03,
        "McCoy-compatible profile expression is most concentrated\n"
        "at the anomaly peak. Fixed-peak sensitivity is supplementary.",
        transform=ax.transAxes, fontsize=8.3, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
    )

    ax = axes[2]
    y = np.arange(len(forest_rows))[::-1]
    for yi, (label, row) in zip(y, forest_rows):
        mean = float(row["median_diff"])
        low, high = (float(value) for value in row["bootstrap_ci"])
        ax.plot([low, high], [yi, yi], color=PURPLE, linewidth=3)
        ax.plot([low, low], [yi - 0.12, yi + 0.12], color=PURPLE, linewidth=1.2)
        ax.plot([high, high], [yi - 0.12, yi + 0.12], color=PURPLE, linewidth=1.2)
        ax.plot(mean, yi, "o", color=INK, markersize=7)
        ax.text(high + 0.004, yi, f"{mean:+.3f}", va="center", fontsize=8.5)
    ax.axvline(0, color=INK, linestyle="--", linewidth=1)
    ax.set_yticks(y, [row[0] for row in forest_rows])
    ax.set_xlabel("Carrier − non-carrier median decay-slope difference")
    ax.set_xlim(-0.06, 0.42)
    ax.set_title("Retention across reference frames")
    ax.text(
        0.03, 0.92,
        "Direction is consistent with a slower decay for carriers,\n"
        "but the tendency is statistically limited (CIs and MW tests remain in the report).",
        transform=ax.transAxes, fontsize=8.2, va="top",
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
    )
    _style_axes(axes[1:])
    fig.suptitle("Rotational organization and finite water-mass retention", fontsize=14)
    return _save(fig, output_dir, "Figure5_rotational_organization.png", dpi)


def render_manuscript_figures(
    tables_root: Path,
    global_thresholds: Path,
    global_sweep: Path,
    quality_summary: Path,
    mccoy_audit: Path,
    transition_audit: Path,
    fixed_water_audit: Path,
    moving_water_audit: Path,
    stage_diagnostics: Path,
    trajectory_ventilation_event: Path,
    output_dir: Path,
    dpi: int = 260,
) -> list[Path]:
    """Render Figures 1–5 from explicit formal result products."""

    tables_root = tables_root.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = pd.read_parquet(global_thresholds.expanduser().resolve())
    sweep = pd.read_parquet(global_sweep.expanduser().resolve())
    quality = pd.read_parquet(quality_summary.expanduser().resolve())
    mccoy_event = _read_table(tables_root, "mccoy_event_summary.parquet")
    lifecycle = _read_table(tables_root, "lifecycle_event_summary.parquet")
    retention = _read_table(tables_root, "transition_retention_comparison.parquet")
    trajectory = _read_table(tables_root, "trajectory3d_population_classification.parquet")
    transition = _read_table(tables_root, "transition_phase_table.parquet")
    ventilation_group = _read_table(tables_root, "ventilation_group_comparison.parquet")
    walong = _read_table(tables_root, "walong_formal_audit.parquet")
    stage = pd.read_parquet(stage_diagnostics.expanduser().resolve())
    ventilation_event = pd.read_parquet(trajectory_ventilation_event.expanduser().resolve())
    mccoy_summary = _read_json(mccoy_audit)
    transition_summary = _read_json(transition_audit)
    fixed_summary = _read_json(fixed_water_audit)
    moving_summary = _read_json(moving_water_audit)
    strict_count = len(mccoy_event)
    return [
        _plot_figure1(thresholds, sweep, output_dir, dpi),
        _plot_figure2(quality, mccoy_event, mccoy_summary, output_dir, dpi),
        _plot_figure3(ventilation_group, trajectory, transition, walong, strict_count, output_dir, dpi),
        _plot_figure4(ventilation_event, output_dir, dpi),
        _plot_figure5(
            lifecycle, stage, retention, transition_summary,
            fixed_summary, moving_summary, output_dir, dpi,
        ),
    ]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables-root", required=True, type=Path)
    parser.add_argument("--global-thresholds", required=True, type=Path)
    parser.add_argument("--global-sweep", required=True, type=Path)
    parser.add_argument("--quality-summary", required=True, type=Path)
    parser.add_argument("--mccoy-audit", required=True, type=Path)
    parser.add_argument("--transition-audit", required=True, type=Path)
    parser.add_argument("--fixed-water-audit", required=True, type=Path)
    parser.add_argument("--moving-water-audit", required=True, type=Path)
    parser.add_argument("--stage-diagnostics", required=True, type=Path)
    parser.add_argument("--trajectory-ventilation-event", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=260)
    args = parser.parse_args(argv)
    for path in render_manuscript_figures(
        args.tables_root,
        args.global_thresholds,
        args.global_sweep,
        args.quality_summary,
        args.mccoy_audit,
        args.transition_audit,
        args.fixed_water_audit,
        args.moving_water_audit,
        args.stage_diagnostics,
        args.trajectory_ventilation_event,
        args.output_dir,
        args.dpi,
    ):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
