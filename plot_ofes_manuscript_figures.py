"""Render the OFES manuscript Figures 1–5 from completed results.

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
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
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
FIGURE_WIDTH_IN = 7.1
PANEL_LABEL_SIZE = 9.5

plt.rcParams.update(
    {
        "font.size": 8.0,
        "axes.titlesize": 8.6,
        "axes.labelsize": 8.0,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "figure.titlesize": 10.2,
        "lines.linewidth": 1.0,
    }
)

# The two frames are deliberately kept distinct in the KE concentration
# panel: the KE frame is the descriptive regional window, while the OFES
# frame is the model delivery domain used by the process diagnostics.
KE_BOUNDS = (140.0, 180.0, 25.0, 45.0)
OFES_BOUNDS = (140.0, 170.0, 25.0, 45.0)


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


def _add_panel_labels(axes: Any) -> None:
    """Add standard manuscript panel labels in axes coordinates."""

    for label, axis in zip("abcdefghijklmnopqrstuvwxyz", np.ravel(axes)):
        axis.text(
            -0.08,
            1.03,
            label,
            transform=axis.transAxes,
            fontsize=PANEL_LABEL_SIZE,
            fontweight="bold",
            ha="left",
            va="bottom",
            clip_on=False,
            zorder=20,
        )


def _save(
    fig: Any,
    output_dir: Path,
    name: str,
    dpi: int,
    bottom: float = 0.22,
    top: float = 0.88,
    left: float = 0.12,
    wspace: float | None = None,
    hspace: float | None = None,
) -> Path:
    path = output_dir / name
    # Keep the requested physical width exact. The margins are reserved in the
    # figure canvas so outside panel labels are not recovered by tight cropping.
    kwargs = {"left": left, "right": 0.985, "top": top, "bottom": bottom}
    if wspace is not None:
        kwargs["wspace"] = wspace
    if hspace is not None:
        kwargs["hspace"] = hspace
    fig.subplots_adjust(**kwargs)
    fig.savefig(path, dpi=dpi, facecolor="white", edgecolor="none")
    plt.close(fig)
    return path


def _coastline_segments() -> list[np.ndarray]:
    """Read a cached Natural Earth coastline without downloading data."""

    try:
        import cartopy
        import shapefile
    except ImportError:
        return []
    coastline = (
        Path(cartopy.config["data_dir"])
        / "shapefiles"
        / "natural_earth"
        / "physical"
        / "ne_110m_coastline.shp"
    )
    if not coastline.exists():
        return []
    segments: list[np.ndarray] = []
    for shape in shapefile.Reader(str(coastline)).shapes():
        points = np.asarray(shape.points, dtype=float)
        starts = list(shape.parts) + [len(points)]
        for start, end in zip(starts[:-1], starts[1:]):
            part = points[start:end]
            if len(part) < 2:
                continue
            # Avoid drawing a spurious line across the map at the dateline.
            breaks = np.flatnonzero(np.abs(np.diff(part[:, 0])) > 180.0) + 1
            segments.extend(chunk for chunk in np.split(part, breaks) if len(chunk) >= 2)
    return segments


def _add_coastlines(axis: Any) -> None:
    """Overlay cached coastlines when the plotting environment provides them."""

    segments = _coastline_segments()
    if segments:
        axis.add_collection(
            LineCollection(
                segments,
                colors="#94a3b8",
                linewidths=0.45,
                alpha=0.9,
                zorder=1,
            )
        )


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
    threshold_order = [50, 35, 20]
    scv = (
        sweep.loc[sweep["group"].eq("DO-evaluable McCoy SCVs")]
        .set_index("threshold_umol_kg")
        .loc[threshold_order]
        .reset_index()
    )
    meta = (
        sweep.loc[sweep["group"].eq("META-matched DO-evaluable Argo")]
        .set_index("threshold_umol_kg")
        .loc[threshold_order]
        .reset_index()
    )
    if len(scv) != 3 or len(meta) != 3:
        raise ValueError("The formal OR table must contain three SCV and META rows")

    colors = {
        "DO50 carriers": RED,
        "DO35–50 carriers": ORANGE,
        "DO20–35 carriers": PURPLE,
        "below DO20": "#c7cfdb",
    }
    fig, axes = plt.subplots(
        1, 3, figsize=(FIGURE_WIDTH_IN, 2.75),
        gridspec_kw={"width_ratios": (1.30, 1.02, 0.98)},
    )

    # Panel a: mutually exclusive global categories.
    _add_coastlines(axes[0])
    for label, color in colors.items():
        rows = categories.loc[categories["category"].eq(label)]
        axes[0].scatter(
            rows["lon"], rows["lat"], s=15 if "below" in label else 24,
            color=color, alpha=0.82 if "below" not in label else 0.62,
            linewidths=0, label=f"{label}, n={len(rows)}",
        )
    axes[0].set_xlim(-180, 180)
    axes[0].set_ylim(-70, 70)
    axes[0].set_xlabel("Longitude (°)")
    axes[0].set_ylabel("Latitude (°)")
    axes[0].set_title("244 DO-evaluable McCoy SCVs")
    axes[0].add_patch(
        Rectangle(
            (KE_BOUNDS[0], KE_BOUNDS[2]),
            KE_BOUNDS[1] - KE_BOUNDS[0],
            KE_BOUNDS[3] - KE_BOUNDS[2],
            fill=False,
            edgecolor=INK,
            linewidth=1.4,
        )
    )
    axes[0].text(141, 46, "KE frame", fontsize=7.5, color=INK)

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
    axes[1].set_title("Threshold-dependent\nodds ratios")

    # Panel c: the KE concentration is descriptive, not a regional probability test.
    ke = categories.loc[
        categories["lon"].between(KE_BOUNDS[0], KE_BOUNDS[1])
        & categories["lat"].between(KE_BOUNDS[2], KE_BOUNDS[3])
    ]
    if int((categories["category"] == "DO50 carriers").sum()) != 17 or int(
        (ke["category"] == "DO50 carriers").sum()
    ) != 16:
        raise ValueError("The frozen KE concentration count changed")
    for label, color in colors.items():
        rows = ke.loc[ke["category"].eq(label)]
        axes[2].scatter(rows["lon"], rows["lat"], s=26, color=color, alpha=0.85, label=label)
    _add_coastlines(axes[2])
    axes[2].add_patch(
        Rectangle(
            (KE_BOUNDS[0], KE_BOUNDS[2]),
            KE_BOUNDS[1] - KE_BOUNDS[0],
            KE_BOUNDS[3] - KE_BOUNDS[2],
            fill=False,
            edgecolor=INK,
            linewidth=1.4,
        )
    )
    axes[2].add_patch(
        Rectangle(
            (OFES_BOUNDS[0], OFES_BOUNDS[2]),
            OFES_BOUNDS[1] - OFES_BOUNDS[0],
            OFES_BOUNDS[3] - OFES_BOUNDS[2],
            fill=False,
            edgecolor=BLUE,
            linewidth=1.4,
            linestyle="--",
        )
    )
    axes[2].set_xlim(138, 181)
    axes[2].set_ylim(23, 47)
    axes[2].set_xlabel("Longitude (°)")
    axes[2].set_ylabel("Latitude (°)")
    axes[2].set_title("KE concentration\nand OFES domain")
    fig.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="none", color=color, label=f"{label}, n={len(categories.loc[categories['category'].eq(label)])}")
            for label, color in colors.items()
        ]
        + [
            Line2D([], [], marker="o", linestyle="none", color=RED, label="McCoy SCV"),
            Line2D([], [], marker="s", linestyle="none", color=TEAL, label="META"),
            Line2D([], [], color=INK, linewidth=1.4, label="KE frame"),
            Line2D([], [], color=BLUE, linewidth=1.4, linestyle="--", label="OFES domain"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.55, 0.02),
        ncol=4,
        frameon=False,
        fontsize=7.5,
        handletextpad=0.3,
        columnspacing=0.8,
    )
    _style_axes(axes)
    _add_panel_labels(axes)
    return _save(fig, output_dir, "Figure1_global_scv_ke.png", dpi, bottom=0.31)


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

    fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH_IN, 5.75))
    ax = axes[0, 0]
    ax.scatter(
        quality_summary.loc[~strict, "peak_water_mass_do_contrast"],
        quality_summary.loc[~strict, "peak_heave_do_contribution"],
        s=17, color=GRAY, alpha=0.45, label="other quality events (n=105)",
    )
    ax.scatter(
        quality_summary.loc[strict, "peak_water_mass_do_contrast"],
        quality_summary.loc[strict, "peak_heave_do_contribution"],
        s=27, color=RED, alpha=0.8, label="strict DO50 events (n=56)",
    )
    x_limits = ax.get_xlim()
    y_limits = ax.get_ylim()
    boundary_x = np.linspace(x_limits[0], x_limits[1], 400)
    ax.fill_between(
        boundary_x,
        -np.abs(boundary_x),
        np.abs(boundary_x),
        color=RED,
        alpha=0.045,
        zorder=-1,
    )
    ax.plot(boundary_x, boundary_x, color=GRAY, linestyle="--", linewidth=0.9)
    ax.plot(boundary_x, -boundary_x, color=GRAY, linestyle="--", linewidth=0.9)
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.axvline(0, color=INK, linewidth=0.8)
    ax.set_xlabel("Peak water-mass DO contrast (µmol kg⁻¹)")
    ax.set_ylabel("Peak heave DO contribution (µmol kg⁻¹)")
    ax.set_title("Water-mass versus heave DO contrast")

    ax = axes[0, 1]
    finite_quality = quality_fraction.dropna()
    finite_strict = strict_fraction.dropna()
    ax.hist(
        finite_quality,
        bins=np.linspace(0, 1, 18),
        color=GRAY,
        alpha=0.55,
        label="quality catalogue (n=161; includes strict)",
    )
    ax.hist(
        finite_strict,
        bins=np.linspace(0, 1, 18),
        color=RED,
        alpha=0.62,
        label="strict DO50 subset (n=56)",
    )
    ax.set_xlabel("Absolute water-mass contribution fraction")
    ax.set_ylabel("Events")
    ax.set_xlim(0, 1)
    ax.set_title("Water-mass contribution fraction")
    ax.axvline(
        float(strict_fraction.median()), color=RED, linestyle="--", linewidth=1.4,
        label=f"strict median ({strict_fraction.median():.1%})",
    )
    ax.axvline(
        float(quality_fraction.median()), color=INK, linestyle=":", linewidth=1.4,
        label=f"quality median ({quality_fraction.median():.1%})",
    )

    ax = axes[1, 0]
    ax.scatter(
        quality_summary.loc[~strict, "peak_same_sigma_theta_contrast"],
        quality_summary.loc[~strict, "peak_same_sigma_salinity_contrast"],
        s=17, color=GRAY, alpha=0.42,
    )
    ax.scatter(
        quality_summary.loc[strict, "peak_same_sigma_theta_contrast"],
        quality_summary.loc[strict, "peak_same_sigma_salinity_contrast"],
        s=27, color=RED, alpha=0.8,
    )
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.axvline(0, color=INK, linewidth=0.8)
    ax.set_xlabel("Same-σ₀ θ contrast (°C)")
    ax.set_ylabel("Same-σ₀ salinity contrast (psu)")
    ax.set_title("Same-σ₀ thermohaline contrasts")

    ax = axes[1, 1]
    order = np.argsort(event_fraction - control_fraction)
    for rank, idx in enumerate(order):
        ax.plot(
            [control_fraction[idx], event_fraction[idx]], [rank, rank],
            color=GRAY, alpha=0.34, linewidth=0.8,
        )
    ax.scatter(control_fraction[order], np.arange(len(order)), color=BLUE, s=16, label="same-event control")
    ax.scatter(event_fraction[order], np.arange(len(order)), color=RED, s=16, label="event core")
    ax.axvline(0, color=INK, linestyle="--", linewidth=0.8)
    ax.set_yticks([])
    ax.set_xlabel("McCoy-compatible profile fraction")
    ax.set_title("McCoy-compatible profile fractions")
    fig.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="none", color=GRAY, label="other quality events (n=105)"),
            Line2D([], [], marker="o", linestyle="none", color=RED, label="strict DO50 events (n=56)"),
            Line2D([], [], color=GRAY, linewidth=7, alpha=0.55, label="quality catalogue (n=161)"),
            Line2D([], [], color=RED, linewidth=7, alpha=0.62, label="strict DO50 subset (n=56)"),
            Line2D([], [], color=RED, linestyle="--", label=f"strict median ({strict_fraction.median():.1%})"),
            Line2D([], [], color=INK, linestyle=":", label=f"quality median ({quality_fraction.median():.1%})"),
            Line2D([], [], marker="o", linestyle="none", color=BLUE, label="same-event control"),
            Line2D([], [], marker="o", linestyle="none", color=RED, label="event core"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.52, 0.98),
        ncol=4,
        frameon=False,
        fontsize=7.5,
        handletextpad=0.3,
        columnspacing=0.8,
    )
    _style_axes(axes)
    _add_panel_labels(axes)
    return _save(fig, output_dir, "Figure2_ofes_water_mass_mccoy.png", dpi, top=0.75, hspace=0.85)


def _plot_figure3(
    ventilation_group: pd.DataFrame,
    trajectory: pd.DataFrame,
    transition: pd.DataFrame,
    walong: pd.DataFrame,
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
        "ever_direct_mld_contact_fraction": "direct MLD",
        "ever_near_mld_contact_fraction": "near MLD",
        "ever_outcrop_opportunity_fraction": "outcrop",
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
    w_along_downward = int(core.gt(0).sum())
    if w_along_downward != 13:
        raise ValueError("Frozen downward w_along count changed")

    fig, axes = plt.subplots(
        1, 2, figsize=(FIGURE_WIDTH_IN, 3.25),
        gridspec_kw={"width_ratios": (1.9, 0.9)},
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
    ax.set_ylim(-0.28, 2.68)
    ax.set_xlabel("Event-equal anomaly − control fraction")
    ax.set_title("30-day event-equal ventilation contrasts")
    fig.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="none", color=BLUE, label="hydrographic control (n=28)"),
            Line2D([], [], marker="o", linestyle="none", color=TEAL, label="kinematic control (n=27)"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.55, 0.98),
        ncol=2,
        frameon=False,
        fontsize=7.5,
        handletextpad=0.3,
        columnspacing=0.8,
    )

    ax = axes[1]
    light_gray = "#d1d5db"
    composition_rows = [
        (
            "Displacement classification\nn=56",
            [("downward", downward), ("upward", upward), ("unclassified", 32)],
        ),
        (
            "Strict resolved pathways\nn=18",
            [("downward", resolved_down), ("upward", resolved_up)],
        ),
        (
            "$w_{\\mathrm{along}}$ subset\nn=19",
            [("downward", w_along_downward), ("negative $w_{\\mathrm{along}}$", 6)],
        ),
    ]
    category_colors = {"downward": RED, "upward": BLUE, "unclassified": light_gray, "negative $w_{\\mathrm{along}}$": light_gray}
    y_positions = np.arange(len(composition_rows))[::-1]
    ax.set_xlim(0, 1.90)
    ax.set_ylim(-0.5, 2.6)
    for y, (row_label, components) in zip(y_positions, composition_rows):
        total = sum(count for _, count in components)
        if total <= 0:
            raise ValueError(f"Invalid direction composition for {row_label}")
        ax.text(0.0, y + 0.30, row_label, ha="left", va="bottom", fontsize=6.8)
        left = 0.0
        for category, count in components:
            width = count / total
            ax.barh(
                y, width, left=left, height=0.42,
                color=category_colors[category], edgecolor="white", linewidth=0.7,
            )
            text_color = "white" if category in {"downward", "upward"} else INK
            ax.text(left + width / 2, y, str(count), ha="center", va="center", fontsize=7.5, color=text_color, weight="bold")
            left += width
        if row_label.startswith("$w_"):
            ax.text(
                1.02, y, f"+{core.mean():.1f} m d⁻¹\npositive =\ndownward",
                ha="left", va="center", fontsize=6.0, color=INK,
            )
    ax.set_yticks([])
    ax.set_xticks([0, 0.5, 1.0], ["0", "50%", "100%"])
    ax.set_xlabel("Complete row")
    ax.grid(axis="x", alpha=0.2, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    fig.legend(
        handles=[
            Line2D([], [], marker="s", linestyle="none", markersize=7, color=RED, label="downward"),
            Line2D([], [], marker="s", linestyle="none", markersize=7, color=BLUE, label="upward"),
            Line2D([], [], marker="s", linestyle="none", markersize=7, color=light_gray, label="unclassified / other"),
        ],
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.72, 0.04),
        ncol=3,
        fontsize=6.8,
        handletextpad=0.25,
        columnspacing=0.55,
    )
    ax.set_title("Direction-resolved pathways", loc="right", pad=4, fontsize=7.8)
    _style_axes(axes[:1])
    _add_panel_labels(axes)
    return _save(fig, output_dir, "Figure3_ventilation_downward.png", dpi, top=0.75, bottom=0.28, left=0.18, wspace=0.30)


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
    if "event_id" in data.columns:
        data = data.loc[data["event_id"].eq("OFES_DO50_E000073")].copy()
    data["date"] = pd.to_datetime(data["date"])
    data = data.sort_values(["particle_group", "particle_id", "date"]).reset_index(drop=True)
    groups = {"anomaly": RED, "hydrographic_control": BLUE}
    group_labels = {"anomaly": "anomaly", "hydrographic_control": "hydrographic control"}
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

    fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH_IN, 5.8))
    peak_date = pd.to_datetime(anomaly_rep["date"]).max()

    # a. Horizontal ensemble trajectories, with one representative particle emphasized.
    ax = axes[0, 0]
    for group, color in groups.items():
        group_data = data.loc[data["particle_group"].eq(group)]
        for _, particle in group_data.groupby("particle_id", sort=False):
            particle = particle.sort_values("date")
            ax.plot(particle["lon"], particle["lat"], color=color, alpha=0.08, linewidth=0.45)
        rep = reps[group]
        ax.plot(
            rep["lon"],
            rep["lat"],
            color=color,
            linewidth=2.2,
            label=group_labels[group],
        )
        ax.scatter(rep["lon"].iloc[0], rep["lat"].iloc[0], color=color, s=30, marker="o")
        ax.scatter(rep["lon"].iloc[-1], rep["lat"].iloc[-1], color=color, s=55, marker="*", zorder=4)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title("Horizontal trajectories")

    # b. Depth and MLD timeline with the exact contact intervals.
    ax = axes[0, 1]
    for group, color in groups.items():
        group_data = data.loc[data["particle_group"].eq(group)]
        summary = group_data.groupby("date")["depth_m"].agg(
            ["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]
        )
        summary.columns = ["median", "q25", "q75"]
        dates = pd.to_datetime(summary.index)
        ax.plot(
            dates,
            summary["median"],
            color=color,
            linewidth=2,
            label=f"{group_labels[group]} depth",
        )
        ax.fill_between(dates, summary["q25"], summary["q75"], color=color, alpha=0.14)
    for group, color, linestyle in (
        ("anomaly", INK, "--"),
        ("hydrographic_control", BLUE, ":"),
    ):
        mld = data.loc[data["particle_group"].eq(group)].groupby("date")["mld_m"].agg(
            ["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]
        )
        mld.columns = ["median", "q25", "q75"]
        mld_dates = pd.to_datetime(mld.index)
        ax.plot(
            mld_dates,
            mld["median"],
            color=color,
            linestyle=linestyle,
            linewidth=1.4,
            label=f"{group_labels[group]}-path MLD",
        )
        ax.fill_between(mld_dates, mld["q25"], mld["q75"], color=color, alpha=0.06)
    for start, end in direct_spans:
        ax.axvspan(start, end + pd.Timedelta(days=1), color=RED, alpha=0.10)
    ax.invert_yaxis()
    ax.set_ylabel("Depth (m; positive downward)")
    ax.set_title("Depth / MLD", loc="right", pad=4)
    ax.text(
        0.02,
        1.09,
        "MLD contact: 13 d initial; 20 d total",
        transform=ax.transAxes,
        fontsize=7.5,
        va="bottom",
        clip_on=False,
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
        ax.plot(
            dates,
            summary["median"],
            color=color,
            linewidth=2,
            label=f"{group_labels[group]} σ₀ median",
        )
        ax.fill_between(dates, summary["q25"], summary["q75"], color=color, alpha=0.14)
        rep = reps[group]
        ax.plot(rep["date"], rep["particle_sigma0"], color=color, linewidth=1.0, alpha=0.75)
    ax.axvline(peak_date, color=INK, linestyle="--", linewidth=1, label="event peak")
    ax.set_ylabel("σ₀ (kg m⁻³)")
    ax.set_title("Density histories")

    # d. Density-depth phase space.
    ax = axes[1, 1]
    for group, color in groups.items():
        group_data = data.loc[data["particle_group"].eq(group)]
        summary = group_data.groupby("date").agg(
            sigma=("particle_sigma0", "median"), depth=("depth_m", "median")
        )
        ax.plot(
            summary["sigma"],
            summary["depth"],
            color=color,
            linewidth=2,
            label=f"{group_labels[group]} ensemble median",
        )
        rep = reps[group]
        ax.plot(rep["particle_sigma0"], rep["depth_m"], color=color, linewidth=1.1, alpha=0.75)
        ax.scatter(rep["particle_sigma0"].iloc[0], rep["depth_m"].iloc[0], color=color, s=30)
        ax.scatter(rep["particle_sigma0"].iloc[-1], rep["depth_m"].iloc[-1], color=color, s=55, marker="*")
    ax.invert_yaxis()
    ax.set_xlabel("σ₀ (kg m⁻³)")
    ax.set_ylabel("Depth (m; positive downward)")
    ax.set_title("Density–depth phase space")
    fig.legend(
        handles=[
            Line2D([], [], color=RED, linewidth=2.2, label="anomaly"),
            Line2D([], [], color=BLUE, linewidth=2.2, label="control"),
            Line2D([], [], marker="o", color=INK, linestyle="none", markersize=4, label="earliest position"),
            Line2D([], [], marker="*", color=INK, linestyle="none", markersize=7, label="event peak"),
            Line2D([], [], color=RED, linewidth=2.0, label="anomaly depth"),
            Line2D([], [], color=BLUE, linewidth=2.0, label="control depth"),
            Line2D([], [], color=INK, linestyle="--", linewidth=1.4, label="anomaly-path MLD"),
            Line2D([], [], color=BLUE, linestyle=":", linewidth=1.4, label="control-path MLD"),
            Line2D([], [], color=RED, linewidth=2.0, label="anomaly σ₀ median"),
            Line2D([], [], color=BLUE, linewidth=2.0, label="control σ₀ median"),
            Line2D([], [], color=RED, linewidth=2.0, label="anomaly ensemble median"),
            Line2D([], [], color=BLUE, linewidth=2.0, label="control ensemble median"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.52, 0.98),
        ncol=4,
        frameon=False,
        fontsize=7.5,
        handletextpad=0.3,
        columnspacing=0.8,
    )
    for axis in (axes[0, 1], axes[1, 0]):
        axis.xaxis.set_major_locator(mdates.AutoDateLocator())
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        axis.tick_params(axis="x", rotation=35)
    _style_axes(axes)
    _add_panel_labels(axes)
    return _save(fig, output_dir, "Figure4_E000073_case.png", dpi, top=0.68, bottom=0.22, hspace=0.55)


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
        (r"$\frac{\Delta\mathrm{DO\ proxy}}{\mathrm{fixed\ site}}$", primary),
        (r"$\frac{\mathrm{water\ mass}}{\mathrm{fixed\ site}}$", fixed),
        (r"$\frac{\mathrm{water\ mass}}{\mathrm{moving\ core}}$", moving),
    ]

    fig, axes = plt.subplots(
        1, 3,
        figsize=(FIGURE_WIDTH_IN, 3.15),
        gridspec_kw={"width_ratios": (1.0, 0.95, 1.15)},
    )
    ax = axes[0]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        (0.05, 0.08, 0.90, 0.84, "56 strict events", INK),
        (0.16, 0.20, 0.68, 0.60, f"{carrier} persistent carriers", RED),
        (0.28, 0.32, 0.44, 0.36, f"{scv} SCV-compatible", PURPLE),
        (0.38, 0.40, 0.24, 0.20, f"{obscured} obscured SCV", TEAL),
    ]
    for index, (x, y, width, height, label, color) in enumerate(boxes):
        ax.add_patch(Rectangle((x, y), width, height, fill=False, edgecolor=color, linewidth=2.0))
        ax.text(
            x + 0.015,
            y + height - 0.025,
            label,
            ha="left",
            va="top",
            fontsize=7.5,
            color=color,
            clip_on=False,
        )
    ax.set_ylim(-0.12, 1)
    ax.set_title("Hierarchy", loc="right")

    ax = axes[1]
    stages = ["start", "peak", "last"]
    counts = [stage_counts[s] for s in stages]
    bars = ax.bar(np.arange(3), counts, color=[BLUE, RED, BLUE], width=0.55)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, count + 0.5, f"{count}/56", ha="center", fontsize=10)
    ax.set_xticks(np.arange(3), ["start", "peak", "last"])
    ax.set_ylim(0, 24)
    ax.set_ylabel("Events with ≥1 compatible profile")
    ax.set_title("McCoy stage counts", loc="right")

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
    ax.set_xlabel("Median decay-slope difference")
    ax.set_xlim(-0.06, 0.50)
    ax.set_title("Decay-slope contrasts")
    ax.tick_params(axis="y", labelsize=8.5, pad=6)
    _style_axes(axes[1:])
    _add_panel_labels(axes)
    return _save(fig, output_dir, "Figure5_rotational_organization.png", dpi, left=0.06, wspace=0.48)


def render_manuscript_figures_from_results(
    thresholds: pd.DataFrame,
    sweep: pd.DataFrame,
    quality: pd.DataFrame,
    mccoy_event: pd.DataFrame,
    mccoy_summary: Mapping[str, Any],
    ventilation_group: pd.DataFrame,
    trajectory: pd.DataFrame,
    transition: pd.DataFrame,
    walong: pd.DataFrame,
    ventilation_event: pd.DataFrame,
    lifecycle: pd.DataFrame,
    stage: pd.DataFrame,
    retention: pd.DataFrame,
    transition_summary: Mapping[str, Any],
    fixed_water_summary: Mapping[str, Any],
    moving_water_summary: Mapping[str, Any],
    output_dir: Path,
    dpi: int = 600,
) -> list[Path]:
    """Render Figures 1–5 from loaded formal producer results."""

    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    return [
        _plot_figure1(thresholds, sweep, output_dir, dpi),
        _plot_figure2(quality, mccoy_event, mccoy_summary, output_dir, dpi),
        _plot_figure3(
            ventilation_group, trajectory, transition, walong,
            output_dir, dpi,
        ),
        _plot_figure4(ventilation_event, output_dir, dpi),
        _plot_figure5(
            lifecycle, stage, retention, transition_summary,
            fixed_water_summary, moving_water_summary, output_dir, dpi,
        ),
    ]


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
    dpi: int = 600,
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
    return render_manuscript_figures_from_results(
        thresholds,
        sweep,
        quality,
        mccoy_event,
        mccoy_summary,
        ventilation_group,
        trajectory,
        transition,
        walong,
        ventilation_event,
        lifecycle,
        stage,
        retention,
        transition_summary,
        fixed_summary,
        moving_summary,
        output_dir,
        dpi,
    )


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
    parser.add_argument("--dpi", type=int, default=600)
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
