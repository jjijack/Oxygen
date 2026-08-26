"""Render the final OFES surface-expression manuscript figure.

The script consumes the completed surface-eddy parquet products and the
audited summary JSON. It does not read OFES NetCDF files or rerun the PET
detector. All input and output locations are explicit command-line arguments
so a result directory can be regenerated in place.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd

BLUE = "#4c78a8"
TEAL = "#72b7b2"
PURPLE = "#7b61a8"
RED = "#d62728"
INK = "#111827"


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise KeyError(f"{label} lacks required columns: {missing}")


def _true_mask(values: pd.Series) -> pd.Series:
    return values.eq(True).fillna(False)  # noqa: E712 - explicit nullable-bool match


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


def _load_outputs(
    surface_root: Path,
    summary_path: Path | None,
) -> tuple[Mapping[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read only the completed PET products required by Figure 6."""

    summary_file = summary_path or (surface_root / "surface_eddy_summary.json")
    summary = json.loads(summary_file.read_text())
    association = pd.read_parquet(surface_root / "surface_eddy_event_association.parquet")
    quality = pd.read_parquet(surface_root / "surface_eddy_quality_eligible_161.parquet")
    rotation = pd.read_parquet(surface_root / "surface_eddy_surface_ro_crosstab.parquet")
    return summary, association, quality, rotation


def _plot_case_panel(
    axis: Any,
    case: pd.Series,
    association: pd.Series,
    case_map: Path | None = None,
    case_map_extent: Sequence[float] | None = None,
    case_map_crop: Sequence[int] | None = None,
) -> None:
    """Draw a data-faithful local schematic for the E000002 surface case.

    The formal tables retain the event core, core-weighted surface Rossby
    number, and nearest PET distance, but not gridded SSH values or contour
    vertices. If an already-rendered audited case map is supplied, it is used
    as a visual backdrop only; no numerical values are inferred from it.
    The surface-Ro mark is a core point with a colored halo, not an inferred
    spatial footprint. Numerical annotations always come from parquet rows.
    """

    lon = float(case["peak_lon"])
    lat = float(case["peak_lat"])
    ro = float(case["surface_core_weighted_rossby_number"])
    distance = float(association["nearest_pet_center_distance_km"])
    radius_ratio = float(association["nearest_pet_center_distance_over_effective_radius"])

    if case_map is not None:
        if case_map_extent is None or len(case_map_extent) != 4:
            raise ValueError("--case-map requires --case-map-extent lon_min lon_max lat_min lat_max")
        image = plt.imread(case_map.expanduser().resolve())
        if case_map_crop is not None:
            if len(case_map_crop) != 4:
                raise ValueError("--case-map-crop requires x0 y0 x1 y1")
            x0, y0, x1, y1 = (int(value) for value in case_map_crop)
            image = image[y0:y1, x0:x1]
        axis.imshow(
            image,
            extent=tuple(float(value) for value in case_map_extent),
            origin="upper",
            aspect="auto",
        )
        axis.set_xlim(float(case_map_extent[0]), float(case_map_extent[1]))
        axis.set_ylim(float(case_map_extent[2]), float(case_map_extent[3]))
        axis.set_facecolor("#f4f7fb")
        axis.set_title("E000002: audited PET contour and surface Ro response")
    else:
        span = 4.0
        axis.set_xlim(lon - span, lon + span)
        axis.set_ylim(lat - span * 0.75, lat + span * 0.75)
        axis.set_facecolor("#f4f7fb")
        axis.set_title("E000002: surface Ro response without closed-SSH containment")
    axis.set_xlabel("Longitude (°E)")
    axis.set_ylabel("Latitude (°N)")

    patch_color = BLUE if ro < 0 else RED
    axis.scatter(
        [lon],
        [lat],
        s=330,
        facecolors="none",
        edgecolors=patch_color,
        linewidths=2.2,
        alpha=0.9,
        zorder=4,
    )
    axis.scatter(
        [lon],
        [lat],
        s=82,
        color=patch_color,
        edgecolor=INK,
        linewidth=0.8,
        zorder=5,
        label="core-weighted surface Ro",
    )
    axis.scatter(
        [lon],
        [lat],
        marker="*",
        s=180,
        color=INK,
        zorder=5,
        label="deep rotational core",
    )

    if case_map is None:
        # The table has no contour vertices; use the audited nearest-centre
        # distance only as a locator for the formal no-containment result.
        outline_lon = lon + min(3.0, max(1.5, distance / 111.0))
        outline_lat = lat + 1.0
        axis.add_patch(
            Ellipse(
                (outline_lon, outline_lat),
                width=1.35,
                height=0.9,
                fill=False,
                edgecolor=PURPLE,
                linewidth=2.0,
                linestyle="--",
                label="closed PET/SSH contour (nearest audited object)",
            )
        )
        axis.annotate(
            "core outside closed contour",
            xy=(lon, lat),
            xytext=(lon - 3.5, lat - 2.1),
            arrowprops={"arrowstyle": "->", "color": INK, "lw": 0.9},
            fontsize=8.5,
        )
    else:
        axis.annotate(
            "deep rotational core\n(core-weighted surface Ro)",
            xy=(lon, lat),
            xytext=(lon - 2.6, lat - 2.0),
            arrowprops={"arrowstyle": "->", "color": INK, "lw": 0.9},
            fontsize=8.5,
        )
    axis.text(
        0.02,
        0.03,
        f"core-weighted surface Ro = {ro:+.3f}\nPET containment: False\n"
        f"nearest centre = {distance:.0f} km ({radius_ratio:.2f}R)",
        transform=axis.transAxes,
        fontsize=8.5,
        va="bottom",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    if case_map is None:
        axis.legend(frameon=False, fontsize=7.5, loc="upper left")


def _plot_figure6(
    summary: Mapping[str, Any],
    association: pd.DataFrame,
    quality: pd.DataFrame,
    rotation: pd.DataFrame,
    output_dir: Path,
    dpi: int,
    case_map: Path | None = None,
    case_map_extent: Sequence[float] | None = None,
    case_map_crop: Sequence[int] | None = None,
) -> Path:
    """Render Figure 6: surface Ro correspondence versus closed-SSH capture."""

    _require_columns(
        association,
        (
            "event_id",
            "peak_core_analysis_eligible",
            "peak_core_contained_by_actual_pet_effective_contour",
            "core_minus_local_ring_occupancy",
            "core_minus_annual_stratified_occupancy",
            "nearest_pet_center_distance_km",
            "nearest_pet_center_distance_over_effective_radius",
        ),
        "strict surface association",
    )
    _require_columns(quality, ("pet_analysis_eligible",), "quality surface association")
    _require_columns(
        rotation,
        (
            "event_id",
            "surface_core_rotation_polarity_match",
            "surface_core_weighted_rossby_number",
            "peak_core_analysis_eligible",
            "peak_core_contained_by_actual_pet_effective_contour",
        ),
        "rotation surface association",
    )

    eligible = _true_mask(rotation["peak_core_analysis_eligible"])
    contained = _true_mask(rotation["peak_core_contained_by_actual_pet_effective_contour"])
    same_sign = _true_mask(rotation["surface_core_rotation_polarity_match"])
    if len(rotation) != 29 or int(eligible.sum()) != 16:
        raise ValueError(
            f"Unexpected rotation denominator: n={len(rotation)}, "
            f"eligible={int(eligible.sum())}"
        )
    if int(same_sign[eligible].sum()) != 15 or int(contained[eligible].sum()) != 0:
        raise ValueError("Frozen 16-event surface/deep correspondence counts changed")

    strict_eligible = _true_mask(association["peak_core_analysis_eligible"])
    quality_eligible = _true_mask(quality["pet_analysis_eligible"])
    if int(strict_eligible.sum()) != 31 or int(quality_eligible.sum()) != 58:
        raise ValueError("Frozen PET eligibility counts changed")

    case = rotation.loc[rotation["event_id"].eq("OFES_DO50_E000002")].iloc[0]
    case_assoc = association.loc[association["event_id"].eq("OFES_DO50_E000002")].iloc[0]

    local_audit = summary["local_ring_null"]
    annual_audit = summary["annual_month_latitude_null"]
    local_values = association["core_minus_local_ring_occupancy"].dropna()
    annual_values = association["core_minus_annual_stratified_occupancy"].dropna()
    if len(local_values) != int(local_audit["event_counts"]["n"]):
        raise ValueError("Local-ring parquet and summary counts disagree")
    if len(annual_values) != int(annual_audit["event_counts"]["n"]):
        raise ValueError("Annual-null parquet and summary counts disagree")

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(16.2, 5.5),
        gridspec_kw={"width_ratios": (1.35, 0.95, 1.05)},
        constrained_layout=True,
    )
    _plot_case_panel(
        axes[0], case, case_assoc, case_map, case_map_extent, case_map_crop,
    )

    values = [15 / 16, 0 / 16]
    labels = [
        "same-sign surface–deep\nRo correspondence",
        "closed-SSH\ncontainment",
    ]
    bars = axes[1].bar(np.arange(2), values, color=[BLUE, PURPLE], width=0.58)
    for bar, count in zip(bars, (15, 0)):
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            max(bar.get_height() + 0.04, 0.05),
            f"{count}/16",
            ha="center",
            fontsize=11,
        )
    axes[1].set_xticks(np.arange(2), labels, fontsize=8.5)
    axes[1].set_ylim(0, 1.15)
    axes[1].set_ylabel("Fraction of rotation PET-eligible events")
    axes[1].set_title("One denominator: 16 rotation events")
    axes[1].text(
        0.02,
        0.02,
        "All rotation events: 26/29\n(strict PET eligibility: 31/56)",
        transform=axes[1].transAxes,
        fontsize=8.2,
        va="bottom",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )

    null_rows = [
        ("local ring", local_audit, BLUE),
        ("annual month × latitude", annual_audit, TEAL),
    ]
    y = np.arange(len(null_rows))[::-1]
    for yi, (label, audit, color) in zip(y, null_rows):
        mean = float(
            audit["core_minus_ring_mean"]
            if "core_minus_ring_mean" in audit
            else audit["core_minus_annual_mean"]
        )
        low, high = (float(value) for value in audit["bootstrap_95ci"])
        axes[2].plot([low, high], [yi, yi], color=color, linewidth=3)
        axes[2].plot([low, low], [yi - 0.12, yi + 0.12], color=color, linewidth=1.2)
        axes[2].plot([high, high], [yi - 0.12, yi + 0.12], color=color, linewidth=1.2)
        axes[2].plot(mean, yi, "o", color=INK, markersize=7)
        axes[2].text(high + 0.0006, yi, f"{mean:+.4f}", va="center", fontsize=8.5)
    axes[2].axvline(0, color=INK, linestyle="--", linewidth=1)
    axes[2].set_yticks(y, [row[0] for row in null_rows])
    axes[2].set_xlabel("Event-equal core − null occupancy")
    axes[2].set_xlim(-0.033, 0.004)
    axes[2].set_title(
        "Formal null contrasts (n=31)\n"
        f"local CI [{local_audit['bootstrap_95ci'][0]:+.4f}, "
        f"{local_audit['bootstrap_95ci'][1]:+.4f}]"
    )
    _style_axes(axes)
    fig.suptitle(
        "Surface Ro correspondence persists while closed-SSH containment remains absent",
        fontsize=14,
    )
    return _save(fig, output_dir, "Figure6_surface_ssh_decoupling.png", dpi)


def render_paper_figures(
    surface_root: Path,
    output_dir: Path,
    dpi: int = 260,
    summary_path: Path | None = None,
    case_map: Path | None = None,
    case_map_extent: Sequence[float] | None = None,
    case_map_crop: Sequence[int] | None = None,
) -> list[Path]:
    """Render Figure 6 from completed PET association products."""

    surface_root = surface_root.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    summary, association, quality, rotation = _load_outputs(surface_root, summary_path)
    return [
        _plot_figure6(
            summary, association, quality, rotation, output_dir, dpi,
            case_map, case_map_extent, case_map_crop,
        )
    ]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--summary", type=Path, default=None)
    parser.add_argument("--case-map", type=Path, default=None)
    parser.add_argument("--case-map-extent", type=float, nargs=4, default=None)
    parser.add_argument("--case-map-crop", type=int, nargs=4, default=None)
    parser.add_argument("--dpi", type=int, default=260)
    args = parser.parse_args(argv)
    for path in render_paper_figures(
        args.surface_root, args.output_dir, args.dpi, args.summary,
        args.case_map, args.case_map_extent, args.case_map_crop,
    ):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
