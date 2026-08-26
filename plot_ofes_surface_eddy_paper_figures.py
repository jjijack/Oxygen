"""Render the final OFES surface-expression manuscript Figure 6.

The figure consumes completed PET association products and reconstructs the
representative map from the source tables used by
``run_ofes_surface_eddy_postprocess._write_plot_outputs``. In particular, it
does not place an old rendered case-map raster behind the figure: contours,
centers, and the DO peak pixels are plotted directly from their parquet
products, so no baked-in axes or annotations can survive into Figure 6a.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd

BLUE = "#4c78a8"
TEAL = "#72b7b2"
PURPLE = "#7b61a8"
RED = "#d62728"
INK = "#111827"
GRAY = "#94a3b8"
OCHRE = "#a16207"
FIGURE_WIDTH_IN = 7.1
PANEL_LABEL_SIZE = 9.5
MAP_REGION = (140.0, 170.0, 25.0, 45.0)

plt.rcParams.update(
    {
        "font.size": 8.0,
        "axes.titlesize": 8.6,
        "axes.labelsize": 8.0,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "lines.linewidth": 1.0,
    }
)


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


def _add_panel_labels(axes: Any) -> None:
    """Place panel labels outside the upper-left axes corner."""

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
    wspace: float | None = None,
) -> Path:
    path = output_dir / name
    # Reserve margins on the fixed canvas so the PNG remains exactly 7.1 in
    # wide at the requested dpi while outside labels stay inside the canvas.
    kwargs = {"left": 0.12, "right": 0.985, "top": 0.88, "bottom": bottom}
    if wspace is not None:
        kwargs["wspace"] = wspace
    fig.subplots_adjust(**kwargs)
    fig.savefig(path, dpi=dpi, facecolor="white", edgecolor="none")
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


def _load_case_map_data(
    surface_root: Path,
    do_catalog_root: Path,
    population_path: Path,
    event_id: str,
) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame]:
    """Load the source rows needed for a pure, data-generated case map."""

    objects = pd.read_parquet(surface_root / "surface_eddy_daily_objects.parquet")
    _require_columns(
        objects,
        (
            "date",
            "object_id",
            "center_lon",
            "center_lat",
            "effective_contour_lon",
            "effective_contour_lat",
            "speed_contour_lon",
            "speed_contour_lat",
            "is_virtual",
            "filter_valid",
            "boundary_censored",
        ),
        "surface daily objects",
    )
    population = pd.read_parquet(population_path)
    _require_columns(
        population,
        ("event_id", "peak_date", "peak_lon", "peak_lat", "daily_object_key"),
        "event population",
    )
    population_row = population.loc[population["event_id"].astype(str).eq(event_id)]
    if len(population_row) != 1:
        raise ValueError(f"Expected one population row for {event_id}, found {len(population_row)}")
    population_row = population_row.iloc[0]
    peak_date = pd.Timestamp(population_row["peak_date"]).normalize()
    objects["date"] = pd.to_datetime(objects["date"]).dt.normalize()
    rows = objects.loc[
        objects["date"].eq(peak_date)
        & ~_true_mask(objects["is_virtual"])
        & _true_mask(objects["filter_valid"])
        & ~_true_mask(objects["boundary_censored"])
    ].copy()
    if rows.empty:
        raise ValueError(f"No valid non-virtual surface objects for {event_id} on {peak_date:%Y-%m-%d}")

    key = str(population_row["daily_object_key"])
    try:
        object_id = int(key.rsplit("_", 1)[-1])
    except ValueError as error:
        raise ValueError(f"Cannot parse daily_object_key {key}") from error
    pixel_path = do_catalog_root / "days" / f"peak_pixels_{peak_date:%Y%m%d}.parquet"
    if not pixel_path.is_file():
        raise FileNotFoundError(pixel_path)
    pixels = pd.read_parquet(pixel_path)
    _require_columns(pixels, ("lon", "lat", "object_id_do50"), "peak-pixel table")
    pixels = pixels.loc[pixels["object_id_do50"].astype(int).eq(object_id)].copy()
    if pixels.empty:
        raise ValueError(f"No DO50 peak pixels for {event_id} / {key}")
    return population_row, rows, pixels


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
            breaks = np.flatnonzero(np.abs(np.diff(part[:, 0])) > 180.0) + 1
            segments.extend(chunk for chunk in np.split(part, breaks) if len(chunk) >= 2)
    return segments


def _add_coastlines(axis: Any) -> None:
    segments = _coastline_segments()
    if segments:
        axis.add_collection(
            LineCollection(
                segments,
                colors="#b6c0cc",
                linewidths=0.45,
                alpha=0.9,
                zorder=1,
            )
        )


def _plot_contour(axis: Any, row: Mapping[str, Any], lon_key: str, lat_key: str, **kwargs: Any) -> None:
    lon = np.asarray(row[lon_key], dtype=float)
    lat = np.asarray(row[lat_key], dtype=float)
    if lon.size < 2 or lat.size != lon.size:
        return
    axis.plot(np.r_[lon, lon[0]], np.r_[lat, lat[0]], **kwargs)


def _plot_case_panel(
    axis: Any,
    case: pd.Series,
    case_rows: pd.DataFrame,
    pixels: pd.DataFrame,
) -> None:
    """Draw Figure 6a directly from daily contours and DO peak pixels."""

    peak_lon = float(case["peak_lon"])
    peak_lat = float(case["peak_lat"])
    peak_date = pd.Timestamp(case["peak_date"]).normalize()
    xlim = (
        max(MAP_REGION[0], peak_lon - 4.0),
        min(MAP_REGION[1], peak_lon + 4.0),
    )
    ylim = (
        max(MAP_REGION[2], peak_lat - 4.0),
        min(MAP_REGION[3], peak_lat + 4.0),
    )
    axis.set_facecolor("#f7f9fc")
    _add_coastlines(axis)
    contour_handles: dict[str, Any] = {}
    for row in case_rows.to_dict("records"):
        before = len(axis.lines)
        _plot_contour(
            axis,
            row,
            "effective_contour_lon",
            "effective_contour_lat",
            color=GRAY,
            linewidth=1.2,
            zorder=2,
        )
        if len(axis.lines) > before:
            contour_handles.setdefault("PET effective contour", axis.lines[-1])
        before = len(axis.lines)
        _plot_contour(
            axis,
            row,
            "speed_contour_lon",
            "speed_contour_lat",
            color=PURPLE,
            linewidth=1.0,
            linestyle="--",
            zorder=2,
        )
        if len(axis.lines) > before:
            contour_handles.setdefault("PET speed contour", axis.lines[-1])
        axis.plot(row["center_lon"], row["center_lat"], "o", ms=3.2, color=INK, zorder=3)

    do_handle = axis.scatter(
        pixels["lon"],
        pixels["lat"],
        s=8,
        color=OCHRE,
        alpha=0.48,
        linewidths=0,
        label="DO peak pixels",
        zorder=2,
    )
    core_handle = axis.scatter(
        peak_lon,
        peak_lat,
        marker="*",
        s=105,
        color=INK,
        label="DO peak core",
        zorder=5,
    )
    distances = (case_rows["center_lon"].to_numpy(float) - peak_lon) ** 2 + (
        case_rows["center_lat"].to_numpy(float) - peak_lat
    ) ** 2
    nearest = case_rows.iloc[int(np.argmin(distances))]
    nearest_handle = axis.scatter(
        [float(nearest["center_lon"])],
        [float(nearest["center_lat"])],
        s=82,
        facecolors="none",
        edgecolors=RED,
        linewidths=1.5,
        label="nearest PET center",
        zorder=4,
    )
    ro_handle = axis.scatter(
        peak_lon,
        peak_lat,
        s=150,
        facecolors="none",
        edgecolors=BLUE,
        linewidths=2.0,
        label="core-weighted surface Ro",
        zorder=4,
    )
    axis.set_xlim(*xlim)
    axis.set_ylim(*ylim)
    axis.set_xlabel("Longitude (°E)")
    axis.set_ylabel("Latitude (°N)")
    axis.set_title("E000002", loc="right")
    handles = [
        do_handle,
        core_handle,
        contour_handles["PET effective contour"],
        contour_handles["PET speed contour"],
        nearest_handle,
        ro_handle,
    ]
    labels = [
        "DO peak pixels",
        "DO peak core",
        "PET effective contour",
        "PET speed contour",
        "nearest PET center",
        "core-weighted surface Ro",
    ]
    axis.legend(
        handles,
        labels,
        frameon=False,
        fontsize=7.5,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.22),
        ncol=2,
        handlelength=1.5,
        columnspacing=0.8,
    )


def _plot_figure6(
    summary: Mapping[str, Any],
    association: pd.DataFrame,
    quality: pd.DataFrame,
    rotation: pd.DataFrame,
    case: pd.Series,
    case_rows: pd.DataFrame,
    pixels: pd.DataFrame,
    output_dir: Path,
    dpi: int,
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
            f"Unexpected rotation denominator: n={len(rotation)}, eligible={int(eligible.sum())}"
        )
    if int(same_sign[eligible].sum()) != 15 or int(contained[eligible].sum()) != 0:
        raise ValueError("Frozen 16-event surface/deep correspondence counts changed")

    strict_eligible = _true_mask(association["peak_core_analysis_eligible"])
    quality_eligible = _true_mask(quality["pet_analysis_eligible"])
    if int(strict_eligible.sum()) != 31 or int(quality_eligible.sum()) != 58:
        raise ValueError("Frozen PET eligibility counts changed")

    case_assoc = association.loc[association["event_id"].eq("OFES_DO50_E000002")].iloc[0]
    local_audit = summary["local_ring_null"]
    annual_audit = summary["annual_month_latitude_null"]
    local_values = association["core_minus_local_ring_occupancy"].dropna()
    annual_values = association["core_minus_annual_stratified_occupancy"].dropna()
    if len(local_values) != int(local_audit["event_counts"]["n"]):
        raise ValueError("Local-ring parquet and summary counts disagree")
    if len(annual_values) != int(annual_audit["event_counts"]["n"]):
        raise ValueError("Annual-null parquet and summary counts disagree")
    if not np.isfinite(float(case_assoc["nearest_pet_center_distance_km"])):
        raise ValueError("Representative case lacks a nearest PET center")

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(FIGURE_WIDTH_IN, 3.15),
        gridspec_kw={"width_ratios": (1.38, 0.95, 1.12)},
    )
    _plot_case_panel(axes[0], case, case_rows, pixels)

    values = [15 / 16, 0 / 16]
    labels = ["same-sign\nsurface Ro", "closed SSH\ncontained"]
    bars = axes[1].bar(np.arange(2), values, color=[BLUE, PURPLE], width=0.58)
    for bar, count in zip(bars, (15, 0)):
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            max(bar.get_height() + 0.04, 0.05),
            f"{count}/16",
            ha="center",
            fontsize=8.0,
        )
    axes[1].set_xticks(np.arange(2), labels, fontsize=7.5)
    axes[1].set_ylim(0, 1.15)
    axes[1].set_ylabel("Fraction of eligible events")
    axes[1].set_title("16-event subset", loc="right")

    null_rows = [
        ("local ring", local_audit, BLUE),
        ("annual stratified", annual_audit, TEAL),
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
        axes[2].text(high + 0.0006, yi, f"{mean:+.4f}", va="center", fontsize=7.5)
    axes[2].axvline(0, color=INK, linestyle="--", linewidth=1)
    axes[2].set_yticks(y, [row[0] for row in null_rows])
    axes[2].set_xlabel("Core − null occupancy")
    axes[2].set_xlim(-0.033, 0.004)
    axes[2].set_title("Occupancy nulls", loc="right")
    _style_axes(axes)
    _add_panel_labels(axes)
    return _save(fig, output_dir, "Figure6_surface_ssh_decoupling.png", dpi, bottom=0.30, wspace=0.45)


def render_paper_figures(
    surface_root: Path,
    output_dir: Path,
    dpi: int = 600,
    summary_path: Path | None = None,
    do_catalog_root: Path | None = None,
    population_path: Path | None = None,
    case_event_id: str = "OFES_DO50_E000002",
) -> list[Path]:
    """Render Figure 6 from completed PET association and map-source products."""

    surface_root = surface_root.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if do_catalog_root is None or population_path is None:
        raise ValueError("Figure 6 requires --do-catalog-root and --population for a pure case map")
    summary, association, quality, rotation = _load_outputs(surface_root, summary_path)
    case, case_rows, pixels = _load_case_map_data(
        surface_root,
        do_catalog_root.expanduser().resolve(),
        population_path.expanduser().resolve(),
        str(case_event_id),
    )
    return [
        _plot_figure6(
            summary,
            association,
            quality,
            rotation,
            case,
            case_rows,
            pixels,
            output_dir,
            dpi,
        )
    ]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--summary", type=Path, default=None)
    parser.add_argument("--do-catalog-root", required=True, type=Path)
    parser.add_argument("--population", required=True, type=Path)
    parser.add_argument("--case-event-id", default="OFES_DO50_E000002")
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args(argv)
    for path in render_paper_figures(
        args.surface_root,
        args.output_dir,
        args.dpi,
        args.summary,
        args.do_catalog_root,
        args.population,
        args.case_event_id,
    ):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
