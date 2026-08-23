"""Render the paper-facing OFES PET figures from audited association outputs.

The plotting entry point reads denominators from the parquet products rather
than hard-coding the strict, rotation, or sensitivity population sizes.  PET
domain failures remain ``unassessable`` and are never counted as no-PET
events.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


def _style_axes(axes: Any) -> None:
    for axis in np.ravel(axes):
        axis.grid(alpha=0.2)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)


def _save(fig: Any, output_dir: Path, name: str, dpi: int) -> Path:
    path = output_dir / name
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def _load_outputs(surface_root: Path) -> tuple[dict[str, Any], ...]:
    summary = json.loads((surface_root / "surface_eddy_summary.json").read_text())
    association = pd.read_parquet(surface_root / "surface_eddy_event_association.parquet")
    sensitivity = pd.read_parquet(surface_root / "surface_eddy_quality_eligible_161.parquet")
    rotation = pd.read_parquet(surface_root / "surface_eddy_surface_ro_crosstab.parquet")
    return summary, association, sensitivity, rotation


def _plot_f14(
    summary: Mapping[str, Any],
    association: pd.DataFrame,
    sensitivity: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> Path:
    del summary
    _require_columns(
        association,
        (
            "peak_core_analysis_eligible",
            "peak_core_contained_by_actual_pet_effective_contour",
            "peak_core_contained_by_actual_pet_speed_contour",
            "nearest_pet_center_distance_over_effective_radius",
        ),
        "strict association",
    )
    _require_columns(
        sensitivity,
        (
            "pet_analysis_eligible",
            "pet_effective_contained",
            "pet_speed_contained",
        ),
        "quality sensitivity association",
    )
    strict_eligible = _true_mask(association["peak_core_analysis_eligible"])
    quality_eligible = _true_mask(sensitivity["pet_analysis_eligible"])
    eligible_counts = [
        int(strict_eligible.sum()),
        int(strict_eligible.sum()),
        int(quality_eligible.sum()),
        int(quality_eligible.sum()),
    ]
    contained_counts = [
        int(_true_mask(association["peak_core_contained_by_actual_pet_effective_contour"]).sum()),
        int(_true_mask(association["peak_core_contained_by_actual_pet_speed_contour"]).sum()),
        int(_true_mask(sensitivity["pet_effective_contained"]).sum()),
        int(_true_mask(sensitivity["pet_speed_contained"]).sum()),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(15, 4.8), constrained_layout=True)
    groups = [
        f"Strict {len(association)}\neffective",
        f"Strict {len(association)}\nspeed",
        f"Quality {len(sensitivity)}\neffective",
        f"Quality {len(sensitivity)}\nspeed",
    ]
    x = np.arange(len(groups))
    axes[0].bar(x - 0.18, eligible_counts, width=0.36, color=BLUE, label="analysis-eligible")
    axes[0].bar(x + 0.18, contained_counts, width=0.36, color=RED, label="contained by PET")
    offset = max(eligible_counts) * 0.025
    for index, (eligible, contained) in enumerate(zip(eligible_counts, contained_counts)):
        axes[0].text(index - 0.18, eligible + offset, str(eligible), ha="center", fontsize=9)
        axes[0].text(index + 0.18, contained + offset, str(contained), ha="center", fontsize=9, color=RED)
    axes[0].set_xticks(x, groups, fontsize=9)
    axes[0].set_ylabel("Events")
    axes[0].set_ylim(0, max(eligible_counts) * 1.18)
    axes[0].set_title(
        "Peak-core containment among PET-domain-eligible events\n"
        f"(strict {eligible_counts[0]}/{len(association)}; quality "
        f"{eligible_counts[2]}/{len(sensitivity)} eligible)"
    )
    axes[0].legend(frameon=False, fontsize=8.5, loc="upper right")

    distance = association.loc[
        strict_eligible, "nearest_pet_center_distance_over_effective_radius"
    ].dropna()
    if distance.empty:
        raise ValueError("No eligible strict-event nearest-PET distances are available")
    median_distance = float(distance.median())
    axes[1].hist(distance, bins=16, color=PURPLE, alpha=0.8)
    axes[1].axvline(
        median_distance,
        color=ORANGE,
        linewidth=1.8,
        label=f"eligible median {median_distance:.1f}× effective radius",
    )
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Nearest PET center distance / effective radius")
    axes[1].set_ylabel("Eligible strict events")
    axes[1].set_title("Eligible event cores remain far from catalogued surface eddies")
    axes[1].legend(frameon=False, fontsize=9, loc="upper right")
    _style_axes(axes)
    fig.suptitle(
        "OFES DO50 cores are not captured by catalogued closed-SSH surface eddies "
        "within the PET-analysis domain"
    )
    return _save(fig, output_dir, "F14_surface_containment.png", dpi)


def _plot_f15(
    summary: Mapping[str, Any],
    association: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> Path:
    _require_columns(
        association,
        ("core_minus_local_ring_occupancy", "core_minus_annual_stratified_occupancy"),
        "strict association",
    )
    rows = [
        (
            "Annual month × 1° latitude",
            association["core_minus_annual_stratified_occupancy"].dropna(),
            summary["annual_month_latitude_null"],
            TEAL,
        ),
        (
            "Same-day 120–240 km ring",
            association["core_minus_local_ring_occupancy"].dropna(),
            summary["local_ring_null"],
            BLUE,
        ),
    ]
    for label, values, audit, _ in rows:
        if len(values) != int(audit["event_counts"]["n"]):
            raise ValueError(f"{label} count differs between parquet and summary")

    fig, axes = plt.subplots(1, 2, figsize=(15, 4.8), constrained_layout=True)
    for index, (label, _, audit, _) in enumerate(rows[::-1]):
        mean = float(
            audit.get("core_minus_ring_mean", audit.get("core_minus_annual_mean"))
        )
        low, high = [float(value) for value in audit["bootstrap_95ci"]]
        axes[0].plot([low, high], [index, index], color=RED, linewidth=2.5, zorder=2)
        axes[0].plot([mean], [index], "o", color=RED, markersize=8, zorder=3)
        axes[0].plot([low, low], [index - 0.15, index + 0.15], color=RED, linewidth=1.4)
        axes[0].plot([high, high], [index - 0.15, index + 0.15], color=RED, linewidth=1.4)
        axes[0].text(low - 0.0008, index, f"{mean:+.4f}", ha="right", va="center", fontsize=9, color=RED)
    axes[0].axvline(0, color=INK, linewidth=1, linestyle="--")
    axes[0].set_yticks([0, 1], [row[0] for row in rows][::-1])
    axes[0].set_xlabel("Event-equal mean difference (core − null)")
    axes[0].set_xlim(-0.032, 0.004)
    annual_p = float(summary["annual_month_latitude_null"]["wilcoxon_two_sided"][1])
    local_p = float(summary["local_ring_null"]["wilcoxon_two_sided"][1])
    axes[0].set_title(
        "Two-sided audit: both 95% CIs lie below zero\n"
        f"(annual p={annual_p:.2g}; local ring p={local_p:.2g}; paired n={len(rows[0][1])})"
    )

    for index, (label, values, audit, color) in enumerate(rows[::-1]):
        short_label = "ring null" if "ring" in label.lower() else "annual null"
        axes[1].hist(values, bins=12, alpha=0.55, color=color, label=short_label, density=True)
        counts = audit["event_counts"]
        axes[1].text(
            0.02,
            0.78 if index == 0 else 0.88,
            f"{counts['negative']}/{counts['n']} negative; "
            f"{counts['zero']} zero; {counts['positive']} positive",
            transform=axes[1].transAxes,
            fontsize=8.5,
            color=color,
        )
    axes[1].axvline(0, color=INK, linewidth=1, linestyle="--")
    axes[1].set_xlabel("Per-event core − null")
    axes[1].set_ylabel("Density")
    axes[1].set_title("Paired event distribution within the PET-analysis domain")
    axes[1].legend(frameon=False, fontsize=9, loc="upper right")
    _style_axes(axes)
    fig.suptitle(
        "PET-domain-eligible event cores have lower closed-SSH-eddy occupancy than both backgrounds"
    )
    return _save(fig, output_dir, "F15_surface_null_two_sided.png", dpi)


def _plot_f16(rotation: pd.DataFrame, output_dir: Path, dpi: int) -> Path:
    _require_columns(
        rotation,
        (
            "peak_core_analysis_eligible",
            "peak_core_contained_by_actual_pet_effective_contour",
            "surface_core_rotation_polarity_match",
        ),
        "rotation association",
    )
    eligible = _true_mask(rotation["peak_core_analysis_eligible"])
    contained = _true_mask(rotation["peak_core_contained_by_actual_pet_effective_contour"])
    eligible_no_pet = int((eligible & ~contained).sum())
    pet_contained = int(contained.sum())
    unassessable = int((~eligible).sum())

    fig, axes = plt.subplots(1, 2, figsize=(15, 4.8), constrained_layout=True)
    labels = ["eligible\nno PET", "PET-contained", "unassessable"]
    counts = [eligible_no_pet, pet_contained, unassessable]
    axes[0].bar(np.arange(3), counts, color=[PURPLE, RED, GRAY], width=0.55)
    for index, value in enumerate(counts):
        axes[0].text(index, value + 0.45, str(value), ha="center", fontsize=11)
    axes[0].set_xticks(np.arange(3), labels)
    axes[0].set_ylabel("Rotation-dominated events")
    axes[0].set_ylim(0, max(counts) + 4)
    axes[0].set_title(
        f"No eligible rotation event contains a PET eddy "
        f"(0/{int(eligible.sum())}); {unassessable} unassessable"
    )

    match = rotation["surface_core_rotation_polarity_match"]
    match_counts = [int(_true_mask(match).sum()), int(match.eq(False).fillna(False).sum())]  # noqa: E712
    axes[1].bar([0, 1], match_counts, color=[BLUE, GRAY], width=0.45)
    for index, value in enumerate(match_counts):
        axes[1].text(index, value + 0.5, str(value), ha="center", fontsize=11)
    axes[1].set_xticks([0, 1], ["surface Ro matches\ndeep polarity", "mismatch"])
    axes[1].set_ylabel("Events")
    axes[1].set_ylim(0, len(rotation) + 1)
    axes[1].set_title(
        f"Surface and deep core-weighted Ro share polarity "
        f"({match_counts[0]}/{len(rotation)})"
    )
    _style_axes(axes)
    fig.suptitle(
        "Rotation expression can share polarity without an assessable catalogued closed-SSH surface eddy"
    )
    return _save(fig, output_dir, "F16_surface_rotation_pet.png", dpi)


def render_paper_figures(surface_root: Path, output_dir: Path, dpi: int = 180) -> list[Path]:
    """Render F14--F16 and return their output paths."""

    surface_root = surface_root.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    summary, association, sensitivity, rotation = _load_outputs(surface_root)
    return [
        _plot_f14(summary, association, sensitivity, output_dir, dpi),
        _plot_f15(summary, association, output_dir, dpi),
        _plot_f16(rotation, output_dir, dpi),
    ]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args(argv)
    for path in render_paper_figures(args.surface_root, args.output_dir, args.dpi):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
