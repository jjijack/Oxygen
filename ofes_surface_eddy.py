"""Standalone PET detector and tracker for the OFES NP30 surface product.

This module deliberately does not import :mod:`track`.  It is intended to be
run from the dedicated ``ofes-pet`` environment, where the PET API and its
NetCDF dependencies are available without the broader Oxygen plotting stack.
"""

from __future__ import annotations

import datetime as dt
import logging
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml
from netCDF4 import Dataset


def _configure_writable_runtime_caches() -> None:
    """Keep optional PET/Numba caches out of read-only home trees when used."""
    cache_root = Path(tempfile.gettempdir()) / "ofes-pet-runtime-cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("NUMBA_CACHE_DIR", str(cache_root / "numba"))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root / "xdg"))


_configure_writable_runtime_caches()


# PET is intentionally imported only inside producer helpers below.

LOGGER = logging.getLogger("ofes_surface_eddy")
PET_EPOCH = dt.datetime(1950, 1, 1)
POLARITY_NAME = {1: "anticyclonic", -1: "cyclonic"}
POLARITY_CODE = {"anticyclonic": 1, "cyclonic": -1}

CONTOUR_COLUMNS = (
    "effective_contour_lon",
    "effective_contour_lat",
    "speed_contour_lon",
    "speed_contour_lat",
)
LIST_COLUMNS = CONTOUR_COLUMNS
NUMERIC_INTERPOLATION_COLUMNS = (
    "center_lon",
    "center_lat",
    "ssh_extremum_lon",
    "ssh_extremum_lat",
    "effective_radius_km",
    "speed_radius_km",
    "amplitude_m",
    "shape_error_effective_pct",
    "shape_error_speed_pct",
    "effective_area_m2",
    "speed_area_m2",
)
DAILY_COLUMNS = (
    "date",
    "object_id",
    "polarity",
    "polarity_code",
    "center_lon",
    "center_lat",
    "ssh_extremum_lon",
    "ssh_extremum_lat",
    "effective_radius_km",
    "speed_radius_km",
    "amplitude_m",
    "shape_error_effective_pct",
    "shape_error_speed_pct",
    "effective_area_m2",
    "speed_area_m2",
    *LIST_COLUMNS,
    "is_virtual",
    "boundary_censored",
    "filter_valid",
    "coastal",
    "source_eta_path",
    "track_id",
    "track_observed_days",
    "track_virtual_days",
    "track_duration_days",
    "track_duration_class",
)


@dataclass(frozen=True)
class EtaSnapshot:
    """One validated OFES ETA snapshot in metres and ``(lat, lon)`` order."""

    date: dt.date
    path: Path
    lon: np.ndarray
    lat: np.ndarray
    eta_m: np.ma.MaskedArray
    raw_mask: np.ndarray
    raw_units: str


@dataclass(frozen=True)
class PreparedPetGrid:
    """PET grid and validity masks after the frozen high-pass operation."""

    grid: RegularGridDataset
    filter_valid_mask: np.ndarray
    ocean_valid_mask: np.ndarray
    raw_mask_xy: np.ndarray
    lon: np.ndarray
    lat: np.ndarray


def load_yaml_config(path: str | Path) -> dict[str, Any]:
    """Load and minimally validate the standalone detector YAML."""

    config_path = Path(path).resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration is not a mapping: {config_path}")
    for section in ("input", "filter", "detection", "tracking", "output"):
        if section not in config:
            raise KeyError(f"Missing configuration section: {section}")
    return config


def parse_date(value: str | dt.date | dt.datetime) -> dt.date:
    """Return a date from an ISO string or datetime-like value."""

    if isinstance(value, dt.datetime):
        return value.date()
    if isinstance(value, dt.date):
        return value
    text = str(value).strip()
    # cftime objects returned by netCDF4 commonly stringify as
    # ``YYYY-MM-DD HH:MM:SS`` rather than a date-only ISO value.
    if "T" in text:
        text = text.split("T", 1)[0]
    elif " " in text:
        text = text.split(" ", 1)[0]
    return dt.date.fromisoformat(text)


def date_range(start: str | dt.date, end: str | dt.date) -> list[dt.date]:
    """Return an inclusive daily date list."""

    start_date, end_date = parse_date(start), parse_date(end)
    if end_date < start_date:
        raise ValueError("end date precedes start date")
    return [
        start_date + dt.timedelta(days=offset)
        for offset in range((end_date - start_date).days + 1)
    ]


def eta_path_for_date(eta_root: str | Path, date_value: str | dt.date) -> Path:
    """Resolve the task-book ETA filename for one date."""

    day = parse_date(date_value)
    return Path(eta_root) / f"eta.{day:%m.%d.%Y}.nc"


def _validate_coordinates(
    lon: np.ndarray, lat: np.ndarray, eta_lat_lon: np.ndarray | np.ma.MaskedArray
) -> None:
    """Reject ambiguous coordinate order and malformed regular-grid input."""

    if lon.ndim != 1 or lat.ndim != 1:
        raise ValueError("OFES lon/lat must both be one-dimensional")
    if lon.size < 3 or lat.size < 3:
        raise ValueError("OFES lon/lat arrays are too short")
    if not np.isfinite(lon).all() or not np.isfinite(lat).all():
        raise ValueError("OFES lon/lat contain non-finite values")
    if not np.all(np.diff(lon) > 0) or not np.all(np.diff(lat) > 0):
        raise ValueError("OFES lon/lat must be strictly increasing")
    expected = (lat.size, lon.size)
    if eta_lat_lon.shape == (lon.size, lat.size) and expected != eta_lat_lon.shape:
        raise ValueError(
            "ETA appears to be in (lon, lat) order; the loader requires (lat, lon)"
        )
    if eta_lat_lon.shape != expected:
        raise ValueError(
            f"ETA shape {eta_lat_lon.shape} does not match (lat, lon)={expected}"
        )


def convert_eta_to_m(
    eta: np.ndarray | np.ma.MaskedArray, raw_units: str
) -> np.ma.MaskedArray:
    """Convert raw ETA to metres while preserving the source mask."""

    value = np.ma.asarray(eta)
    mask = np.ma.getmaskarray(value).copy()
    units = str(raw_units).strip().lower()
    if units in {"cm", "centimeter", "centimeters"}:
        converted = value.astype("f8") * 0.01
    elif units in {"m", "meter", "meters"}:
        converted = value.astype("f8")
    else:
        raise ValueError(f"Unsupported ETA units: {raw_units!r}")
    return np.ma.array(np.asarray(converted), mask=mask, copy=False)


def load_eta_snapshot(
    path: str | Path,
    expected_date: str | dt.date | None = None,
    raw_units: str = "cm",
) -> EtaSnapshot:
    """Read and validate one OFES ETA file without filling land values."""

    source_path = Path(path).resolve()
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    expected = parse_date(expected_date) if expected_date is not None else None
    with Dataset(source_path, "r") as dataset:
        required = {"lon", "lat", "eta"}
        missing = required.difference(dataset.variables)
        if missing:
            raise KeyError(f"Missing ETA variables {sorted(missing)} in {source_path}")
        lon = np.asarray(dataset.variables["lon"][:], dtype="f8")
        lat = np.asarray(dataset.variables["lat"][:], dtype="f8")
        eta_variable = dataset.variables["eta"]
        dimensions = tuple(eta_variable.dimensions)
        if dimensions != ("time", "lat", "lon"):
            raise ValueError(
                "ETA dimensions must be exactly ('time', 'lat', 'lon'), "
                f"got {dimensions!r}"
            )
        if len(eta_variable.shape) != 3 or eta_variable.shape[0] != 1:
            raise ValueError(f"ETA must have one time record, got {eta_variable.shape}")
        eta_raw = np.ma.asarray(eta_variable[0])
        _validate_coordinates(lon, lat, eta_raw)
        eta_m = convert_eta_to_m(eta_raw, raw_units)
        raw_mask = np.ma.getmaskarray(eta_m).copy()
        time_value = None
        if "time" in dataset.variables:
            time_variable = dataset.variables["time"]
            if time_variable.size:
                time_value = float(np.asarray(time_variable[:]).ravel()[0])
                units = getattr(time_variable, "units", None)
                if units and "since" in units:
                    try:
                        import netCDF4

                        time_date = netCDF4.num2date(
                            time_value,
                            units,
                            calendar=getattr(time_variable, "calendar", "standard"),
                        )
                        time_date = parse_date(time_date)
                        if expected is not None and time_date != expected:
                            raise ValueError(
                                f"ETA time {time_date} disagrees with requested {expected}"
                            )
                    except (TypeError, ValueError, OverflowError):
                        raise
        if expected is None:
            match = source_path.name.split(".")
            if len(match) >= 4:
                expected = dt.date(int(match[3]), int(match[1]), int(match[2]))
            else:
                raise ValueError(f"Cannot infer date from ETA filename: {source_path}")
    return EtaSnapshot(
        date=expected,
        path=source_path,
        lon=lon,
        lat=lat,
        eta_m=eta_m,
        raw_mask=raw_mask,
        raw_units=raw_units,
    )


def _kernel_valid_mask(
    grid: RegularGridDataset,
    wave_length_km: float,
    order: int,
    lat_max: float,
) -> np.ndarray:
    """Return the complete-kernel interior mask in PET ``(lon, lat)`` order."""

    valid = np.zeros((grid.x_c.size, grid.y_c.size), dtype=bool)
    x_indices = np.arange(grid.x_c.size)
    for j, latitude in enumerate(np.asarray(grid.y_c)):
        if abs(float(latitude)) > lat_max:
            continue
        kernel = grid.kernel_bessel(float(latitude), wave_length_km, order=order)
        half_x = (kernel.shape[0] - 1) // 2
        half_y = (kernel.shape[1] - 1) // 2
        valid[:, j] = (x_indices >= half_x) & (
            x_indices < grid.x_c.size - half_x
        )
        if j < half_y or j >= grid.y_c.size - half_y:
            valid[:, j] = False
    return valid


def prepare_pet_grid(snapshot: EtaSnapshot, config: Mapping[str, Any]) -> PreparedPetGrid:
    """Build a PET grid, apply the configured filter, and add PET geostrophic speed.

    PET is loaded lazily here so ordinary Oxygen imports and loaders work in the
    main environment without the dedicated detector installation.
    """
    from py_eddy_tracker.dataset.grid import RegularGridDataset

    filter_config = config["filter"]
    data_xy = np.ma.array(snapshot.eta_m, copy=True).T
    raw_mask_xy = np.ma.getmaskarray(data_xy).copy()
    grid = RegularGridDataset.with_array(
        ("lon", "lat"),
        {
            "lon": snapshot.lon,
            "lat": snapshot.lat,
            "eta": data_xy,
        },
        variables_description={"eta": {"units": "m", "long_name": "model SSH"}},
        centered=True,
    )
    grid.bessel_high_filter(
        "eta",
        float(filter_config["wavelength_km"]),
        order=int(filter_config["order"]),
        lat_max=float(filter_config["lat_max"]),
        extend=bool(filter_config["extend"]),
    )
    filter_valid = _kernel_valid_mask(
        grid,
        float(filter_config["wavelength_km"]),
        int(filter_config["order"]),
        float(filter_config["lat_max"]),
    )
    filtered = np.ma.array(grid.grid("eta"), copy=True)
    # A complete kernel footprint is necessary but not sufficient: PET may
    # still return a masked high-pass value when the source mask intersects
    # the convolution.  Keep that distinction in the audit mask instead of
    # treating a numerically present contour as scientifically valid.
    filtered_mask = np.ma.getmaskarray(filtered)
    filter_valid &= ~filtered_mask
    filtered.mask = filtered_mask | ~filter_valid
    grid.vars["eta"] = filtered
    # PET's documented order is high-pass ETA first, then stencil u/v.
    grid.add_uv(
        "eta",
        "u",
        "v",
        stencil_halfwidth=int(config["detection"]["stencil_halfwidth"]),
    )
    return PreparedPetGrid(
        grid=grid,
        filter_valid_mask=filter_valid,
        ocean_valid_mask=~raw_mask_xy,
        raw_mask_xy=raw_mask_xy,
        lon=snapshot.lon,
        lat=snapshot.lat,
    )


def _finite_contour(lon: Any, lat: Any) -> tuple[list[float], list[float]]:
    """Convert a PET fixed-size contour array to finite Python lists."""

    lon_array = np.asarray(lon, dtype="f8").ravel()
    lat_array = np.asarray(lat, dtype="f8").ravel()
    mask = np.isfinite(lon_array) & np.isfinite(lat_array)
    return lon_array[mask].tolist(), lat_array[mask].tolist()


def _contour_indices(
    lon: Sequence[float],
    lat: Sequence[float],
    grid_lon: np.ndarray,
    grid_lat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map contour vertices to nearest grid cells, retaining out-of-domain flags."""

    lon_array, lat_array = np.asarray(lon, dtype="f8"), np.asarray(lat, dtype="f8")
    x = np.searchsorted(grid_lon, lon_array).clip(0, grid_lon.size - 1)
    y = np.searchsorted(grid_lat, lat_array).clip(0, grid_lat.size - 1)
    return x.astype(int), y.astype(int)


def _contour_boundary_flags(
    contour_lon: Sequence[float],
    contour_lat: Sequence[float],
    prepared: PreparedPetGrid,
) -> tuple[bool, bool]:
    """Return ``(boundary_censored, filter_valid)`` for one contour."""

    if not contour_lon or not contour_lat:
        return True, False
    lon_array = np.asarray(contour_lon, dtype="f8")
    lat_array = np.asarray(contour_lat, dtype="f8")
    finite = np.isfinite(lon_array) & np.isfinite(lat_array)
    if not finite.all():
        return True, False
    x, y = _contour_indices(lon_array, lat_array, prepared.lon, prepared.lat)
    outside = (
        (lon_array < prepared.lon[0])
        | (lon_array > prepared.lon[-1])
        | (lat_array < prepared.lat[0])
        | (lat_array > prepared.lat[-1])
    )
    filter_ok = bool((~outside).all() and prepared.filter_valid_mask[x, y].all())
    outer = (x == 0) | (x == prepared.lon.size - 1) | (y == 0) | (
        y == prepared.lat.size - 1
    )
    return bool(outside.any() or outer.any() or not filter_ok), filter_ok


def _contour_is_coastal(
    contour_lon: Sequence[float],
    contour_lat: Sequence[float],
    prepared: PreparedPetGrid,
) -> bool:
    """Flag a contour whose vertices lie within one cell of the source mask."""

    if not contour_lon or not contour_lat:
        return False
    x, y = _contour_indices(contour_lon, contour_lat, prepared.lon, prepared.lat)
    for x_index, y_index in zip(x, y):
        x0, x1 = max(0, x_index - 1), min(prepared.lon.size, x_index + 2)
        y0, y1 = max(0, y_index - 1), min(prepared.lat.size, y_index + 2)
        if prepared.raw_mask_xy[x0:x1, y0:y1].any():
            return True
    return False


def _field_value(observations: Any, name: str, index: int, default: float = np.nan) -> Any:
    """Read an optional PET observation field without hiding real errors."""

    try:
        return getattr(observations, name)[index]
    except (AttributeError, KeyError, IndexError):
        return default


def _observation_records(
    observations: Any,
    prepared: PreparedPetGrid,
    date_value: dt.date,
    source_eta_path: Path,
) -> list[dict[str, Any]]:
    """Flatten PET observations to the stable project table schema."""

    records: list[dict[str, Any]] = []
    sign_type = int(observations.sign_type)
    polarity = POLARITY_NAME.get(sign_type, str(sign_type))
    for index in range(len(observations)):
        effective_lon, effective_lat = _finite_contour(
            observations.contour_lon_e[index], observations.contour_lat_e[index]
        )
        speed_lon, speed_lat = _finite_contour(
            observations.contour_lon_s[index], observations.contour_lat_s[index]
        )
        effective_boundary, effective_valid = _contour_boundary_flags(
            effective_lon, effective_lat, prepared
        )
        speed_boundary, speed_valid = _contour_boundary_flags(
            speed_lon, speed_lat, prepared
        )
        object_id = f"{date_value:%Y%m%d}_{polarity[0].upper()}_{index:04d}"
        records.append(
            {
                "date": date_value.isoformat(),
                "object_id": object_id,
                "polarity": polarity,
                "polarity_code": sign_type,
                "center_lon": float(_field_value(observations, "lon", index)),
                "center_lat": float(_field_value(observations, "lat", index)),
                "ssh_extremum_lon": float(
                    _field_value(observations, "lon_max", index)
                ),
                "ssh_extremum_lat": float(
                    _field_value(observations, "lat_max", index)
                ),
                "effective_radius_km": float(
                    _field_value(observations, "radius_e", index) / 1000.0
                ),
                "speed_radius_km": float(
                    _field_value(observations, "radius_s", index) / 1000.0
                ),
                "amplitude_m": float(_field_value(observations, "amplitude", index)),
                "shape_error_effective_pct": float(
                    _field_value(observations, "shape_error_e", index)
                ),
                "shape_error_speed_pct": float(
                    _field_value(observations, "shape_error_s", index)
                ),
                "effective_area_m2": float(
                    _field_value(observations, "effective_area", index)
                ),
                "speed_area_m2": float(
                    _field_value(observations, "speed_area", index)
                ),
                "effective_contour_lon": effective_lon,
                "effective_contour_lat": effective_lat,
                "speed_contour_lon": speed_lon,
                "speed_contour_lat": speed_lat,
                "is_virtual": False,
                "boundary_censored": bool(effective_boundary or speed_boundary),
                "filter_valid": bool(effective_valid and speed_valid),
                "coastal": bool(
                    _contour_is_coastal(effective_lon, effective_lat, prepared)
                ),
                "source_eta_path": str(source_eta_path),
                "track_id": None,
                "track_observed_days": np.nan,
                "track_virtual_days": np.nan,
                "track_duration_days": np.nan,
                "track_duration_class": None,
            }
        )
    return records


def empty_daily_frame() -> pd.DataFrame:
    """Return an empty daily frame with stable column order."""

    frame = pd.DataFrame({column: pd.Series(dtype="object") for column in DAILY_COLUMNS})
    for column in (
        "polarity_code",
        "center_lon",
        "center_lat",
        "ssh_extremum_lon",
        "ssh_extremum_lat",
        "effective_radius_km",
        "speed_radius_km",
        "amplitude_m",
        "shape_error_effective_pct",
        "shape_error_speed_pct",
        "effective_area_m2",
        "speed_area_m2",
        "track_observed_days",
        "track_virtual_days",
        "track_duration_days",
    ):
        frame[column] = pd.Series(dtype="float64")
    for column in ("is_virtual", "boundary_censored", "filter_valid", "coastal"):
        frame[column] = pd.Series(dtype="bool")
    return frame


def detect_ofes_eddies_day(
    date_value: str | dt.date,
    config: Mapping[str, Any],
    eta_root: str | Path,
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any]]:
    """Detect one date with PET and return objects, native observations, and audit."""

    day = parse_date(date_value)
    source_path = eta_path_for_date(eta_root, day)
    snapshot = load_eta_snapshot(
        source_path,
        expected_date=day,
        raw_units=str(config["input"]["raw_eta_units"]),
    )
    region = config["input"].get("region", {})
    if region:
        region_bounds = {
            "lon": (float(region["lon_min"]), float(region["lon_max"])),
            "lat": (float(region["lat_min"]), float(region["lat_max"])),
        }
        for name, values in (("lon", snapshot.lon), ("lat", snapshot.lat)):
            lower, upper = region_bounds[name]
            if float(values.min()) < lower or float(values.max()) > upper:
                raise ValueError(
                    f"OFES {name} coordinates exceed locked region "
                    f"[{lower}, {upper}]: "
                    f"[{float(values.min())}, {float(values.max())}]"
                )
    prepared = prepare_pet_grid(snapshot, config)
    detection = config["detection"]
    pet_date = dt.datetime.combine(day, dt.time())
    observations = prepared.grid.eddy_identification(
        "eta",
        "u",
        "v",
        pet_date,
        step=float(detection["contour_interval_m"]),
        shape_error=float(detection["shape_error_pct"]),
        pixel_limit=(
            int(detection["pixel_limit_min"]),
            int(detection["pixel_limit_max"]),
        ),
        presampling_multiplier=int(detection["presampling_multiplier"]),
        sampling=int(detection["sampling"]),
        sampling_method=str(detection["sampling_method"]),
        mle=int(detection["maximum_local_extrema"]),
        nb_step_min=int(detection["nb_step_min"]),
        nb_step_to_be_mle=int(detection["nb_step_to_be_mle"]),
        force_height_unit=str(detection["force_height_unit"]),
        force_speed_unit=str(detection["force_speed_unit"]),
    )
    records: list[dict[str, Any]] = []
    native: dict[str, Any] = {}
    for obs in observations:
        records.extend(_observation_records(obs, prepared, day, source_path))
        native[POLARITY_NAME[int(obs.sign_type)]] = obs
    frame = pd.DataFrame.from_records(records, columns=DAILY_COLUMNS)
    if frame.empty:
        frame = empty_daily_frame()
    else:
        frame = frame.sort_values(["polarity", "object_id"]).reset_index(drop=True)
    audit = {
        "date": day.isoformat(),
        "source_eta_path": str(source_path),
        "source_eta_shape_lat_lon": list(snapshot.eta_m.shape),
        "lon_size": int(snapshot.lon.size),
        "lat_size": int(snapshot.lat.size),
        "eta_masked_cells": int(snapshot.raw_mask.sum()),
        "eta_cells": int(snapshot.eta_m.size),
        "filter_valid_cells": int(prepared.filter_valid_mask.sum()),
        "filter_valid_fraction": float(prepared.filter_valid_mask.mean()),
        "objects_total": int(len(frame)),
        "objects_anticyclonic": int((frame["polarity_code"] == 1).sum())
        if len(frame)
        else 0,
        "objects_cyclonic": int((frame["polarity_code"] == -1).sum())
        if len(frame)
        else 0,
    }
    return frame, native, audit


def _contour_overlap(record_a: Mapping[str, Any], record_b: Mapping[str, Any]) -> float:
    """Compute PET polygon IoU for two effective contours."""

    lon_a = np.asarray(record_a["effective_contour_lon"], dtype="f8")
    lat_a = np.asarray(record_a["effective_contour_lat"], dtype="f8")
    lon_b = np.asarray(record_b["effective_contour_lon"], dtype="f8")
    lat_b = np.asarray(record_b["effective_contour_lat"], dtype="f8")
    if min(lon_a.size, lat_a.size, lon_b.size, lat_b.size) < 3:
        return 0.0
    # Cheap bounding-box rejection before PET's polygon operation.
    lon_b_wrapped = lon_b.copy()
    if abs(lon_a[0] - lon_b_wrapped[0]) > 180:
        lon_b_wrapped = (lon_b_wrapped - (lon_a[0] - 180)) % 360 + lon_a[0] - 180
    if (
        lon_a.max() < lon_b_wrapped.min()
        or lon_a.min() > lon_b_wrapped.max()
        or lat_a.max() < lat_b.min()
        or lat_a.min() > lat_b.max()
    ):
        return 0.0
    from py_eddy_tracker.poly import vertice_overlap

    score = vertice_overlap(
        np.asarray([lon_a]),
        np.asarray([lat_a]),
        np.asarray([lon_b]),
        np.asarray([lat_b]),
        min_overlap=0.0,
    )
    return float(score[0])


def _interpolate_list(a: Any, b: Any, fraction: float) -> list[float]:
    """Linearly interpolate fixed-size contour lists when possible."""

    first, second = np.asarray(a, dtype="f8"), np.asarray(b, dtype="f8")
    if first.size == second.size and first.size:
        return ((1.0 - fraction) * first + fraction * second).tolist()
    return (first if fraction < 0.5 else second).tolist()


def _interpolate_record(
    first: Mapping[str, Any], second: Mapping[str, Any], day: dt.date, fraction: float
) -> dict[str, Any]:
    """Construct one explicitly flagged virtual observation."""

    result = dict(first)
    result["date"] = day.isoformat()
    result["object_id"] = f"{first['track_id']}_V_{day:%Y%m%d}"
    result["is_virtual"] = True
    for column in NUMERIC_INTERPOLATION_COLUMNS:
        first_value, second_value = first.get(column), second.get(column)
        if pd.notna(first_value) and pd.notna(second_value):
            result[column] = (1.0 - fraction) * float(first_value) + fraction * float(
                second_value
            )
    for column in CONTOUR_COLUMNS:
        result[column] = _interpolate_list(first[column], second[column], fraction)
    result["filter_valid"] = bool(first["filter_valid"] and second["filter_valid"])
    result["boundary_censored"] = bool(
        first["boundary_censored"] or second["boundary_censored"]
    )
    result["coastal"] = bool(first["coastal"] or second["coastal"])
    result["source_eta_path"] = None
    return result


def _duration_class(duration_days: int) -> str:
    if duration_days <= 1:
        return "untracked"
    if duration_days < 10:
        return "short"
    if duration_days < 30:
        return "long"
    return "persistent"


def track_surface_eddies(
    daily_frames: Sequence[pd.DataFrame], config: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Track daily objects with PET effective-contour IoU and add virtual rows."""

    tracking = config["tracking"]
    overlap_min = float(tracking["minimum_overlap"])
    search_days = int(tracking["maximum_search_days"])
    max_virtual = int(tracking["maximum_consecutive_virtual"])
    actual = pd.concat(
        [frame.loc[~frame["is_virtual"].astype(bool)].copy() for frame in daily_frames],
        ignore_index=True,
    ) if daily_frames else empty_daily_frame()
    if actual.empty:
        return actual, empty_daily_frame(), pd.DataFrame(
            columns=[
                "track_id",
                "polarity",
                "polarity_code",
                "first_date",
                "last_date",
                "duration_days",
                "observed_days",
                "virtual_days",
                "duration_class",
                "boundary_censored_any",
                "coastal_any",
                "start_censored",
                "end_censored",
            ]
        )
    actual["date_dt"] = pd.to_datetime(actual["date"]).dt.date
    actual = actual.sort_values(["date_dt", "polarity", "object_id"]).reset_index(drop=True)
    next_track = {"anticyclonic": 1, "cyclonic": 1}
    states: dict[str, dict[str, Any]] = {}
    assignments: dict[int, str] = {}
    dates = sorted(actual["date_dt"].unique())
    for current_date in dates:
        current_indices = actual.index[actual["date_dt"] == current_date].tolist()
        active = {
            track_id: state
            for track_id, state in states.items()
            if (current_date - state["last_date"]).days <= search_days
        }
        candidates: list[tuple[float, str, int]] = []
        for index in current_indices:
            row = actual.loc[index]
            for track_id, state in active.items():
                if row["polarity"] != state["polarity"]:
                    continue
                gap = (current_date - state["last_date"]).days
                if gap < 1 or gap > search_days:
                    continue
                score = _contour_overlap(state["row"], row)
                if score >= overlap_min:
                    candidates.append((score, str(row["object_id"]), track_id))
        candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
        used_current: set[int] = set()
        used_tracks: set[str] = set()
        object_to_index = {
            str(actual.loc[index, "object_id"]): index for index in current_indices
        }
        for _, object_id, track_id in candidates:
            index = object_to_index[object_id]
            if index in used_current or track_id in used_tracks:
                continue
            assignments[index] = track_id
            used_current.add(index)
            used_tracks.add(track_id)
            states[track_id] = {
                "polarity": states[track_id]["polarity"],
                "last_date": current_date,
                "row": actual.loc[index].to_dict(),
            }
        for index in current_indices:
            if index in assignments:
                continue
            row = actual.loc[index]
            polarity = str(row["polarity"])
            track_id = f"{polarity[0].upper()}{next_track[polarity]:06d}"
            next_track[polarity] += 1
            assignments[index] = track_id
            states[track_id] = {
                "polarity": polarity,
                "last_date": current_date,
                "row": row.to_dict(),
            }
    actual["track_id"] = [assignments[index] for index in actual.index]
    actual = actual.drop(columns=["date_dt"])

    virtual_records: list[dict[str, Any]] = []
    for track_id, group in actual.groupby("track_id", sort=True):
        ordered = group.sort_values("date")
        rows = ordered.to_dict("records")
        for first, second in zip(rows[:-1], rows[1:]):
            first_date, second_date = parse_date(first["date"]), parse_date(second["date"])
            gap = (second_date - first_date).days
            if 1 < gap <= max_virtual + 1:
                for offset in range(1, gap):
                    virtual_records.append(
                        _interpolate_record(
                            first,
                            second,
                            first_date + dt.timedelta(days=offset),
                            offset / gap,
                        )
                    )
    virtual = pd.DataFrame.from_records(virtual_records, columns=DAILY_COLUMNS)
    if virtual.empty:
        virtual = empty_daily_frame()
    all_observations = pd.concat([actual, virtual], ignore_index=True, sort=False)
    all_observations = all_observations.sort_values(["track_id", "date", "is_virtual"])
    track_rows = []
    analysis_start = parse_date(config["input"]["start_date"])
    analysis_end = parse_date(config["input"]["end_date"])
    for track_id, group in all_observations.groupby("track_id", sort=True):
        group_dates = [parse_date(value) for value in group["date"]]
        observed = group.loc[~group["is_virtual"].astype(bool)]
        first_date, last_date = min(group_dates), max(group_dates)
        duration = (last_date - first_date).days + 1
        track_rows.append(
            {
                "track_id": track_id,
                "polarity": str(group["polarity"].iloc[0]),
                "polarity_code": int(group["polarity_code"].iloc[0]),
                "first_date": first_date.isoformat(),
                "last_date": last_date.isoformat(),
                "duration_days": duration,
                "observed_days": int((~group["is_virtual"].astype(bool)).sum()),
                "virtual_days": int(group["is_virtual"].astype(bool).sum()),
                "duration_class": _duration_class(duration),
                "boundary_censored_any": bool(group["boundary_censored"].astype(bool).any()),
                "coastal_any": bool(group["coastal"].astype(bool).any()),
                "start_censored": first_date == analysis_start,
                "end_censored": last_date == analysis_end,
            }
        )
    tracks = pd.DataFrame(track_rows)
    metadata_by_track = tracks.set_index("track_id").to_dict("index")
    for frame in (actual, all_observations):
        frame["track_observed_days"] = frame["track_id"].map(
            lambda value: metadata_by_track[value]["observed_days"]
        )
        frame["track_virtual_days"] = frame["track_id"].map(
            lambda value: metadata_by_track[value]["virtual_days"]
        )
        frame["track_duration_days"] = frame["track_id"].map(
            lambda value: metadata_by_track[value]["duration_days"]
        )
        frame["track_duration_class"] = frame["track_id"].map(
            lambda value: metadata_by_track[value]["duration_class"]
        )
        frame.drop(columns=[column for column in ("date_dt",) if column in frame], inplace=True)
    return (
        actual[list(DAILY_COLUMNS)].sort_values(["date", "polarity", "object_id"]).reset_index(drop=True),
        all_observations[list(DAILY_COLUMNS)]
        .sort_values(["track_id", "date", "is_virtual"])
        .reset_index(drop=True),
        tracks.sort_values("track_id").reset_index(drop=True),
    )


def validate_basic_helpers() -> dict[str, Any]:
    """Run dependency-free parser/unit/mask/translation checks."""

    source = np.ma.array([[0.0, 100.0], [200.0, 300.0]], mask=[[False, True], [False, False]])
    converted = convert_eta_to_m(source, "cm")
    if not np.isclose(converted[1, 0], 2.0) or not converted.mask[0, 1]:
        raise AssertionError("cm-to-m conversion or mask preservation failed")
    lon = np.arange(5.0)
    lat = np.arange(4.0)
    eta = np.zeros((lat.size, lon.size))
    _validate_coordinates(lon, lat, eta)
    try:
        _validate_coordinates(lon, lat, eta.T)
    except ValueError:
        swapped_rejected = True
    else:
        swapped_rejected = False
    if not swapped_rejected:
        raise AssertionError("dimension-swapped ETA was not rejected")

    x = np.arange(0.0, 20.0, 0.2)
    y = np.arange(20.0, 40.0, 0.2)

    def center(field: np.ndarray) -> tuple[float, float]:
        i, j = np.unravel_index(np.argmin(field), field.shape)
        return float(x[i]), float(y[j])

    first = -np.exp(-(((x[:, None] - 10.0) / 1.0) ** 2 + ((y[None, :] - 30.0) / 1.0) ** 2))
    second = -np.exp(-(((x[:, None] - 11.0) / 1.0) ** 2 + ((y[None, :] - 30.4) / 1.0) ** 2))
    first_center, second_center = center(first), center(second)
    if abs((second_center[0] - first_center[0]) - 1.0) > 0.21 or abs(
        (second_center[1] - first_center[1]) - 0.4
    ) > 0.21:
        raise AssertionError("analytic translated center did not translate as expected")
    return {
        "cm_to_m": True,
        "mask_preserved": True,
        "dimension_swap_rejected": True,
        "analytic_center_translation": True,
    }


def validate_tracking_helpers() -> dict[str, Any]:
    """Exercise PET IoU matching and the four-virtual-observation rule."""

    def make_row(day: str, object_id: str, shift: float) -> dict[str, Any]:
        contour_lon = [
            9.0 + shift,
            11.0 + shift,
            11.0 + shift,
            9.0 + shift,
            9.0 + shift,
        ]
        contour_lat = [29.0, 29.0, 31.0, 31.0, 29.0]
        return {
            "date": day,
            "object_id": object_id,
            "polarity": "anticyclonic",
            "polarity_code": 1,
            "center_lon": 10.0 + shift,
            "center_lat": 30.0,
            "ssh_extremum_lon": 10.0 + shift,
            "ssh_extremum_lat": 30.0,
            "effective_radius_km": 100.0,
            "speed_radius_km": 80.0,
            "amplitude_m": 0.01,
            "shape_error_effective_pct": 10.0,
            "shape_error_speed_pct": 10.0,
            "effective_area_m2": 1.0,
            "speed_area_m2": 1.0,
            "effective_contour_lon": contour_lon,
            "effective_contour_lat": contour_lat,
            "speed_contour_lon": contour_lon,
            "speed_contour_lat": contour_lat,
            "is_virtual": False,
            "boundary_censored": False,
            "filter_valid": True,
            "coastal": False,
            "source_eta_path": "synthetic",
            "track_id": None,
            "track_observed_days": np.nan,
            "track_virtual_days": np.nan,
            "track_duration_days": np.nan,
            "track_duration_class": None,
        }

    first = pd.DataFrame(
        [make_row("2003-01-01", "a", 0.0)], columns=DAILY_COLUMNS
    )
    second = pd.DataFrame(
        [make_row("2003-01-06", "b", 0.2)], columns=DAILY_COLUMNS
    )
    actual, observations, tracks = track_surface_eddies(
        [first, second],
        {
            "input": {
                "start_date": "2003-01-01",
                "end_date": "2003-01-06",
            },
            "tracking": {
                "minimum_overlap": 0.05,
                "maximum_search_days": 5,
                "maximum_consecutive_virtual": 4,
            }
        },
    )
    if tracks.shape[0] != 1 or tracks.iloc[0]["observed_days"] != 2:
        raise AssertionError("PET IoU synthetic objects did not form one track")
    if tracks.iloc[0]["virtual_days"] != 4 or len(observations) != 6:
        raise AssertionError("four-day virtual-observation rule failed")
    if observations["is_virtual"].astype(bool).sum() != 4:
        raise AssertionError("virtual observations were not flagged separately")
    return {
        "iou_match": True,
        "virtual_days": 4,
        "virtual_flag_separate": True,
    }
