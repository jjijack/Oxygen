import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.path import Path as MplPath
import geopandas as gpd
import inspect
from netCDF4 import Dataset
import re
import os
from pathlib import Path
import pickle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import glob
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
import copy
import gsw
from collections import defaultdict
import multiprocessing
import h5py
import time as tm
import shutil
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster, as_completed
from dask import delayed, compute
from dask.diagnostics import ProgressBar
from tqdm.auto import tqdm
import yaml
import json, pyarrow as pa, pyarrow.parquet as pq
import math
import zarr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from adjustText import adjust_text
# -------------------- 区域配置加载 --------------------
def _load_region_config(config_path: str | Path = 'config/regions.yml', region: str | None = None):
    """加载区域配置文件并返回指定区域字典。

    若未安装 PyYAML 或文件缺失，则回退到黑潮延伸体默认范围。
    """
    fallback = {
        'lon_min': 140 - 2.5,
        'lon_max': 180 + 2.5,
        'lat_min': 28 - 2.5,
        'lat_max': 40 + 2.5,
        'crosses_dateline': True,
        '_fallback': True
    }
    cfg_path = Path(config_path)
    if yaml is None or not cfg_path.exists():
        return fallback
    try:
        with open(cfg_path, 'r') as f:
            cfg = yaml.safe_load(f)
        if region is None:
            region = cfg.get('default_region')
        region_dict = cfg['regions'][region]
        # 统一键名，缺失时使用 fallback
        for k, v in fallback.items():
            region_dict.setdefault(k, v)
        region_dict['_fallback'] = False
        return region_dict
    except Exception:
        return fallback

_REGION_CFG = _load_region_config()
lonmin, lonmax = _REGION_CFG['lon_min'], _REGION_CFG['lon_max']
latmin, latmax = _REGION_CFG['lat_min'], _REGION_CFG['lat_max']

# -------------------- 数据路径与处理参数配置加载（Paths & Processing Config） --------------------
def _load_yaml(path: str | Path) -> dict:
    if yaml is None:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    try:
        with open(p, 'r') as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}

_PATHS_CFG = _load_yaml('config/paths.yml')
_PROC_CFG = _load_yaml('config/processing.yml')

argo_origin_path = Path(_PATHS_CFG.get('paths', {}).get('argo_origin', './Argo_origin'))
tmp_parquet_path = Path(_PATHS_CFG.get('paths', {}).get('argo_intermediate', './Argo_data_tmp'))
argo_path = Path(_PATHS_CFG.get('paths', {}).get('argo_parquet', './Argo_data'))
argo_mat_input_path = Path(_PATHS_CFG.get('paths', {}).get('argo_mat_input', './Argo_addFloat'))
Glorys_path = _PATHS_CFG.get('paths', {}).get('glorys_root', '../copernicus/GLORYS')
META_root_path = Path(_PATHS_CFG.get('paths', {}).get('meta_root', '../META3.2_DT_allsat'))
plots_output_root = Path(_PATHS_CFG.get('paths', {}).get('plots_output', './plot_outputs'))

# 用下划线隐藏内部配置值，提供 getter 避免随处写死名称
circle_enlargement_factor = float(
    _PROC_CFG.get('processing', {}).get('circle_enlargement_factor', 1.2)
)
# 已弃用的平均米/度常量删除，统一使用 approximate_degree_length 计算局地尺度
_default_delta_do_threshold = float(
    _PROC_CFG.get('processing', {}).get('default_delta_do_threshold', 50.0)
)
_default_salinity_threshold = float(
    _PROC_CFG.get('processing', {}).get('default_salinity_threshold', 0.0)
)
_default_temperature_threshold = float(
    _PROC_CFG.get('processing', {}).get('default_temperature_threshold', 0.0)
)
_default_depth_interval = float(
    _PROC_CFG.get('processing', {}).get('depth_interval', 100.0)
)
_default_depth_merge_tolerance = float(
    _PROC_CFG.get('processing', {}).get('depth_merge_tolerance', 10.0)
)
_default_duplicate_depth_strategy = _PROC_CFG.get('processing', {}).get('duplicate_depth_strategy', 'best_qc')
_cfg_anomaly_min_depth = float(
    _PROC_CFG.get('processing', {}).get('anomaly_min_depth', 300.0)
)
_default_adaptive_lat_threshold = float(
    _PROC_CFG.get('processing', {}).get('adaptive_lat_threshold', 70.0)
)
_default_adaptive_distance_threshold_km = float(
    _PROC_CFG.get('processing', {}).get('adaptive_distance_threshold_km', 300.0)
)

def _load_basemap_colors() -> dict:
    """加载底图颜色配置。

    未配置的键使用默认值；缺失 plot.basemap_colors 则全部取默认。
    """
    defaults = {
        'ocean': '#f8fafc',
        'land': '#d9d9d9',
        'coastline': '#787d85',
        'border': '#9aa0a8',
        'grid': '#b6bdc6',
    }
    if not isinstance(_PROC_CFG, dict):
        return defaults
    bm = _PROC_CFG.get('plot', {}).get('basemap_colors', {})
    merged = defaults.copy()
    if isinstance(bm, dict):
        merged.update({k: str(v) for k, v in bm.items() if k in defaults})
    return merged

_BASEMAP_COLORS = _load_basemap_colors()

# -------------------- 自带底图加载（使用本地 Natural Earth, 简洁版） --------------------
def _load_world_geodataframe():
    """从配置读取底图路径，默认 external/natural_earth/ne_110m_admin_0_countries.shp。"""
    cfg_paths = _PATHS_CFG.get('paths', {}) if isinstance(_PATHS_CFG, dict) else {}
    shp_path = Path(cfg_paths.get('gpd_world_shp', './external/natural_earth/ne_110m_admin_0_countries.shp'))
    return gpd.read_file(str(shp_path))

# -------------------------------------------------------------------------------

def print_current_processing_defaults():
        """打印当前生效的处理参数全局默认值（从 processing.yml 读取/回退）。

        包含：
            circle_enlargement_factor,
            distance_deg_per_meter,
            delta_do_threshold / salinity_threshold / temperature_threshold,
            depth_interval / depth_merge_tolerance / duplicate_depth_strategy。

        目的：调试与运行时确认当前配置，无返回值。
        """
        print("[Processing Defaults]")
        print(f"  circle_enlargement_factor : {circle_enlargement_factor}")
        print(f"  delta_do_threshold        : {_default_delta_do_threshold}")
        print(f"  salinity_threshold        : {_default_salinity_threshold}")
        print(f"  temperature_threshold     : {_default_temperature_threshold}")
        print(f"  depth_interval            : {_default_depth_interval}")
        print(f"  depth_merge_tolerance     : {_default_depth_merge_tolerance}")
        print(f"  duplicate_depth_strategy  : {_default_duplicate_depth_strategy}")
        print(f"  anomaly_min_depth         : {_cfg_anomaly_min_depth}")

def approximate_degree_length(lat: float | np.ndarray, lon: float | np.ndarray | None = None) -> dict:
    """计算指定纬度（可选经度）处经纬度与距离的近似换算关系。

    采用 WGS84 椭球常用的近似级数公式 (单位: 米/度)，适用于绝大多数海洋学分析精度需求。

    公式来源（展开项保留到 cos(5φ)/cos(6φ)）：
        meters_per_degree_lat ≈ 111132.92 - 559.82*cos(2φ) + 1.175*cos(4φ) - 0.0023*cos(6φ)
        meters_per_degree_lon ≈ 111412.84*cos(φ) - 93.5*cos(3φ) + 0.118*cos(5φ)

    参数:
        lat (float | np.ndarray): 纬度（度）。可为标量或 numpy 数组。
        lon (float | np.ndarray | None): 经度（度）。对当前计算不影响，只是为了接口对称；
            可传入与 lat 同形状数组（将被忽略）。保留此参数便于未来扩展（比如考虑地形加权等）。

    返回:
        dict: 包含以下键：
            meters_per_degree_lat: 指定纬度上一度纬差对应的米数
            meters_per_degree_lon: 指定纬度上一度经差对应的米数
            degrees_per_meter_lat: 上述量的倒数（度/米）
            degrees_per_meter_lon: 上述量的倒数（度/米）

    备注:
        1. 若输入为数组，则返回值各字段为同形状 numpy 数组。
        2. 该函数提供纬度依赖的更精细米/度估计，替代旧的单一平均值。
        3. 用于将“米单位半径”换算到“角度半径”时，推荐： radius_deg_lat = radius_m / meters_per_degree_lat。
    """
    # 转换为 numpy 数组以统一处理
    lat_arr = np.asarray(lat, dtype=float)
    phi = np.deg2rad(lat_arr)

    # 分别计算纬度与经度方向一度的米数（WGS84 近似）
    meters_per_degree_lat = (
        111132.92
        - 559.82 * np.cos(2 * phi)
        + 1.175 * np.cos(4 * phi)
        - 0.0023 * np.cos(6 * phi)
    )
    meters_per_degree_lon = (
        111412.84 * np.cos(phi)
        - 93.5 * np.cos(3 * phi)
        + 0.118 * np.cos(5 * phi)
    )

    # 倒数（度/米）。对 0 做防护：若 cos(phi) ~ 0（极点），经向一度长度 → 0，避免除零设为 np.inf
    with np.errstate(divide='ignore', invalid='ignore'):
        degrees_per_meter_lat = 1.0 / meters_per_degree_lat
        degrees_per_meter_lon = np.where(meters_per_degree_lon != 0, 1.0 / meters_per_degree_lon, np.inf)

    result = {
        'meters_per_degree_lat': meters_per_degree_lat,
        'meters_per_degree_lon': meters_per_degree_lon,
        'degrees_per_meter_lat': degrees_per_meter_lat,
        'degrees_per_meter_lon': degrees_per_meter_lon,
    }

    # 若输入是标量（非数组或0维），将结果中对应项转换为原生 float，避免下游出现 0-d ndarray 序列化/打印差异
    if np.asarray(lat).shape == ():
        for k, v in result.items():
            # v 可能是 ndarray 标量或 Python float
            if isinstance(v, np.ndarray) and v.shape == ():
                result[k] = v.item()
    return result

def _minimal_lon_diff_deg(lon: float | np.ndarray, lon0: float) -> np.ndarray:
    """计算经度差并映射到 (-180, 180] 区间，支持标量或数组。

    用于跨日界线区域，保证 179.8° 与 -179.7° 的差为 0.5° 而不是 -359.5°。
    """
    d = np.asarray(lon) - lon0
    return (d + 180.0) % 360.0 - 180.0

def _normalize_lon_array(val: np.ndarray | float) -> np.ndarray:
    """将任意经度值映射到 (-180, 180] 区间，返回 np.ndarray。"""
    arr = np.asarray(val, dtype=float)
    return (arr + 180.0) % 360.0 - 180.0

def _region_lon_mask(lon_vals: np.ndarray, lon_min_cfg: float, lon_max_cfg: float) -> np.ndarray:
    """根据区域经度范围（允许跨日界线）生成布尔掩码。"""
    lon_arr = np.asarray(lon_vals, dtype=float)
    if lon_arr.size == 0:
        return np.zeros(0, dtype=bool)

    lon_min_norm = float(_normalize_lon_array(lon_min_cfg))
    lon_max_norm = float(_normalize_lon_array(lon_max_cfg))

    raw_span = abs(float(lon_max_cfg) - float(lon_min_cfg))
    eff_span = (lon_max_norm - lon_min_norm) % 360.0
    is_global_lon = (
        (raw_span >= 359.5)
        or (eff_span >= 359.5)
        or np.isclose(eff_span, 0.0, atol=1e-6)
    )
    if is_global_lon:
        return np.ones(lon_arr.shape, dtype=bool)

    lon_norm = _normalize_lon_array(lon_arr)
    if lon_min_norm <= lon_max_norm:
        return (lon_norm >= lon_min_norm) & (lon_norm <= lon_max_norm)
    return (lon_norm >= lon_min_norm) | (lon_norm <= lon_max_norm)

def local_xy_distance_m(lon: float | np.ndarray, lat: float | np.ndarray,
                        lon0: float, lat0: float, wrap_dateline: bool = True) -> np.ndarray:
    """计算点 (lon,lat) 到参考点 (lon0,lat0) 的局地平面近似距离（米）。

    参数:
        lon, lat: 点的经纬度，可为标量或数组（广播到与 lon 相同形状）。
        lon0, lat0: 参考中心（标量）。
        wrap_dateline: 是否对经度差进行跨日界线 (±180°) 最短差处理。

    说明:
        1. 使用纬度依赖的经/纬一度长度（WGS84 近似）。在中低纬、距离 <~500 km 下平面近似足够。
        2. 若距离很大或靠近极区，平面近似误差增大，可考虑改用大圆距离（后续可扩展）。
        3. wrap_dateline=True 时能正确处理 179.9° 与 -179.9° 仅 0.2° 之差的情况。
    """
    scale = approximate_degree_length(lat0)
    m_per_deg_lat = scale['meters_per_degree_lat']
    m_per_deg_lon = scale['meters_per_degree_lon']
    dlon = _minimal_lon_diff_deg(lon, lon0) if wrap_dateline else (np.asarray(lon) - lon0)
    dlat = np.asarray(lat) - lat0
    dx_m = dlon * m_per_deg_lon
    dy_m = dlat * m_per_deg_lat
    return np.hypot(dx_m, dy_m)

def great_circle_distance_m(lon: float | np.ndarray, lat: float | np.ndarray,
                            lon0: float, lat0: float, wrap_dateline: bool = True,
                            radius_earth_m: float = 6371000.0) -> np.ndarray:
    """计算球面大圆距离（Haversine 公式），单位: 米。

    设计目标: 简单、稳定、无自动切换逻辑；与 local_xy_distance_m 并存，供需要更精确/大尺度/高纬场景手动调用。

    参数:
        lon, lat : 目标点经纬度（标量或一维数组）。
        lon0, lat0 : 中心点经纬度（标量）。
        wrap_dateline : True 时对经度差做跨日界线最短差归一（±180°）。
        radius_earth_m : 地球半径（可根据需要调整为更精确椭球平均半径）。

    返回:
        与输入 (lon, lat) 形状一致的距离（米）。标量输入 → 标量输出。

    说明:
        Haversine 公式: 
            a = sin²(Δφ/2) + cos φ1 * cos φ2 * sin²(Δλ/2)
            c = 2 * asin( sqrt(a) )
            d = R * c
        在中短距离下与球面真值非常接近；对极区与大尺度优于局地平面近似。
    """
    lon_arr = np.asarray(lon, dtype=float)
    lat_arr = np.asarray(lat, dtype=float)
    dlon_deg = _minimal_lon_diff_deg(lon_arr, lon0) if wrap_dateline else (lon_arr - lon0)
    dlat_deg = lat_arr - lat0
    dlon = np.deg2rad(dlon_deg)
    dlat = np.deg2rad(dlat_deg)
    lat1 = np.deg2rad(lat0)
    lat2 = np.deg2rad(lat_arr)
    sin_dlat = np.sin(dlat / 2.0)
    sin_dlon = np.sin(dlon / 2.0)
    a = sin_dlat**2 + np.cos(lat1) * np.cos(lat2) * sin_dlon**2
    # 数值稳定保护：浮点误差可能导致 a 略超 1
    c = 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    dist = radius_earth_m * c
    if np.ndim(dist) == 0:
        return float(dist)
    return dist

def adaptive_distance_m(
    lon: float | np.ndarray,
    lat: float | np.ndarray,
    lon0: float | np.ndarray,
    lat0: float | np.ndarray,
    wrap_dateline: bool = True,
    gc_lat_threshold: float | None = None,
    gc_distance_threshold_km: float | None = None,
    force_great_circle: bool = False,
    radius_earth_m: float = 6371000.0
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """自适应距离 (米)。在保持平面估算速度的前提下，自动在高纬或大尺度条件改用大圆距离。

    策略:
        1. force_great_circle=True 时直接使用大圆（适合全球/大范围批处理）。
        2. 若 |lat0| >= gc_lat_threshold → 直接大圆，跳过平面计算。
        3. 否则先计算一次平面近似 (planar)；若最大平面距离 > gc_distance_threshold_km → 改用大圆；否则保留平面结果。

    使用建议:
        - 全球或大量远距点：直接 force_great_circle=True 减少双重计算。
        - 中低纬涡旋附近（大多数点在半径数倍内）：保持默认，可获得接近平面距离的速度与足够精度。
        - 极区或大半径场景对精度敏感：force_great_circle=True。

    参数:
        lon, lat : 目标点（标量或数组）。
        lon0, lat0 : 中心经纬度（标量或数组，可与目标点一一对应广播）。
        gc_lat_threshold : 高纬触发大圆的绝对纬度阈值；None 时使用配置文件 adaptive_lat_threshold。
        gc_distance_threshold_km : 平面最大距离超过该值 (km) 触发大圆；None 时使用配置文件 adaptive_distance_threshold_km。
        force_great_circle : 强制使用大圆。
        radius_earth_m : 地球半径。

    返回:
        与 (lon, lat) 形状一致的距离；标量输入返回 float。
    """
    if gc_lat_threshold is None:
        gc_lat_threshold = _default_adaptive_lat_threshold
    if gc_distance_threshold_km is None:
        gc_distance_threshold_km = _default_adaptive_distance_threshold_km

    lon_arr = np.asarray(lon, dtype=float)
    lat_arr = np.asarray(lat, dtype=float)
    lon0_arr = np.asarray(lon0, dtype=float)
    lat0_arr = np.asarray(lat0, dtype=float)

    # 广播到统一形状，便于处理逐点中心（例如逐日涡旋中心）
    try:
        broadcast_shape = np.broadcast(lon_arr, lat_arr, lon0_arr, lat0_arr).shape
    except ValueError as exc:
        raise ValueError("adaptive_distance_m: inputs cannot be broadcast to a common shape") from exc
    lon_arr = np.broadcast_to(lon_arr, broadcast_shape)
    lat_arr = np.broadcast_to(lat_arr, broadcast_shape)
    lon0_arr = np.broadcast_to(lon0_arr, broadcast_shape)
    lat0_arr = np.broadcast_to(lat0_arr, broadcast_shape)

    # (1) 强制模式：直接大圆
    if force_great_circle:
        gc = great_circle_distance_m(lon_arr, lat_arr, lon0_arr, lat0_arr, wrap_dateline=wrap_dateline, radius_earth_m=radius_earth_m)
        return gc

    # (2) 高纬提前判定：无需先算平面
    if gc_lat_threshold is not None and np.any(np.abs(lat0_arr) >= gc_lat_threshold):
        gc = great_circle_distance_m(lon_arr, lat_arr, lon0_arr, lat0_arr, wrap_dateline=wrap_dateline, radius_earth_m=radius_earth_m)
        return gc

    # (3) 需要距离阈值判定，才计算平面距离
    scale = approximate_degree_length(lat0_arr)
    dlon = _minimal_lon_diff_deg(lon_arr, lon0_arr) if wrap_dateline else (lon_arr - lon0_arr)
    dlat = lat_arr - lat0_arr
    planar = np.hypot(
        dlon * scale['meters_per_degree_lon'],
        dlat * scale['meters_per_degree_lat']
    )

    if gc_distance_threshold_km is not None and np.max(planar) > gc_distance_threshold_km * 1000.0:
        gc = great_circle_distance_m(lon_arr, lat_arr, lon0_arr, lat0_arr, wrap_dateline=wrap_dateline, radius_earth_m=radius_earth_m)
        return gc

    if np.ndim(planar) == 0:
        return float(planar)
    return planar

def inside_radius_m_adaptive(
    lon: float | np.ndarray,
    lat: float | np.ndarray,
    lon0: float,
    lat0: float,
    radius_m: float | np.ndarray,
    enlarge: float = 1.0,
    wrap_dateline: bool = True,
    gc_lat_threshold: float | None = None,
    gc_distance_threshold_km: float | None = None,
    force_great_circle: bool = False
) -> np.ndarray | bool:
    """自适应包含判定：与 adaptive_distance_m 策略一致。"""
    dist = adaptive_distance_m(
        lon, lat, lon0, lat0,
        wrap_dateline=wrap_dateline,
        gc_lat_threshold=gc_lat_threshold,
        gc_distance_threshold_km=gc_distance_threshold_km,
        force_great_circle=force_great_circle,
    )
    r_arr = np.asarray(radius_m, dtype=float)
    if r_arr.shape not in ((), np.shape(dist)):
        r_arr = np.broadcast_to(r_arr, np.shape(dist))
    inside = dist <= r_arr * enlarge
    if np.shape(inside) == ():
        return bool(inside)
    return inside

def ellipse_patch_for_eddy(lon0: float, lat0: float, radius_m: float, enlarge: float = 1.0,
                           **kwargs):
    """构造在经纬度坐标下表示真实米尺度涡旋半径的椭圆补丁。

    在 lon-lat 图上，同一米半径在经向与纬向角度跨度不同（lon 方向按 cos(lat) 收缩），
    因此绘制为椭圆 (width=Δlon, height=Δlat)。
    """
    scale = approximate_degree_length(lat0)
    m_per_deg_lat = scale['meters_per_degree_lat']
    m_per_deg_lon = scale['meters_per_degree_lon']
    radius_eff = radius_m * enlarge
    dlat = radius_eff / m_per_deg_lat
    dlon = radius_eff / m_per_deg_lon
    default_kwargs = dict(edgecolor='gray', facecolor='none', linestyle='--', linewidth=1.0, zorder=3)
    default_kwargs.update(kwargs)
    return Ellipse((lon0, lat0), width=2*dlon, height=2*dlat, **default_kwargs)

def switch_region(region_name: str, config_path: str | Path = 'config/regions.yml'):
    """在运行时切换默认区域（无需改 YAML），并刷新全局经纬度。

    更新内容: `lonmin, lonmax, latmin, latmax, _REGION_CFG`。

    推荐用法 (确保后续访问得到最新值):
        import track
        track.switch_region('global')
        print(track.lonmin, track.lonmax)

    不推荐:
        from track import lonmin  # 后续 switch_region 不会自动更新这个已绑定的数值副本

    参数:
        region_name: 在 regions.yml 中定义的区域 key。
        config_path: 配置文件路径，默认 'config/regions.yml'。

    异常:
        KeyError: 区域未找到或加载失败（进入 fallback）。
    """
    global _REGION_CFG, lonmin, lonmax, latmin, latmax
    new_cfg = _load_region_config(config_path=config_path, region=region_name)
    if new_cfg.get('_fallback'):
        raise KeyError(f"Region '{region_name}' not found or config load failed; fallback config in use.")
    _REGION_CFG = new_cfg
    lonmin, lonmax = _REGION_CFG['lon_min'], _REGION_CFG['lon_max']
    latmin, latmax = _REGION_CFG['lat_min'], _REGION_CFG['lat_max']
    if _REGION_CFG.get('crosses_dateline') and (lonmax < lonmin):
        print(f"[RegionConfig] Region '{region_name}' crosses dateline; implement split-range filtering if needed.")
    else:
        print(f"[RegionConfig] Switched to region '{region_name}': lon[{lonmin}, {lonmax}], lat[{latmin}, {latmax}]")

def load_meta_data(path: str | os.PathLike | None = None, version: float = 3.2):
    '''
    加载 META 涡旋数据，返回 (ACS, ACL, CS, CL)。

    参数:
        path : META 数据根目录；None 时使用配置文件 `paths.meta_root` (默认 '../META3.2_DT_allsat')。
        version : 3.1 或 3.2。

    目录期望包含相应命名模式的 NetCDF 文件，例如 (version=3.2):
        META3.2_DT_allsat_Anticyclonic_short_*.nc
        META3.2_DT_allsat_Anticyclonic_long_*.nc
        META3.2_DT_allsat_Cyclonic_short_*.nc
        META3.2_DT_allsat_Cyclonic_long_*.nc
    '''
    if path is None:
        path = META_root_path
    else:
        path = Path(path)
    if version == 3.2:
        ACS= Dataset(path / 'META3.2_DT_allsat_Anticyclonic_short_19930101_20220209.nc')
        ACL= Dataset(path / 'META3.2_DT_allsat_Anticyclonic_long_19930101_20220209.nc')
        CS= Dataset(path / 'META3.2_DT_allsat_Cyclonic_short_19930101_20220209.nc')
        CL= Dataset(path / 'META3.2_DT_allsat_Cyclonic_long_19930101_20220209.nc')
    elif version == 3.1:
        ACS=Dataset(path / 'META3.1exp_DT_allsat_Anticyclonic_short_19930101_20200307.nc')
        ACL=Dataset(path / 'META3.1exp_DT_allsat_Anticyclonic_long_19930101_20200307.nc')
        CS=Dataset(path / 'META3.1exp_DT_allsat_Cyclonic_short_19930101_20200307.nc')
        CL=Dataset(path / 'META3.1exp_DT_allsat_Cyclonic_long_19930101_20200307.nc')
    else:
        raise ValueError("Unsupported version. Please use 3.1 or 3.2.")

    return ACS, ACL, CS, CL


def area_limit(DS, latmin, latmax, lonmin, lonmax):
    '''
    限制区域范围
    
    输出：涡旋序号，时间，中心点经度，中心点纬度，最值点经度，最值点纬度，边界经度，边界纬度，半径，速度边界经度，速度边界纬度
    '''
    time = DS.variables['time'][:].data
    if np.issubdtype(time.dtype, np.floating):
        # 如果时间是浮点数，转换为整数
        time = np.round(time).astype(np.uint32)
    # time = convert_date(time)
    center_lon = DS.variables['longitude'][:].data
    center_lat = DS.variables['latitude'][:].data
    max_lon = DS.variables['longitude_max'][:].data
    max_lat = DS.variables['latitude_max'][:].data
    contour_lon = DS.variables['effective_contour_longitude'][:].data
    contour_lat = DS.variables['effective_contour_latitude'][:].data
    radius = DS.variables['effective_radius'][:].data
    track = DS.variables['track'][:].data
    speed_contour_lon = DS.variables['speed_contour_longitude'][:].data
    speed_contour_lat = DS.variables['speed_contour_latitude'][:].data

    mask = (center_lat >= latmin) & (center_lat <= latmax) & \
           (center_lon >= lonmin) & (center_lon <= lonmax)
    indices = np.nonzero(mask)[0]
    ds = []
    current_track = None
    current_list = []
    for i in indices:
       if current_track is None:
            current_track = track[i]
       # 当track变化时，保存之前的列表，并新建一个列表
       if track[i] != current_track:
            ds.append(current_list)
            current_list = []
            current_track = track[i]
       current_list.append([i,time[i], center_lon[i], center_lat[i], 
              max_lon[i], max_lat[i], contour_lon[i], contour_lat[i], 
              radius[i], speed_contour_lon[i], speed_contour_lat[i]])
    if current_list:
              ds.append(current_list)
    return ds

def _ensure_meta_tracks_root(root: str | Path | None = None) -> Path:
    """获取/创建 META_tracks 根目录。优先读取配置 `paths.META_tracks_root`，兼容 `paths.meta_tracks_root`；默认 `./META_tracks`。"""
    default_root = Path('./META_tracks')
    paths_cfg = _PATHS_CFG.get('paths', {})
    cfg_val = paths_cfg.get('META_tracks_root', default_root)
    cfg_root = Path(cfg_val)
    if root is not None:
        cfg_root = Path(root)
    cfg_root.mkdir(parents=True, exist_ok=True)
    return cfg_root

def _current_region_key() -> str:
    # regions.yml 中的 key 并未直接保存在 _REGION_CFG；此处基于 fallback 标志无法反推 key。
    # 方案：若 _REGION_CFG 包含 name 且不是 fallback，取 name 的 slug；否则 generic。
    name = _REGION_CFG.get('name') or 'region'
    slug = name.lower().replace(' ', '_')
    return slug

def export_meta_tracks(ds: Dataset,
                       kind: str,
                       region_key: str | None = None,
                       output_root: str | Path | None = None,
                       chunk_size: int = 100_000,
                       write_contours: bool = False,
                       build_track_summary: bool = True,
                       keep_legacy_pickle: bool = False,
                       use_dask: bool = False,
                       dask_num_workers: int | None = None,
                       compact_after: bool = True) -> dict:
    """流式导出 META 涡旋数据为标准化拆表格式（Parquet + Zarr）。

    输出：
        - <kind>_daily.parquet   : 每个涡旋每日一行（标量 & 计数）
        - <kind>_contours.zarr   : 轮廓顶点按 1D 扁平数组存储，含 prefix index（N+1）
        - <kind>_contours.parquet: (可选) 展开轮廓顶点逐行记录（用于 SQL/分析）
        - <kind>_tracks.parquet  : (可选) 汇总轨迹层（聚合统计）
        - <kind>_metadata.json   : 基本元信息（区域、生成时间、列说明、Zarr 路径）

    参数:
        ds : netCDF4.Dataset (已打开的 META 文件)
        kind : 标识（如 'acs','acl','cs','cl'）用于文件前缀
        region_key : 区域 key / slug（默认基于当前 _REGION_CFG['name'] 生成）
        output_root : 根输出目录（默认 ./META_tracks 或配置 paths.META_tracks_root/meta_tracks_root）
        chunk_size : indices 分块大小（按满足区域筛选的记录索引数量）
        write_contours : 是否同时写 contours 的 Parquet 拆表（顶点逐行），Zarr 始终会生成
        build_track_summary : 是否生成轨迹汇总表
        keep_legacy_pickle : 是否额外写出旧嵌套结构 pickle（用于临时兼容调试）
        use_dask : 是否用 Dask 并行按块处理
        dask_num_workers : Dask 并行进程数
        compact_after : 是否在完成后直接将目录形式的 Parquet 压实为单文件，并删除目录

    返回:
        dict: 包含写出的关键文件路径（如 daily_file/daily_dir、tracks_path、contours_zarr、metadata 等）。
    """
    # 记录总长度用于分块；避免一次性加载变量
    try:
        total_len = ds.dimensions['time'].size
    except Exception:
        total_len = ds.variables['time'].shape[0]

    region_slug = region_key or _current_region_key()
    root = _ensure_meta_tracks_root(output_root)
    region_dir = root / region_slug
    region_dir.mkdir(parents=True, exist_ok=True)

    # 输出路径：分片写入使用临时目录，压实后生成同名 .parquet 单文件，避免目录/文件冲突
    daily_stage = region_dir / f"{kind}_daily_tmp"
    contours_stage = region_dir / f"{kind}_contours_tmp"
    contours_zarr_path = region_dir / f"{kind}_contours.zarr"
    contours_parts_dir = region_dir / f"{kind}_contours_zarr_parts"
    tracks_path = region_dir / f"{kind}_tracks.parquet"
    meta_path = region_dir / f"{kind}_metadata.json"
    legacy_pickle_path = region_dir / f"{kind}_legacy_nested.pkl"

    # 清理旧的分片临时目录（若存在）
    if daily_stage.exists():
        if daily_stage.is_dir():
            shutil.rmtree(daily_stage)
        else:
            daily_stage.unlink()
    if write_contours and contours_stage.exists():
        if contours_stage.is_dir():
            shutil.rmtree(contours_stage)
        else:
            contours_stage.unlink()
    # 清理旧 Zarr 目标和分片目录
    if contours_zarr_path.exists():
        shutil.rmtree(contours_zarr_path)
    if contours_parts_dir.exists():
        shutil.rmtree(contours_parts_dir)

    daily_stage.mkdir(parents=True, exist_ok=True)
    if write_contours:
        contours_stage.mkdir(parents=True, exist_ok=True)
    contours_parts_dir.mkdir(parents=True, exist_ok=True)

    daily_schema = pa.schema([
        pa.field('track_id', pa.int64()),
        pa.field('time', pa.timestamp('ns')),
        pa.field('center_lon', pa.float64()),
        pa.field('center_lat', pa.float64()),
        pa.field('max_lon', pa.float64()),
        pa.field('max_lat', pa.float64()),
        pa.field('radius', pa.float64()),
        pa.field('contour_vertex_count', pa.int32()),
        pa.field('speed_contour_vertex_count', pa.int32()),
        pa.field('orig_index', pa.int64()),
    ])

    contour_schema = pa.schema([
        pa.field('track_id', pa.int64()),
        pa.field('time', pa.timestamp('ns')),
        pa.field('kind', pa.string()),  # 'effective' or 'speed'
        pa.field('vertex_index', pa.int32()),
        pa.field('lon', pa.float64()),
        pa.field('lat', pa.float64()),
    ]) if write_contours else None

    # 构造 chunk 写入
    def write_daily_chunk(tbl: pa.Table, idx: int):
        pq.write_table(tbl, daily_stage / f"part_{idx:05d}.parquet")

    def write_contour_chunk(tbl: pa.Table, idx: int, suffix: str):
        pq.write_table(tbl, contours_stage / f"part_{suffix}_{idx:05d}.parquet")

    # 旧结构临时兼容：按轨迹分组的嵌套列表（外层按 track，内层按日）
    # 为减少内存占用，使用 dict 累积，每个 track_id -> list[day_item]
    # day_item 结构与 area_limit 输出一致：[orig_i, time_YYYYMMDD, center_lon, center_lat, max_lon, max_lat, contour_lon[], contour_lat[], radius, speed_contour_lon[], speed_contour_lat[]]
    legacy_nested = {} if (keep_legacy_pickle and not use_dask) else None

    # 自动决定是否跳过区域筛选：当配置等同于全球范围时跳过
    is_global_lon = (lonmin <= -179.999) and (lonmax >= 179.999)
    is_global_lat = (latmin <= -89.999) and (latmax >= 89.999)
    auto_skip = bool(is_global_lon and is_global_lat)

    # 单块处理函数：可被串行调用，或用 dask.delayed 并发调用
    def _process_chunk(file_path: str, part_idx: int, start: int, end: int) -> int:
        ds_local = Dataset(file_path)
        try:
            # 读取该块变量
            time_var = ds_local.variables['time'][start:end]
            center_lon = ds_local.variables['longitude'][start:end]
            center_lat = ds_local.variables['latitude'][start:end]
            n_chunk = len(center_lon)
            # 区域过滤（单块）；若跳过则全选
            if auto_skip:
                sel = np.arange(n_chunk, dtype=int)
            else:
                mask = (center_lat >= latmin) & (center_lat <= latmax) & \
                       (center_lon >= lonmin) & (center_lon <= lonmax)
                if not np.any(mask):
                    return 0
                sel = np.nonzero(mask)[0]
            # 仅在需要时读取其它变量
            max_lon = ds_local.variables['longitude_max'][start:end]
            max_lat = ds_local.variables['latitude_max'][start:end]
            radius = ds_local.variables['effective_radius'][start:end]
            track_id = ds_local.variables['track'][start:end]
            eff_lon = ds_local.variables['effective_contour_longitude'][start:end]
            eff_lat = ds_local.variables['effective_contour_latitude'][start:end]
            spd_lon = ds_local.variables['speed_contour_longitude'][start:end]
            spd_lat = ds_local.variables['speed_contour_latitude'][start:end]

            # 时间转换到 datetime64[ns]
            time_dt = convert_date(pd.Series(time_var)).to_numpy().astype('datetime64[ns]')

            # 组装并写出
            rows = []
            eff_rows = []
            spd_rows = []
            # Zarr 分片缓存
            e_counts_list = []
            s_counts_list = []
            eff_concat_lon = []
            eff_concat_lat = []
            spd_concat_lon = []
            spd_concat_lat = []
            for off in sel:
                orig_i = start + int(off)
                v_track = int(track_id[off])
                v_time = time_dt[off]
                c_lon = float(center_lon[off])
                c_lat = float(center_lat[off])
                mx_lon = float(max_lon[off])
                mx_lat = float(max_lat[off])
                rad = float(radius[off])
                e_lon = eff_lon[off]
                e_lat = eff_lat[off]
                s_lon = spd_lon[off]
                s_lat = spd_lat[off]
                e_cnt = len(e_lon) if hasattr(e_lon, '__len__') else 0
                s_cnt = len(s_lon) if hasattr(s_lon, '__len__') else 0
                rows.append((v_track, v_time, c_lon, c_lat, mx_lon, mx_lat, rad, e_cnt, s_cnt, int(orig_i)))

                if write_contours and e_cnt:
                    for vi in range(e_cnt):
                        eff_rows.append((v_track, v_time, 'effective', vi, float(e_lon[vi]), float(e_lat[vi])))
                if write_contours and s_cnt:
                    for vi in range(s_cnt):
                        spd_rows.append((v_track, v_time, 'speed', vi, float(s_lon[vi]), float(s_lat[vi])))

                # 记录 Zarr 分片内容（每日计数与扁平顶点）
                e_counts_list.append(int(e_cnt))
                s_counts_list.append(int(s_cnt))
                if e_cnt:
                    eff_concat_lon.append(np.asarray(e_lon, dtype='f8'))
                    eff_concat_lat.append(np.asarray(e_lat, dtype='f8'))
                if s_cnt:
                    spd_concat_lon.append(np.asarray(s_lon, dtype='f8'))
                    spd_concat_lat.append(np.asarray(s_lat, dtype='f8'))

                if legacy_nested is not None:
                    time_ymd = int(str(v_time)[:10].replace('-', ''))
                    item = [
                        int(orig_i), time_ymd, c_lon, c_lat,
                        mx_lon, mx_lat, e_lon, e_lat, rad, s_lon, s_lat
                    ]
                    legacy_nested.setdefault(v_track, []).append(item)

            if rows:
                daily_tbl = pa.Table.from_pylist([
                    {
                        'track_id': r[0], 'time': r[1], 'center_lon': r[2], 'center_lat': r[3],
                        'max_lon': r[4], 'max_lat': r[5], 'radius': r[6],
                        'contour_vertex_count': r[7], 'speed_contour_vertex_count': r[8], 'orig_index': r[9]
                    } for r in rows
                ], schema=daily_schema)
                write_daily_chunk(daily_tbl, part_idx)

            if write_contours:
                if eff_rows:
                    eff_tbl = pa.Table.from_pylist([
                        {
                            'track_id': r[0], 'time': r[1], 'kind': r[2], 'vertex_index': r[3], 'lon': r[4], 'lat': r[5]
                        } for r in eff_rows
                    ], schema=contour_schema)
                    write_contour_chunk(eff_tbl, part_idx, 'eff')
                if spd_rows:
                    spd_tbl = pa.Table.from_pylist([
                        {
                            'track_id': r[0], 'time': r[1], 'kind': r[2], 'vertex_index': r[3], 'lon': r[4], 'lat': r[5]
                        } for r in spd_rows
                    ], schema=contour_schema)
                    write_contour_chunk(spd_tbl, part_idx, 'spd')

            # 将本分块的 Zarr 数据写为 NPZ 分片，供后续合并
            if rows:
                e_counts = np.asarray(e_counts_list, dtype='i8')
                s_counts = np.asarray(s_counts_list, dtype='i8')
                eff_lon_concat = np.concatenate(eff_concat_lon) if eff_concat_lon else np.empty((0,), dtype='f8')
                eff_lat_concat = np.concatenate(eff_concat_lat) if eff_concat_lat else np.empty((0,), dtype='f8')
                spd_lon_concat = np.concatenate(spd_concat_lon) if spd_concat_lon else np.empty((0,), dtype='f8')
                spd_lat_concat = np.concatenate(spd_concat_lat) if spd_concat_lat else np.empty((0,), dtype='f8')
                np.savez_compressed(
                    contours_parts_dir / f"part_{part_idx:05d}.npz",
                    eff_lon=eff_lon_concat,
                    eff_lat=eff_lat_concat,
                    eff_counts=e_counts,
                    spd_lon=spd_lon_concat,
                    spd_lat=spd_lat_concat,
                    spd_counts=s_counts,
                )
            return len(rows)
        finally:
            try:
                ds_local.close()
            except Exception:
                pass

    # 生成块边界
    file_path = ds.filepath() if hasattr(ds, 'filepath') else None
    if file_path is None:
        # netCDF4.Dataset 在某些场景没有 filepath 属性，尝试从 repr 中提取失败则报错
        raise RuntimeError("Dataset filepath unavailable; Dask/streaming requires file path.")

    # 在执行前，根据数据量与目标并发自动调整 chunk_size，避免只有 1 个分块导致只能单核
    expected_workers = 1
    if use_dask:
        expected_workers = int(dask_num_workers) if dask_num_workers else max(1, (os.cpu_count() or 1))
    parts = max(1, (total_len + chunk_size - 1) // chunk_size)
    if use_dask and parts < expected_workers:
        target_parts = max(expected_workers * 3, 2)
        new_chunk = max(25_000, math.ceil(total_len / target_parts))
        if new_chunk != chunk_size:
            print(f"[export_meta_tracks] Auto-tune chunk_size: {chunk_size} -> {new_chunk} (total={total_len}, workers={expected_workers}, target_parts≈{target_parts})")
            chunk_size = new_chunk
            parts = max(1, (total_len + chunk_size - 1) // chunk_size)

    print(f"[export_meta_tracks] Plan kind={kind}: total={total_len}, chunk_size={chunk_size}, parts={parts}, parallel={use_dask}, workers={expected_workers}")

    # dask 并行或串行执行
    wrote_any = False
    if use_dask:
        if keep_legacy_pickle:
            print("[export_meta_tracks] keep_legacy_pickle=True 在 Dask 模式下会占用大量内存，已强制关闭。")
            legacy_nested = None

    if use_dask:
        tasks = []
        part_idx = 0
        for start in range(0, total_len, chunk_size):
            end = min(total_len, start + chunk_size)
            task = delayed(_process_chunk)(file_path, part_idx, start, end)
            tasks.append(task)
            part_idx += 1
        scheduler_kwargs = {}
        if dask_num_workers is not None:
            scheduler_kwargs['num_workers'] = dask_num_workers
        with ProgressBar():
            results = compute(*tasks, scheduler='processes', **scheduler_kwargs)
        total_rows = int(np.sum(results))
        wrote_any = total_rows > 0
        print(f"[export_meta_tracks] kind={kind} dask-complete parts={len(tasks)} rows={total_rows}")
    else:
        parts = (total_len + chunk_size - 1) // chunk_size
        for part_idx, start in enumerate(tqdm(range(0, total_len, chunk_size), total=parts, desc=f"export {kind} serial", unit="part")):
            end = min(total_len, start + chunk_size)
            nrows = _process_chunk(file_path, part_idx, start, end)
            wrote_any = wrote_any or (nrows > 0)
            # 串行时不需要额外计数，最终通过分片合并生成 index

    # 轨迹汇总
    written = {
        'daily_dir': str(daily_stage),
        'contours_dir': str(contours_stage) if write_contours else None,
    }

    if build_track_summary and wrote_any:
        # 读回 daily parquet dataset 聚合（数据量相比 contours 小得多）
        daily_dataset = pq.ParquetDataset(daily_stage)
        daily_tbl = daily_dataset.read()
        df = daily_tbl.to_pandas()
        grp = df.groupby('track_id', sort=False)
        summary = grp.agg(
            n_points=('time', 'count'),
            time_start=('time', 'min'),
            time_end=('time', 'max'),
            center_lon_mean=('center_lon', 'mean'),
            center_lat_mean=('center_lat', 'mean'),
            radius_mean=('radius', 'mean'),
            radius_max=('radius', 'max'),
        ).reset_index()
        pq.write_table(pa.Table.from_pandas(summary), tracks_path)
        written['tracks_path'] = str(tracks_path)

    # 合并 NPZ 分片为单一 Zarr 存储（支持串行和 Dask）
    def _merge_zarr_parts(parts_dir: Path, out_path: Path, meta_info: dict):
        part_files = sorted(parts_dir.glob('part_*.npz'))
        if not part_files:
            return None
        # 第一遍：统计总长度
        eff_rows_total = 0
        spd_rows_total = 0
        eff_vert_total = 0
        spd_vert_total = 0
        for pf in part_files:
            with np.load(pf) as npz:
                ec = npz['eff_counts']
                sc = npz['spd_counts']
                eff_rows_total += int(ec.size)
                spd_rows_total += int(sc.size)
                eff_vert_total += int(npz['eff_lon'].size)
                spd_vert_total += int(npz['spd_lon'].size)
        # 创建 Zarr 存储
        if out_path.exists():
            shutil.rmtree(out_path)
        root = zarr.open_group(out_path, mode='w')
        z_eff = root.create_group('effective')
        z_spd = root.create_group('speed')
        z_meta = root.create_group('meta')
        eff_lon_arr = z_eff.create_array('lon', shape=(eff_vert_total,), chunks=(262144,), dtype='f8')
        eff_lat_arr = z_eff.create_array('lat', shape=(eff_vert_total,), chunks=(262144,), dtype='f8')
        spd_lon_arr = z_spd.create_array('lon', shape=(spd_vert_total,), chunks=(262144,), dtype='f8')
        spd_lat_arr = z_spd.create_array('lat', shape=(spd_vert_total,), chunks=(262144,), dtype='f8')
        eff_idx_arr = z_eff.create_array('index', shape=(eff_rows_total + 1,), chunks=(262144,), dtype='i8')
        spd_idx_arr = z_spd.create_array('index', shape=(spd_rows_total + 1,), chunks=(262144,), dtype='i8')
        eff_idx_arr[0] = 0
        spd_idx_arr[0] = 0
        z_meta.attrs.update(meta_info)
        # 第二遍：按顺序写入
        e_vert_pos = 0
        s_vert_pos = 0
        e_row_pos = 0
        s_row_pos = 0
        for pf in part_files:
            with np.load(pf) as npz:
                ev_lon = npz['eff_lon']; ev_lat = npz['eff_lat']; ec = npz['eff_counts']
                sv_lon = npz['spd_lon']; sv_lat = npz['spd_lat']; sc = npz['spd_counts']
                if ev_lon.size:
                    eff_lon_arr[e_vert_pos:e_vert_pos+ev_lon.size] = ev_lon
                    eff_lat_arr[e_vert_pos:e_vert_pos+ev_lat.size] = ev_lat
                    e_vert_pos += int(ev_lon.size)
                if sv_lon.size:
                    spd_lon_arr[s_vert_pos:s_vert_pos+sv_lon.size] = sv_lon
                    spd_lat_arr[s_vert_pos:s_vert_pos+sv_lat.size] = sv_lat
                    s_vert_pos += int(sv_lon.size)
                if ec.size:
                    csum = np.cumsum(ec, dtype=np.int64)
                    eff_idx_arr[e_row_pos+1:e_row_pos+1+csum.size] = eff_idx_arr[e_row_pos] + csum
                    e_row_pos += int(ec.size)
                if sc.size:
                    csum = np.cumsum(sc, dtype=np.int64)
                    spd_idx_arr[s_row_pos+1:s_row_pos+1+csum.size] = spd_idx_arr[s_row_pos] + csum
                    s_row_pos += int(sc.size)
        return str(out_path)

    contours_zarr_written = None
    if True:
        meta_info = {
            'kind': kind,
            'region': region_slug,
            'lon_min': float(lonmin), 'lon_max': float(lonmax),
            'lat_min': float(latmin), 'lat_max': float(latmax),
        }
        try:
            contours_zarr_written = _merge_zarr_parts(contours_parts_dir, contours_zarr_path, meta_info)
        finally:
            if contours_parts_dir.exists():
                try:
                    shutil.rmtree(contours_parts_dir)
                except Exception:
                    pass

    # 元信息写出
    meta = {
        'kind': kind,
        'region': region_slug,
        'lon_min': float(lonmin), 'lon_max': float(lonmax),
        'lat_min': float(latmin), 'lat_max': float(latmax),
        'generated_at': str(np.datetime64('now')),
        'columns_daily': [f.name for f in daily_schema],
        # 兼容旧键：Parquet 顶点拆表
        'has_contours': write_contours,
        'contour_schema': [f.name for f in contour_schema] if contour_schema else None,
        # 新增：Zarr 轮廓信息
        'has_contours_zarr': bool(contours_zarr_written),
        'contours_format': 'zarr' if contours_zarr_written else None,
        'contours_zarr': contours_zarr_written,
        'track_summary': bool(build_track_summary),
        'legacy_nested_pickle': bool(keep_legacy_pickle),
    }
    with open(meta_path, 'w') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    written['metadata'] = str(meta_path)

    if keep_legacy_pickle and legacy_nested is not None:
        # 把 dict 转为 list[list]，保持“首次出现”的轨迹顺序，以及每轨迹内的原始追加顺序
        nested_list = list(legacy_nested.values())
        with open(legacy_pickle_path, 'wb') as f:
            pickle.dump(nested_list, f)
        written['legacy_pickle'] = str(legacy_pickle_path)
        print(f"[export_meta_tracks] Legacy nested pickle written: {legacy_pickle_path}")

    # 压实：将目录 Dataset 合并为单文件，并删除目录
    if wrote_any and compact_after:
        try:
            daily_file = compact_parquet_dataset(kind, 'daily', region_slug, root, delete_source=True)
            written['daily_file'] = daily_file
            written['daily_dir'] = None
        except Exception as e:
            print(f"[export_meta_tracks] compact daily failed: {e}")
        if write_contours:
            try:
                contours_file = compact_parquet_dataset(kind, 'contours', region_slug, root, delete_source=True)
                written['contours_file'] = contours_file
                written['contours_dir'] = None
            except Exception as e:
                print(f"[export_meta_tracks] compact contours failed: {e}")

    # 简单 Zarr 校验
    if contours_zarr_written:
        try:
            zroot = zarr.open_group(contours_zarr_path, mode='r')
            eff_idx_len = int(zroot['effective']['index'].shape[0])
            spd_idx_len = int(zroot['speed']['index'].shape[0])
            if eff_idx_len != spd_idx_len:
                print(f"[export_meta_tracks] Warning: effective/speed index length mismatch: {eff_idx_len} vs {spd_idx_len}")
        except Exception as e:
            print(f"[export_meta_tracks] Zarr open failed for validation: {e}")

    if not wrote_any:
        print(f"[export_meta_tracks] No records inside region bounds for kind={kind}.")
    out_hint = written.get('daily_file', str(daily_stage))
    if contours_zarr_written:
        written['contours_zarr'] = contours_zarr_written
    print(f"[export_meta_tracks] Completed kind={kind}. Daily: {out_hint}; Contours(zarr): {contours_zarr_written}")
    return written

def convert_mat_to_parquet(year: int, input_dir: str | Path = None, output_dir: str | Path = None):
    """将某年份的逐月 Argo .mat 原始文件合并为单一 Parquet。

    参数:
        year (int): 年份 (如 2008)。
        input_dir (str | Path | None): 每月 .mat 文件所在目录；None → 配置 paths.yml 的 argo_mat_input。
        output_dir (str | Path | None): 输出目录；None → 配置 paths.yml 的 argo_parquet。

    说明:
        该函数当前仅示例性地提取 'do' 数据集（如果存在），并按固定列顺序写出。
        若未来需要更多变量，可在此扩展。
    """
    year_str = str(year)
    if input_dir is None:
        input_dir = argo_mat_input_path
    else:
        input_dir = Path(input_dir)
    if output_dir is None:
        output_dir = argo_path
    else:
        output_dir = Path(output_dir)

    columns = [
        "Year", "Month", "Day", "Longitude", "Latitude", "Depth_m", "DO_mol_kg", "Salinity_psu",
        "Temperature_degC", "Oxygen_flag", "Oxygen_flag2", "Profile_number", "Datasets_number",
        "Platform_number", "Cycle_number", "Float_serial_no"
    ]

    paths = [input_dir / f'Argo{year_str}_{m}.mat' for m in range(1, 13)]
    print(f"[convert_mat_to_parquet] Year {year_str}: expecting {len(paths)} monthly files under {input_dir}")

    all_monthly = []
    for p in paths:
        if not p.exists():
            print(f"  - Missing file: {p}")
            continue
        try:
            with h5py.File(p, 'r') as f:
                if 'do' not in f:
                    print(f"  - 'do' dataset absent in {p}, skipping")
                    continue
                do_data = f['do'][:].T  # 原数据转置为 (n_rows, n_cols)
                df = pd.DataFrame(do_data, columns=columns[:do_data.shape[1]])  # 防御性截断
                all_monthly.append(df)
                print(f"  + Loaded {p}")
        except Exception as e:
            print(f"  ! Error reading {p}: {e}")

    if not all_monthly:
        print(f"[convert_mat_to_parquet] No monthly data loaded for {year_str}; abort.")
        return

    final_df = pd.concat(all_monthly, ignore_index=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f'Argo{year_str}.parquet'
    try:
        final_df.to_parquet(out_path, index=False)
        print(f"[convert_mat_to_parquet] Saved merged parquet: {out_path}")
    except Exception as e:
        print(f"[convert_mat_to_parquet] Failed to save {out_path}: {e}")

def compact_parquet_dataset(kind: str,
                            which: str,
                            region_key: str | None = None,
                            output_root: str | Path | None = None,
                            delete_source: bool = True) -> str:
    """将目录形式的 Parquet Dataset (many part_*.parquet) 压实为单一文件，并可选删除源目录。

    参数:
        kind: 数据类型前缀，如 'acs'、'acl'、'cs'、'cl'。
        which: 'daily' 或 'contours'。
        region_key: 区域 slug；None 时基于当前区域生成。
    output_root: 根输出目录；None 则读取配置或使用默认 './META_tracks'。
        delete_source: 压实完成后是否删除源目录。

    返回:
        目标单一 Parquet 文件路径（字符串）。
    """
    region_slug = region_key or _current_region_key()
    root = _ensure_meta_tracks_root(output_root)
    region_dir = root / region_slug
    if which not in { 'daily', 'contours' }:
        raise ValueError("which 必须是 'daily' 或 'contours'")

    # 优先使用 *_tmp 作为分片目录；若不存在，则兼容旧的 *_<which>.parquet 目录
    candidate_dirs = [
        region_dir / f"{kind}_{which}_tmp",
        region_dir / f"{kind}_{which}.parquet",
    ]
    dir_path = None
    for cand in candidate_dirs:
        if cand.exists() and cand.is_dir():
            dir_path = cand
            break
    if dir_path is None:
        raise FileNotFoundError(f"未找到分片目录: {candidate_dirs[0]} 或 {candidate_dirs[1]}")

    # 目标文件路径（最终与目录同名，但先写 tmp，再替换）
    tmp_path = region_dir / f"{kind}_{which}.parquet.tmp"
    final_file = region_dir / f"{kind}_{which}.parquet"

    # 收集分片文件
    part_files = sorted([p for p in dir_path.glob('*.parquet') if p.is_file()])
    if not part_files:
        raise RuntimeError(f"未发现任何分片文件于 {dir_path}")

    # 流式合并
    first_tbl = pq.read_table(part_files[0])
    schema = first_tbl.schema
    # 若存在旧 tmp，先移除
    if tmp_path.exists():
        tmp_path.unlink()
    rows_total = 0
    with pq.ParquetWriter(tmp_path, schema) as writer:
        for pf in part_files:
            tbl = pq.read_table(pf)
            writer.write_table(tbl)
            rows_total += tbl.num_rows

    # 删除源目录，重命名 tmp 为最终文件（若已存在则先删除）
    if delete_source:
        shutil.rmtree(dir_path)
    if final_file.exists():
        final_file.unlink()
    tmp_path.rename(final_file)
    print(f"[compact] {which} -> {final_file} (rows≈{rows_total})")
    return str(final_file)

def parse_argo_txt_file(file_path: Path) -> pd.DataFrame | None:
    """
    解析单个扁平的、由制表符分隔的Argo表格文件。

    功能:
        此函数假设输入文件是一个标准的表格文件，其中第一行是列标题，
        后续行是数据，所有字段均由制表符 ('\t') 分隔。
        它会自动将指定的占位符（如 -999）转换成标准的 NaN。

    参数:
        file_path (Path): 单个Argo txt文件的路径。

    返回:
        pd.DataFrame | None: 如果文件成功解析且内容不为空，返回DataFrame；否则返回None。
    """
    try:
        # 定义需要被视作NaN的占位符列表
        missing_values = [-999, -999.0, -999.00, -999.000, -999.0000, -999.000000]
        
        # 在读取时使用 na_values 参数直接完成替换
        df = pd.read_csv(
            file_path,
            sep='\t',
            na_values=missing_values
        )
        
        # 如果文件只有表头或完全为空，df.empty会是True
        if df.empty:
            return None
        return df
    except Exception as e:
        # print(f"  [Warning] 读取文件 {file_path.name} 失败: {e}")
        return None

def worker_process_file(task_args: tuple[Path, Path]):
    """
    (Worker函数) 这是Dask的工作单元，负责处理单个文件。

    功能:
        接收一个包含输入和输出路径的元组，调用解析函数，
        并将结果保存为中间Parquet文件。

    参数:
        task_args (tuple[Path, Path]): 一个元组，包含:
            - 原始txt文件路径 (input_path)
            - 中间Parquet文件输出路径 (output_path)
    """
    input_path, output_path = task_args
    try:
        df = parse_argo_txt_file(input_path)
        if df is not None:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(output_path, index=False)
    except Exception as e:
        print(f"\n[Error] 在工作进程中处理文件 {input_path.name} 时出错: {e}")


def process_argo_txt_to_yearly_parquet_dask(
    origin_dir: Path | str | None = None,
    temp_dir: Path | str | None = None,
    final_dir: Path | str | None = None,
    cleanup_temp_dir: bool = True
):
    """
    (主流程函数) 使用Dask统一并行处理Argo txt文件，高效处理大规模数据。

    功能:
        这是一个完整的、适用于超大数据集的ETL流程，完全由Dask驱动：
        1. (并行) [Dask] 将所有原始 .txt 文件转换为临时的 .parquet 文件。
        2. (并行) [Dask] 对所有临时文件进行并行地合并、清洗、排序，并输出为按年份分区的临时目录。
        3. (串行) 将分区目录合并为最终的单个年度文件。
        此版本利用Dask的HPC环境自适应能力，避免内存瓶颈，全程并行执行。
        最后，根据参数选择是否清理所有临时文件。

    参数:
        origin_dir (Path | str | None): 存放原始 Argo .txt 文件的目录。
            - None 时自动读取 paths.yml 中的 paths.argo_txt_input；若缺省，则回退到 './Argo_origin'。
        temp_dir (Path | str | None): 存放中间产物（初始 Parquet、映射表、分区数据）的临时目录。
            - None 时优先读取 paths.yml 中的 paths.tmp_parquet_path；若缺省，则回退到 final_dir / '_tmp_txt2parquet_dask'。
        final_dir (Path | str | None): 保存最终年度 Parquet 文件（ArgoYYYY.parquet）的目录。
            - None 时使用配置 paths.yml 中的 paths.argo_parquet（即全局 argo_path）。
        cleanup_temp_dir (bool, optional): 是否在任务结束后删除临时目录。默认为 True。
    """
    # --- 0. 解析/回退目录参数（与 META 一致的配置驱动） ---
    start_total_time = tm.time()
    origin_dir = Path(origin_dir) if origin_dir is not None else Path(_PATHS_CFG.get('paths', {}).get('argo_txt_input', './Argo_origin'))
    final_dir = Path(final_dir) if final_dir is not None else Path(argo_path)
    if temp_dir is None:
        # 使用已在模块顶层解析好的全局 tmp_parquet_path（paths.yml: argo_intermediate），
        # 若用户未配置则回退到 final_dir/_tmp_txt2parquet_dask
        temp_dir = Path(tmp_parquet_path) if tmp_parquet_path else (final_dir / '_tmp_txt2parquet_dask')
    else:
        temp_dir = Path(temp_dir)
    print("[*] Using directories for Argo TXT → Parquet:")
    print(f"    - origin_dir = {origin_dir}")
    print(f"    - temp_dir   = {temp_dir}")
    print(f"    - final_dir  = {final_dir}")
    
    # --- 准备工作：初始化Dask客户端并创建目录 ---
    client = Client()
    print(f"[*] Dask客户端已启动，仪表盘链接: {client.dashboard_link}")
    
    scheduler_info = client.scheduler_info()
    workers_count = len(scheduler_info['workers'])
    total_threads = sum(worker['nthreads'] for worker in scheduler_info['workers'].values())
    print(f"[*] Dask自动检测到 {workers_count} 个 worker，共 {total_threads} 个工作核心。")
    
    print("[*] 准备工作环境...")
    if temp_dir.exists():
        print(f"  - 发现旧的临时目录，正在清理: {temp_dir}")
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    final_dir.mkdir(parents=True, exist_ok=True)
    
    # --- 阶段一: 使用Dask并行解析 -> 中间文件 ---
    print("\n--- 阶段一: 开始并行解析原始 .txt 文件 ---")
    txt_files = list(origin_dir.glob('*.txt'))
    if not txt_files:
        print("[Error] 输入目录中未找到任何 .txt 文件。程序终止。")
        client.close()
        return
        
    tasks = [(p, temp_dir / f"{p.stem}.parquet") for p in txt_files]
    print(f"[*] 找到 {len(txt_files)} 个文件，将任务分发给Dask处理...")
    
    start_map_time = tm.time()
    futures = client.map(worker_process_file, tasks)
    client.gather(futures)
    print(f"--- 阶段一完成 --- (耗时: {tm.time() - start_map_time:.2f} 秒)")
    
    # --- 阶段二: Dask并行ETL并输出分区文件 ---
    print("\n--- 阶段二: Dask开始并行处理中间文件 ---")
    start_reduce_time = tm.time()

    print(f"[*] Dask正在逻辑读取 {temp_dir} 中的Parquet文件...")
    ddf = dd.read_parquet(temp_dir / '*.parquet')
    record_count = len(ddf)
    print(f"[*] 数据逻辑加载完成，总计 {record_count:,} 条记录。")

    print("[*] 正在重命名列...")
    rename_map = {'WMO': 'Platform_number', 'Cycle': 'Profile_number', 'Lon': 'Longitude', 'Lat': 'Latitude'}
    ddf = ddf.rename(columns=rename_map)

    print("[*] 正在并行生成全局唯一的 Profile_number...")
    unique_profiles = ddf[['Platform_number', 'Profile_number']].drop_duplicates().compute()
    unique_profiles_sorted = unique_profiles.sort_values(by=['Platform_number', 'Profile_number'])
    unique_profiles_sorted['Global_Profile_Number'] = range(1, len(unique_profiles_sorted) + 1)
    mapping_file_path = temp_dir / 'profile_mapping.parquet'
    unique_profiles_sorted.to_parquet(mapping_file_path, index=False)
    mapping_ddf = dd.read_parquet(mapping_file_path)
    ddf = ddf.merge(mapping_ddf, on=['Platform_number', 'Profile_number'], how='left')
    ddf = ddf.drop(columns=['Profile_number']).rename(columns={'Global_Profile_Number': 'Profile_number'})

    print("[*] 正在优化并转换数据类型...")
    float64_cols = list(ddf.select_dtypes(include='float64').columns)
    if float64_cols:
        for col in float64_cols:
            ddf[col] = ddf[col].astype('float32')
    int_conversion_map = {'Year': 'Int16', 'Month': 'Int8', 'Day': 'Int8', 'Platform_number': 'Int64', 'Profile_number': 'Int64'}
    for col, dtype in int_conversion_map.items():
        if col in ddf.columns:
            ddf[col] = ddf[col].astype(dtype)
    
    print("[*] Dask开始执行所有计算并按年份分区写入临时目录...")
    ddf = ddf.dropna(subset=['Year'])
    partitioned_output_dir = temp_dir / 'dask_partitioned_output'
    ddf.to_parquet(partitioned_output_dir, write_index=False, partition_on=['Year'])
    
    end_reduce_time = tm.time()
    print(f"--- 阶段二完成 --- (耗时: {end_reduce_time - start_reduce_time:.2f} 秒)")
    
    # --- 阶段三: 将分区文件合并为最终产物 ---
    print("\n--- 阶段三: 开始合并分区文件并重排定序 ---")
    start_consolidation_time = tm.time()
    
    year_dirs = sorted([d for d in partitioned_output_dir.iterdir() if d.is_dir() and d.name.startswith('Year=')])
    for year_dir in year_dirs:
        year_str = year_dir.name.split('=')[1]
        year_int = int(year_str)
        final_single_file_path = final_dir / f"Argo{year_int}.parquet"
        print(f"  -> 合并年份 {year_int} -> {final_single_file_path}")
        
        year_df = pd.read_parquet(year_dir)
        year_df['Year'] = year_int
        year_df.sort_values(by=['Year', 'Month', 'Day', 'Platform_number'], inplace=True)
        
        start_cols = ['Year', 'Month', 'Day']
        end_cols = ['Profile_number', 'Platform_number']
        middle_cols = [col for col in year_df.columns if col not in start_cols + end_cols]
        final_column_order = start_cols + middle_cols + end_cols
        year_df = year_df[final_column_order]
        
        year_df.to_parquet(final_single_file_path, index=False)
        
    end_consolidation_time = tm.time()
    print(f"--- 阶段三完成 --- (耗时: {end_consolidation_time - start_consolidation_time:.2f} 秒)")

    # --- 清理工作 ---
    if cleanup_temp_dir:
        print("\n--- 清理阶段: 删除临时文件 ---")
        try:
            shutil.rmtree(temp_dir)
            print(f"[*] 临时目录 {temp_dir} 已成功删除。")
        except Exception as e:
            print(f"[Error] 删除临时目录失败: {e}")
    else:
        print("\n--- 清理阶段: 已跳过 ---")
        print(f"[*] 临时文件已保留在: {temp_dir}")

    client.close()
    end_total_time = tm.time()
    print("\n==================================================")
    print(f"[Success] 所有任务完成！总耗时: {(end_total_time - start_total_time)/60:.2f} 分钟。")
    print("==================================================")

def load_argo_data(year: int, data_dir: str | Path = None,
                   variable_selection: dict | None = None,
                   verbose: bool = False) -> pd.DataFrame:
    """
    加载指定年份的 Argo Parquet 数据文件，并进行列名规范化和变量选择。

    功能:
        1. 自动将旧版列名 (如 'Depth_m') 转换为新版标准名 ('Depth')。
        2. 通过 `variable_selection` 参数，允许用户灵活选择最终输出的标准变量
           (Temperature, DO, Salinity) 分别来源于文件中的哪一列数据。

    参数:
        year (int): 需要加载的数据年份 (例如 2014)。
        data_dir (str | Path | None): Argo Parquet 所在目录；None → paths.yml: argo_parquet。
        variable_selection (dict | None): 覆盖默认变量来源映射，例如 {'Salinity':'PSAL_WOA'}。
            默认: {'Temperature': 'Temp_Adjusted', 'DO': 'DOXY_Adjusted', 'Salinity': 'PSAL_Adjusted'}
        verbose (bool): 是否输出详细日志，默认 False（安静模式）。

    返回:
        pd.DataFrame: 一个包含处理后 Argo 数据的 pandas DataFrame，其列名和数据源
                      均已根据参数进行了标准化。
    """
    if data_dir is None:
        data_dir = argo_path
    else:
        data_dir = Path(data_dir)

    # --- 1. 定义并合并变量选择 ---
    # 定义默认选择
    default_selection = {
        'Temperature': 'Temp_Adjusted',
        'Temperature_Flag': 'Temp_Adjusted_Flag',
        'DO': 'DOXY_Adjusted',
        'DO_Flag': 'DOXY_Adjusted_Flag',
        'Salinity': 'PSAL_Adjusted',
        'Salinity_Flag': 'PSAL_Adjusted_Flag'
    }
    # 如果用户提供了自定义选择，则用它来更新（覆盖）默认值
    if variable_selection:
        default_selection.update(variable_selection)
    
    # 最终生效的选择
    final_selection = default_selection
    
    # --- 2. 构建路径并加载文件 ---
    file_path = Path(data_dir) / f'Argo{year}.parquet'
    if verbose:
        print(f"Attempting to load Argo data from: {file_path}")
    if not file_path.exists():
        raise FileNotFoundError(f"Error: The file '{file_path}' was not found.")
    try:
        argo_df = pd.read_parquet(file_path)
    except Exception as e:
        raise RuntimeError(f"Failed to read Argo parquet file: {file_path}") from e

    # --- 3. 列名规范化：将所有已知的旧列名统一为新版中对应的名称 ---
    normalization_map = {
        'Depth_m': 'Depth',
        'Temperature_degC': 'Temp_Adjusted',
        'DO_mol_kg': 'DOXY_Adjusted',
        'Salinity_psu': 'PSAL_Adjusted'
    }
    argo_df.rename(columns=normalization_map, inplace=True)

    # --- 4. 变量选择与最终DataFrame构建 ---
    base_columns = [
        'Year', 'Month', 'Day', 'Longitude', 'Latitude', 'Depth',
        'Profile_number', 'Platform_number'
    ]
    existing_base_columns = [col for col in base_columns if col in argo_df.columns]
    final_df = argo_df[existing_base_columns].copy()

    for standard_name, source_col in final_selection.items():
        if source_col in argo_df.columns:
            final_df[standard_name] = argo_df[source_col]
        else:
            if verbose:
                print(f"Warning: Column '{source_col}' not found in {file_path}. Creating empty column '{standard_name}'.")
            final_df[standard_name] = pd.NA

    if verbose:
        print("Argo data loaded and processed successfully.")
    return final_df

def find_track(
    DS_or_kind: list | str,
    num: int | list | tuple | set | np.ndarray,
    *,
    region: str | None = None,
    output_root: str | Path | None = None,
    include_contours: bool = True,
    return_list: bool = False
) -> pd.DataFrame | list | dict:
    """查找一个或多个涡旋编号的整条轨迹（兼容老/新两种数据源），并统一输出约定。

    支持的数据源:
    - DS_or_kind 为旧的嵌套列表（例如从 pickle 读取的 legacy `acl/acs/...` 列表）。
    - DS_or_kind 为字符串 kind（'acs'|'acl'|'cs'|'cl'）：从 META_tracks/<region>
      读取 <kind>_daily.parquet 与（可选）<kind>_contours.zarr。

    返回契约（默认 return_list=False):
    - 单 ID: 返回 DataFrame（列依数据源而异）
            • 新结构（DS_or_kind 为字符串 kind）：首列为 'track_id'，列为
                ['track_id','time','center_lon','center_lat','max_lon','max_lat',
                    'contour_lon','contour_lat','radius','speed_contour_lon','speed_contour_lat','date']
                说明：'date' 由 'time' 规范化得到；不附加末尾重复的 'track_id' 列。
            • 旧结构（DS_or_kind 为嵌套列表 legacy）：首列为 'index_org'，并在末尾追加 'track_id'：
                ['index_org','time','center_lon','center_lat','max_lon','max_lat',
                    'contour_lon','contour_lat','radius','speed_contour_lon','speed_contour_lat','date','track_id']
    - 多 ID: 返回合并后的 DataFrame（含 'track_id' 列），并按 ['track_id','date'] 排序。

    当 return_list=True 时:
    - 单 ID: 返回旧版 list[list] 结构，每个日项的首元素为真实 track_id：
        [track_id, YYYYMMDD, center_lon, center_lat, max_lon, max_lat,
         eff_contour_lon[], eff_contour_lat[], radius, speed_contour_lon[], speed_contour_lat[]]
    - 多 ID: 返回 { track_id: list[list] } 的字典。

    参数:
    - DS_or_kind: 旧结构列表或新结构 kind 字符串。
    - num: 单个 track_id 或可迭代的多个 track_id。
    - region: 区域 slug，None 使用当前默认区域。
    - output_root: META_tracks 根目录，None 读取配置或使用默认 './META_tracks'。
    - include_contours: True 时加载等值线（新结构路径）；False 时 'contour_*' 为空数组。

    异常:
    - TypeError: DS_or_kind 类型非法。
    - ValueError: 指定 track_id 未找到。
    """
    # 规范化 num 输入，判断是否批量
    def _is_scalar_id(x) -> bool:
        return isinstance(x, (int, np.integer))

    if _is_scalar_id(num):
        id_set = {int(num)}
        multi = False
    else:
        try:
            id_set = set(int(x) for x in list(num))
        except Exception:
            id_set = {int(num)}
            multi = False
        else:
            multi = len(id_set) > 1
    # 1) 旧数据结构兼容：DS_or_kind 为 list 时沿用旧逻辑
    if isinstance(DS_or_kind, list):
        DS = DS_or_kind
        id_to_track = {}
        for track in DS:
            if not track:
                continue
            try:
                tid = int(track[0][0])
            except Exception:
                continue
            if tid in id_set:
                id_to_track[tid] = track

        missing = sorted(list(id_set - set(id_to_track.keys())))
        if missing:
            raise ValueError(f"Track not found for id(s): {missing}")

        if not multi:
            tid = next(iter(id_set))
            track = id_to_track[tid]
            if return_list:
                return track
            df = pd.DataFrame(track, columns=[
                'index_org','time','center_lon','center_lat','max_lon','max_lat',
                'contour_lon','contour_lat','radius','speed_contour_lon','speed_contour_lat'
            ])
            try:
                df['date'] = convert_date(df['time'])
            except Exception:
                df['date'] = pd.NaT
            
            ACS, ACL, CS, CL = load_meta_data()
            # 匹配传入列表变量名对应的 META 数据集（ACS/ACL/CS/CL）
            matched_ds_name = 'UNKNOWN'
            try:
                caller_locals = inspect.currentframe().f_back.f_locals
                for var_name, var_val in caller_locals.items():
                    if var_val is DS_or_kind:
                        matched_ds_name = var_name.upper()
                        break
            except Exception:
                matched_ds_name = 'UNKNOWN'

            meta_map = {'ACS': ACS, 'ACL': ACL, 'CS': CS, 'CL': CL}
            matched_meta = meta_map.get(matched_ds_name)

            df['track_id'] = matched_meta['track'][df['index_org']] if matched_meta is not None else int(tid)
            return df
        else:
            if return_list:
                return {tid: id_to_track[tid] for tid in sorted(id_set)}
            frames = []
            for tid, track in id_to_track.items():
                df = pd.DataFrame(track, columns=[
                    'index_org','time','center_lon','center_lat','max_lon','max_lat',
                    'contour_lon','contour_lat','radius','speed_contour_lon','speed_contour_lat'
                ])
                try:
                    df['date'] = convert_date(df['time'])
                except Exception:
                    df['date'] = pd.NaT

                ACS, ACL, CS, CL = load_meta_data()
                # 匹配传入列表变量名对应的 META 数据集（ACS/ACL/CS/CL）
                matched_ds_name = 'UNKNOWN'
                try:
                    caller_locals = inspect.currentframe().f_back.f_locals
                    for var_name, var_val in caller_locals.items():
                        if var_val is DS_or_kind:
                            matched_ds_name = var_name.upper()
                            break
                except Exception:
                    matched_ds_name = 'UNKNOWN'

                meta_map = {'ACS': ACS, 'ACL': ACL, 'CS': CS, 'CL': CL}
                matched_meta = meta_map.get(matched_ds_name)

                df['track_id'] = matched_meta['track'][df['index_org']] if matched_meta is not None else int(tid)

                frames.append(df)
            out = pd.concat(frames, ignore_index=True)
            out.sort_values(['track_id', 'date'], inplace=True)
            out.reset_index(drop=True, inplace=True)
            return out

    # 2) 新数据结构：按 kind 从 Parquet + Zarr 中提取
    if not isinstance(DS_or_kind, str):
        raise TypeError("find_track 新用法：请传入 kind 字符串（'acs'|'acl'|'cs'|'cl'）作为第一个参数，或传旧的 DS 列表。")
    kind = DS_or_kind.lower()
    if kind not in {'acs', 'acl', 'cs', 'cl'}:
        raise ValueError(f"未知 kind='{DS_or_kind}', 期望 'acs'|'acl'|'cs'|'cl'.")

    region_slug = region or _current_region_key()
    root = _ensure_meta_tracks_root(output_root)
    region_dir = Path(root) / region_slug
    if not region_dir.exists():
        raise FileNotFoundError(f"区域目录不存在：{region_dir}")

    # 2.1 定位 daily 源：优先单文件，其次目录（_tmp 或 .parquet 目录）
    daily_file = region_dir / f"{kind}_daily.parquet"
    daily_tmp_dir = region_dir / f"{kind}_daily_tmp"
    daily_dir = region_dir / f"{kind}_daily.parquet"  # 目录形式
    daily_source = None
    source_type = None  # 'file' | 'dir'
    if daily_file.exists() and daily_file.is_file():
        daily_source = daily_file
        source_type = 'file'
    elif daily_tmp_dir.exists() and daily_tmp_dir.is_dir():
        daily_source = daily_tmp_dir
        source_type = 'dir'
    elif daily_dir.exists() and daily_dir.is_dir():
        daily_source = daily_dir
        source_type = 'dir'
    else:
        raise FileNotFoundError(f"未找到 daily 数据：{daily_file} 或 {daily_tmp_dir} / {daily_dir}")

    # 2.2 读取元信息以定位 Zarr（若需要）
    zarr_path = region_dir / f"{kind}_contours.zarr"
    meta_path = region_dir / f"{kind}_metadata.json"
    if include_contours and meta_path.exists():
        try:
            with open(meta_path, 'r') as f:
                meta = json.load(f)
            if meta.get('contours_zarr'):
                zarr_path = Path(meta['contours_zarr'])
        except Exception:
            pass

    # 2.3 流式扫描 Parquet（按行组/分片顺序）收集该 track 的全部行及其“全局行号”
    def _scan_daily_rows_file(parquet_file: Path, target_ids: set[int]):
        pf = pq.ParquetFile(parquet_file)
        # 找到 track_id 列的索引（用于读取统计信息）
        try:
            track_idx = [i for i, n in enumerate(pf.schema.names) if n == 'track_id'][0]
        except Exception:
            track_idx = None
        cum = 0
        min_tid = min(target_ids)
        max_tid = max(target_ids)
        for rg in range(pf.num_row_groups):
            # 先用统计信息快速判断是否可能包含 num
            if track_idx is not None and pf.metadata is not None:
                try:
                    col_meta = pf.metadata.row_group(rg).column(track_idx)
                    stats = col_meta.statistics
                    if stats is not None and stats.has_min_max:
                        min_v, max_v = stats.min, stats.max
                        # 某些版本可能返回 bytes，需要安全转换
                        try:
                            min_v = int(min_v); max_v = int(max_v)
                        except Exception:
                            pass
                        if isinstance(min_v, (int, np.integer)) and isinstance(max_v, (int, np.integer)):
                            if (max_tid < min_v) or (min_tid > max_v):
                                # 不可能命中，跳过读取该行组
                                cum += pf.metadata.row_group(rg).num_rows
                                continue
                except Exception:
                    pass

            tbl = pf.read_row_group(rg, columns=[
                'track_id', 'time', 'center_lon', 'center_lat', 'max_lon', 'max_lat', 'radius', 'orig_index'
            ])
            t_id = tbl.column('track_id').to_numpy()
            mask = np.isin(t_id, list(target_ids))
            idxs = np.nonzero(mask)[0]
            if idxs.size:
                df = tbl.to_pandas()
                for j in idxs.tolist():
                    yield cum + j, df.iloc[j]
            cum += tbl.num_rows

    def _scan_daily_rows_dir(parquet_dir: Path, target_ids: set[int]):
        parts = sorted([p for p in parquet_dir.glob('*.parquet') if p.is_file()])
        cum = 0
        min_tid = min(target_ids)
        max_tid = max(target_ids)
        for p in parts:
            pf = pq.ParquetFile(p)
            # 找到 track_id 列的索引
            try:
                track_idx = [i for i, n in enumerate(pf.schema.names) if n == 'track_id'][0]
            except Exception:
                track_idx = None
            for rg in range(pf.num_row_groups):
                # 统计信息快速跳过
                if track_idx is not None and pf.metadata is not None:
                    try:
                        col_meta = pf.metadata.row_group(rg).column(track_idx)
                        stats = col_meta.statistics
                        if stats is not None and stats.has_min_max:
                            min_v, max_v = stats.min, stats.max
                            try:
                                min_v = int(min_v); max_v = int(max_v)
                            except Exception:
                                pass
                            if isinstance(min_v, (int, np.integer)) and isinstance(max_v, (int, np.integer)):
                                if (max_tid < min_v) or (min_tid > max_v):
                                    cum += pf.metadata.row_group(rg).num_rows
                                    continue
                    except Exception:
                        pass

                tbl = pf.read_row_group(rg, columns=[
                    'track_id', 'time', 'center_lon', 'center_lat', 'max_lon', 'max_lat', 'radius', 'orig_index'
                ])
                t_id = tbl.column('track_id').to_numpy()
                mask = np.isin(t_id, list(target_ids))
                idxs = np.nonzero(mask)[0]
                if idxs.size:
                    df = tbl.to_pandas()
                    for j in idxs.tolist():
                        yield cum + j, df.iloc[j]
                cum += tbl.num_rows

    scanner = _scan_daily_rows_file if source_type == 'file' else _scan_daily_rows_dir
    rows = list(scanner(daily_source, id_set))
    if not rows:
        raise ValueError(f"Track id(s) {sorted(list(id_set))} 未在 {daily_source} 中找到")

    # 2.4 可选：打开 Zarr 并准备索引数组
    eff_lon_arr = eff_lat_arr = spd_lon_arr = spd_lat_arr = None
    eff_idx = spd_idx = None
    if include_contours:
        if not zarr_path.exists():
            print(f"[find_track] 警告：未找到 Zarr 存储 {zarr_path}，将仅返回无轮廓的轨迹基础信息。")
            include_contours = False
        else:
            zroot = zarr.open_group(zarr_path, mode='r')
            eff_grp = zroot['effective']
            spd_grp = zroot['speed']
            # 注意：保留 Zarr 数组对象本身，按需切片，避免整数组载入内存
            eff_lon_arr = eff_grp['lon']; eff_lat_arr = eff_grp['lat']
            spd_lon_arr = spd_grp['lon']; spd_lat_arr = spd_grp['lat']
            eff_idx = eff_grp['index'];  spd_idx = spd_grp['index']

    # 2.5 组装为结果（按时间升序）
    # rows: list of (global_row_idx, pandas_series)
    def _ymd_i8(ts: pd.Timestamp) -> int:
        ts = pd.Timestamp(ts)
        return ts.year * 10000 + ts.month * 100 + ts.day

    # 排序（多数情况下已是有序，这里稳妥起见按 time 排）
    rows_sorted = sorted(rows, key=lambda kv: pd.Timestamp(kv[1]['time']))

    # 将结果分配到各个 track_id 下
    rows_per_id: dict[int, list] = {}
    for global_idx, r in rows_sorted:
        tid = int(r['track_id']) if 'track_id' in r else None
        if tid is None:
            continue
        # 基础字段
        t_ymd = _ymd_i8(r['time'])
        center_lon = float(r['center_lon'])
        center_lat = float(r['center_lat'])
        max_lon = float(r['max_lon']) if 'max_lon' in r else float('nan')
        max_lat = float(r['max_lat']) if 'max_lat' in r else float('nan')
        radius = float(r['radius']) if 'radius' in r else float('nan')

        # 轮廓
        eff_lon = np.array([], dtype='f8')
        eff_lat = np.array([], dtype='f8')
        spd_lon = np.array([], dtype='f8')
        spd_lat = np.array([], dtype='f8')
        if include_contours and eff_idx is not None:
            if 0 <= global_idx < (eff_idx.shape[0] - 1):
                s_e = eff_idx[global_idx:global_idx+2]
                s = int(s_e[0]); e = int(s_e[1])
                if e > s:
                    eff_lon = np.asarray(eff_lon_arr[s:e])
                    eff_lat = np.asarray(eff_lat_arr[s:e])
            if 0 <= global_idx < (spd_idx.shape[0] - 1):
                s_e = spd_idx[global_idx:global_idx+2]
                s = int(s_e[0]); e = int(s_e[1])
                if e > s:
                    spd_lon = np.asarray(spd_lon_arr[s:e])
                    spd_lat = np.asarray(spd_lat_arr[s:e])

        # 注意：为简化后续标注/检索，首字段使用 track_id（tid）而非原始行索引
        rows_per_id.setdefault(tid, []).append([
            int(tid), t_ymd, center_lon, center_lat, max_lon, max_lat,
            eff_lon, eff_lat, radius, spd_lon, spd_lat
        ])

    # 检查是否全部命中
    hit_ids = set(rows_per_id.keys())
    missing = sorted(list(id_set - hit_ids))
    if missing:
        raise ValueError(f"Track not found for id(s): {missing}")

    if not multi:
        tid = next(iter(id_set))
        rows_list = rows_per_id[tid]
        if return_list:
            return rows_list
        df = pd.DataFrame(rows_list, columns=[
            'track_id','time','center_lon','center_lat','max_lon','max_lat',
            'contour_lon','contour_lat','radius','speed_contour_lon','speed_contour_lat'
        ])
        try:
            df['date'] = convert_date(df['time'])
        except Exception:
            df['date'] = pd.NaT
        return df
    else:
        if return_list:
            return {tid: rows_per_id[tid] for tid in sorted(rows_per_id.keys())}
        frames = []
        for tid, rows_list in rows_per_id.items():
            df = pd.DataFrame(rows_list, columns=[
                'track_id','time','center_lon','center_lat','max_lon','max_lat',
                'contour_lon','contour_lat','radius','speed_contour_lon','speed_contour_lat'
            ])
            try:
                df['date'] = convert_date(df['time'])
            except Exception:
                df['date'] = pd.NaT
            frames.append(df)
        out = pd.concat(frames, ignore_index=True)
        out.sort_values(['track_id', 'date'], inplace=True)
        out.reset_index(drop=True, inplace=True)
        return out

def is_point_in_contour(row: pd.Series) -> bool:
    """
    检查一个数据行中的点坐标是否在其涡旋轮廓坐标内部。

    功能:
        这是一个辅助函数，通常与DataFrame的.apply()方法配合使用。
        它从一行数据中提取'Longitude', 'Latitude', 'contour_lon', 'contour_lat'，
        并判断点是否在多边形内部。

    参数:
        row (pd.Series): 
            一个Pandas DataFrame的数据行，必须包含经纬度和轮廓坐标。

    返回:
        bool: 如果点在轮廓内部，则返回True，否则返回False。
    """
    try:
        contour_lon = np.asarray(row['contour_lon'], dtype=float)
        contour_lat = np.asarray(row['contour_lat'], dtype=float)
        if contour_lon.size < 3 or contour_lat.size < 3:
            return False
        if contour_lon.shape != contour_lat.shape:
            return False

        center_lon = row.get('center_lon', np.nan)
        if pd.isna(center_lon):
            if contour_lon.size:
                center_lon = float(contour_lon[0])
            else:
                center_lon = float(row['Longitude'])

        contour_lon_norm = center_lon + _minimal_lon_diff_deg(contour_lon, center_lon)
        point_lon_norm = center_lon + _minimal_lon_diff_deg(row['Longitude'], center_lon)

        contour_coords = list(zip(contour_lon_norm, contour_lat))
        if len(contour_coords) < 3:
            return False

        path = MplPath(contour_coords)
        return path.contains_point((point_lon_norm, row['Latitude']))
    except Exception:
        return False

def filtered_float_data(
    DS: list | str,
    no: int,
    argo_data_dir: str | Path = None,
    circle_enlargement_factor: float | None = None,
    use_adaptive_circle: bool = True,
    adaptive_lat_threshold: float = 70.0,
    adaptive_distance_threshold_km: float = 300.0,
    force_great_circle_circle: bool = False,
    track: list | pd.DataFrame | None = None
) -> pd.DataFrame:
    """
    根据涡旋轨迹，动态加载并筛选匹配的Argo浮标剖面数据。

    功能:
        1. 分析涡旋轨迹覆盖的年份，并自动加载对应年份的Argo数据。
        2. 使用向量化的方式高效匹配在同一天出现在涡旋区域内的Argo浮标。
        3. 筛选标准为：浮标位置处于涡旋的有效轮廓内，或处于扩大一定倍数后的有效半径内。

    参数:
        DS (list | str):
            • 旧结构：legacy 轨迹列表（如从 pickle 读取的 acl/acs/... 数据）。
            • 新结构：字符串 kind（'acs'|'acl'|'cs'|'cl'），内部会调用 find_track 动态装载 META_tracks/<region> 数据。
        no (int): 需要筛选的涡旋的唯一编号。
        argo_data_dir (str | Path, optional): 存放Argo Parquet文件的目录。None 时使用配置文件 paths.yml 中的 argo_parquet。
        circle_enlargement_factor (float | None, optional): None 时回退配置值。
        use_adaptive_circle (bool): 若为 True，则半径匹配的距离用 `adaptive_distance_m`（高纬或大距离自动切换大圆），否则使用局地平面近似。
        adaptive_lat_threshold (float): |lat| 高于该阈值触发自适应大圆距离计算。
        adaptive_distance_threshold_km (float): 平面近似距离超过该阈值(km)触发大圆距离计算。
        force_great_circle_circle (bool): 强制半径距离全部使用大圆（忽略阈值条件）。

    返回:
        pd.DataFrame: 匹配的Argo剖面完整数据（所有深度层级）。若无匹配返回空 DataFrame。

    额外优化:
        track 参数可传入事先加载好的涡旋轨迹（list 或 DataFrame），以避免函数内部再次调用 find_track 触发重复磁盘 I/O。
    """
    # 参数回退
    if argo_data_dir is None:
        argo_data_dir = argo_path
    if circle_enlargement_factor is None:
        circle_enlargement_factor = globals().get('circle_enlargement_factor', 1.2)

    # --- 1. 准备涡旋数据 ---
    # print(f"[*] Preparing track data for eddy ID {no}...")
    if track is None:
        wanted_track = find_track(DS, no)
        # 兼容 find_track 返回 DataFrame 的新默认行为
        if wanted_track is None or (isinstance(wanted_track, pd.DataFrame) and wanted_track.empty):
            print(f"  - Track for eddy {no} not found, returning empty result.")
            return pd.DataFrame()
        if isinstance(wanted_track, pd.DataFrame):
            track_df = wanted_track.copy()
        else:
            # 将涡旋轨迹列表转换为DataFrame，方便后续处理（旧结构）
            track_df = pd.DataFrame(
                wanted_track,
                columns=['track_id', 'time', 'center_lon', 'center_lat', 'max_lon', 'max_lat', 
                         'contour_lon', 'contour_lat', 'radius', 'speed_contour_lon', 'speed_contour_lat']
            )
    elif isinstance(track, pd.DataFrame):
        track_df = track.copy()
        # 尝试确保必要列存在
        needed_cols = {'track_id','time','center_lon','center_lat','max_lon','max_lat','contour_lon','contour_lat','radius','speed_contour_lon','speed_contour_lat'}
        missing = needed_cols - set(track_df.columns)
        if missing:
            raise ValueError(f"track DataFrame is missing required columns: {sorted(list(missing))}")
    else:
        # 视为 list[list | tuple]
        track_df = pd.DataFrame(
            track,
            columns=['track_id', 'time', 'center_lon', 'center_lat', 'max_lon', 'max_lat', 
                     'contour_lon', 'contour_lat', 'radius', 'speed_contour_lon', 'speed_contour_lat']
        )
    # 使用convert_date函数转换日期（如果没有或类型不对）
    if 'date' not in track_df.columns:
        track_df['date'] = convert_date(track_df['time'])

    # --- 2. 动态加载所需的Argo数据 ---
    min_year, max_year = track_df['date'].min().year, track_df['date'].max().year
    print(f"[*] Eddy track spans years: {min_year} - {max_year}. Loading corresponding Argo data...")
    
    all_argo_data = []
    for year in range(min_year, max_year + 1):
        try:
            yearly_argo_data = load_argo_data(year, data_dir=argo_data_dir)
            if not yearly_argo_data.empty:
                all_argo_data.append(yearly_argo_data)
        except FileNotFoundError:
            print(f"  - Warning: Argo data file for year {year} not found, skipping.")

    if not all_argo_data:
        print(f"  - No Argo data was loaded, returning empty result.")
        return pd.DataFrame()
        
    argo_data = pd.concat(all_argo_data, ignore_index=True)
    # print(f"[*] All relevant yearly Argo data loaded, total records: {len(argo_data):,}.")

    # --- 3. 准备Argo数据用于合并 ---
    # 预处理：删除没有有效坐标的记录，并创建日期列
    argo_data.dropna(subset=['Longitude', 'Latitude'], inplace=True)
    argo_data['date'] = pd.to_datetime(argo_data[['Year', 'Month', 'Day']])
    
    # 每个剖面只需要一个位置点来做匹配，因此我们按Profile_number去重，只保留第一个出现的点
    argo_positions = argo_data.drop_duplicates(subset='Profile_number', keep='first').copy()
    
    # --- 4. 核心步骤：将Argo与涡旋数据按日期合并 ---
    # pd.merge会高效地找出在同一天既有涡旋轨迹、又有Argo浮标的记录
    merged_df = pd.merge(argo_positions, track_df, on='date')
    if merged_df.empty:
        print(f"  - No Argo floats found on the same day as the eddy track.")
        return pd.DataFrame()
    # print(f"[*] Found {len(merged_df)} potential 'float-eddy' pairs based on date matching.")

    # --- 5. 向量化地理筛选 ---
    # print(f"[*] Performing geographic location matching...")
    # 5.1 检查是否在多边形内部
    inside_poly_mask = merged_df.apply(is_point_in_contour, axis=1)

    # 5.2 圆内判定
    if use_adaptive_circle:
        dist_m = adaptive_distance_m(
            merged_df['Longitude'].values,
            merged_df['Latitude'].values,
            merged_df['center_lon'].values,
            merged_df['center_lat'].values,
            wrap_dateline=True,
            gc_lat_threshold=adaptive_lat_threshold,
            gc_distance_threshold_km=adaptive_distance_threshold_km,
            force_great_circle=force_great_circle_circle
        )
    else:
        scale_all = approximate_degree_length(merged_df['center_lat'].values)
        scale_lat = scale_all['meters_per_degree_lat']
        scale_lon = scale_all['meters_per_degree_lon']
        dlon_deg = _minimal_lon_diff_deg(
            merged_df['Longitude'].values,
            merged_df['center_lon'].values
        )
        dx_m = dlon_deg * scale_lon
        dy_m = (merged_df['Latitude'].values - merged_df['center_lat'].values) * scale_lat
        dist_m = np.hypot(dx_m, dy_m)
    inside_circle_mask = dist_m <= (merged_df['radius'].values * circle_enlargement_factor)
    
    # 5.3 合并两种筛选条件
    final_mask = inside_poly_mask | inside_circle_mask
    
    # 筛选出真正匹配的记录
    matched_profiles = merged_df[final_mask]
    
    if matched_profiles.empty:
        print(f"  - After geographic filtering, no matching Argo floats were found.")
        return pd.DataFrame()

    # --- 6. 根据筛选结果获取完整的剖面数据 ---
    # 从原始的、包含所有深度数据的argo_data中，筛选出所有匹配上的Profile_number
    matching_profile_numbers = matched_profiles['Profile_number'].unique()
    final_argo_data = argo_data[argo_data['Profile_number'].isin(matching_profile_numbers)].copy()
    
    # 在返回结果前，删除临时的'date'辅助列
    if 'date' in final_argo_data.columns:
        final_argo_data = final_argo_data.drop(columns=['date'])
    
    print(f"[*] Filtering complete! Found {len(matching_profile_numbers)} matching Argo profiles, with a total of {len(final_argo_data):,} data records.")
    
    return final_argo_data

def _resolve_track_context(
    DS_input: list | str | tuple | dict,
    track_id: int,
    *,
    include_contours: bool = True,
) -> tuple[pd.DataFrame, str, list | str | tuple | dict]:
    """Normalize various dataset inputs into a track DataFrame and metadata."""

    def _infer_dataset_name_from_caller(obj) -> str:
        frame = inspect.currentframe()
        plot_frame = None
        caller_frame = None
        try:
            plot_frame = frame.f_back
            caller_frame = plot_frame.f_back if plot_frame else None
            search_frame = caller_frame or plot_frame
            while search_frame:
                for var_name, var_val in search_frame.f_locals.items():
                    if var_val is obj:
                        return var_name.upper()
                search_frame = search_frame.f_back if search_frame is caller_frame else None
        except Exception:
            pass
        finally:
            del frame
            if plot_frame:
                del plot_frame
            if caller_frame:
                del caller_frame
        return "UNKNOWN"

    def _is_kind_sequence(value) -> bool:
        return isinstance(value, (list, tuple)) and len(value) > 0 and all(isinstance(item, str) for item in value)

    def _build_track_df(track_payload):
        required_cols = [
            'track_id', 'time', 'center_lon', 'center_lat', 'max_lon', 'max_lat',
            'contour_lon', 'contour_lat', 'radius', 'speed_contour_lon', 'speed_contour_lat'
        ]
        if isinstance(track_payload, pd.DataFrame):
            missing_cols = [col for col in required_cols if col not in track_payload.columns]
            if missing_cols:
                raise ValueError(f"Track data missing required columns: {missing_cols}")
            df = track_payload[required_cols].copy()
        else:
            df = pd.DataFrame(track_payload, columns=required_cols)
        if 'date' not in df.columns:
            df['date'] = convert_date(df['time'])
        return df

    ds_name = "UNKNOWN"
    ds_source_for_filter: list | str | tuple | dict = DS_input
    wanted_track = None
    candidate_kind_names: list[str] = []

    if isinstance(DS_input, str):
        candidate_kind_names = [DS_input]
        ds_name = DS_input.upper()
    elif _is_kind_sequence(DS_input):
        candidate_kind_names = [str(k) for k in DS_input]
    elif isinstance(DS_input, dict):
        for key, value in DS_input.items():
            if isinstance(value, list) and not _is_kind_sequence(value):
                try:
                    wanted_track = find_track(value, track_id)
                except ValueError:
                    continue
                ds_name = str(key).upper()
                ds_source_for_filter = value
                break
        if wanted_track is None:
            candidate_kind_names = [str(k) for k in DS_input.keys()]
    elif isinstance(DS_input, list):
        wanted_track = find_track(DS_input, track_id)
        ds_name = _infer_dataset_name_from_caller(DS_input)
    elif isinstance(DS_input, tuple):
        ds_list = list(DS_input)
        wanted_track = find_track(ds_list, track_id)
        ds_source_for_filter = ds_list
        ds_name = _infer_dataset_name_from_caller(ds_list)
    else:
        raise TypeError("Unsupported DS type. Expected legacy list, kind string, tuple, or dataset dict.")

    if wanted_track is None and candidate_kind_names:
        last_error: Exception | None = None
        for kind in candidate_kind_names:
            kind_str = str(kind).strip()
            if not kind_str:
                continue
            try:
                wanted_track = find_track(kind_str.lower(), track_id, include_contours=include_contours)
                ds_source_for_filter = kind_str.lower()
                ds_name = kind_str.upper()
                break
            except Exception as exc:
                last_error = exc
        if wanted_track is None:
            if last_error is not None:
                raise ValueError(f"Track for eddy {track_id} not found in datasets {candidate_kind_names}: {last_error}")
            raise ValueError(f"Track for eddy {track_id} not found in datasets {candidate_kind_names}.")

    if wanted_track is None:
        raise ValueError(f"Track for eddy {track_id} not found.")

    try:
        track_df = _build_track_df(wanted_track)
    except Exception as exc:
        raise ValueError(f"Failed to normalize track data for eddy {track_id}: {exc}") from exc

    if track_df.empty:
        raise ValueError(f"Track for eddy {track_id} is empty.")

    return track_df, ds_name, ds_source_for_filter

def plot_track(
    DS: list | str | tuple | dict, 
    no: int,
    save_fig: bool = False,
    show_fig: bool = True,
    plot_radius: bool = False,
    connection_threshold_days: int = 5,
    do_threshold: float | None = None,
    salinity_threshold: float | None = None,
    temperature_threshold: float | None = None,
    depth_interval: float | None = None,
    depth_merge_tolerance: float | None = None,
    duplicate_depth_strategy: str | None = None,
    anomaly_min_depth: float | None = None,
    plot_unrelated_argo: bool = True,
    fix_delta_do_colorbar: bool = True,
    delta_do_cbar_min: float = 50.0,
    delta_do_cbar_max: float = 100.0,
    delta_do_cbar_ticks: list | None = None,
    min_anomaly_count: int = 0
):
    """
    绘制指定编号涡旋的详细轨迹，并智能高亮显示与 Argo 剖面的 ΔDO 异常交互情况。

    功能:
        1. 自动加载并筛选与指定涡旋匹配的Argo数据。
        2. 绘制涡旋的完整轨迹（虚线）。
        3. 智能高亮Argo浮标存在的时期：当浮标存在日的间隔小于阈值时，
           会将这段完整的涡旋轨迹绘制为连续实线（孤立的单日则标记为点）。
        4. 使用 calculate_delta_do 识别 ΔDO 异常（取每个剖面最大 ΔDO 一条，支持按最小深度过滤），并按 ΔDO 着色。
        5. 可选地绘制涡旋在交互日的有效半径和轮廓。

    参数:
        DS (list | str | sequence[str] | dict): 
            legacy 模式传入已加载的数据列表（如 ACL/ACS）；
            新模式可直接传入字符串 kind（'acs' 等）或字符串列表/元组，
            函数会自动从本地 META_tracks 中检索对应轨迹；
            亦支持传入 {"ACS": acs, ...} 字典以兼容旧流程。
        no (int): 
            需要绘制的涡旋的唯一编号。
        save_fig (bool, optional): 
            是否将生成的图像保存到文件。默认为 False。
        show_fig (bool, optional): 
            是否在交互式环境中显示生成的图像。默认为 True。
        plot_radius (bool, optional): 
            是否以圆的形式绘制涡旋在交互日的有效半径。默认为 False。
        connection_threshold_days (int, optional):
            连接 Argo 交互点的最大天数阈值。默认为 5 天。
        do_threshold / salinity_threshold / temperature_threshold / depth_interval / depth_merge_tolerance / duplicate_depth_strategy: 
            传递给 calculate_delta_do；若为 None 则回退到全局配置默认。
        anomaly_min_depth (float | None): 
            ΔDO 异常最小深度限制；≤0 表示不限制；None 表示使用 processing.yml 的 anomaly_min_depth 配置值。
        plot_unrelated_argo (bool): 
            是否绘制被 ΔDO 筛选掉的所有匹配 Argo 剖面基准位置（空心灰圈）。
        fix_delta_do_colorbar (bool): 
            是否固定 ΔDO 色标范围。
        delta_do_cbar_min / delta_do_cbar_max / delta_do_cbar_ticks: 
            色标范围与刻度设置。
        min_anomaly_count (int): 
            若 >0：要求 ΔDO 异常数量 ≥ 该值才绘图；=0 表示不做数量阈值过滤（默认 0）。
    """
    # --- 1. 准备涡旋和Argo数据 ---
    print(f"[*] Preparing data for eddy ID {no}...")

    # 参数回退
    if do_threshold is None:
        do_threshold = _default_delta_do_threshold
    if salinity_threshold is None:
        salinity_threshold = _default_salinity_threshold
    if temperature_threshold is None:
        temperature_threshold = _default_temperature_threshold
    if depth_interval is None:
        depth_interval = _default_depth_interval
    if depth_merge_tolerance is None:
        depth_merge_tolerance = _default_depth_merge_tolerance
    if duplicate_depth_strategy is None:
        duplicate_depth_strategy = _default_duplicate_depth_strategy
    if anomaly_min_depth is None:
        anomaly_min_depth = _cfg_anomaly_min_depth
    
    try:
        track_df, ds_name, ds_source_for_filter = _resolve_track_context(DS, no, include_contours=True)
    except Exception as exc:
        print(f"  - Error: {exc}")
        return
    num = track_df['track_id'].iloc[0]

    # 调用筛选函数，获取所有匹配的 Argo 数据（包含所有深度）
    argo_data_filtered = filtered_float_data(ds_source_for_filter, no, track=track_df)

    # 预筛选：若匹配到的剖面数不足 min_anomaly_count，直接跳过
    if argo_data_filtered.empty:
        print(f"  - Skip plotting {ds_name}{no}: no matched Argo profiles.")
        return
    matched_profile_count = argo_data_filtered['Profile_number'].nunique()
    if min_anomaly_count > 0 and matched_profile_count < min_anomaly_count:
        print(f"\033[31m  - Skip plotting {ds_name}{no}: matched profiles ({matched_profile_count}) < min_anomaly_count ({min_anomaly_count}).\033[0m")
        return

    # 使用 ΔDO 异常检测
    anomalies = pd.DataFrame()
    base_argo_positions = pd.DataFrame()
    if not argo_data_filtered.empty:
        # 基准剖面位置（去除深度重复）
        base_argo_positions = (
            argo_data_filtered.sort_values(['Profile_number','Depth'])
            .groupby('Profile_number', as_index=False)
            .first()[['Profile_number','Longitude','Latitude','Year','Month','Day']]
        )
        anomalies = calculate_delta_do(
            argo_data_filtered,
            depth_interval=depth_interval,
            do_threshold=do_threshold,
            salinity_threshold=salinity_threshold,
            temperature_threshold=temperature_threshold,
            anomaly_min_depth=anomaly_min_depth,
            depth_merge_tolerance=depth_merge_tolerance,
            duplicate_depth_strategy=duplicate_depth_strategy,
            remove_outliers=True,
            verbose=False
        )
        if not anomalies.empty:
            # 每个剖面取最大 ΔDO
            anomalies = (
                anomalies.sort_values('delta_do', ascending=False)
                .drop_duplicates(subset='Profile_number', keep='first')
            )

    # 若异常数量不足阈值，直接跳过
    anomaly_count = 0 if anomalies.empty else len(anomalies)
    if min_anomaly_count > 0 and anomaly_count < min_anomaly_count:
        print(f"\033[31m  - Skip plotting {ds_name}{no}: anomalies ({anomaly_count}) < min_anomaly_count ({min_anomaly_count}).\033[0m")
        return

    # --- 2. 准备绘图 ---
    prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    eddy_color = prop_colors[1] if 'AC' in ds_name else prop_colors[0]
    
    # 加载地理底图
    world = _load_world_geodataframe()

    # 初始化画布与坐标轴
    fig, ax = plt.subplots(figsize=(30, 20))
    ax.set_title(f'Track Analysis for Eddy {ds_name}{num}', fontsize=20)
    ax.set_xlabel('Longitude', fontsize=20); ax.set_ylabel('Latitude', fontsize=20)
    world.plot(color='lightgrey', edgecolor='white', ax=ax)

    # --- 3. 绘制涡旋轨迹 (背景和高亮) ---
    # 绘制完整的背景轨迹 (虚线)
    ax.plot(track_df['center_lon'], track_df['center_lat'], color=eddy_color, linestyle='--', alpha=0.5, label='Full Track Path')
    
    # 找出有Argo交互的日期
    interaction_track_df = pd.DataFrame(columns=track_df.columns)
    if not anomalies.empty:
        interaction_dates = pd.to_datetime(anomalies[['Year', 'Month', 'Day']]).unique()
        interaction_track_df = track_df[track_df['date'].isin(interaction_dates)].copy()
        
        # 识别并分别绘制不连续的实线段
        if not interaction_track_df.empty:
            date_diffs = interaction_track_df['date'].diff()
            break_points = date_diffs > pd.Timedelta(f'{connection_threshold_days} days')
            segment_ids = break_points.cumsum()
            
            labeled_interaction_path = False
            labeled_interaction_point = False
            # 按段落ID分组并分别绘图
            for _, segment_df in interaction_track_df.groupby(segment_ids):
                start_segment_date = segment_df['date'].min()
                end_segment_date = segment_df['date'].max()
                full_path_segment = track_df[
                    (track_df['date'] >= start_segment_date) & 
                    (track_df['date'] <= end_segment_date)
                ]

                # 如果段内只有一个交互点，画成散点
                if len(segment_df) == 1:
                    label = 'Interaction Point' if not labeled_interaction_point else None
                    ax.scatter(full_path_segment['center_lon'], full_path_segment['center_lat'], 
                               facecolors='none', edgecolors=eddy_color, 
                                linewidths=1.5, s=100, zorder=5, label=label)
                    labeled_interaction_point = True
                # 否则，将这段时间的完整轨迹画成实线
                else:
                    label = 'Interaction Path' if not labeled_interaction_path else None
                    ax.plot(full_path_segment['center_lon'], full_path_segment['center_lat'], 
                            color=eddy_color, linestyle='-', linewidth=2.5, label=label)
                    labeled_interaction_path = True
    
    # 标记轨迹的起点和终点
    ax.plot(track_df['center_lon'].iloc[0], track_df['center_lat'].iloc[0], marker='o', color=eddy_color, markersize=10, label='Start')
    ax.plot(track_df['center_lon'].iloc[-1], track_df['center_lat'].iloc[-1], marker='x', color=eddy_color, markersize=12, mew=2.5, label='End')

    # --- 4. 绘制交互日的轮廓、半径和日期 ---
    labeled_contour = False
    labeled_radius = False
    for _, eddy_day in interaction_track_df.iterrows():
        # 绘制轮廓线
        if not labeled_contour:
            ax.plot(eddy_day['contour_lon'], eddy_day['contour_lat'], color=eddy_color, linewidth=1, alpha=0.6, label='Effective Contour')
            labeled_contour = True
        else:
            ax.plot(eddy_day['contour_lon'], eddy_day['contour_lat'], color=eddy_color, linewidth=1, alpha=0.6)
        
        # 绘制有效半径
        if plot_radius:
            radius_color = 'r' if 'AC' in ds_name else 'purple'
            circle_label = 'Effective Radius' if not labeled_radius else None
            scale = approximate_degree_length(eddy_day['center_lat'])
            deg_height = (eddy_day['radius'] * circle_enlargement_factor) / scale['meters_per_degree_lat']
            deg_width = (eddy_day['radius'] * circle_enlargement_factor) / scale['meters_per_degree_lon']
            ell = Ellipse((eddy_day['center_lon'], eddy_day['center_lat']), width=2*deg_width, height=2*deg_height,
                          edgecolor=radius_color, facecolor='none', linestyle='--', alpha=0.4, linewidth=1.5, label=circle_label)
            ax.add_patch(ell)
            labeled_radius = True
            
        # 标记交互的起始和结束日期
        if eddy_day['date'] == interaction_track_df['date'].min() or eddy_day['date'] == interaction_track_df['date'].max():
             ax.text(eddy_day['center_lon'], eddy_day['center_lat'] + 0.1, eddy_day['date'].strftime('%Y-%m-%d'), 
                     fontsize=16, color='black', ha='center', zorder=11)

    # --- 5. 绘制 ΔDO 异常与基准剖面 ---
    if plot_unrelated_argo and not base_argo_positions.empty:
        ax.scatter(
            base_argo_positions['Longitude'], base_argo_positions['Latitude'],
            facecolors='none', edgecolors='gray', linewidths=0.8, s=50,
            label='All Matched Argo Profiles', zorder=5
        )

    if not anomalies.empty:
        scatter_kwargs = {}
        if fix_delta_do_colorbar:
            scatter_kwargs.update(dict(vmin=delta_do_cbar_min, vmax=delta_do_cbar_max))
        depth_label = (
            f' @ depth ≥ {anomaly_min_depth} m'
            if anomaly_min_depth is not None and anomaly_min_depth > 0
            else ''
        )
        sc = ax.scatter(
            anomalies['Longitude'], anomalies['Latitude'],
            c=anomalies['delta_do'], cmap='Reds', s=90,
            edgecolors='black', linewidths=0.6,
            label=f'ΔDO ≥ {do_threshold} μmol kg⁻¹{depth_label}',
            zorder=10,
            **scatter_kwargs
        )
        cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.08)
        cbar.set_label('ΔDO / μmol·kg⁻¹', fontsize=18)
        cbar.ax.tick_params(labelsize=14)
        if fix_delta_do_colorbar:
            if delta_do_cbar_ticks is not None:
                cbar.set_ticks(delta_do_cbar_ticks)
            else:
                rng = delta_do_cbar_max - delta_do_cbar_min
                if rng > 30:
                    mid = (delta_do_cbar_min + delta_do_cbar_max) / 2
                    cbar.set_ticks([delta_do_cbar_min, mid, delta_do_cbar_max])
                else:
                    cbar.set_ticks([delta_do_cbar_min, delta_do_cbar_max])

    # --- 6. 最终化绘图设置 ---
    # 设定边界时排除META中错误的contour数据
    valid_contours_lon = [lon for day_lons in track_df['contour_lon'] for lon in day_lons if lon != 180.0]
    valid_contours_lat = [lat for day_lats in track_df['contour_lat'] for lat in day_lats if lat != 0.0]
    if valid_contours_lon and valid_contours_lat:
        ax.set_xlim(min(valid_contours_lon) - 0.5, max(valid_contours_lon) + 0.5)
        ax.set_ylim(min(valid_contours_lat) - 0.5, max(valid_contours_lat) + 0.5)
        
    ax.set_aspect('equal')
    ax.legend(fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # --- 7. 输出控制 ---
    if save_fig:
        region_slug = _current_region_key()
        output_dir = Path(plots_output_root) / region_slug / "plot_track"
        output_dir.mkdir(exist_ok=True, parents=True)
        base_filename = f"Track_Analysis_{ds_name}{num}.png"
        save_path = output_dir / base_filename
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
    
    if show_fig:
        plt.show()
    
    plt.close(fig)

def plot_track_variable_timeseries(
    DS: list | str | tuple | dict,
    no: int,
    variable: str = 'DO',
    threshold: float | None = None,
    only_above_threshold: bool = True,
    depth_col: str = 'Depth',
    max_depth: float | None = 1000.0,
    depth_bin_size: float = 25.0,
    cmap: str = 'RdYlBu_r',
    show_profile_hist: bool = True,
    start_date: str | int | float | pd.Timestamp | None = None,
    end_date: str | int | float | pd.Timestamp | None = None,
    platform_number: int | str | list | tuple | set | None = None,
    salinity_threshold: float | None = None,
    temperature_threshold: float | None = None,
    depth_interval: float | None = None,
    depth_merge_tolerance: float | None = None,
    duplicate_depth_strategy: str | None = None,
    anomaly_min_depth: float | None = None,
    save_fig: bool = False,
    show_fig: bool = True,
):
    """绘制位于涡旋内的 Argo 剖面变量在时间-深度平面上的连续等值分布。
    

    参数:
        DS: legacy轨迹列表、kind字符串、字符串序列或数据集字典。
        no: 涡旋编号。
        variable: 需要统计的变量列名，默认 'DO'。
        threshold: 变量阈值；only_above_threshold=True 时用于筛选，None 回退到 `_default_delta_do_threshold`。
        only_above_threshold: True 时仅统计变量值 ≥ threshold 的剖面；False 使用全部可用剖面。
        depth_col: 深度列名，默认 'Depth'。
        max_depth: 最大深度（单位与 depth_col 一致），默认 1000 dbar；None 表示使用观测最大值。
        depth_bin_size: 深度分箱大小（dbar），默认 25。
        cmap: 色标名称，默认 'RdYlBu_r'。
        show_profile_hist: 是否在底部显示每日剖面数柱状条，默认 True。
        start_date / end_date: 限制横轴开始/结束日期，可传 pandas.Timestamp、YYYYMMDD 整数、days-since-1950、ISO 字符串。
        platform_number: 可选，仅绘制指定 Argo 浮标（单个或列表/集合），先于 ΔDO 筛选过滤。
        salinity_threshold / temperature_threshold: 传递给 calculate_delta_do 的 ΔS/ΔT 条件，None 时使用 processing.yml 默认。
        depth_interval / depth_merge_tolerance / duplicate_depth_strategy: 传递给 calculate_delta_do 的参数，控制 ΔDO 计算窗口与合并方式。
        anomaly_min_depth: ΔDO 最小深度阈值（同 plot_all_tracks_in_range）；None 回退配置。
        save_fig / show_fig: 控制图像输出。
    """
    print(f"[*] Building {variable} time series for eddy {no}...")

    if threshold is None:
        threshold = _default_delta_do_threshold

    def _normalize_date_input(value, label):
        if value is None:
            return None
        if isinstance(value, pd.Timestamp):
            return value.normalize()
        try:
            converted = convert_date(value)
        except Exception as exc:
            raise ValueError(f"Invalid {label}: {exc}")
        if isinstance(converted, pd.Series):
            converted = converted.iloc[0]
        if pd.isna(converted):
            raise ValueError(f"Invalid {label}: could not parse {value!r}")
        return pd.to_datetime(converted).normalize()

    try:
        start_dt = _normalize_date_input(start_date, 'start_date')
        end_dt = _normalize_date_input(end_date, 'end_date')
    except ValueError as exc:
        print(f"  - {exc}")
        return

    if start_dt is not None and end_dt is not None and start_dt > end_dt:
        print("  - start_date must be earlier than or equal to end_date.")
        return

    try:
        track_df, ds_name, ds_source_for_filter = _resolve_track_context(DS, no, include_contours=False)
    except Exception as exc:
        print(f"  - Error: {exc}")
        return

    argo_data = filtered_float_data(ds_source_for_filter, no, track=track_df)
    if argo_data.empty:
        print(f"  - Skip plotting {ds_name}{no}: no matched Argo profiles.")
        return

    if variable not in argo_data.columns:
        print(f"  - Column '{variable}' not found in matched Argo data.")
        return

    if depth_col not in argo_data.columns:
        print(f"  - Column '{depth_col}' not found in matched Argo data.")
        return

    argo_data = argo_data.copy()
    if platform_number is not None:
        if 'Platform_number' not in argo_data.columns:
            print("  - Column 'Platform_number' not available; cannot filter by platform_number.")
            return
        pn_list = platform_number
        if not isinstance(pn_list, (list, tuple, set)):
            pn_list = [pn_list]
        argo_data = argo_data[argo_data['Platform_number'].isin(pn_list)]
        if argo_data.empty:
            print("  - No data after filtering by platform_number.")
            return
    argo_data['date'] = pd.to_datetime(argo_data[['Year', 'Month', 'Day']], errors='coerce')
    if argo_data['date'].isna().all():
        print("  - No valid timestamps/depths available after parsing Year/Month/Day.")
        return
    if start_dt is not None:
        argo_data = argo_data[argo_data['date'] >= start_dt]
    if end_dt is not None:
        argo_data = argo_data[argo_data['date'] <= end_dt]
    if argo_data.empty:
        print("  - No data within the requested date window.")
        return

    if 'Profile_number' not in argo_data.columns:
        print("  - Column 'Profile_number' is required to align with ΔDO anomaly filtering.")
        return

    anomaly_filtered = False
    if only_above_threshold:
        anomalies = calculate_delta_do(
            argo_data,
            depth_col=depth_col,
            depth_interval=depth_interval,
            do_threshold=threshold,
            salinity_threshold=salinity_threshold,
            temperature_threshold=temperature_threshold,
            anomaly_min_depth=anomaly_min_depth,
            depth_merge_tolerance=depth_merge_tolerance,
            duplicate_depth_strategy=duplicate_depth_strategy,
            remove_outliers=True,
            verbose=False,
        )
        if anomalies.empty or 'Profile_number' not in anomalies.columns:
            print("  - No profiles satisfy the ΔDO anomaly criteria.")
            return
        anomalies_sorted = anomalies.sort_values(by='delta_do', ascending=False)
        anomalies_unique = anomalies_sorted.drop_duplicates(subset='Profile_number', keep='first')
        qualifying_profiles = anomalies_unique['Profile_number'].dropna().unique()
        if qualifying_profiles.size == 0:
            print("  - No profiles satisfy the ΔDO anomaly criteria.")
            return
        argo_data = argo_data[argo_data['Profile_number'].isin(qualifying_profiles)].copy()
        if argo_data.empty:
            print("  - No profiles remain after applying ΔDO anomaly filter.")
            return
        anomaly_filtered = True

    argo_data['_var'] = pd.to_numeric(argo_data[variable], errors='coerce')
    argo_data['_depth'] = pd.to_numeric(argo_data[depth_col], errors='coerce')
    argo_data.dropna(subset=['date', '_depth'], inplace=True)
    if argo_data.empty:
        print("  - No valid timestamps/depths available after parsing Year/Month/Day.")
        return

    mask = ~argo_data['_var'].isna()
    if not anomaly_filtered and threshold is not None:
        mask &= argo_data['_var'] >= threshold
    if max_depth is not None:
        mask &= argo_data['_depth'] <= max_depth
    selected = argo_data.loc[mask]
    if selected.empty:
        print("  - No profiles satisfy the selected threshold/variable conditions.")
        return

    if depth_bin_size <= 0:
        raise ValueError("depth_bin_size must be positive.")

    if max_depth is None:
        max_depth_val = float(selected['_depth'].max())
    else:
        max_depth_val = float(max_depth)
    depth_min_val = float(max(selected['_depth'].min(), 0.0))
    depth_bins = np.arange(depth_min_val, max_depth_val + depth_bin_size, depth_bin_size)
    if depth_bins.size < 2:
        depth_bins = np.array([depth_min_val, depth_min_val + depth_bin_size])
    selected = selected.copy()
    selected['depth_bin'] = pd.cut(selected['_depth'], bins=depth_bins, include_lowest=True, right=False)
    selected.dropna(subset=['depth_bin'], inplace=True)
    selected['depth_mid'] = selected['depth_bin'].apply(lambda iv: iv.left + depth_bin_size / 2 if pd.notna(iv) else np.nan)
    selected.dropna(subset=['depth_mid'], inplace=True)

    grouped = (
        selected.groupby(['depth_mid', 'date'], observed=False)['_var']
        .mean()
        .reset_index()
    )
    if grouped.empty:
        print("  - No averages could be computed for the selected bins.")
        return

    pivot = grouped.pivot(index='depth_mid', columns='date', values='_var')
    pivot.sort_index(inplace=True)
    pivot = pivot.sort_index(axis=1)
    if pivot.isna().all().all():
        print("  - All grid cells are NaN; cannot plot.")
        return

    depth_vals = pivot.index.to_numpy()
    date_vals = pivot.columns.to_pydatetime()
    if len(date_vals) == 0:
        print("  - No valid dates after grouping.")
        return

    Z = np.ma.masked_invalid(pivot.to_numpy())
    if np.ma.count_masked(Z) == Z.size:
        print("  - No valid data points remain after masking.")
        return

    vmin = np.nanmin(pivot.values)
    vmax = np.nanmax(pivot.values)
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        print("  - Unable to determine color scale (all values NaN/inf).")
        return
    if np.isclose(vmin, vmax):
        levels = 10
    else:
        levels = 20

    date_nums = mdates.date2num(date_vals)
    if show_profile_hist:
        hist_counts = selected.groupby('date').size().reindex(pivot.columns, fill_value=0)

    height_ratios = [5, 1] if show_profile_hist else [1]
    nrows = 2 if show_profile_hist else 1
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=1,
        sharex=True,
        figsize=(20, 10 if show_profile_hist else 8),
        gridspec_kw={'height_ratios': height_ratios},
        constrained_layout=True
    )
    if show_profile_hist:
        ax, ax_hist = axes
    else:
        ax = axes

    cf = ax.contourf(date_vals, depth_vals, Z, levels=levels, cmap=cmap)
    ax.invert_yaxis()
    ax.set_ylabel('Depth (dbar)', fontsize=14)
    title_suffix = ' (thresholded)' if only_above_threshold else ' (all profiles)'
    ax.set_title(f"{ds_name}{no} {variable} mean inside eddy{title_suffix}", fontsize=16)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    if start_dt is not None or end_dt is not None:
        xmin = start_dt if start_dt is not None else date_vals[0]
        xmax = end_dt if end_dt is not None else date_vals[-1]
        ax.set_xlim(xmin, xmax)

    cbar_axes = [ax, ax_hist] if show_profile_hist else ax
    cbar = fig.colorbar(
        cf,
        ax=cbar_axes,
        orientation='vertical',
        pad=0.12,
        fraction=0.035,
        location='right'
    )
    cbar.set_label(f'{variable} mean', fontsize=14)

    if show_profile_hist:
        ax_hist.bar(hist_counts.index, hist_counts.values, width=0.8, color='gray')
        ax_hist.set_ylabel('Profiles', fontsize=12)
        ax_hist.set_xlabel('Date', fontsize=14)
        ax_hist.set_ylim(0, max(hist_counts.max() * 1.2, 1))
        ax_hist.xaxis_date()
        ax_hist.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        if start_dt is not None or end_dt is not None:
            xmin = start_dt if start_dt is not None else hist_counts.index.min()
            xmax = end_dt if end_dt is not None else hist_counts.index.max()
            ax_hist.set_xlim(xmin, xmax)
    else:
        ax.set_xlabel('Date', fontsize=14)

    fig.autofmt_xdate()

    if save_fig:
        region_slug = _current_region_key()
        output_dir = Path(plots_output_root) / region_slug / "plot_track_timeseries"
        output_dir.mkdir(parents=True, exist_ok=True)
        thr_suffix = ''
        if only_above_threshold and threshold is not None:
            thr_suffix = f"_thr{str(threshold).replace('.', 'p')}"
        depth_suffix = f"_depth{str(max_depth).replace('.', 'p')}" if max_depth is not None else ''
        base_filename = f"{ds_name}{no}_{variable}_timeseries{thr_suffix}{depth_suffix}.png"
        save_path = output_dir / base_filename
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")

    if show_fig:
        plt.show()
    plt.close(fig)
    
def convert_date(values: pd.Series | np.ndarray | list | str | int | float) -> pd.Series | pd.Timestamp:
    """统一将整数/字符串编码日期转为 pandas datetime（日精度）。

    支持两种编码：
      1) 自 1950-01-01 起的天数 (CF 常见 time 轴)
      2) YYYYMMDD 8 位整数 (如 20220131)

    返回：
      - 标量输入（如单个字符串/数字）→ 返回单个 `pd.Timestamp`
      - 序列输入（Series/ndarray/list）→ 返回 `pd.Series` (dtype=datetime64[ns])，未能解析的元素为 NaT。
    """
    is_scalar_input = isinstance(values, (str, bytes, int, float, np.integer, np.floating, np.str_))
    if isinstance(values, pd.Series):
        ser = values
    elif is_scalar_input:
        ser = pd.Series([values])
    else:
        # 注意：若 values 是字符串且不想被逐字符拆分，需封装为列表
        if isinstance(values, (str, bytes, np.str_)):
            ser = pd.Series([values])
        else:
            ser = pd.Series(values)
    # 结果容器，避免在原 float Series 上就地赋值导致 dtype 警告
    result = pd.Series(pd.NaT, index=ser.index, dtype='datetime64[ns]')

    non_na = ser.dropna()
    if non_na.empty:
        return result.iloc[0] if is_scalar_input else result

    # 若输入已是 datetime64 类型，直接标准化返回
    if pd.api.types.is_datetime64_any_dtype(non_na):
        dt = pd.to_datetime(non_na).dt.normalize()
        result.loc[dt.index] = dt.values
        return result.iloc[0] if is_scalar_input else result

    # 判定是否全为 8 位 YYYYMMDD
    try:
        strs = non_na.astype('int64').astype(str)
    except Exception:
        strs = pd.Series([str(x) for x in non_na], index=non_na.index)
    is_all_yyyymmdd = strs.map(len).eq(8).all()

    if is_all_yyyymmdd:
        parsed = []
        for s in strs:
            try:
                y = int(s[0:4]); m = int(s[4:6]); d = int(s[6:8])
                parsed.append(pd.Timestamp(year=y, month=m, day=d))
            except Exception as e:
                raise ValueError(f"Invalid YYYYMMDD integer {s}: {e}")
        dt = pd.Series(parsed, index=strs.index)
        result.loc[dt.index] = dt.values
        return result.iloc[0] if is_scalar_input else result

    # 默认按 days since 1950-01-01 解析（对浮点会取整）
    t0 = pd.Timestamp('1950-01-01')
    try:
        deltas = pd.to_timedelta(non_na.astype('int64'), unit='D')
        dt = t0 + deltas
        result.loc[dt.index] = dt.values
        return result.iloc[0] if is_scalar_input else result
    except Exception as e:
        raise ValueError(f"Values neither valid YYYYMMDD nor days-since-1950: {e}")

def convert_eddy_number(
    kind: str,
    value: str | int | float | list | tuple | np.ndarray,
    order: str,
    *,
    meta_root: str | os.PathLike | None = None,
    version: float = 3.2
) -> int | list[int]:
    """
    在 legacy orig_index 与标准 track_id 之间转换涡旋编号。

    参数:
        kind: 涡旋类型，'acs'|'acl'|'cs'|'cl'。必须指定，以保证索引映射到正确数据集。
        value: 单个编号或编号可迭代。可为旧编号（orig_index）或新编号（track_id）。
        order: 转换方向，支持别名：
            • legacy->new: 'old_to_new', 'orig_to_track', 'legacy_to_track', 'to_track'.
            • new->legacy: 'new_to_old', 'track_to_orig', 'track_to_legacy', 'to_orig'.
        meta_root: META 原始 NetCDF 目录，默认读取配置 paths.meta_root。
        version: META 版本（3.1 或 3.2）。

    返回:
        按输入形状返回转换后的编号：
            • 标量输入 -> int
            • 可迭代输入 -> list[int]

    说明:
        - 旧编号对应 NetCDF 中的全局行索引（orig_index）。
        - 新编号为轨迹唯一的 track_id。
        - new->old 方向返回该轨迹首次出现的行索引（最早时间点）。
    """

    direction = order.lower().strip()
    if direction in {'old_to_new', 'orig_to_track', 'legacy_to_track', 'to_track', 'orig_to_new'}:
        mode = 'orig_to_track'
    elif direction in {'new_to_old', 'track_to_orig', 'track_to_legacy', 'to_orig'}:
        mode = 'track_to_orig'
    else:
        raise ValueError("order must be one of: old_to_new/orig_to_track/legacy_to_track/to_track or new_to_old/track_to_orig/track_to_legacy/to_orig")

    kind_norm = kind.lower().strip()
    if kind_norm not in {'acs', 'acl', 'cs', 'cl'}:
        raise ValueError("kind must be one of 'acs', 'acl', 'cs', 'cl'")

    ACS, ACL, CS, CL = load_meta_data(path=meta_root, version=version)
    ds_map = {'acs': ACS, 'acl': ACL, 'cs': CS, 'cl': CL}
    ds = ds_map[kind_norm]

    # 将 track 数组载入内存，便于快速索引或反查
    track_arr = np.asarray(ds.variables['track'][:], dtype=np.int64)

    def _convert_one(val) -> int:
        try:
            n = int(val)
        except Exception as e:
            raise ValueError(f"Cannot convert value '{val}' to int") from e

        if mode == 'orig_to_track':
            if n < 0 or n >= track_arr.shape[0]:
                raise ValueError(f"orig_index {n} out of range for {kind_norm}")
            return int(track_arr[n])

        # track_id -> 最早出现的 orig_index
        hits = np.nonzero(track_arr == n)[0]
        if hits.size == 0:
            raise ValueError(f"track_id {n} not found in {kind_norm}")
        return int(hits[0])

    is_scalar = np.isscalar(value) or isinstance(value, (str, bytes, np.str_))
    if is_scalar:
        return _convert_one(value)

    try:
        values = list(value)
    except Exception as e:
        raise ValueError("value must be scalar or iterable of scalars") from e

    return [_convert_one(v) for v in values]

def plot_vertical(
    DS: list,
    no: int,
    show_fig: bool = False,
    save_fig: bool = False,
    color_mode: str = 'distance',
    variables: list = ['DO', 'Temp', 'Salinity'],
    show_colorbar: bool = False,
    remove_outliers: bool = True,
    aggregated: bool = False,
    argo_required: list | None = None,
    year_required: list | None = None,
    month_required: list | None = None,
    day_required: list | None = None
):
    '''
    根据涡旋轨迹与匹配到的 Argo 剖面，绘制变量-深度的垂直剖面。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        show_fig (bool): 是否显示图片，默认 False。
        save_fig (bool): 是否保存图片，默认 False。
        color_mode (str): 颜色模式，'distance' 或 'time'，默认 'distance'。
        variables (list): 需要绘制的变量名称，默认 ['DO', 'Temp', 'Salinity']。
        show_colorbar (bool): 是否显示颜色条，默认 False。
        remove_outliers (bool): 是否执行 QC 过滤与规则法去极值，默认 True。
        aggregated (bool): 是否进行跨平台聚合绘制，默认 False。
        argo_required (list | None): 平台过滤；None 表示不过滤；传入平台编号列表时仅保留指定平台。
        year_required (list | None): 年份过滤；None 表示不过滤；传入年份列表时仅保留指定年份。
        month_required (list | None): 月份过滤；聚合模式 None 表示使用所有可用月份；逐平台模式 None 表示不过滤。
        day_required (list | None): 日期过滤（按日1-31）；None 表示不过滤；传入日期列表时仅保留指定日期。

    功能:
        - 为 variables 中的每个变量创建一个子图，按剖面绘制变量随深度变化的曲线。
        - 曲线颜色可根据与涡旋中心的相对距离（distance）或采样时间（time）变化。
        - 可选显示颜色条；支持图片保存与显示。
        - remove_outliers=True 时执行基础质量控制（QC 仅保留 {1,2,5,8}；DO<=1 置为 NaN）。
        - 可用 month_required 和 argo_required 对数据进行月份与平台筛选。

    模式差异:
        - aggregated=False（逐平台）：
            • 为每个浮标平台单独出图；每图包含 variables 中各变量的一个子图。
        - aggregated=True（聚合）：
            • 所有平台剖面在同一张图上聚合绘制（每个变量一个子图）。
    '''
    wanted_track = find_track(DS, no)
    # 复用已读取的 wanted_track，避免 filtered_float_data 内部重复 find_track
    track_df = pd.DataFrame(
        wanted_track,
        columns=['track_id', 'time', 'center_lon', 'center_lat', 'max_lon', 'max_lat', 
                 'contour_lon', 'contour_lat', 'radius', 'speed_contour_lon', 'speed_contour_lat']
    )
    track_df['date'] = convert_date(track_df['time'])
    argo_data_filtered = filtered_float_data(DS, no, track=track_df)

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE: Could not automatically determine the dataset name.")

    if argo_data_filtered.empty:
        msg = "plot vertical profiles." if not aggregated else "to plot aggregated vertical profiles."
        print(f"No Argo data found for eddy {ds_names}{no} {msg}")
        return

    # 统一的月份和平台过滤（在两种模式下共享）
    # 先处理平台过滤（argo_required）
    if argo_required is not None:
        try:
            argo_required_set = set(int(x) for x in argo_required)
        except Exception:
            argo_required_set = set(argo_required)
        before_cnt = len(argo_data_filtered)
        argo_data_filtered = argo_data_filtered[argo_data_filtered['Platform_number'].isin(argo_required_set)]
        if argo_data_filtered.empty:
            print(f"No Argo data left after filtering by platforms {sorted(list(argo_required_set))}.")
            return

    # 处理年月日过滤（year/month/day）
    # 注意：aggregated=True 时，按原逻辑在未提供月份时默认使用所有可用月份；
    #       aggregated=False 时，仅在提供了 month_required 时进行过滤。
    #       year_required/day_required 在两种模式下仅当提供时进行过滤。
    
    # Year
    argo_data_filtered['Year'] = pd.to_numeric(argo_data_filtered['Year'], errors='coerce')
    argo_data_filtered.dropna(subset=['Year'], inplace=True)
    argo_data_filtered['Year'] = argo_data_filtered['Year'].astype(int)

    # Month
    argo_data_filtered['Month'] = pd.to_numeric(argo_data_filtered['Month'], errors='coerce')
    argo_data_filtered.dropna(subset=['Month'], inplace=True)
    argo_data_filtered['Month'] = argo_data_filtered['Month'].astype(int)

    # Day
    argo_data_filtered['Day'] = pd.to_numeric(argo_data_filtered['Day'], errors='coerce')
    argo_data_filtered.dropna(subset=['Day'], inplace=True)
    argo_data_filtered['Day'] = argo_data_filtered['Day'].astype(int)

    if aggregated:
        # 若未指定月份，采用所有可用月份（保持原 plot_vertical_monthly 行为）
        if not month_required:
            all_available_months = sorted(argo_data_filtered['Month'].unique().tolist())
            if not all_available_months:
                print(f"No valid months found in data for eddy {ds_names}{no}.")
                return
            print(f"No months specified. Defaulting to all available months: {all_available_months}")
            month_required = all_available_months
        # 根据年月日过滤
        if month_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Month'].isin(month_required)].copy()
        if year_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Year'].isin(year_required)].copy()
        if day_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Day'].isin(day_required)].copy()
        if argo_data_filtered.empty:
            print(f"No data found for eddy {ds_names}{no} after applying filters: year={year_required}, month={month_required}, day={day_required}.")
            return
    else:
        # 非聚合模式：仅当提供了对应过滤条件才过滤
        if year_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Year'].isin(year_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in years {year_required} (per-platform mode).")
                return
        if month_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Month'].isin(month_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in months {month_required} (per-platform mode).")
                return
        if day_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Day'].isin(day_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in days {day_required} (per-platform mode).")
                return

    # =========================
    # 分支一：非聚合（原 plot_vertical 行为，新增月份/平台筛选）
    # =========================
    if not aggregated:
        for platform_id_val, platform_data in argo_data_filtered.groupby("Platform_number"):
            profile_num_agg = platform_data['Profile_number'].agg(['min', 'max'])
            min_profile_num = profile_num_agg['min']
            max_profile_num = profile_num_agg['max']

            num_variables = len(variables)
            fig, axes = plt.subplots(1, num_variables, figsize=(10 * num_variables, 20), sharey=True)
            if num_variables == 1:
                axes = [axes]

            cmap = plt.cm.coolwarm
            profile_dates_for_title = []

            for i, var_name in enumerate(variables):
                ax = axes[i]
                original_variable_name = var_name
                db_variable_name = 'Temperature' if var_name == 'Temp' else var_name

                if db_variable_name not in platform_data.columns:
                    ax.text(0.5, 0.5, f"Variable '{db_variable_name}'\nnot found in data.",
                            ha='center', va='center', transform=ax.transAxes, fontsize=16)
                    ax.set_title(f"Variable: {original_variable_name}", fontsize=20)
                    continue

                for profile_num, rows in platform_data.groupby("Profile_number"):
                    if rows.empty:
                        continue

                    rows_to_plot = rows.dropna(subset=[db_variable_name, "Depth"])
                    if rows_to_plot.empty:
                        continue

                    if remove_outliers:
                        qc_column_name = f"{db_variable_name}_Flag"
                        if qc_column_name in rows_to_plot.columns:
                            good_qc_flags = ['1', '2', '5', '8', 1, 2, 5, 8]
                            bad_qc_mask = ~rows_to_plot[qc_column_name].isin(good_qc_flags)
                            rows_to_plot.loc[bad_qc_mask, db_variable_name] = np.nan
                        if db_variable_name == 'DO':
                            bad_value_mask = rows_to_plot[db_variable_name] <= 1.0
                            rows_to_plot.loc[bad_value_mask, db_variable_name] = np.nan

                    if rows_to_plot.empty:
                        continue

                    try:
                        current_profile_date = pd.Timestamp(year=int(rows.iloc[0]['Year']),
                                                            month=int(rows.iloc[0]['Month']),
                                                            day=int(rows.iloc[0]['Day']))
                        if i == 0:
                            profile_dates_for_title.append(current_profile_date)
                    except (ValueError, TypeError):
                        continue

                    color_value_normalized = 0.5

                    if color_mode == 'distance':
                        if 'Longitude' in rows.iloc[0] and 'Latitude' in rows.iloc[0] and wanted_track is not None and len(wanted_track) > 0:
                            track_dates_converted = wanted_track['date'] if 'date' in wanted_track.columns else convert_date(wanted_track['time'])
                            idx_track_list = [j for j, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_profile_date.date()]
                            if idx_track_list:
                                idx_track = idx_track_list[0]
                                center_lon = float(wanted_track.iloc[idx_track]['center_lon'])
                                center_lat = float(wanted_track.iloc[idx_track]['center_lat'])
                                radius = float(wanted_track.iloc[idx_track]['radius'])
                                if radius > 1e-6:
                                    dist_m = adaptive_distance_m(
                                        rows.iloc[0]['Longitude'], rows.iloc[0]['Latitude'],
                                        center_lon, center_lat
                                    )
                                    distance = dist_m / radius
                                    color_value_normalized = 1.0 - np.clip(distance, 0.0, 1.0)
                    elif color_mode == 'time':
                        if max_profile_num > min_profile_num:
                            color_value_normalized = (profile_num - min_profile_num) / (max_profile_num - min_profile_num)
                        else:
                            color_value_normalized = 0.0

                    color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
                    ax.plot(rows_to_plot[db_variable_name], rows_to_plot["Depth"], color=color, alpha=0.7)

                # 子图属性
                ax.set_ylim(-50, 2050)
                if db_variable_name == 'DO':
                    ax.set_xlim(10, 350)
                elif db_variable_name == 'Temperature':
                    ax.set_xlim(1, 32)
                elif db_variable_name == 'Salinity':
                    ax.set_xlim(32.5, 35.5)

                ax.set_xlabel(original_variable_name, fontsize=20)
                ax.tick_params(axis='x', labelsize=16)
                ax.grid(True)

                if show_colorbar:
                    norm_for_cbar = Normalize(vmin=0, vmax=1)
                    scalar_mappable = ScalarMappable(cmap=cmap, norm=norm_for_cbar)
                    scalar_mappable.set_array([])
                    cbar = plt.colorbar(scalar_mappable, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
                    if color_mode == 'distance':
                        current_ticks = cbar.get_ticks()
                        new_tick_labels = [f"{1.0 - t:.1f}" for t in current_ticks]
                        cbar.set_ticks(current_ticks)
                        cbar.set_ticklabels(new_tick_labels)
                        cbar.set_label('Normalized Distance from Eddy Center (0=center, 1=edge)', fontsize=14)
                    elif color_mode == 'time':
                        cbar.set_label(f'Normalized Profile Sequence (Range: {min_profile_num} to {max_profile_num})', fontsize=14)

            # 整体属性
            if not profile_dates_for_title:
                plt.close(fig)
                continue

            date_start_platform = min(profile_dates_for_title)
            date_end_platform = max(profile_dates_for_title)

            axes[0].set_ylabel("Depth/m", fontsize=20)
            axes[0].tick_params(axis='y', labelsize=16)
            axes[0].invert_yaxis()

            fig.suptitle(f"{ds_names}{no}, Platform: {int(platform_id_val)}, {date_start_platform.date()}~{date_end_platform.date()}", fontsize=24, y=0.95)
            plt.tight_layout(rect=[0, 0, 1, 0.93])

            if save_fig:
                output_dir = "plot_vertical_profiles"
                os.makedirs(output_dir, exist_ok=True)
                plt.savefig(os.path.join(output_dir, f"{ds_names}{no}_Platform_{int(platform_id_val)}.png"), dpi=300, bbox_inches='tight')

            if show_fig:
                plt.show()
            plt.close(fig)
        return

    # =========================
    # 分支二：聚合（等价于原 plot_vertical_monthly）
    # =========================
    # 收集各剖面的基本信息
    profiles_to_plot = []
    try:
        argo_data_filtered['date_ts'] = pd.to_datetime(argo_data_filtered[['Year', 'Month', 'Day']])
    except (ValueError, TypeError) as e:
        print(f"Could not create timestamps for all rows due to invalid date components: {e}.")
        argo_data_filtered = argo_data_filtered.dropna(subset=['Year', 'Month', 'Day'])
        argo_data_filtered['date_ts'] = pd.to_datetime(argo_data_filtered[['Year', 'Month', 'Day']])

    for _, profile_rows in argo_data_filtered.groupby(['Platform_number', 'Profile_number']):
        if not profile_rows.empty:
            profiles_to_plot.append({
                'rows': profile_rows,
                'date': profile_rows.iloc[0]['date_ts'],
                'lon': profile_rows.iloc[0]['Longitude'],
                'lat': profile_rows.iloc[0]['Latitude']
            })

    if not profiles_to_plot:
        print(f"No data found for eddy {ds_names}{no} in months {month_required}.")
        return

    all_dates = [p['date'] for p in profiles_to_plot]
    min_time_for_norm, max_time_for_norm = (min(all_dates), max(all_dates)) if all_dates else (None, None)

    num_variables = len(variables)
    fig, axes = plt.subplots(1, num_variables, figsize=(10 * num_variables, 20), sharey=True)
    if num_variables == 1:
        axes = [axes]

    cmap = plt.cm.coolwarm

    for i, var_name in enumerate(variables):
        ax = axes[i]
        original_variable_name = var_name
        db_variable_name = 'Temperature' if var_name == 'Temp' else var_name

        if db_variable_name not in argo_data_filtered.columns:
            ax.text(0.5, 0.5, f"Variable '{db_variable_name}'\nnot found in data.",
                    ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_xlabel(original_variable_name, fontsize=20)
            ax.grid(True)
            continue

        for profile_info in profiles_to_plot:
            rows = profile_info['rows']
            rows_to_plot = rows.dropna(subset=[db_variable_name, "Depth"])
            if rows_to_plot.empty:
                continue

            if remove_outliers:
                qc_column_name = f"{db_variable_name}_Flag"
                if qc_column_name in rows_to_plot.columns:
                    good_qc_flags = ['1', '2', '5', '8', 1, 2, 5, 8]
                    bad_qc_mask = ~rows_to_plot[qc_column_name].isin(good_qc_flags)
                    rows_to_plot.loc[bad_qc_mask, db_variable_name] = np.nan
                if db_variable_name == 'DO':
                    bad_value_mask = rows_to_plot[db_variable_name] <= 1.0
                    rows_to_plot.loc[bad_value_mask, db_variable_name] = np.nan

            if rows_to_plot.empty:
                continue

            current_date = profile_info['date']
            color_value_normalized = 0.5

            if color_mode == 'distance':
                if wanted_track is not None and len(wanted_track) > 0 and 'lon' in profile_info and 'lat' in profile_info:
                    track_dates_converted = wanted_track['date'] if 'date' in wanted_track.columns else convert_date(wanted_track['time'])
                    idx_track_list = [j for j, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_date.date()]
                    if idx_track_list:
                        idx_track = idx_track_list[0]
                        center_lon = float(wanted_track.iloc[idx_track]['center_lon'])
                        center_lat = float(wanted_track.iloc[idx_track]['center_lat'])
                        radius = float(wanted_track.iloc[idx_track]['radius'])
                        if radius > 1e-6:
                            dist_m = adaptive_distance_m(
                                profile_info['lon'], profile_info['lat'],
                                center_lon, center_lat
                            )
                            distance = dist_m / radius
                            color_value_normalized = 1.0 - np.clip(distance, 0.0, 1.0)
            elif color_mode == 'time':
                if min_time_for_norm and max_time_for_norm and max_time_for_norm > min_time_for_norm:
                    total_delta = (max_time_for_norm - min_time_for_norm).total_seconds()
                    current_delta = (current_date - min_time_for_norm).total_seconds()
                    color_value_normalized = current_delta / total_delta if total_delta > 0 else 0.0
                else:
                    color_value_normalized = 0.0

            color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
            ax.plot(rows_to_plot[db_variable_name], rows_to_plot["Depth"], color=color, alpha=0.7)

        # 子图属性
        ax.set_ylim(-50, 2050)
        if db_variable_name == 'DO':
            ax.set_xlim(10, 350)
        elif db_variable_name == 'Temperature':
            ax.set_xlim(1, 32)
        elif db_variable_name == 'Salinity':
            ax.set_xlim(32.5, 35.5)

        ax.set_xlabel(original_variable_name, fontsize=20)
        ax.tick_params(axis='x', labelsize=16)
        ax.grid(True)

        if show_colorbar:
            norm = Normalize(vmin=0, vmax=1)
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
            if color_mode == 'distance':
                cbar.set_ticks([0, 0.5, 1])
                cbar.set_ticklabels(['1.0', '0.5', '0.0'])
                cbar.set_label('Normalized Distance (0=center, 1=edge)', fontsize=14)
            elif color_mode == 'time' and min_time_for_norm and max_time_for_norm:
                cbar.set_label(f'Normalized Time\n({min_time_for_norm.strftime("%Y-%m-%d")} to {max_time_for_norm.strftime("%Y-%m-%d")})', fontsize=12)

    # 整体属性
    axes[0].set_ylabel("Depth/m", fontsize=20)
    axes[0].tick_params(axis='y', labelsize=16)
    axes[0].invert_yaxis()

    month_str = "All" if month_required and len(month_required) > 6 else (", ".join(map(str, month_required)) if month_required else "All")
    date_range_str = f"{min_time_for_norm.date()}~{max_time_for_norm.date()}" if all_dates else "No date range"
    fig.suptitle(f"{ds_names}{no}, Months: {month_str}, Data: {date_range_str}", fontsize=24, y=0.95)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    if save_fig:
        output_dir = "plot_vertical_monthly_aggregated"
        os.makedirs(output_dir, exist_ok=True)
        month_suffix = "all" if not month_required or (month_required and len(month_required) > 6) else "_".join(map(str, month_required))
        filename = f"{ds_names}{no}_months_{month_suffix}_aggregated.png"
        plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, filename)}")

    if show_fig:
        plt.show()
    plt.close(fig)
        
def plot_relative_position(
    DS: list,
    no: int,
    show_fig: bool = False,
    save_fig: bool = False,
    color_mode: str = 'distance',
    show_colorbar: bool = False,
    aggregated: bool = False,
    argo_required: list | None = None,
    year_required: list | None = None,
    month_required: list | None = None,
    day_required: list | None = None
):
    '''
    根据涡旋轨迹和浮标数据，绘制浮标在单位圆涡旋中的相对位置分布图。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        show_fig (bool): 是否显示图片，默认 False。
        save_fig (bool): 是否保存图片，默认 False。
        color_mode (str): 'distance' 或 'time'，默认 'distance'。
        show_colorbar (bool): 是否显示颜色条，默认 False。
        aggregated (bool): 是否聚合所有平台于一图，默认 False。
        argo_required (list | None): 平台筛选；None 表示不过滤；传入平台编号列表时仅保留指定平台。
        year_required (list | None): 年份筛选；None 表示不过滤。
        month_required (list | None): 月份筛选；聚合模式 None 表示使用所有可用月份；逐平台模式 None 表示不过滤。
        day_required (list | None): 日期筛选（按日1-31）；None 表示不过滤。

        功能:
        - 计算剖面代表点在单位圆涡旋中的相对位置，并以散点标注（中心×，单位圆圈）。
        - 颜色模式: 'distance'（距离中心归一化，0=中心，1=边缘）或 'time'（按时间顺序归一化）。
        - 横纵轴刻度包含真实经纬度与相对坐标。可选显示颜色条；支持图片保存与显示。
        - 支持按照月份（month_required）与平台编号（argo_required）进行筛选。

    模式差异:
        - aggregated=False（逐平台）：
            • 对每个浮标平台分别出图，点内数字代表该平台内部的剖面时序，从1开始递增。
        - aggregated=True（聚合）：
            • 所有平台的代表点聚合到同一张图，点内数字代表相对于所选月份范围起始日的累积天数。例如，若数据从7月29日开始，则该日所有剖面的数字为29，7月30日为30，8月1日则为32。
    '''
    wanted_track = find_track(DS, no)
    argo_data_filtered = filtered_float_data(DS, no, track=wanted_track)

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE: Could not automatically determine the dataset name.")

    if argo_data_filtered.empty:
        print(f"No Argo data found for eddy {ds_names}{no} to plot relative positions.")
        return

    # 平台筛选（与 plot_vertical 一致）
    if argo_required is not None:
        try:
            argo_required_set = set(int(x) for x in argo_required)
        except Exception:
            argo_required_set = set(argo_required)
        argo_data_filtered = argo_data_filtered[argo_data_filtered['Platform_number'].isin(argo_required_set)]
        if argo_data_filtered.empty:
            print(f"No Argo data left after filtering by platforms {sorted(list(argo_required_set))}.")
            return

    # 年月日筛选（与 plot_vertical 一致）
    argo_data_filtered['Year'] = pd.to_numeric(argo_data_filtered['Year'], errors='coerce')
    argo_data_filtered.dropna(subset=['Year'], inplace=True)
    argo_data_filtered['Year'] = argo_data_filtered['Year'].astype(int)

    argo_data_filtered['Month'] = pd.to_numeric(argo_data_filtered['Month'], errors='coerce')
    argo_data_filtered.dropna(subset=['Month'], inplace=True)
    argo_data_filtered['Month'] = argo_data_filtered['Month'].astype(int)

    argo_data_filtered['Day'] = pd.to_numeric(argo_data_filtered['Day'], errors='coerce')
    argo_data_filtered.dropna(subset=['Day'], inplace=True)
    argo_data_filtered['Day'] = argo_data_filtered['Day'].astype(int)

    if aggregated:
        if not month_required:
            all_months = sorted(argo_data_filtered['Month'].unique().tolist())
            if not all_months:
                print(f"No valid months found in data for eddy {ds_names}{no}.")
                return
            print(f"No months specified for relative position plot. Defaulting to all available months: {all_months}")
            month_required = all_months
        # 按提供的年/月/日过滤（月份在聚合模式下为必有列表）
        if month_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Month'].isin(month_required)].copy()
        if year_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Year'].isin(year_required)].copy()
        if day_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Day'].isin(day_required)].copy()
        if argo_data_filtered.empty:
            print(f"No data found for eddy {ds_names}{no} after applying filters: year={year_required}, month={month_required}, day={day_required}.")
            return
    else:
        if year_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Year'].isin(year_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in years {year_required} (per-platform mode).")
                return
        if month_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Month'].isin(month_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in months {month_required} (per-platform mode).")
                return
        if day_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Day'].isin(day_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in days {day_required} (per-platform mode).")
                return

    # =========================
    # 分支一：逐平台模式（aggregated=False）
    # =========================
    if not aggregated:
        for platform_id_val, platform_data in argo_data_filtered.groupby("Platform_number"):
        # 获取每个剖面的第一行数据作为代表点
        # 原文: needed_data = platform.groupby("Profile_number").apply(lambda group: group.iloc[0])
        #       needed_data.index.name = None
        # 使用 .first() 更简洁，并且可以直接用 Profile_number 作为索引或后续处理
            profile_first_rows = platform_data.groupby("Profile_number").first().reset_index() # reset_index 使 Profile_number 成为列

            if profile_first_rows.empty:
            # print(f"No profile data for platform {platform_id_val}.")
                continue

            # 准备每个剖面点对应的涡旋轨迹数据
            points_for_this_platform = []
            track_info_for_this_platform = [] # 用于计算该平台下的平均涡旋参数
        
            track_dates_converted = wanted_track['date'] if 'date' in wanted_track.columns else convert_date(wanted_track['time'])

            for i, p_row in profile_first_rows.iterrows(): # i 将用作顺序编号
                try:
                    current_date_profile = pd.Timestamp(year=int(p_row['Year']),
                                                        month=int(p_row['Month']),
                                                        day=int(p_row['Day']))
                except (ValueError, TypeError):
                    # print(f"Skipping profile {p_row.get('Profile_number')} for platform {platform_id_val} due to invalid date.")
                    continue

                center_lon, center_lat, radius = None, None, None
                if wanted_track is not None and len(wanted_track) > 0:
                    matches = [k for k, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_date_profile.date()]
                    if matches:
                        idx_track = matches[0]
                        center_lon = float(wanted_track.iloc[idx_track]['center_lon'])
                        center_lat = float(wanted_track.iloc[idx_track]['center_lat'])
                        radius = float(wanted_track.iloc[idx_track]['radius'])
            
                if center_lon is not None and radius is not None and radius > 1e-6:
                    if 'Longitude' not in p_row or 'Latitude' not in p_row:
                        # print(f"Skipping point on {current_date_profile.date()} due to missing Longitude/Latitude.")
                        continue
                    scale = approximate_degree_length(center_lat)
                    dlon_deg = _minimal_lon_diff_deg(p_row['Longitude'], center_lon)
                    dx_m = dlon_deg * scale['meters_per_degree_lon']
                    dy_m = (p_row['Latitude'] - center_lat) * scale['meters_per_degree_lat']
                    rel_x = dx_m / radius
                    rel_y = dy_m / radius
                    points_for_this_platform.append({
                        'rel_x': rel_x, 'rel_y': rel_y, 'date': current_date_profile,
                        'profile_num_original_idx': p_row['Profile_number'], # 用于 'time' mode
                        'sequence_label': i + 1 # 平台内时序标签
                    })
                    track_info_for_this_platform.append([center_lon, center_lat, radius])
            # else: 涡旋数据不匹配或半径过小，则跳过该点
        
            if not points_for_this_platform:
            # print(f"No valid points with track data for platform {platform_id_val}.")
                continue

            fig, ax = plt.subplots(figsize=(30, 20)) # 原为 (30,20) 但相对位置图常用正方形，可考虑 (20,20) 或 (15,15)
            cmap = plt.cm.coolwarm

            # 准备 'time' 模式颜色归一化 (基于 Profile_number)
            min_prof_num_platform = min(p['profile_num_original_idx'] for p in points_for_this_platform)
            max_prof_num_platform = max(p['profile_num_original_idx'] for p in points_for_this_platform)

            for point_data in points_for_this_platform:
                rel_x = point_data['rel_x']
                rel_y = point_data['rel_y']
            
                color_value_normalized = 0.5
                if color_mode == 'distance':
                    distance = np.sqrt(rel_x**2 + rel_y**2)
                    color_value_normalized = 1.0 - np.clip(distance, 0.0, 1.0)
                elif color_mode == 'time':
                    if max_prof_num_platform > min_prof_num_platform:
                        color_value_normalized = (point_data['profile_num_original_idx'] - min_prof_num_platform) / \
                                                 (max_prof_num_platform - min_prof_num_platform)
                    else:
                        color_value_normalized = 0.0
            
                color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
                ax.scatter(rel_x, rel_y, color=color, s=300, zorder=5)
                ax.text(rel_x, rel_y, str(point_data['sequence_label']), weight='bold', fontsize=9, color='black', ha='center', va='center', zorder=6)
            
            ax.plot(0, 0, marker='x', color='black', markersize=16, markeredgewidth=3, label='Eddy Center (Relative)', zorder=3)
            circle = plt.Circle((0, 0), 1, color='gray', fill=False, linestyle='--', linewidth=2, label='Unit Eddy Boundary', zorder=2)
            ax.add_patch(circle)
            ax.set_aspect('equal')

            # 标题和坐标轴标签
            platform_dates = [p['date'] for p in points_for_this_platform]
            date_start_platform = min(platform_dates)
            date_end_platform = max(platform_dates)
            ax.set_title(f"{ds_names}{no}, Platform: {int(platform_id_val)}, {date_start_platform.date()}~{date_end_platform.date()}, Points: {len(points_for_this_platform)}", fontsize=20)
            ax.set_xlabel('Relative X (Eddy Radii)', fontsize=20)
            ax.set_ylabel('Relative Y (Eddy Radii)', fontsize=20)
            ax.tick_params(axis='both', which='major', labelsize=16) # 使用tick_params统一设置

            # 设置坐标轴刻度 (与聚合版一致)
            if track_info_for_this_platform:
                mean_center_lon = np.mean([info[0] for info in track_info_for_this_platform])
                mean_center_lat = np.mean([info[1] for info in track_info_for_this_platform])
                mean_radius = np.mean([info[2] for info in track_info_for_this_platform])

                if not np.isnan(mean_center_lon) and not np.isnan(mean_center_lat) and not np.isnan(mean_radius) and mean_radius > 1e-6:
                    scale_mean = approximate_degree_length(mean_center_lat)
                    mean_deg_x = mean_radius / scale_mean['meters_per_degree_lon']
                    mean_deg_y = mean_radius / scale_mean['meters_per_degree_lat']
                
                    tick_locs = [-1, -0.5, 0, 0.5, 1] # 使用更详细的刻度
                
                    x_tick_labels = [f"{(mean_center_lon + tick_loc * mean_deg_x):.2f}°\n({tick_loc:.1f})" for tick_loc in tick_locs]
                    ax.set_xticks(tick_locs)
                    ax.set_xticklabels(x_tick_labels)

                    y_tick_labels = [f"{(mean_center_lat + tick_loc * mean_deg_y):.2f}°\n({tick_loc:.1f})" for tick_loc in tick_locs]
                    ax.set_yticks(tick_locs)
                    ax.set_yticklabels(y_tick_labels)
        
            ax.set_xlim([-1.25, 1.25])
            ax.set_ylim([-1.25, 1.25])
            # ax.legend(fontsize=14) # 原代码legend注释掉了，保持一致

            # 添加颜色条
            if show_colorbar:
                norm_for_cbar = Normalize(vmin=0, vmax=1)
                scalar_mappable = ScalarMappable(cmap=cmap, norm=norm_for_cbar)
                scalar_mappable.set_array([])
                cbar = plt.colorbar(scalar_mappable, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)

                if color_mode == 'distance':
                    current_ticks = cbar.get_ticks()
                    new_tick_labels = [f"{1.0 - t:.1f}" for t in current_ticks]
                    cbar.set_ticks(current_ticks)
                    cbar.set_ticklabels(new_tick_labels)
                    cbar.set_label('Normalized Distance from Eddy Center (0=center, 1=edge)', fontsize=14)
                elif color_mode == 'time':
                    cbar.set_label(f'Normalized Profile Sequence (Platform Range: {min_prof_num_platform} to {max_prof_num_platform})', fontsize=14)
        
            if save_fig:
                output_dir = "plot_relative_position"
                os.makedirs(output_dir, exist_ok=True)
                plt.savefig(os.path.join(output_dir, f"{ds_names}{no}RP{int(platform_id_val)}.png"), dpi=300, bbox_inches='tight')
        
            if show_fig:
                plt.show()
            plt.close(fig)
        return

    # =========================
    # 分支二：聚合模式（aggregated=True）
    # =========================
    # 1) 收集聚合绘图所需点
    points_to_process = []  # {date, data_row}
    monthly_filtered_data = argo_data_filtered.copy()
    if month_required:
        monthly_filtered_data = monthly_filtered_data[monthly_filtered_data['Month'].isin(month_required)].copy()

    if not monthly_filtered_data.empty:
        # 获取每个剖面的第一行数据作为代表点
        profile_first_rows = monthly_filtered_data.groupby(["Platform_number", "Profile_number"]).first().reset_index()
        for _, p_row in profile_first_rows.iterrows():
            try:
                current_date_profile = pd.Timestamp(year=int(p_row['Year']), month=int(p_row['Month']), day=int(p_row['Day']))
            except (ValueError, TypeError):
                continue
            if not month_required or current_date_profile.month in month_required:
                points_to_process.append({'date': current_date_profile, 'data_row': p_row})

    if not points_to_process:
        print(f"No data found for eddy {ds_names}{no} in months {month_required}.")
        return

    # 2) 计算日期标签的参考起点（按所选月份的最小月的一号）
    min_plot_date_overall = min(p['date'] for p in points_to_process)
    ref_month = min(month_required) if month_required else min_plot_date_overall.month
    reference_start_date_for_labels = pd.Timestamp(year=min_plot_date_overall.year, month=ref_month, day=1)

    points_to_plot = []
    all_track_info_for_overall_mean = []
    all_profile_dates_for_title = []
    all_profile_timestamps_for_time_mode = []

    track_dates_converted = wanted_track['date'] if 'date' in wanted_track.columns else convert_date(wanted_track['time'])

    for point_info in points_to_process:
        current_date = point_info['date']
        p_row = point_info['data_row']

        day_label = (current_date - reference_start_date_for_labels).days + 1

        center_lon, center_lat, radius = None, None, None
        if wanted_track is not None and len(wanted_track) > 0:
            matches = [i for i, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_date.date()]
            if matches:
                idx_track = matches[0]
                center_lon = float(wanted_track.iloc[idx_track]['center_lon'])
                center_lat = float(wanted_track.iloc[idx_track]['center_lat'])
                radius = float(wanted_track.iloc[idx_track]['radius'])

        if center_lon is not None and radius is not None and radius > 1e-6:
            if 'Longitude' not in p_row or 'Latitude' not in p_row:
                continue

            scale = approximate_degree_length(center_lat)
            dlon_deg = _minimal_lon_diff_deg(p_row['Longitude'], center_lon)
            dx_m = dlon_deg * scale['meters_per_degree_lon']
            dy_m = (p_row['Latitude'] - center_lat) * scale['meters_per_degree_lat']
            rel_x = dx_m / radius
            rel_y = dy_m / radius

            points_to_plot.append({'rel_x': rel_x, 'rel_y': rel_y, 'date': current_date, 'day_label': day_label})
            all_track_info_for_overall_mean.append([center_lon, center_lat, radius])
            all_profile_dates_for_title.append(current_date)
            if color_mode == 'time':
                all_profile_timestamps_for_time_mode.append(current_date)

    if not points_to_plot:
        print(f"No valid points with track data found for eddy {ds_names}{no} in months {month_required}.")
        return

    # 3) 确定时间归一化范围（仅用于 color_mode='time'）
    min_time_for_norm, max_time_for_norm = (None, None)
    if color_mode == 'time' and all_profile_timestamps_for_time_mode:
        min_time_for_norm = min(all_profile_timestamps_for_time_mode)
        max_time_for_norm = max(all_profile_timestamps_for_time_mode)
        if min_time_for_norm == max_time_for_norm and len(all_profile_timestamps_for_time_mode) > 1:
            max_time_for_norm = min_time_for_norm + pd.Timedelta(days=1)

    date_start_overall = min(all_profile_dates_for_title) if all_profile_dates_for_title else None
    date_end_overall = max(all_profile_dates_for_title) if all_profile_dates_for_title else None

    # 4) 绘图
    fig, ax = plt.subplots(figsize=(30, 20))
    cmap = plt.cm.coolwarm

    for point in points_to_plot:
        rel_x, rel_y = point['rel_x'], point['rel_y']
        current_date, day_label = point['date'], point['day_label']

        color_value_normalized = 0.5
        if color_mode == 'distance':
            distance_from_center = np.sqrt(rel_x**2 + rel_y**2)
            color_value_normalized = 1.0 - np.clip(distance_from_center, 0.0, 1.0)
        elif color_mode == 'time':
            if min_time_for_norm and max_time_for_norm and min_time_for_norm < max_time_for_norm:
                total_delta = (max_time_for_norm - min_time_for_norm).total_seconds()
                current_delta = (current_date - min_time_for_norm).total_seconds()
                color_value_normalized = (current_delta / total_delta) if total_delta > 0 else 0.0
            elif min_time_for_norm and max_time_for_norm and min_time_for_norm == max_time_for_norm:
                color_value_normalized = 0.0

        color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
        ax.scatter(rel_x, rel_y, color=color, s=300, zorder=5)
        ax.text(rel_x, rel_y, str(day_label), weight='bold', fontsize=9, color='black', ha='center', va='center', zorder=6)

    # 中心与单位圆
    ax.plot(0, 0, marker='x', color='black', markersize=16, markeredgewidth=3, label='Eddy Center (Relative)', zorder=3)
    circle = plt.Circle((0, 0), 1, color='gray', fill=False, linestyle='--', linewidth=2, label='Unit Eddy Boundary', zorder=2)
    ax.add_patch(circle)
    ax.set_aspect('equal')

    # 标题与坐标轴
    month_str = "All" if (month_required and len(month_required) > 6) else ", ".join(map(str, month_required)) if month_required else "All"
    title_str = f"{ds_names}{no}, Months: {month_str}, Relative Positions"
    if date_start_overall and date_end_overall:
        title_str += f"\nData: {date_start_overall.date()}~{date_end_overall.date()}, Total Points: {len(points_to_plot)}"
    ax.set_title(title_str, fontsize=20)
    ax.set_xlabel('Relative X (Eddy Radii)', fontsize=20)
    ax.set_ylabel('Relative Y (Eddy Radii)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # 地理刻度（均值）
    if all_track_info_for_overall_mean:
        mean_center_lon = np.mean([info[0] for info in all_track_info_for_overall_mean])
        mean_center_lat = np.mean([info[1] for info in all_track_info_for_overall_mean])
        mean_radius = np.mean([info[2] for info in all_track_info_for_overall_mean])
        if not np.isnan(mean_center_lon) and not np.isnan(mean_center_lat) and not np.isnan(mean_radius) and mean_radius > 1e-6:
            scale_mean = approximate_degree_length(mean_center_lat)
            mean_deg_x = mean_radius / scale_mean['meters_per_degree_lon']
            mean_deg_y = mean_radius / scale_mean['meters_per_degree_lat']
            tick_locs = [-1, -0.5, 0, 0.5, 1]
            x_tick_labels = [f"{(mean_center_lon + t * mean_deg_x):.2f}°\n({t})" for t in tick_locs]
            y_tick_labels = [f"{(mean_center_lat + t * mean_deg_y):.2f}°\n({t})" for t in tick_locs]
            ax.set_xticks(tick_locs)
            ax.set_xticklabels(x_tick_labels)
            ax.set_yticks(tick_locs)
            ax.set_yticklabels(y_tick_labels)

    ax.set_xlim([-1.25, 1.25])
    ax.set_ylim([-1.25, 1.25])

    # 颜色条
    if show_colorbar:
        norm_for_cbar = Normalize(vmin=0, vmax=1)
        scalar_mappable = ScalarMappable(cmap=cmap, norm=norm_for_cbar)
        scalar_mappable.set_array([])
        cbar = plt.colorbar(scalar_mappable, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
        if color_mode == 'distance':
            current_ticks = cbar.get_ticks()
            new_tick_labels = [f"{1.0 - t:.1f}" for t in current_ticks]
            cbar.set_ticks(current_ticks)
            cbar.set_ticklabels(new_tick_labels)
            cbar.set_label('Normalized Distance from Eddy Center (0=center, 1=edge)', fontsize=14)
        elif color_mode == 'time':
            if min_time_for_norm and max_time_for_norm:
                cbar.set_label(f'Normalized Time ({min_time_for_norm.strftime("%Y-%m-%d")} to {max_time_for_norm.strftime("%Y-%m-%d")})', fontsize=14)
            else:
                cbar.set_label('Normalized Time', fontsize=14)

    # 保存/显示
    if save_fig:
        output_dir = "plot_relative_position_monthly_aggregated"
        os.makedirs(output_dir, exist_ok=True)
        month_suffix = "all" if not month_required or (month_required and len(month_required) > 6) else "_".join(map(str, month_required))
        base_filename = f"{ds_names}{no}_RP_months_{month_suffix}_aggregated.png"
        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

    if show_fig:
        plt.show()
    plt.close(fig)

def plot_vertical_monthly(
    DS: list,
    no: int,
    month_required: list = None,
    show_fig: bool = False,
    save_fig: bool = False,
    color_mode: str = 'distance',
    variables: list = ['DO', 'Temp', 'Salinity'],
    show_colorbar: bool = False,
    remove_outliers: bool = True,
    argo_required: list | None = None
):
    '''
    兼容封装：保持旧接口，内部调用统一的 plot_vertical。
    '''
    return plot_vertical(
        DS=DS,
        no=no,
        show_fig=show_fig,
        save_fig=save_fig,
        color_mode=color_mode,
        variables=variables,
        show_colorbar=show_colorbar,
        remove_outliers=remove_outliers,
        aggregated=True,
        month_required=month_required,
        argo_required=argo_required
    )
    
def plot_relative_position_monthly(
    DS: list,
    no: int,
    show_fig: bool = False,
    save_fig: bool = False,
    color_mode: str = 'distance',
    show_colorbar: bool = False,
    month_required: list | None = None,
    argo_required: list | None = None
):
    '''
    兼容封装：保持旧接口，内部调用统一的 plot_relative_position（aggregated=True）。
    '''
    return plot_relative_position(
        DS=DS,
        no=no,
        show_fig=show_fig,
        save_fig=save_fig,
        color_mode=color_mode,
        show_colorbar=show_colorbar,
        aggregated=True,
        month_required=month_required,
        argo_required=argo_required,
    )

def get_glorys_filepath(date)-> str:
    '''
    根据给定日期返回对应的GLORYS NetCDF文件路径。

    构造目标日期的文件夹路径和文件名模式，查找匹配的文件。如果找到且仅有一个匹配文件，则返回其完整路径；
    如果没有找到或找到多个匹配文件，则抛出异常。

    参数:
        date: 由convert_date得出的需要查找的日期。

    返回:
        str: 匹配的GLORYS NetCDF文件的完整路径。

    异常:
        RuntimeError: 如果找到多个匹配文件。
        FileNotFoundError: 如果没有找到匹配文件。
    '''
    nc_path = os.path.join(Glorys_path, date.strftime("%Y"), date.strftime("%m"), 'mercatorglorys12v1_gl12_mean_')
    nc_path += date.strftime("%Y%m%d")
    matches = glob.glob(nc_path + '_R*')
    if matches:
        if len(matches) > 1:
            raise RuntimeError(f"Multiple matching files found for {date.strftime('%Y-%m-%d')}: {matches}")
        nc_path = matches[0]
        #print("Matching file:", nc_path)
        return nc_path
    else:
        raise FileNotFoundError(f"No matching file found for {date.strftime('%Y-%m-%d')}.")
        
def print_glorys_variable(nc_path: str):
    '''
    打印指定路径对应NetCDF文件中的所有变量名称、标准名称、维度和形状信息。
    '''
    nc_file = Dataset(nc_path, 'r')
    
    with Dataset(nc_path, 'r') as nc_file:
        for var in nc_file.variables:
            print(var, nc_file.variables[var].standard_name)
            print(nc_file.variables[var].dimensions, nc_file.variables[var].shape)

def find_track_glorys_filepath(DS: list, no: int) -> dict:
    '''
    根据涡旋编号在GLORYS数据集中查找对应的轨迹文件路径。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。

    返回:
        dict: 涡旋轨迹数据文件路径的字典，格式为 {date: glorys_filepath}。
              如果未找到对应的轨迹数据或文件路径，则返回空字典。
    '''
    wanted_track = find_track(DS, no)
    if wanted_track is None or len(wanted_track) == 0:
        print(f"未找到涡旋 {no} 的轨迹数据。")
        return {}
    
    glorys_filepaths_dict = {}
    for _, track_point in wanted_track.iterrows():
        try:
            date = convert_date(track_point['time'])  
            glorys_filepath = get_glorys_filepath(date)
            glorys_filepaths_dict[date] = glorys_filepath          
        except (RuntimeError, FileNotFoundError) as e:
            print(f"为涡旋 {no} 在日期 {track_point.get('time', '未知')} (转换后: {date if 'date' in locals() else '未知'}) 查找 GLORYS 文件时出错: {e}")
        except IndexError:
            print(f"涡旋 {no} 的轨迹点数据格式不正确: {track_point}")
        except Exception as e: # 捕获其他可能的 convert_date 或 get_glorys_filepath 异常
            print(f"处理涡旋 {no} 在日期 {track_point.get('time', '未知')} 时发生未知错误: {e}")

    if not glorys_filepaths_dict:
        print(f"未找到涡旋 {no} 的 GLORYS 文件。")
        return {}
        
    return glorys_filepaths_dict

class LineDrawer:
    """
    一个用于在matplotlib图像上通过点击两点来交互式绘制直线的辅助类。
    新增了撤回功能：按 'z' 键可以撤回上一个点击的点或已绘制的直线。
    """
    def __init__(self, ax):
        self.ax = ax
        self.points = []  # 用于存储用户点击的坐标
        self.lines = []   # 用于存储已绘制的直线对象
        self.markers = [] # 用于存储用户点击时产生的临时视觉标记
        # 将键盘事件与 onkey 方法连接
        self.cid_key = self.ax.figure.canvas.mpl_connect('key_press_event', self.onkey)
        
        print("\n--- 交互模式已激活 (稳定版) ---")
        print("请在图像上点击两个点来定义一条新的直线。")
        print("提示：按 'z' 键可以撤回上一个点或刚绘制的直线。")

    def onclick(self, event):
        """鼠标点击事件的处理函数。"""
        # 确保点击事件发生在指定的坐标轴内
        if event.inaxes != self.ax: return
        
        # 记录点击的坐标点
        self.points.append((event.xdata, event.ydata))
        # 在点击位置添加一个临时的"+"标记，给予用户视觉反馈
        marker = self.ax.plot(event.xdata, event.ydata, 'm+', markersize=12, markeredgewidth=2)
        self.markers.extend(marker)
        # 触发一次完整的重绘，确保标记点被永久画上
        self.ax.figure.canvas.draw()

        # 当记录的点达到两个时，开始绘制直线
        if len(self.points) == 2:
            self.draw_line()

    def onkey(self, event):
        """键盘按键事件的处理函数，用于实现撤回功能。"""
        if event.key == 'z':
            # 如果当前正在定义一条线（已点击了一个点）
            if self.markers:
                # 移除记录的点和对应的标记
                self.points.pop()
                last_marker = self.markers.pop()
                last_marker.remove()
                self.ax.figure.canvas.draw()
                print("上一个点已撤回。")
            # 如果要撤回一条已完成的线
            elif self.lines:
                # 移除上一条绘制的直线
                last_line = self.lines.pop()
                last_line.remove()
                self.ax.legend() # 更新图例
                self.ax.figure.canvas.draw()
                print("上一条线已撤回。")
            else:
                print("没有可撤回的操作。")

    def draw_line(self):
        """根据记录的两个点计算并绘制直线。"""
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]

        # 计算直线方程 y = kx + b
        if abs(x2 - x1) < 1e-6: # 处理斜率不存在的垂直线情况
            k, b, eq_text = np.inf, np.nan, f"x = {x1:.2f}"
        else:
            k = (y2 - y1) / (x2 - x1)
            b = y1 - k * x1
            eq_text = f"y = {k:.4f}x + {b:.4f}"

        # 根据坐标轴的当前范围，计算直线的端点并绘制
        x_vals = np.array(self.ax.get_xlim())
        if np.isinf(k): # 绘制垂直线
             y_vals = np.array(self.ax.get_ylim())
             line = self.ax.plot([x1, x1], y_vals, 'purple', linestyle='--', linewidth=2, label=f'Interactive: {eq_text}')
        else: # 绘制普通斜率的直线
            y_vals = k * x_vals + b
            line = self.ax.plot(x_vals, y_vals, 'purple', linestyle='--', linewidth=2, label=f'Interactive: {eq_text}')
        
        self.lines.extend(line)
        
        # 清理本次操作的临时标记点和坐标
        for marker in self.markers:
            marker.remove()
        self.markers.clear()
        self.points.clear()

        # 更新图例并重绘图像
        self.ax.legend()
        self.ax.figure.canvas.draw()
        print("\n准备就绪，可继续点击绘制下一条直线。")

def plot_track_area_horizontal_glorys(DS: list, no: int, needed_idx: int, variable: str = 'vorticity',
                                   show_fig: bool = False, save_fig: bool = False, deep_argo: bool = False,
                                   k: float | list[float] | None = None, b: float | list[float] | None = None, 
                                   needed_depth: float | int = 0, inline_mode: bool = True):
    '''
    绘制指定涡旋在特定时刻的表层物理场快照及相关的Argo浮标数据。

    该函数支持两种模式：
    1. inline_mode=True (默认): 适用于静态图表生成和高分辨率保存，行为与原始版本完全一致。
    2. inline_mode=False: 适用于在Jupyter Notebook中使用 %matplotlib widget 进行交互式分析。

    参数:
        DS (list): 包含所有涡旋轨迹信息的数据集。
        no (int): 需要绘制的涡旋的唯一编号。
        needed_idx (int): 涡旋轨迹的时间点索引，用于确定绘图的具体日期。
        variable (str): 作为背景场绘制的GLORYS物理变量。默认为 'vorticity'。
        show_fig (bool): 是否在运行时显示生成的图像。默认为 False。
        save_fig (bool): 是否将生成的图像保存为文件。默认为 False。
        deep_argo (bool): 是否使用深层Argo数据模式。默认为 False。
        k (float | list[float], optional): 直线方程 y=kx+b 的斜率或斜率列表。默认为 None。
        b (float | list[float], optional): 直线方程 y=kx+b 的截距或截距列表。默认为 None。
        needed_depth (float | int): 需要绘制的GLORYS数据深度，默认为0（表层）。
        inline_mode (bool): 是否为静态内联模式。默认为True。设为False以启用交互式widget模式的优化。
    '''
    # 根据模式定义一套协调的尺寸参数
    if inline_mode:
        # 适用于高分辨率保存的静态模式尺寸，与原始版本完全一致
        figsize = (30, 25)
        title_fs, label_fs, tick_fs, legend_fs = 20, 20, 16, 18
        cbar_label_fs, cbar_tick_fs = 20, 14
        argo_text_fs = 7
        track_lw, contour_lw, circle_lw, line_lw = 1.0, 1.0, 1.0, 2.0
        cbar_pad = 0.046 # 原始cbar间距
    else:
        # 适用于交互式widget的屏幕友好尺寸
        figsize = (12, 10) 
        title_fs, label_fs, tick_fs, legend_fs = 16, 14, 12, 10
        cbar_label_fs, cbar_tick_fs = 12, 10
        argo_text_fs = 6
        track_lw, contour_lw, circle_lw, line_lw = 1.0, 1.0, 1.0, 2.0
        cbar_pad = 0.12 # 交互模式下增大cbar间距

    # 加载轨迹（DataFrame）并准备常用列
    wanted_track = find_track(DS, no)
    dates = wanted_track['date'] if 'date' in wanted_track.columns else convert_date(wanted_track['time'])
    center_lon_arr = wanted_track['center_lon'].to_numpy()
    center_lat_arr = wanted_track['center_lat'].to_numpy()
    radius_arr = wanted_track['radius'].to_numpy()
    # 当前时刻（needed_idx）的轮廓坐标
    curr_contour_lon = wanted_track.iloc[needed_idx]['contour_lon']
    curr_contour_lat = wanted_track.iloc[needed_idx]['contour_lat']

    # 获取Argo浮标数据（复用已读取轨迹，避免重复 I/O）
    argo_data_filtered = filtered_float_data(DS, no, track=wanted_track)
    argo_data_filtered = argo_data_filtered[pd.to_datetime(argo_data_filtered[['Year', 'Month', 'Day']])==dates[needed_idx]]
    
    if deep_argo:
        filtered_by_depth = argo_data_filtered[argo_data_filtered['Depth'] >= 500].copy()
        if filtered_by_depth.empty:
            print("Warning: No data found with Depth >= 500.")
            needed_data = pd.DataFrame(columns=argo_data_filtered.columns)
        else:
            idx_max_do = filtered_by_depth.groupby('Profile_number')['DO'].idxmax()
            needed_data = filtered_by_depth.loc[idx_max_do]
            needed_data.index.name = None
    else:
        if argo_data_filtered.empty:
            needed_data = pd.DataFrame(columns=argo_data_filtered.columns)
        else:
            needed_data = argo_data_filtered.groupby('Profile_number').apply(lambda group: group.iloc[0])
            needed_data.index.name = None

    # 获取区域边界（基于当前时刻轮廓）
    contour_lon_filtered = np.ma.masked_equal(curr_contour_lon, 180.0)
    contour_lat_filtered = np.ma.masked_equal(curr_contour_lat, 0.0)
    glorys_lon_min = np.min(contour_lon_filtered) - 0.5
    glorys_lon_max = np.max(contour_lon_filtered) + 0.5
    glorys_lat_min = np.min(contour_lat_filtered) - 0.5
    glorys_lat_max = np.max(contour_lat_filtered) + 0.5

    #获取背景场数据
    if variable == 'vorticity':
        glorys_lon_filtered, glorys_lat_filtered, glorys_depth_filtered, glorys_variables_filtered = get_track_area_glorys(DS, no, needed_idx, variables=['u', 'v'], depth=needed_depth)
        zeta, f = calculate_vorticity(glorys_lon_filtered, glorys_lat_filtered, glorys_variables_filtered['u'], glorys_variables_filtered['v'])
        glorys_variable_filtered = zeta/f
    else:
        glorys_lon_filtered, glorys_lat_filtered, glorys_depth_filtered, glorys_variables_filtered = get_track_area_glorys(DS, no, needed_idx, variables=[variable], depth=needed_depth)
        glorys_variable_filtered = glorys_variables_filtered[variable]

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE")

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    if ds_names == 'ACS' or ds_names == 'ACL':
        colors =colors[1]
    else:
        colors =colors[0]
    world = _load_world_geodataframe()

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(f'Track {ds_names}{no} at {glorys_depth_filtered[0]:.2f}m, {dates[needed_idx].strftime("%Y-%m-%d")}', fontsize=title_fs)
    ax.set_xlabel('Longitude', fontsize=label_fs)
    ax.set_ylabel('Latitude', fontsize=label_fs)
    world.plot(color='green', ax=ax)

    ax.tick_params(axis='both', which='major', labelsize=tick_fs)

    # 绘制涡旋轨迹
    ax.plot(center_lon_arr, center_lat_arr, color=colors, linewidth=track_lw, label='Center Track')
    ax.plot(center_lon_arr[0], center_lat_arr[0], marker='o', color=colors, markersize=10)
    ax.plot(center_lon_arr[-1], center_lat_arr[-1], marker='x', color=colors, markersize=10)

    # 绘制背景场
    pc = ax.pcolormesh(glorys_lon_filtered, glorys_lat_filtered, glorys_variable_filtered, cmap='seismic', shading='auto', alpha=0.5)
    cbar = plt.colorbar(pc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.046)
    if variable == 'vorticity':
        cbar.set_label(r'$\zeta/f$', fontsize=cbar_label_fs)
        pc.set_clim(-0.7,0.7)
    elif variable == 'thetao':
        cbar.set_label('Temperature (°C)', fontsize=cbar_label_fs)
    elif variable == 'so':
        cbar.set_label('Salinity (psu)', fontsize=cbar_label_fs)
    elif variable == 'u':
        cbar.set_label('Zonal Velocity (m/s)', fontsize=cbar_label_fs)
    elif variable == 'v':
        cbar.set_label('Meridional Velocity (m/s)', fontsize=cbar_label_fs)
    elif variable == 'ssh':
        cbar.set_label('Sea Surface Height (m)', fontsize=cbar_label_fs)
    else:
        cbar.set_label(variable, fontsize=cbar_label_fs)
    cbar.ax.tick_params(labelsize=cbar_tick_fs)

    # 绘制Argo浮标数据
    if not needed_data.empty:
        if deep_argo:
            sc = ax.scatter(needed_data['Longitude'], needed_data['Latitude'], c=needed_data['DO'], cmap = 'bwr', s=180,
                            vmin=150, vmax=240, edgecolors='black', linewidths=0.5, label='Argo with max DO under 500m', zorder=5)
            cbar2 = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=cbar_pad)
            cbar2.set_label('DO/μmol·kg⁻¹', fontsize=cbar_label_fs)
            cbar2.ax.tick_params(labelsize=cbar_tick_fs)
        else:
            ax.scatter(needed_data['Longitude'], needed_data['Latitude'], color='blue', s=180, label='Argo', zorder=5)
        for idx, row in needed_data.iterrows():
            ax.text(row['Longitude'], row['Latitude'], f"{int(row['Depth'])}", fontsize=argo_text_fs, fontweight='bold', ha='center', va='center', color='black', zorder=6)
    else:
        print(f"No Argo data available for eddy {ds_names}{no} at the specified index {needed_idx}.")

    # 绘制当前时刻涡旋
    scale_now = approximate_degree_length(center_lat_arr[needed_idx])
    deg_h = radius_arr[needed_idx] / scale_now['meters_per_degree_lat']
    deg_w = radius_arr[needed_idx] / scale_now['meters_per_degree_lon']
    ell_now = Ellipse((center_lon_arr[needed_idx], center_lat_arr[needed_idx]), width=2*deg_w, height=2*deg_h,
                      edgecolor='r', facecolor='none', linestyle='--', alpha=0.2, linewidth=circle_lw, label='Effective Radius')
    ax.add_patch(ell_now)
    ax.scatter(center_lon_arr[needed_idx], center_lat_arr[needed_idx], color='black', s=20, label='Eddy Center', zorder=5)
    ax.plot(curr_contour_lon, curr_contour_lat, color=colors, linewidth=contour_lw, alpha=0.5, label='Effective Contour')

    # 绘制 y = kx + b 直线
    if k is not None and b is not None:
        k_list = [k] if isinstance(k, (int, float)) else k
        b_list = [b] if isinstance(b, (int, float)) else b

        if len(k_list) != len(b_list):
            raise ValueError("The lists for k and b must have the same length.")

        line_x = np.array([glorys_lon_min, glorys_lon_max])

        for i, (k_val, b_val) in enumerate(zip(k_list, b_list)):
            line_y = k_val * line_x + b_val
            ax.plot(line_x, line_y, color='purple', linestyle='-', linewidth=line_lw, label=f'Profile Line {i+1}: y={k_val:.2f}x+{b_val:.2f}')

    ax.legend(fontsize=legend_fs)
    ax.set_xlim(glorys_lon_min, glorys_lon_max)
    ax.set_ylim(glorys_lat_min, glorys_lat_max)
    ax.set_aspect('equal')

    # 紧凑布局，消除多余空白
    plt.tight_layout()

    # 保存图片
    if save_fig:
        output_dir = "plot_track_area_horizontal_glorys"
        os.makedirs(output_dir, exist_ok=True)
        base_filename = f"{ds_names}{no}_{glorys_depth_filtered[0]:.2f}m_{variable}_{dates[needed_idx].strftime('%Y%m%d')}.png"
        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"\nFigure saved to: {os.path.join(output_dir, base_filename)}")

    # 显示图像
    if show_fig:
        # 仅在非内联模式（即交互模式）下，激活直线绘制器
        if not inline_mode:
            line_drawer = LineDrawer(ax)
            fig.canvas.mpl_connect('button_press_event', line_drawer.onclick)
        
        plt.show()

    # 只有在静态内联模式下，才在函数结束时关闭图像以释放内存
    if inline_mode:
        plt.close(fig)

def get_track_area_glorys(DS: list, no: int, needed_idx: int | pd.Timestamp, variables: list = ['thetao'], depth: float | int | None = None):
    '''
    获取指定涡旋在特定时间点周围的 GLORYS 数据。

    该函数会根据涡旋轮廓确定一个矩形区域，并从相应的 GLORYS 文件中
    提取此区域内的一个或多个物理变量。

    参数:
        DS (list): 包含涡旋轨迹信息的数据集。
        no (int): 涡旋的唯一编号。
        needed_idx (int | pd.Timestamp): 需要提取数据的时间点索引或时间戳。
        variables (list): 需要提取的变量列表，默认为 ['thetao']，可选'salinity', 'u', 'v', 'ssh', 'mlt'。
        depth (float | int | None): 如果指定，提取该深度的 GLORYS 数据；如果为 None，则提取2000米以内的所有深度数据。

    返回:
        一个元组，包含筛选后的经度、纬度、深度数组，以及一个存储了所有请求变量数据的字典。
    '''
    wanted_track = find_track(DS, no)
    contour_lon = wanted_track['contour_lon'].values
    contour_lat = wanted_track['contour_lat'].values

    if type(needed_idx) is int:
        glorys_filepaths_dict = find_track_glorys_filepath(DS, no)
        needed_glorys_data = Dataset(list(glorys_filepaths_dict.values())[needed_idx], 'r')
    elif isinstance(needed_idx, pd.Timestamp):
        glorys_filepaths_dict = {needed_idx: get_glorys_filepath(needed_idx)}
        needed_glorys_data = Dataset(glorys_filepaths_dict[needed_idx], 'r')
    else:
        raise ValueError("needed_idx must be an integer index or a pd.Timestamp.")

    contour_lon_filtered = np.ma.masked_equal(contour_lon, 180.0)
    contour_lat_filtered = np.ma.masked_equal(contour_lat, 0.0)

    glorys_lon_min = np.min(contour_lon_filtered) - 0.5
    glorys_lon_max = np.max(contour_lon_filtered) + 0.5
    glorys_lat_min = np.min(contour_lat_filtered) - 0.5
    glorys_lat_max = np.max(contour_lat_filtered) + 0.5

    glorys_lon = needed_glorys_data.variables['longitude'][:]
    glorys_lat = needed_glorys_data.variables['latitude'][:]
    glorys_depth = needed_glorys_data.variables['depth'][:]

    glorys_lon_mask = (glorys_lon >= glorys_lon_min) & (glorys_lon <= glorys_lon_max)
    glorys_lat_mask = (glorys_lat >= glorys_lat_min) & (glorys_lat <= glorys_lat_max)
    if depth is not None:
        glorys_depth_mask = np.zeros_like(glorys_depth, dtype=bool)
        if glorys_depth.size > 0:
            k = np.argmin(np.abs(glorys_depth - depth))
            glorys_depth_mask[k] = True
    else:
        glorys_depth_mask = (glorys_depth >= 0) & (glorys_depth <= 2000)
        
    glorys_lon_filtered = glorys_lon[glorys_lon_mask]
    glorys_lat_filtered = glorys_lat[glorys_lat_mask]
    glorys_depth_filtered = glorys_depth[glorys_depth_mask]

    # 存储多个变量的字典
    glorys_variables_filtered = {}
    for var in variables:
        if var == 'thetao':
            glorys_variables_filtered['thetao'] = needed_glorys_data.variables['thetao'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if depth is not None:
                glorys_variables_filtered['thetao'] = glorys_variables_filtered['thetao'][0, :, :]
        elif var == 'salinity' or var == 'so':
            glorys_variables_filtered['salinity'] = needed_glorys_data.variables['so'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if depth is not None:
                glorys_variables_filtered['salinity'] = glorys_variables_filtered['salinity'][0, :, :]
        elif var == 'u' or var == 'uo':
            glorys_variables_filtered['u'] = needed_glorys_data.variables['uo'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if depth is not None:
                glorys_variables_filtered['u'] = glorys_variables_filtered['u'][0, :, :]
        elif var == 'v' or var == 'vo':
            glorys_variables_filtered['v'] = needed_glorys_data.variables['vo'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if depth is not None:
                glorys_variables_filtered['v'] = glorys_variables_filtered['v'][0, :, :]
        elif var == 'ssh' or var == 'zos':
            glorys_variables_filtered['ssh'] = needed_glorys_data.variables['zos'][0, glorys_lat_mask, glorys_lon_mask]
        elif var == 'mlt' or var == 'mlotst':
            glorys_variables_filtered['mlt'] = needed_glorys_data.variables['mlotst'][0, glorys_lat_mask, glorys_lon_mask]
        else:
            raise ValueError(f"Unsupported variable: {var}. Supported variables are: 'thetao', 'so', 'uo', 'vo', 'zos', 'mlotst'.")

    needed_glorys_data.close()
    
    return glorys_lon_filtered, glorys_lat_filtered, glorys_depth_filtered, glorys_variables_filtered

def calculate_vorticity(lon, lat, u, v):
    '''
    计算给定速度场的相对涡度 (zeta) 和科里奥利参数 (f)。
    该函数可智能处理2D/3D以及Masked Array输入。

    输入的速度场 u, v 的维度应为 (latitude, longitude) 或
    (depth, latitude, longitude)。如果输入为Masked Array，则输出也会是
    Masked Array，其中在mask边缘计算不准确的点会被自动mask掉。

    参数:
        lon (np.ndarray): 一维经度数组 (单位: 度)。
        lat (np.ndarray): 一维纬度数组 (单位: 度)。
        u (np.ndarray | np.ma.MaskedArray): Zonal 速度数组。
        v (np.ndarray | np.ma.MaskedArray): Meridional 速度数组。

    返回:
        tuple: 包含两个元素的元组 (zeta, f)。
               zeta (np.ndarray | np.ma.MaskedArray): 计算得到的相对涡度数组。
               f (np.ndarray | np.ma.MaskedArray): 计算得到的科里奥利参数数组。
    '''
    # --- 1. 输入校验和维度处理 ---
    if u.shape[-2:] != (len(lat), len(lon)) or u.shape != v.shape:
        raise ValueError("速度数组的最后两个维度必须与经纬度数组的长度匹配，且u,v数组形状需一致。")
    if u.ndim not in [2, 3]:
        raise ValueError("输入速度场必须是 2D 或 3D 数组。")

    # --- 2. 智能处理 Masked Array ---
    # 检查输入是否为 masked array
    is_masked_input = np.ma.is_masked(u)

    if is_masked_input:
        # 如果是，用 NaN 填充被 mask 的位置。梯度计算会正确传播NaN。
        u_data = u.filled(np.nan)
        v_data = v.filled(np.nan)
    else:
        # 如果不是，直接使用原始数据
        u_data = u
        v_data = v

    # --- 3. 统一升维处理 ---
    original_ndim = u.ndim
    if original_ndim == 2:
        u_proc = u_data[np.newaxis, :, :]
        v_proc = v_data[np.newaxis, :, :]
    else:
        u_proc = u_data
        v_proc = v_data

    # --- 4. 计算物理坐标间距 (dx, dy) 和科里奥利参数 (f) ---
    Omega = 7.2921e-5  # 地球自转角速度 (弧度/秒)

    # lon, lat 是 1D 数组
    scale = approximate_degree_length(lat) # shape: (n_lat,)
    m_per_deg_lat = scale['meters_per_degree_lat']
    m_per_deg_lon = scale['meters_per_degree_lon']
    
    # 广播到 2D 网格 (lat, lon)
    # m_per_deg_lat 只随 lat 变化，沿 lon 轴广播
    # m_per_deg_lon 只随 lat 变化，沿 lon 轴广播
    m_per_deg_lat_2d = m_per_deg_lat[:, np.newaxis]
    m_per_deg_lon_2d = m_per_deg_lon[:, np.newaxis]
    
    # 计算梯度 (单位: 度)
    # np.gradient(lat) 返回 lat 方向的梯度 (dlat)
    # np.gradient(lon) 返回 lon 方向的梯度 (dlon)
    dlat_grid = np.gradient(lat)[:, np.newaxis] # shape (n_lat, 1)
    dlon_grid = np.gradient(lon)[np.newaxis, :] # shape (1, n_lon)
    
    # dy: 沿 axis=0 (lat) 的距离变化
    dy = dlat_grid * m_per_deg_lat_2d
    
    # dx: 沿 axis=1 (lon) 的距离变化
    dx = dlon_grid * m_per_deg_lon_2d
    
    # 计算科里奥利参数 f
    lon_rad, lat_rad = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))
    f_2d = 2 * Omega * np.sin(lat_rad)

    # 广播 f_2d 到 u_proc 和 v_proc 的深度维度，以便最终输出的 f 形状与 zeta 保持一致
    if original_ndim == 3:
        f_proc = np.tile(f_2d[np.newaxis, :, :], (u_proc.shape[0], 1, 1))
    else: # original_ndim == 2
        f_proc = f_2d # f_2d 已经是 (latitude, longitude) 形状

    # --- 5. 核心计算逻辑 ---
    vorticity_list = []
    for k_slice_idx in range(u_proc.shape[0]):
        u_slice = u_proc[k_slice_idx, :, :]
        v_slice = v_proc[k_slice_idx, :, :]
        dvdx = np.gradient(v_slice, axis=1) / dx
        dudy = np.gradient(u_slice, axis=0) / dy
        vorticity_list.append(dvdx - dudy)
    vorticity_result = np.stack(vorticity_list, axis=0)
    
    # --- 6. 根据输入类型，决定最终输出 ---
    if is_masked_input:
        # 如果输入是 masked，将结果中的 NaN 转回为 mask
        zeta_final = np.ma.masked_invalid(vorticity_result, copy=False)
        f_final = np.ma.masked_invalid(f_proc, copy=False) # f 也可能是 masked array
    else:
        zeta_final = vorticity_result
        f_final = f_proc

    # 如果原始输入是 2D，则降维回去
    if original_ndim == 2:
        return zeta_final.squeeze(axis=0), f_final
    else:
        return zeta_final, f_final
    
def get_idx(DS: list, no: int, start_date: str, end_date: str = None) -> int | list | None:
    '''
    获取指定涡旋编号在给定时间或时间范围内的索引。

    - 如果只提供 start_date，则返回该日期对应的单个索引。
    - 如果同时提供 start_date 和 end_date，则返回该时间范围内的索引列表。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        start_date (str): 起始日期，格式为 'YYYY-MM-DD'。
        end_date (str, optional): 结束日期，格式为 'YYYY-MM-DD'。默认为 None。

    返回:
        int | list | None: 
        - 如果 end_date 为 None，返回单个整数索引或 None (如果未找到)。
        - 如果提供了 end_date，返回一个整数索引的列表。
    '''
    wanted_track = find_track(DS, no)
    if wanted_track is None or len(wanted_track) == 0:
        print(f"未找到涡旋 {no} 的轨迹数据。")
        # 根据调用模式返回正确的空值
        return None if end_date is None else []

    start_date_dt = pd.to_datetime(start_date)

    # 情况一：只查找单个日期的索引
    if end_date is None:
        for idx, track_point in wanted_track.iterrows():
            track_date = convert_date(track_point['time'])
            # 精确匹配日期，忽略时间部分
            if track_date.date() == start_date_dt.date():
                return idx  # 找到后立即返回单个索引
        
        # 如果循环结束仍未找到匹配的日期
        print(f"在涡旋 {no} 的轨迹中未找到日期 {start_date}。")
        return None

    # 情况二：查找一个日期范围内的索引列表
    else:
        end_date_dt = pd.to_datetime(end_date)
        idx_list = []
        for idx, track_point in wanted_track.iterrows():
            track_date = convert_date(track_point['time'])
            if start_date_dt <= track_date <= end_date_dt:
                idx_list.append(idx)
        
        return idx_list

def get_vertical_glorys(DS: list, no: int, needed_idx: int,
                        k: float | list[float], b: float | list[float],
                        variables: list = ['vorticity']) -> list[dict]:
    '''
    计算并返回指定涡旋在特定时刻沿一条或多条 y=kx+b 剖面的物理量数据。

    该函数封装了从GLORYS数据场中提取剖面数据的核心插值计算，并以字典列表的形式返回结果。
    它能正确处理三维变量（如温度，返回二维垂直剖面）和二维变量（如混合层深度，返回一维水平剖面）。
    无论输入变量使用何种别名，输出字典中的键都将是标准化的变量名。

    参数:
        DS (list): 包含涡旋轨迹信息的数据集。
        no (int): 涡旋的唯一编号。
        needed_idx (int): 涡旋轨迹的时间点索引。
        k (float | list[float]): 直线方程 y=kx+b 的斜率或斜率列表。
        b (float | list[float]): 直线方程 y=kx+b 的截距或截距列表。
        variables (list): 需要提取的GLORYS物理变量列表。默认为 ['vorticity']。

    返回:
        list[dict]: 一个包含一个或多个剖面结果字典的列表。列表中每个字典的结构如下：
        
          - **'profile_data' (dict)**:
            - *键*: 标准化的物理量名称 (字符串, 如 'salinity', 'mlt')。
            - *值*:
                - 对于三维变量 (如'vorticity', 'thetao'), 值为一个二维 `numpy.ma.MaskedArray` 数组，代表其垂直剖面，维度为 `(深度层数, 剖面水平点数)`。
                - 对于二维变量 (如'mlt', 'ssh'), 值为一个一维 `numpy.ndarray` 数组，代表其沿剖面线的水平分布，长度为 `剖面水平点数`。

          - **'y_coords' (np.ndarray)**:
            - 一个一维NumPy数组，表示剖面的横坐标轴（物理距离）。
            - 数值单位为公里 (km)，`0` 点对应涡旋中心在剖面线上的投影位置。

          - **'z_coords' (np.ndarray)**:
            - 一个一维NumPy数组，表示剖面的纵坐标轴（物理深度）。
            - 数值单位为米 (m)，代表了GLORYS数据中的深度分层。
          
          - **'lon_coords' (np.ndarray)**:
            - 一个一维NumPy数组，包含了剖面线上每个点（对应`y_coords`）的经度。
          
          - **'lat_coords' (np.ndarray)**:
            - 一个一维NumPy数组，包含了剖面线上每个点（对应`y_coords`）的纬度。

          - **'projections' (dict)**:
            - 一个字典，包含了涡旋边界在横坐标 (`y_coords`) 上的投影位置。
            - *键*: 边界类型 (字符串)， `'radius'` 代表有效半径，`'contour'` 代表有效轮廓。
            - *值*: 一个包含交点位置(km)的列表。

          - **'metadata' (dict)**:
            - 一个包含涡旋元数据的字典。
            - *键*:
                - 'eddy_no'` (int): 涡旋的唯一编号。
                - 'date_str'` (str): 该剖面对应日期的字符串，格式为 'YYYY-MM-DD'。
                - 'k'` (float): 该剖面所用直线方程的斜率。
                - 'b'` (float): 该剖面所用直线方程的截距。
    '''
    # --- 0. 准备工作：统一输入格式并获取公共数据 ---
    if k is None or b is None:
        raise ValueError("k 和 b 必须提供以计算垂直剖面。")

    k_list = [k] if isinstance(k, (int, float)) else k
    b_list = [b] if isinstance(b, (int, float)) else b

    if len(k_list) != len(b_list):
        raise ValueError("k 和 b 的列表长度必须一致。")

    # 建立别名到标准名，以及变量维度的映射
    alias_map = {
        'thetao': 'thetao',
        'salinity': 'salinity', 'so': 'salinity',
        'u': 'u', 'uo': 'u',
        'v': 'v', 'vo': 'v',
        'ssh': 'ssh', 'zos': 'ssh',
        'mlt': 'mlt', 'mlotst': 'mlt',
        'vorticity': 'vorticity'
    }
    var_dims = {
        'thetao': 3, 'salinity': 3, 'u': 3, 'v': 3, 'vorticity': 3,
        'ssh': 2, 'mlt': 2
    }

    raw_vars_to_fetch = set()
    for var in variables:
        if var == 'vorticity': raw_vars_to_fetch.update(['u', 'v'])
        else: raw_vars_to_fetch.add(var)

    wanted_track = find_track(DS, no)
    if wanted_track is None or len(wanted_track) == 0:
        print(f"未找到涡旋 {no} 的轨迹数据。")
        return [{} for _ in k_list]
    
    dates = wanted_track['date'] if 'date' in wanted_track.columns else convert_date(wanted_track['time'])
    contour_lon = wanted_track['contour_lon'].values
    contour_lat = wanted_track['contour_lat'].values
    center_lon_arr = wanted_track['center_lon'].values
    center_lat_arr = wanted_track['center_lat'].values
    radius_arr = wanted_track['radius'].values
    
    glorys_lon_raw, glorys_lat_raw, glorys_depth_raw, glorys_data_raw = get_track_area_glorys(
        DS, no, needed_idx, variables=list(raw_vars_to_fetch)
    )
    
    if glorys_depth_raw.size == 0 and not all(var_dims.get(alias_map.get(v, v)) == 2 for v in variables):
        return [{} for _ in k_list]

    all_profiles_data = []

    # --- 开始循环，为每一对 k, b 计算一个剖面 ---
    for k_val, b_val in zip(k_list, b_list):
        # --- 1. 计算水平剖面线的坐标 ---
        contour_lon_filtered = np.ma.masked_equal(contour_lon, 180.0)
        glorys_lon_min, glorys_lon_max = np.min(contour_lon_filtered) - 0.5, np.max(contour_lon_filtered) + 0.5
        contour_lat_filtered = np.ma.masked_equal(contour_lat, 0.0)
        glorys_lat_min, glorys_lat_max = np.min(contour_lat_filtered) - 0.5, np.max(contour_lat_filtered) + 0.5
        
        num_points = 500
        if k_val == 0:
            profile_lons = np.linspace(glorys_lon_min, glorys_lon_max, num_points)
            profile_lats = np.full_like(profile_lons, b_val)
        else:
            lons_temp = np.linspace(glorys_lon_min, glorys_lon_max, num_points)
            lats_temp = k_val * lons_temp + b_val
            mask = (lats_temp >= glorys_lat_min) & (lats_temp <= glorys_lat_max)
            profile_lons, profile_lats = lons_temp[mask], lats_temp[mask]
            if len(profile_lons) < 2: 
                all_profiles_data.append({})
                continue

        dlat_deg = np.diff(profile_lats)
        dlon_deg = np.diff(profile_lons) # 假设不跨日界线，因为是局部剖面
        mid_lats_deg = (profile_lats[:-1] + profile_lats[1:]) / 2
        
        scale_mid = approximate_degree_length(mid_lats_deg)
        dist_segments = np.hypot(
            dlon_deg * scale_mid['meters_per_degree_lon'],
            dlat_deg * scale_mid['meters_per_degree_lat']
        )
        y_coords_raw = np.insert(np.cumsum(dist_segments), 0, 0) / 1000.0

        current_center_lon, current_center_lat = center_lon_arr[needed_idx], center_lat_arr[needed_idx]
        if k_val == 0:
            xp, yp = current_center_lon, b_val
        else:
            xp = (current_center_lon + k_val * current_center_lat - k_val * b_val) / (1 + k_val**2)
            yp = k_val * xp + b_val
        center_idx_on_profile = np.argmin((profile_lons - xp)**2 + (profile_lats - yp)**2)
        y_coords_recenter = y_coords_raw - y_coords_raw[center_idx_on_profile]

        # --- 2. 区分维度，计算剖面数据 ---
        z_coords = glorys_depth_raw
        profile_data_dict = {}
        for var in variables:
            standard_name = alias_map.get(var, var)
            dimension = var_dims.get(standard_name, 3)

            # --- 2.1 处理2D变量 (如 mlt, ssh) ---
            if dimension == 2:
                data_2d = glorys_data_raw.get(standard_name)
                if data_2d is None or data_2d.size == 0 or np.all(np.ma.getmask(data_2d)):
                    profile_data_dict[standard_name] = np.ma.masked_all(len(profile_lons))
                    continue
                
                interp_func_2d = RegularGridInterpolator((glorys_lat_raw, glorys_lon_raw), data_2d.filled(np.nan), 
                                                         method='linear', bounds_error=False, fill_value=np.nan)
                interpolated_1d = interp_func_2d(list(zip(profile_lats, profile_lons)))
                profile_data_dict[standard_name] = np.ma.masked_invalid(interpolated_1d)
                continue

            # --- 2.2 处理3D变量 (如 thetao, vorticity) ---
            glorys_variable_3d = None
            if standard_name == 'vorticity':
                u, v = glorys_data_raw.get('u'), glorys_data_raw.get('v')
                if u is not None and v is not None and u.size > 0 and v.size > 0:
                    if u.ndim == 2: u, v = u[np.newaxis, :, :], v[np.newaxis, :, :]
                    zeta_3d, f_3d = calculate_vorticity(glorys_lon_raw, glorys_lat_raw, u, v)
                    glorys_variable_3d = zeta_3d / f_3d
            else:
                glorys_variable_3d = glorys_data_raw.get(standard_name)

            if glorys_variable_3d is None or glorys_variable_3d.size == 0 or np.all(np.ma.getmask(glorys_variable_3d)):
                profile_data_dict[standard_name] = np.ma.masked_all((len(z_coords), len(profile_lons)))
                continue

            if glorys_variable_3d.ndim == 2: glorys_variable_3d = glorys_variable_3d[np.newaxis, :, :]
            
            query_depths, query_lats = np.meshgrid(z_coords, profile_lats, indexing='ij')
            _, query_lons = np.meshgrid(z_coords, profile_lons, indexing='ij')
            xi_points = np.vstack([query_depths.ravel(), query_lats.ravel(), query_lons.ravel()]).T
            
            interp_func = RegularGridInterpolator((z_coords, glorys_lat_raw, glorys_lon_raw), glorys_variable_3d.filled(np.nan), method='linear', bounds_error=False, fill_value=np.nan)
            interpolated_values_flat = interp_func(xi_points)
            profile_data_dict[standard_name] = np.ma.masked_invalid(interpolated_values_flat.reshape(len(z_coords), len(profile_lons)))

        # --- 3. 计算边界投影 ---
        scale_line = approximate_degree_length(current_center_lat)
        effective_radius_deg = radius_arr[needed_idx] / scale_line['meters_per_degree_lon']  # 使用经向角度跨度近似
        A, B = 1 + k_val**2, 2 * (k_val*b_val - k_val*current_center_lat - current_center_lon)
        C = current_center_lon**2 + (b_val - current_center_lat)**2 - effective_radius_deg**2
        discriminant = B**2 - 4*A*C
        radius_intersections_lon = [(-B + s * np.sqrt(discriminant)) / (2*A) for s in [-1, 1]] if discriminant >= 0 else []
        radius_proj_dists = [y_coords_raw[np.argmin((profile_lons - lon_i)**2 + (profile_lats - (k_val*lon_i + b_val))**2)] - y_coords_raw[center_idx_on_profile] for lon_i in radius_intersections_lon]

        contour_lon_valid = contour_lon[needed_idx][contour_lon[needed_idx] != 180.0]
        contour_lat_valid = contour_lat[needed_idx][contour_lat[needed_idx] != 0.0]
        contour_intersections_xy = find_polygon_line_intersections(contour_lon_valid, contour_lat_valid, profile_lons, profile_lats)
        contour_proj_dists = [y_coords_raw[np.argmin((profile_lons - lon_i)**2 + (profile_lats - lat_i)**2)] - y_coords_raw[center_idx_on_profile] for lon_i, lat_i in contour_intersections_xy]

        projections_dict = {'radius': sorted(radius_proj_dists), 'contour': sorted(contour_proj_dists)}
        
        # --- 4. 整合单个剖面的结果 ---
        single_profile_result = {
            'profile_data': profile_data_dict,
            'y_coords': y_coords_recenter,
            'z_coords': z_coords,
            'lon_coords': profile_lons,
            'lat_coords': profile_lats,
            'projections': projections_dict,
            'metadata': {
                'eddy_no': no,
                'date_str': dates[needed_idx].strftime("%Y-%m-%d"),
                'k': k_val,
                'b': b_val,
            }
        }
        all_profiles_data.append(single_profile_result)

    # --- 5. 返回所有剖面的结果列表 ---
    return all_profiles_data

def plot_vertical_glorys(DS: list, no: int, needed_idx: int, 
                         k: float | list[float], b: float | list[float], 
                         variable: str = 'vorticity',
                         show_fig: bool = False, save_fig: bool = False, 
                         xmin: float = None, xmax: float = None,
                         ymin: float = None, ymax: float = None,
                         plot_mlt: bool = False):
    '''
    绘制指定涡旋在特定时刻，沿一条或多条剖面线 (y = kx + b) 的物理量垂直剖面图。

    该函数调用 get_vertical_glorys 获取指定变量的数据，然后进行可视化。
    当 k 和 b 为列表时，会为每一对 (k, b) 生成一张独立的图表。
    可以选择性地在图上叠加混合层深度曲线。

    参数:
        (同 get_vertical_glorys, 但 variable 为 str，一次只处理一个变量)
        ...
        show_fig (bool): 是否显示图像。
        save_fig (bool): 是否保存图像。
        xmin, xmax, ymin, ymax (float): 坐标轴范围。
        plot_mlt (bool): 是否在图上绘制混合层深度分界线。默认为 False。
    '''
    # --- 1. 获取所有计算好的剖面数据包 ---
    
    vars_to_fetch = {variable}
    if plot_mlt:
        vars_to_fetch.add('mlt')
        
    all_data_packages = get_vertical_glorys(DS, no, needed_idx, k, b, variables=list(vars_to_fetch))

    if not all_data_packages:
        print(f"警告: get_vertical_glorys 未能返回任何数据。绘图已取消。")
        return

    # --- 开始循环，为每一个数据包生成一张图 ---
    for data_package in all_data_packages:
        if not data_package:
            print(f"警告: 收到一个空的数据包，跳过此剖面的绘图。")
            continue
            
        profile_variable_2d = data_package['profile_data'].get(variable)
        if profile_variable_2d is None:
             alias_map = {'so': 'salinity', 'uo': 'u', 'vo': 'v'}
             standard_name = alias_map.get(variable, variable)
             profile_variable_2d = data_package['profile_data'].get(standard_name)

        if profile_variable_2d is None or np.all(getattr(profile_variable_2d, 'mask', True)):
            k_meta, b_meta = data_package.get('metadata', {}).get('k'), data_package.get('metadata', {}).get('b')
            print(f"警告: 变量 '{variable}' 在剖面 k={k_meta}, b={b_meta} 上的数据无效。绘图已取消。")
            continue

        # --- 2. 准备绘图 ---
        y_coords = data_package['y_coords']
        z_coords = data_package['z_coords']
        projections = data_package['projections']
        metadata = data_package['metadata'] 

        callers_local_vars = inspect.currentframe().f_back.f_locals
        ds_name = [var_name for var_name, var_val in callers_local_vars.items() if var_val is DS][0].upper()
        
        prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        eddy_color = prop_colors[1] if 'AC' in ds_name else prop_colors[0]
        
        if profile_variable_2d.shape[1] != len(y_coords):
           y_coords = y_coords[:profile_variable_2d.shape[1]]

        Y_mesh, Z_mesh = np.meshgrid(y_coords, z_coords)
        
        # 设置变量相关的绘图属性
        if variable == 'vorticity': cbar_label, cmap, clim = r'$\zeta/f$', 'seismic', (-0.3, 0.3)
        elif variable in ['thetao']: cbar_label, cmap = 'Temperature (°C)', 'rainbow'
        elif variable in ['salinity', 'so']: cbar_label, cmap = 'Salinity (psu)', 'viridis'
        elif variable in ['u', 'v', 'uo', 'vo']: cbar_label, cmap = 'Velocity (m/s)', 'RdBu_r'
        else: cbar_label, cmap, clim = variable, 'viridis', None
        
        if 'clim' not in locals() or clim is None:
            valid_values = profile_variable_2d[~profile_variable_2d.mask]
            clim = (valid_values.min(), valid_values.max()) if valid_values.size > 0 else (0,1)
            if variable in ['u', 'v', 'uo', 'vo']:
                max_abs = np.max(np.abs(valid_values)) if valid_values.size > 0 else 1
                clim = (-max_abs, max_abs)

        # --- 3. 执行绘图 ---
        fig, ax = plt.subplots(figsize=(20, 15))
        
        date_str = metadata['date_str']
        title = (f"Vertical Profile of {cbar_label} for Track {ds_name}{metadata['eddy_no']} "
                 f"on {date_str}, y={metadata['k']:.2f}x+{metadata['b']:.2f}")
        ax.set_title(title, fontsize=20)
        ax.set_xlabel('Distance from Eddy Center Projection (km)', fontsize=18)
        ax.set_ylabel('Depth (m)', fontsize=18)
        ax.tick_params(labelsize=14)
        
        pc = ax.pcolormesh(Y_mesh, Z_mesh, profile_variable_2d, cmap=cmap, shading='auto', vmin=clim[0], vmax=clim[1])
        cbar = fig.colorbar(pc, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
        cbar.set_label(cbar_label, fontsize=18)
        cbar.ax.tick_params(labelsize=14)

        ax.axvline(0, color='black', linestyle='--', linewidth=2, label='Eddy Center Projection')
        for i, dist in enumerate(projections['radius']): 
            ax.axvline(dist, color='r', linestyle='--', linewidth=2, label='Effective Radius Projection' if i == 0 else "")
        for i, dist in enumerate(projections['contour']): 
            ax.axvline(dist, color=eddy_color, linestyle=':', linewidth=2, label='Effective Contour Projection' if i == 0 else "")
            
        # 绘制混合层深度
        mld_lines = () # 创建一个空元组，用于存放MLD的两条线
        if plot_mlt:
            mlt_data = data_package['profile_data'].get('mlt')
            if mlt_data is not None and not np.all(getattr(mlt_data, 'mask', True)):
                # **捕获两条线的艺术家对象**
                # 注意 plot 返回的是一个列表，所以我们用 l, 来解包
                black_line, = ax.plot(y_coords, mlt_data, color='black', linewidth=2.5, zorder=3)
                white_line, = ax.plot(y_coords, mlt_data, color='white', linewidth=1.5, zorder=4, label='Mixed Layer Depth')
                mld_lines = (black_line, white_line) # 将两条线存入元组
            else:
                print(f"注意: 未能在剖面 k={metadata['k']}, b={metadata['b']} 上找到有效的混合层深度数据。")

        ax.set_ylim(z_coords.max(), z_coords.min())
        if xmin is not None and xmax is not None: ax.set_xlim(xmin, xmax)
        if ymin is not None and ymax is not None: ax.set_ylim(ymax, ymin)
        
        # --- 构建并自定义图例 ---
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles)) # 创建一个去重的标签-句柄字典

        # 如果我们画了MLD线，就用我们创建的复合句柄替换掉自动生成的单个白线句柄
        if plot_mlt and mld_lines:
            by_label['Mixed Layer Depth'] = mld_lines
            
        ax.legend(by_label.values(), by_label.keys(), fontsize=18)

        # --- 4. 保存和显示 ---
        if save_fig:
            output_dir = "plot_vertical_glorys"
            os.makedirs(output_dir, exist_ok=True)
            date_fn = date_str.replace('-', '')
            base_filename = (f"{ds_name}{metadata['eddy_no']}_vertical_{variable}_{date_fn}_"
                             f"k{metadata['k']:.2f}b{metadata['b']:.2f}.png")
            plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

        if show_fig:
            plt.show()

        plt.close(fig)

# 帮助函数：判断三个点 (p, q, r) 的方向（共线，顺时针，逆时针）
def _orientation(p, q, r):
    val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
    if val == 0: return 0  # Collinear
    return 1 if val > 0 else 2 # Clockwise or Counterclockwise (1 for clockwise, 2 for counterclockwise)

# 帮助函数：判断点 q 是否在线段 pr 上
def _on_segment(p, q, r):
    return (q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and
            q[1] <= max(p[1], r[1]) and q[1] >= min(p[1], r[1]))

# 核心帮助函数：计算两条线段的交点
def _line_segment_intersect(p1, q1, p2, q2):
    """
    Finds the intersection point of two line segments (p1, q1) and (p2, q2).
    Returns the intersection point as (x, y) if they intersect, otherwise None.
    Handles collinear cases where segments overlap.
    """
    o1 = _orientation(p1, q1, p2)
    o2 = _orientation(p1, q1, q2)
    o3 = _orientation(p2, q2, p1)
    o4 = _orientation(p2, q2, q1)

    # General case: Non-collinear and non-parallel intersection
    if o1 != o2 and o3 != o4:
        denom = (p1[0] - q1[0]) * (p2[1] - q2[1]) - (p1[1] - q1[1]) * (p2[0] - q2[0])
        if denom == 0:
            return None 

        t = ((p1[0] - p2[0]) * (p2[1] - q2[1]) - (p1[1] - p2[1]) * (p2[0] - q2[0])) / denom
        u = -((p1[0] - q1[0]) * (p1[1] - p2[1]) - (p1[1] - q1[1]) * (p1[0] - p2[0])) / denom

        if (0 <= t <= 1) and (0 <= u <= 1):
            intersection_x = p1[0] + t * (q1[0] - p1[0])
            intersection_y = p1[1] + t * (q1[1] - p1[1])
            return (intersection_x, intersection_y)
        return None

    # Special Cases for collinear segments
    if o1 == 0 and _on_segment(p1, p2, q1): return p2
    if o2 == 0 and _on_segment(p1, q2, q1): return q2
    if o3 == 0 and _on_segment(p2, p1, q2): return p1
    if o4 == 0 and _on_segment(p2, q1, q2): return q1

    return None


# 寻找多边形和线段的交点
def find_polygon_line_intersections(polygon_lon, polygon_lat, line_lons, line_lats, tolerance=1e-6):
    """
    Finds intersection points of a line (defined by line_lons, line_lats)
    with a closed polygon (defined by polygon_lon, polygon_lat).

    Args:
        polygon_lon (array): Longitudes of polygon vertices.
        polygon_lat (array): Latitudes of polygon vertices.
        line_lons (array): Longitudes of line points.
        line_lats (array): Latitudes of line points.
        tolerance (float): Tolerance for checking duplicate intersection points.

    Returns:
        list: A list of (lon, lat) tuples for intersection points.
    """
    intersections = []
    
    if len(polygon_lon) < 3: # 需要至少三个点才能构成多边形
        return []

    polygon_points_closed = list(zip(polygon_lon, polygon_lat))
    if not (np.isclose(polygon_points_closed[0][0], polygon_points_closed[-1][0], atol=tolerance) and
            np.isclose(polygon_points_closed[0][1], polygon_points_closed[-1][1], atol=tolerance)):
        polygon_points_closed.append(polygon_points_closed[0])

    line_points = list(zip(line_lons, line_lats))

    for i in range(len(polygon_points_closed) - 1):
        poly_seg_p1 = polygon_points_closed[i]
        poly_seg_p2 = polygon_points_closed[i+1]
        
        for j in range(len(line_points) - 1):
            line_seg_p1 = line_points[j]
            line_seg_p2 = line_points[j+1]
            
            intersect_pt = _line_segment_intersect(poly_seg_p1, poly_seg_p2, line_seg_p1, line_seg_p2)
            if intersect_pt:
                is_duplicate = False
                for existing_pt in intersections:
                    if np.isclose(existing_pt[0], intersect_pt[0], atol=tolerance) and \
                       np.isclose(existing_pt[1], intersect_pt[1], atol=tolerance):
                        is_duplicate = True
                        break
                if not is_duplicate:
                    intersections.append(intersect_pt)
    
    if len(intersections) > 0:
        if np.abs(line_lons[0] - line_lons[-1]) > np.abs(line_lats[0] - line_lats[-1]):
            intersections.sort(key=lambda p: p[0])
        else:
            intersections.sort(key=lambda p: p[1])
            
    return intersections

def regrid_vertical_slice(data_package: dict, dy: float = None, dz: float = None) -> dict:
    '''
    将垂直剖面数据包插值到新的等间距网格上。

    该函数接收 get_vertical_glorys 的输出，并根据用户指定的水平(dy)或
    垂直(dz)间距，生成一个新的、网格化更规整的数据包。

    参数:
        data_package (dict): 从 get_vertical_glorys 函数获取的原始数据包。
        dy (float, optional): 新的水平(y轴)网格间距，单位为公里(km)。
                              默认为 None，表示不改变水平网格。
        dz (float, optional): 新的垂直(z轴)网格间距，单位为米(m)。
                              默认为 None，表示不改变垂直网格。

    返回:
        dict: 一个结构与输入相同但数据已被插值到新网格上的新数据包。
              字典包含以下键值对：

            - **'profile_data' (dict)**:
                - 其内部结构与输入相同，但每个二维数组的值都是插值后的结果。

            - **'y_coords' (np.ndarray)**:
                - 如果提供了 `dy`，这将是一个新生成的一维等间距数组。
                - 否则，与输入的 'y_coords' 相同。

            - **'z_coords' (np.ndarray)**:
                - 如果提供了 `dz`，这将是一个新生成的一维等间距数组。
                - 否则，与输入的 'z_coords' 相同。
            
            - **'lon_coords' (np.ndarray)**:
                - 如果水平坐标被重采样，这将是新的一维插值经度数组。
                - 否则，与输入的 'lon_coords' 相同。

            - **'lat_coords' (np.ndarray)**:
                - 如果水平坐标被重采样，这将是新的一维插值纬度数组。
                - 否则，与输入的 'lat_coords' 相同。

            - **'projections' (dict)**:
                - 从原始数据包中原样复制而来，其数值仍对应原始的坐标系。
                
            - **'metadata' (dict)**:
                - 从原始数据包中原样复制而来。
    '''
    # 如果用户未指定任何新的网格间距，则直接返回原始数据包的深拷贝
    if dy is None and dz is None:
        return copy.deepcopy(data_package)

    # --- 1. 提取原始数据和坐标 ---
    original_y = data_package['y_coords']
    original_z = data_package['z_coords']
    original_profiles = data_package['profile_data']
    original_lon = data_package['lon_coords']
    original_lat = data_package['lat_coords']
    
    # RegularGridInterpolator 要求坐标必须是严格递增的。
    if np.any(np.diff(original_z) <= 0):
        raise ValueError("原始z坐标（深度）必须是严格递增的才能进行插值。")
    if np.any(np.diff(original_y) <= 0):
        raise ValueError("原始y坐标（距离）必须是严格递增的才能进行插值。")

    # --- 2. 创建新的网格坐标 ---
    new_y = np.arange(original_y.min(), original_y.max(), dy) if dy is not None else original_y
    if dz is not None:
        z_min, z_max = original_z.min(), original_z.max()
        new_z = np.arange(z_min, z_max, dz)
    else:
        new_z = original_z

    # --- 3. 对每个变量的剖面数据及经纬度进行插值 ---
    # 插值经纬度 (1D)，仅当水平坐标被重采样时执行
    if dy is not None:
        new_lon = np.interp(new_y, original_y, original_lon)
        new_lat = np.interp(new_y, original_y, original_lat)
    else:
        new_lon, new_lat = original_lon, original_lat

    # 插值剖面数据 (2D)
    new_profiles = {}
    # 创建查询点网格，用于高效插值, 'ij'索引确保维度顺序正确
    zz_grid, yy_grid = np.meshgrid(new_z, new_y, indexing='ij')
    query_points = np.vstack([zz_grid.ravel(), yy_grid.ravel()]).T

    for var, data_array in original_profiles.items():
        if data_array.mask.all():
            new_profiles[var] = np.ma.masked_all((len(new_z), len(new_y)))
            continue

        # 创建插值器实例
        interp_func = RegularGridInterpolator(
            (original_z, original_y), 
            data_array.filled(np.nan),
            method='linear', 
            bounds_error=False, 
            fill_value=np.nan
        )
        
        # 在新网格上执行插值
        new_data_flat = interp_func(query_points)
        
        # 将一维插值结果重塑为二维网格
        new_data = new_data_flat.reshape(len(new_z), len(new_y))
        
        # 将结果转回 masked array 并存入新字典
        new_profiles[var] = np.ma.masked_invalid(new_data)

    # --- 4. 构建并返回新的数据包 ---
    new_data_package = {
        'profile_data': new_profiles,
        'y_coords': new_y,
        'z_coords': new_z,
        'lon_coords': new_lon,
        'lat_coords': new_lat,
        'projections': copy.deepcopy(data_package['projections']),
        'metadata': copy.deepcopy(data_package['metadata'])
    }

    return new_data_package

def plot_data_package(data_package: dict, DS: list, variable: str,
                      show_fig: bool = False, save_fig: bool = False, xmin: float = None, xmax: float = None,
                      ymin: float = None, ymax: float = None):
    '''
    根据一个数据包和原始涡旋数据集，绘制单一物理量的垂直剖面图。

    该函数是一个灵活的可视化接口，它接收一个数据包，并根据传入的涡旋数据集
    (DS)来确定绘图风格（如颜色和标题），适用于 get_vertical_glorys 或
    regrid_vertical_slice 的输出。

    参数:
        data_package (dict): 从 get_vertical_glorys 或 regrid_vertical_slice 获取的数据包。
        DS (list): 原始的涡旋轨迹信息数据集 (如ACL, CL)，用于确定涡旋类型和名称。
        variable (str): 需要从 data_package 中绘制的变量名。
        show_fig (bool): 是否显示图像。
        save_fig (bool): 是否保存图像。
        xmin, xmax, ymin, ymax (float): 坐标轴范围。
    '''
    # --- 1. 验证输入并解包数据 ---
    required_keys = ['profile_data', 'y_coords', 'z_coords', 'lon_coords', 'lat_coords', 'projections', 'metadata']
    if not data_package or not all(k in data_package for k in required_keys):
        print(f"错误: 输入的 data_package 格式不完整。缺少键。")
        return
        
    profile_variable_2d = data_package['profile_data'].get(variable)
    if profile_variable_2d is None or np.all(profile_variable_2d.mask):
        print(f"警告: 变量 '{variable}' 的剖面数据无效，无法绘图。")
        return

    y_coords = data_package['y_coords']
    z_coords = data_package['z_coords']
    projections = data_package['projections']
    metadata = data_package['metadata']

    # --- 2. 在当前上下文中计算与DS相关的元数据 ---
    callers_local_vars = inspect.currentframe().f_back.f_locals
    ds_name_list = [var_name for var_name, var_val in callers_local_vars.items() if var_val is DS]
    if not ds_name_list:
        raise ValueError("无法在调用者环境中找到数据集变量名。请确保DS参数已正确传入。")
    ds_name = ds_name_list[0].upper()
    
    prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    eddy_color = prop_colors[1] if 'AC' in ds_name else prop_colors[0]

    # --- 3. 准备绘图元素 ---
    if profile_variable_2d.shape[1] != len(y_coords):
       y_coords = y_coords[:profile_variable_2d.shape[1]]

    Y_mesh, Z_mesh = np.meshgrid(y_coords, z_coords)
    
    if variable == 'vorticity': cbar_label, cmap, clim = r'$\zeta/f$', 'seismic', (-0.3, 0.3)
    elif variable in ['thetao']: cbar_label, cmap = 'Temperature (°C)', 'rainbow'
    elif variable in ['salinity', 'so']: cbar_label, cmap = 'Salinity (psu)', 'viridis'
    elif variable in ['u', 'v', 'uo', 'vo']: cbar_label, cmap = 'Velocity (m/s)', 'RdBu_r'
    else: cbar_label, cmap, clim = variable, 'viridis', None
    
    if 'clim' not in locals() or clim is None:
        valid_values = profile_variable_2d[~profile_variable_2d.mask]
        clim = (valid_values.min(), valid_values.max()) if valid_values.size > 0 else (0,1)
        if variable in ['u', 'v', 'uo', 'vo']:
            max_abs = np.max(np.abs(valid_values)) if valid_values.size > 0 else 1
            clim = (-max_abs, max_abs)

    # --- 4. 执行绘图 ---
    fig, ax = plt.subplots(figsize=(20, 15))
    
    k_val, b_val = metadata.get('k'), metadata.get('b')
    if k_val is not None and b_val is not None:
        title = (f"Vertical Profile of {cbar_label} for Track {ds_name}{metadata['eddy_no']} "
                 f"on {metadata['date_str']}, y={k_val:.2f}x+{b_val:.2f}")
    else:
         title = (f"Vertical Profile of {cbar_label} for Track {ds_name}{metadata['eddy_no']} "
                 f"on {metadata['date_str']}")

    ax.set_title(title, fontsize=20)
    ax.set_xlabel('Distance from Eddy Center Projection (km)', fontsize=18)
    ax.set_ylabel('Depth (m)', fontsize=18)
    ax.tick_params(labelsize=14)
    
    pc = ax.pcolormesh(Y_mesh, Z_mesh, profile_variable_2d, cmap=cmap, shading='auto', vmin=clim[0], vmax=clim[1])
    cbar = fig.colorbar(pc, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label, fontsize=18)
    cbar.ax.tick_params(labelsize=14)

    ax.axvline(0, color='black', linestyle='--', linewidth=2, label='Eddy Center Projection')
    for i, dist in enumerate(projections['radius']): 
        ax.axvline(dist, color='r', linestyle='--', linewidth=2, label='Effective Radius Projection' if i == 0 else "")
    for i, dist in enumerate(projections['contour']): 
        ax.axvline(dist, color=eddy_color, linestyle=':', linewidth=2, label='Effective Contour Projection' if i == 0 else "")

    ax.set_ylim(z_coords.max(), z_coords.min())
    if xmin is not None and xmax is not None: ax.set_xlim(xmin, xmax)
    if ymin is not None and ymax is not None: ax.set_ylim(ymax, ymin)
    
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=18)

    # --- 5. 保存和显示 ---
    if save_fig:
        output_dir = "plot_vertical_glorys"
        os.makedirs(output_dir, exist_ok=True)
        date_fn = metadata['date_str'].replace('-', '')
        
        if k_val is not None and b_val is not None:
            k_str = f"k{k_val:.2f}"
            b_str = f"b{b_val:.2f}"
            base_filename = (f"{ds_name}{metadata['eddy_no']}_vertical_{variable}_{date_fn}_{k_str}{b_str}.png")
        else:
            base_filename = (f"{ds_name}{metadata['eddy_no']}_vertical_{variable}_{date_fn}.png")

        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

    if show_fig:
        plt.show()

    plt.close(fig)

def get_raytraceR_inputs(data_package: dict,
                                f_coriolis: float, # 科里奥利频率 f (例如 1.454e-4)。**必须提供**。
                                rho0: float = 1025.0, # 参考密度 (kg/m^3)。
                                g: float = 9.81, # 重力加速度 (m/s^2)。
                                m0: float = 1/100.0, # 初始垂直波数。
                                omega_factor: float = 0.97, # 固有频率 omega = omega_factor * f。
                                v0_amplitude: float = 0.1, # 初始波速振幅。
                                thresh_val: float = 0.1e-2, # 内部反射阈值。
                                chstart_val: int = 1, # 初始特征线 (1或2)。
                                filter_sigma_z: float = 0.5, # 沿Z方向高斯平滑的标准差（用于导数计算，防止噪声放大）。
                                filter_sigma_y: float = 0.5 # 沿Y方向高斯平滑的标准差。
                               ) -> dict:
    """
    将来自 track.regrid_vertical_slice 的数据包转换为 MATLAB raytraceR 函数所需的输入格式。

    假设输入 data_package 的结构如下（来自 track.get_vertical_glorys 和 regrid_vertical_slice 处理后）：
    - z_coords (np.ndarray): 等间距的深度 (Z) 坐标，单位为米(m)。
    - y_coords (np.ndarray): 等间距的横流 (Y) 坐标，单位为公里(km)。
    - lon_coords (np.ndarray): 剖面线上每个点对应的经度。
    - lat_coords (np.ndarray): 剖面线上每个点对应的纬度。
    - profile_data (dict): 包含以下二维 np.ndarray (或 np.ma.MaskedArray) 数组的字典:
        - 'u': 纬向速度 (uo)
        - 'v': 经向速度 (vo)
        - 'salinity': 盐度 (Practical Salinity, SP)
        - 'thetao': 位势温度 (potential temperature, pt)
    
    参数:
        data_package (dict): 包含 GLORYS 垂直切片数据的字典，通常是
                             track.regrid_vertical_slice 的输出。
        f_coriolis (float): 科里奥利频率 f (例如 1.454e-4 rad/s)。**必须提供**。
        rho0 (float): 参考密度 (kg/m^3)。
        g (float): 重力加速度 (m/s^2)。
        m0 (float): 初始垂直波数。
        omega_factor (float): 固有频率的倍数，omega = omega_factor * f。
        v0_amplitude (float): 初始波速振幅。
        thresh_val (float): 内部反射阈值。
        chstart_val (int): 初始特征线 (1或2)。
        filter_sigma_z (float): 沿深度方向高斯平滑的标准差（用于导数计算）。
        filter_sigma_y (float): 沿横流方向高斯平滑的标准差。

    返回:
        dict: 包含 raytraceR 所需所有输入参数的字典。
    """

    # 1. 提取和映射输入数据
    # raytraceR 的 Z 是垂直/深度，Y 是水平横流。
    z_g = data_package['z_coords'] # 深度坐标，单位：米 (m)
    y_g_km = data_package['y_coords'] # 横流坐标，单位：公里 (km)
    lon_coords_1d = data_package['lon_coords'] # 剖面经度 (1D array)
    lat_coords_1d = data_package['lat_coords'] # 剖面纬度 (1D array)

    # 将横流坐标从公里转换为米，以保持 raytraceR 内部单位一致性
    y_g = y_g_km * 1000.0 # 横流坐标，单位：米 (m)

    # 从 profile_data 中获取核心物理量
    # 它们已经是 (深度点数, 横流点数) 的形状
    # .filled(np.nan) 用于处理 MaskedArray，将其掩码值填充为 NaN
    u_data = data_package['profile_data']['u'].filled(np.nan)
    v_data = data_package['profile_data']['v'].filled(np.nan)
    salinity_data = data_package['profile_data']['salinity'].filled(np.nan) # Practical Salinity (SP)
    thetao_data = data_package['profile_data']['thetao'].filled(np.nan) # Potential Temperature (pt)

    # 检查维度和数据有效性
    expected_shape_2d = (len(z_g), len(y_g_km))
    if not (u_data.shape == salinity_data.shape == thetao_data.shape == expected_shape_2d):
        raise ValueError(f"Profile data shapes are inconsistent or do not match expected ({len(z_g)}, {len(y_g_km)}).")
    if not np.isfinite(f_coriolis):
        raise ValueError("Coriolis frequency 'f_coriolis' must be a finite number.")

    # 2. 计算网格间距
    dz = np.mean(np.diff(z_g)) # regrid_vertical_slice 保证了等间距
    dy = np.mean(np.diff(y_g)) # 横流坐标已经转换为米，这里计算米为单位的间距

    # 3. 对物理量进行平滑处理 (减少噪声对导数计算的影响)
    # 确保 sigma 值与数据尺寸兼容
    u_smoothed = gaussian_filter(u_data, sigma=(filter_sigma_z, filter_sigma_y))
    v_smoothed = gaussian_filter(v_data, sigma=(filter_sigma_z, filter_sigma_y))
    salinity_smoothed = gaussian_filter(salinity_data, sigma=(filter_sigma_z, filter_sigma_y))
    thetao_smoothed = gaussian_filter(thetao_data, sigma=(filter_sigma_z, filter_sigma_y))
    
    # 处理平滑后可能出现的NaN（高斯滤波的边界效应），尤其是在mask区域外推的NaN
    u_smoothed[np.isnan(u_smoothed)] = 0.0 # 假设用0填充，或选择其他填充值
    v_smoothed[np.isnan(v_smoothed)] = 0.0
    salinity_smoothed[np.isnan(salinity_smoothed)] = np.nanmean(salinity_smoothed) if np.all(np.isnan(salinity_smoothed)) else salinity_smoothed[np.isnan(salinity_smoothed)].mean() # 用平均值填充，避免gsw出错
    thetao_smoothed[np.isnan(thetao_smoothed)] = np.nanmean(thetao_smoothed) if np.all(np.isnan(thetao_smoothed)) else thetao_smoothed[np.isnan(thetao_smoothed)].mean() # 同上

    # 4. 计算核心物理量: 密度, 浮力 (bg), N2, ug, F2, S2, s_M, omegamin
    
    # 4.1. 转换为 gsw 所需的 Absolute Salinity (SA) 和 Conservative Temperature (CT)
    # gsw 需要经度和纬度网格与数据形状匹配
    Z_mesh_m, Y_mesh_m = np.meshgrid(z_g, y_g, indexing='ij') # Z_mesh_m (深度, 米), Y_mesh_m (横流，米)
    
    # 将 1D lon/lat_coords 扩展为 2D 网格，匹配 (深度点数, 横流点数) 形状
    # 假设 lon_coords_1d 和 lat_coords_1d 是与 y_g_km 对应的，即沿横流方向的经纬度
    # 那么它们需要沿着深度方向广播
    lon_2d = np.tile(lon_coords_1d[np.newaxis, :], (len(z_g), 1))
    lat_2d = np.tile(lat_coords_1d[np.newaxis, :], (len(z_g), 1))

    # 计算压力（p_from_z 需要负深度）
    pressure_2d = gsw.p_from_z(-Z_mesh_m, lat_2d)

    # 实用盐度 -> 绝对盐度
    SA_processed = gsw.SA_from_SP(salinity_smoothed, pressure_2d, lon_2d, lat_2d)
    
    # 位温 -> 保守温度
    CT_processed = gsw.CT_from_pt(SA_processed, thetao_smoothed)

    # 计算位势密度 (referenced to 0 dbar)
    rho_potential_processed = gsw.rho(SA_processed, CT_processed, 0)
    
    # 计算浮力 bg = -g * (rho_potential - rho0) / rho0
    bg_processed = -g * (rho_potential_processed - rho0) / rho0

    # 4.2. 计算 N2 (浮力频率的平方)
    # gsw.Nsquared(SA, CT, p, lat=None, axis=0) 返回 (N2, p_mid)
    # axis=0 表示沿深度（第一个）维度计算梯度
    N2_output, p_mid_N2 = gsw.Nsquared(SA_processed, CT_processed, pressure_2d, lat=lat_2d, axis=0)
    
    # N2_output 的形状比原始数据在深度维度上少一个点 (因为是 mid-point 导数)
    # 需要将其插值或映射回原始深度网格的尺寸
    # 常用方法是插值回原始深度点或在第一个点进行外推/复制
    # 这里我们采用在第一个深度点复制 N2_output 的第一行，使其维度匹配
    N2_processed = np.vstack([N2_output[0:1, :], N2_output])
    # 再次裁剪以确保与原始深度维度完全一致（防止因浮点数问题多/少一行）
    N2_processed = N2_processed[:u_data.shape[0], :] 
    
    # 处理 N2 中可能出现的 NaN 或负值（物理上 N2 应为正或零）
    N2_processed[np.isnan(N2_processed)] = 0.0 # 用 0 填充 NaN
    N2_processed[N2_processed < 0] = 0.0 # 强制非负

    # 4.3. 计算 ug (背景地转流场)
    # 论文中 ug 是背景地转流的 x 分量。这里直接使用平滑后的 u 分量作为 ug。
    ug_processed = u_smoothed

    # 4.4. 计算 F2, S2 (根据 Whitt and Thomas 2013 论文定义)
    
    # zeta_g = -du_g/dy
    # np.gradient 自动处理 MaskedArray，但我们已填充为 NaN
    # 注意 dy 已经是米
    du_dy_smoothed = np.gradient(ug_processed, dy, axis=1) # du_g/dy
    zeta_g_processed = -du_dy_smoothed
    
    # F2 = f * (f + zeta_g)
    F2_processed = f_coriolis * (f_coriolis + zeta_g_processed)

    # S2 = f * du_g/dz
    # 注意 dz 已经是米
    du_dz_smoothed = np.gradient(ug_processed, dz, axis=0) # du_g/dz
    S2_processed = f_coriolis * du_dz_smoothed
    
    # 4.5. 计算 s_M = F2 / S2
    # 避免除以零或 NaN
    s_M_processed = np.divide(F2_processed, S2_processed, out=np.zeros_like(F2_processed), where=S2_processed!=0)
    s_M_processed[np.isnan(s_M_processed)] = 0.0 # 填充 NaN
    s_M_processed[np.isinf(s_M_processed)] = 0.0 # 填充 Inf

    # 4.6. 计算 omegamin
    # omegamin = sqrt(q/N2)
    # q = F2 * N2 - S2^4
    PV_q = F2_processed * N2_processed - S2_processed**4
    
    # 避免根号下出现负值或除以零
    omegamin_squared = np.divide(PV_q, N2_processed, out=np.zeros_like(PV_q), where=N2_processed!=0)
    omegamin_squared = np.maximum(0, omegamin_squared) # 确保非负

    omegamin_processed = np.sqrt(omegamin_squared)
    omegamin_processed = np.real(omegamin_processed) # 确保是实数
    omegamin_processed[np.isnan(omegamin_processed)] = 0.0 # 填充 NaN
    omegamin_processed[np.isinf(omegamin_processed)] = 0.0 # 填充 Inf

    # 5. 二阶导数 (根据论文示例，默认可以设为零，或者在这里计算)
    # 鉴于计算二阶导数对噪声敏感且需要额外代码，我们这里保持为零，
    # 除非用户明确要求并提供计算方法。
    d2udy2_processed = np.zeros_like(u_data)
    d2udz2_processed = np.zeros_like(u_data)
    d2bdy2_processed = np.zeros_like(u_data)
    d2bdz2_processed = np.zeros_like(u_data)
    d2bdzdy_processed = np.zeros_like(u_data)
    d2udzdy_processed = np.zeros_like(u_data)

    # 6. 计算领域长度和波参数
    omega = omega_factor * f_coriolis
    Lz = z_g.max() - z_g.min()
    Ly = y_g.max() - y_g.min() # y_g 已经是米

    # m0 = 1/minv，所以 minv = 1/m0
    minv_val = 1.0 / m0

    # 准备输出字典，包含 raytraceR 所需所有参数
    raytraceR_inputs = {
        'Z': -Z_mesh_m, # 2D 深度网格 (深度, 横流)，单位：米 (m)
        'dz': dz, # 深度间距 (米)
        'Y': Y_mesh_m, # 2D 横流网格 (深度, 横流)，单位：米 (m)
        'dy': dy, # 横流间距 (米)
        'F2': F2_processed, # 2D 矩阵 (深度, 横流)
        'S2': S2_processed, # 2D 矩阵 (深度, 横流)
        'N2': N2_processed, # 2D 矩阵 (深度, 横流)
        'ug': ug_processed, # 2D 矩阵 (深度, 横流)
        'bg': bg_processed, # 2D 矩阵 (深度, 横流)
        's_M': s_M_processed, # 2D 矩阵 (深度, 横流)
        'omegamin': omegamin_processed, # 2D 矩阵 (深度, 横流)
        'f': f_coriolis, # 标量
        'g': g, # 标量
        'rho0': rho0, # 标量
        'omega': omega, # 标量
        'minv': minv_val, # 标量
        'v0': v0_amplitude, # 标量
        
        # 二阶导数 (目前设为零)
        'd2udy2': d2udy2_processed,
        'd2udz2': d2udz2_processed,
        'd2bdy2': d2bdy2_processed,
        'd2bdz2': d2bdz2_processed,
        'd2bdzdy': d2bdzdy_processed,
        'd2udzdy': d2udzdy_processed,
        
        # 领域长度 (根据 Z, Y 范围计算)
        'Lz': Lz,
        'Ly': Ly,
        
        # 其他 raytraceR 需要但通常在调用时设定的参数
        'thresh': thresh_val,
        'm0': m0, # 与 minv 互补，此处也返回
        'chstart': chstart_val,
        # projections数据用于Matlab中绘图
        'eddy_radius': data_package['projections'].get('radius', []), # 半径投影
        'eddy_contour': data_package['projections'].get('contour', []), # 等值线投影
    }
    
    print("\n--- Data Preparation Summary for raytraceR ---")
    print(f"Z (Depth) meshgrid shape: {Z_mesh_m.shape}, min/max: {Z_mesh_m.min():.1f}/{Z_mesh_m.max():.1f} m")
    print(f"Y (Cross-stream) meshgrid shape: {Y_mesh_m.shape}, min/max: {Y_mesh_m.min()/1000:.1f}/{Y_mesh_m.max()/1000:.1f} km")
    print(f"dz (Depth spacing): {dz:.3f} m")
    print(f"dy (Cross-stream spacing): {dy:.3f} m")
    print(f"Calculated F2 shape: {F2_processed.shape}, min/max: {np.nanmin(F2_processed):.2e}/{np.nanmax(F2_processed):.2e}")
    print(f"Calculated S2 shape: {S2_processed.shape}, min/max: {np.nanmin(S2_processed):.2e}/{np.nanmax(S2_processed):.2e}")
    print(f"Calculated N2 shape: {N2_processed.shape}, min/max: {np.nanmin(N2_processed):.2e}/{np.nanmax(N2_processed):.2e}")
    print(f"Calculated ug shape: {ug_processed.shape}, min/max: {np.nanmin(ug_processed):.2e}/{np.nanmax(ug_processed):.2e}")
    print(f"Calculated bg shape: {bg_processed.shape}, min/max: {np.nanmin(bg_processed):.2e}/{np.nanmax(bg_processed):.2e}")
    print(f"Calculated omegamin shape: {omegamin_processed.shape}, min/max: {np.nanmin(omegamin_processed):.2e}/{np.nanmax(omegamin_processed):.2e}")
    print(f"Coriolis frequency (f): {f_coriolis:.4e} rad/s")
    print(f"Wave frequency (omega): {omega:.4e} rad/s")
    print("--------------------------------------------\n")

    return raytraceR_inputs

def init_worker(eddy_data_shared: dict):
    """
    (Initializer函数) 为multiprocessing的子进程初始化共享的、只读的数据。

    功能:
        这个函数在每个工作进程启动时仅被调用一次。它接收大的数据集
        (如涡旋数据字典)并将其设置为该进程的全局变量。这可以极大地避免
        在每个任务间重复传输大数据的开销，是性能优化的关键。
    """
    # 限制 OpenMP/MKL 线程数，避免多进程 x 多线程导致 CPU 争抢
    # 注意：这只影响子进程环境，不影响主进程
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    # 强制切换 Matplotlib 后端为非交互式 (Agg)，避免多进程绘图死锁
    # 注意：init_worker 仅在子进程中运行，此处修改不会影响 Jupyter 主进程的交互式绘图
    try:
        plt.switch_backend('Agg')
    except Exception:
        pass
    
    global worker_eddy_datasets
    worker_eddy_datasets = eddy_data_shared

def check_single_track(
    track_data,
    argo_by_date,
    start_date,
    end_date,
    ds_name,
    circle_enlargement_factor: float | None = None,
    use_adaptive_circle: bool = False,
    adaptive_lat_threshold: float = 70.0,
    adaptive_distance_threshold_km: float = 300.0,
    force_great_circle_circle: bool = False,
    save_interacting_argo: bool = False,
):
    """
    (内部辅助函数) 检查单个涡旋轨迹是否与Argo数据有交集。

    功能:
        这是一个纯计算函数，会为所有在时间范围内的涡旋返回结果，
        并附带一个布尔标志来说明其是否与Argo浮标有交集。

    参数:
        track_data (list): 单条涡旋的轨迹数据 (list of lists)。
        argo_by_date (dict): 按日期组织的 Argo 明细，形如 {date: list[dict]}，
            每个 dict 至少包含 'Longitude','Latitude'，可附带 Profile_number、Year/Month/Day、
            delta_do、do_value/DO、Anomaly_depth 等元数据。
        start_date (pd.Timestamp): 检查的开始日期。
        end_date (pd.Timestamp): 检查的结束日期。
        ds_name (str): 数据集名称 (如 'ACS')。
        circle_enlargement_factor (float | None): 半径放大因子，None 回退全局配置。
        use_adaptive_circle (bool): True 时半径距离用 `adaptive_distance_m` 自适应大圆。
        adaptive_lat_threshold (float): |lat| 超过此值触发大圆距离。
        adaptive_distance_threshold_km (float): 平面距离超过此值(km)触发大圆距离。
        force_great_circle_circle (bool): 强制使用大圆距离（忽略阈值）。
        save_interacting_argo (bool): True 时收集并返回每个命中的 Argo 点明细（含 method/track 等），
            False 时为性能考虑在首个命中即停止当日迭代且不返回点明细。

    返回:
        dict | None: 如果涡旋在时间范围内，返回包含绘图/判定信息的字典，否则返回 None。
        主要键：
          - 'track_id': 轨迹编号
          - 'has_interaction': 是否与 Argo 交互
          - 'in_range_segments': 连续片段用于绘图
          - 'contours_to_plot': 命中多边形时用于绘图的等值线
          - 'candidate_dates_for_contour': 圆命中待进一步二次判定的日期
          - 'dates_in_range': 轨迹在时间窗内的日期
          - 'text_info': 绘图标签信息
          - 'is_ace': 是否反气旋
          - 'interacting_argo': 当 save_interacting_argo=True 时返回 list[dict]，
             每个 dict 至少含 {'date','lon','lat','method'(poly/circle),'ds_name'}，
                 并按 argo_by_date 附带 Profile_number/指标等元数据；否则为空列表。
    """
    num, time, center_lon, center_lat, _, _, contour_lon, contour_lat, radius, _, _ = zip(*track_data)
    dates = convert_date(time)
    
    indices_in_range = np.where((dates >= start_date) & (dates <= end_date))[0]
    if indices_in_range.size == 0:
        return None

    if circle_enlargement_factor is None:
        circle_enlargement_factor = globals().get('circle_enlargement_factor', 1.2)

    has_interaction = False
    contours_to_plot = []
    candidate_dates_for_contour = set()
    interacting_points: list[dict] = []
    for i in indices_in_range:
        current_date = dates[i].normalize()
        if current_date in argo_by_date:
            center_lon_i = float(center_lon[i])
            center_lat_i = float(center_lat[i])
            radius_i = float(radius[i])

            # 预计算等值线多边形（若存在）
            contour_poly = None
            contour_lon_norm = None
            try:
                contour_lon_i = contour_lon[i]
                contour_lat_i = contour_lat[i]
                if contour_lon_i is not None and contour_lat_i is not None:
                    contour_lon_arr = np.asarray(contour_lon_i, dtype=float)
                    contour_lat_arr = np.asarray(contour_lat_i, dtype=float)
                    if contour_lon_arr.size >= 3 and contour_lat_arr.size >= 3 and contour_lon_arr.shape == contour_lat_arr.shape:
                        contour_lon_norm = center_lon_i + _minimal_lon_diff_deg(contour_lon_arr, center_lon_i)
                        contour_poly = MplPath(list(zip(contour_lon_norm, contour_lat_arr)))
            except Exception:
                contour_poly = None

            points_today = argo_by_date[current_date]

            any_point_inside_poly = False
            added_contour_for_today = False
            for idx_pt, argo_point in enumerate(points_today):
                # argo_point 为 dict，至少含 Longitude/Latitude；兼容 lon/lat 键
                point_lon = float(argo_point.get('Longitude', argo_point.get('lon')))
                point_lat = float(argo_point.get('Latitude', argo_point.get('lat')))
                if use_adaptive_circle:
                    distance_m = adaptive_distance_m(
                        point_lon,
                        point_lat,
                        center_lon_i,
                        center_lat_i,
                        wrap_dateline=True,
                        gc_lat_threshold=adaptive_lat_threshold,
                        gc_distance_threshold_km=adaptive_distance_threshold_km,
                        force_great_circle=force_great_circle_circle
                    )
                else:
                    scale_ci = approximate_degree_length(center_lat_i)
                    dlon_deg = _minimal_lon_diff_deg(point_lon, center_lon_i)
                    dx_m = dlon_deg * scale_ci['meters_per_degree_lon']
                    dy_m = (point_lat - center_lat_i) * scale_ci['meters_per_degree_lat']
                    distance_m = np.hypot(dx_m, dy_m)
                inside_circle = distance_m <= radius_i * circle_enlargement_factor

                inside_poly_point = False
                if contour_poly is not None:
                    point_lon_norm = center_lon_i + _minimal_lon_diff_deg(point_lon, center_lon_i)
                    try:
                        inside_poly_point = contour_poly.contains_point((point_lon_norm, point_lat))
                    except Exception:
                        inside_poly_point = False
                if inside_poly_point:
                    any_point_inside_poly = True

                if inside_poly_point or inside_circle:
                    has_interaction = True
                    if inside_circle:
                        candidate_dates_for_contour.add(current_date)
                    # 仅在首次发现多边形命中时加入对应等值线以供绘制
                    if inside_poly_point and (not added_contour_for_today) and contour_poly is not None and contour_lon[i] is not None and contour_lat[i] is not None:
                        contours_to_plot.append((contour_lon[i], contour_lat[i]))
                        added_contour_for_today = True
                    if save_interacting_argo:
                        rec = {
                            'method': 'poly' if inside_poly_point else 'circle',
                            'ds_name': ds_name,
                        }
                        if isinstance(argo_point, dict):
                            rec.update(argo_point)
                        # 若缺失日期拆分列，使用当前日期补齐
                        if 'Year' not in rec or pd.isna(rec.get('Year')):
                            rec['Year'] = int(current_date.year)
                        if 'Month' not in rec or pd.isna(rec.get('Month')):
                            rec['Month'] = int(current_date.month)
                        if 'Day' not in rec or pd.isna(rec.get('Day')):
                            rec['Day'] = int(current_date.day)
                        # 移除与原始字段重复的派生列
                        rec.pop('lon', None)
                        rec.pop('lat', None)
                        rec.pop('date', None)
                        interacting_points.append(rec)
                    # 若不需要收集全部点，则命中一个即可跳出当天循环
                    if not save_interacting_argo:
                        break
    
    in_range_segments = []
    splits = np.where(np.diff(indices_in_range) != 1)[0] + 1
    contiguous_blocks = np.split(indices_in_range, splits)
    for block in contiguous_blocks:
        if block.size > 0:
            in_range_segments.append((np.array(center_lon)[block], np.array(center_lat)[block]))

    start_idx = indices_in_range[0]
    # 统一使用旧版结构首字段（已在新路径中改为 track_id）
    _tid_label = int(num[0])
    text_info = {"text": f"{ds_name}{_tid_label}", "lon": center_lon[start_idx], "lat": center_lat[start_idx]}

    return {
        "ds_name": ds_name,
        "track_id": int(_tid_label) if len(num) else None,
        "center_lon": center_lon,
        "center_lat": center_lat,
        "in_range_segments": in_range_segments,
        "contours_to_plot": contours_to_plot,
        "candidate_dates_for_contour": sorted(candidate_dates_for_contour),
        "dates_in_range": [pd.Timestamp(d).normalize() for d in dates[indices_in_range]],
        "text_info": text_info,
        "is_ace": 'AC' in ds_name.upper(),
        "has_interaction": has_interaction,
        "interacting_argo": interacting_points if save_interacting_argo else [],
    }

def plot_all_tracks_in_range(
    start_date_str: str,
    end_date_str: str,
    eddy_datasets: dict | list[str] | tuple[str, ...] | None = None,
    plot_unrelated_eddies: bool = False,
    plot_unrelated_argo: bool = True,
    save_fig: bool = False,
    skip_save_if_empty: bool = False,
    show_labels: bool = True,
    show_fig: bool = True,
    circle_enlargement_factor: float | None = None,
    do_threshold: float | None = None,
    salinity_threshold: float | None = None,
    temperature_threshold: float | None = None,
    depth_interval: float | None = None,
    depth_merge_tolerance: float | None = None,
    duplicate_depth_strategy: str | None = None,
    anomaly_min_depth: float | None = None,
    anomaly_color_by: str = 'delta_do',
    fix_delta_do_colorbar: bool = True,
    delta_do_cbar_min: float = 50.0,
    delta_do_cbar_max: float = 100.0,
    delta_do_cbar_ticks: list | None = None,
    meta_output_root: str | Path | None = None,
    save_interacted_eddies: bool = False,
    save_interacting_argo: bool = False,
):
    """(核心绘图) 指定时间段内涡旋轨迹 + Argo ΔDO 异常代表点（仅采用 ΔDO 方法）。

    依赖 Cartopy 进行制图，自动处理跨国际日期变更线的轨迹连线。

    工作流程：
      1. 装载时间范围内 Argo 数据 → 过滤地理范围 → 计算 ΔDO 异常。
      2. 若 anomaly_min_depth > 0，则按该阈值过滤异常深度。
      3. 每个剖面保留 delta_do（或 do_value）最大的一条。
      4. 按 anomaly_color_by 着色：'delta_do' (默认) 或 'do_value'。

    参数:
        start_date_str, end_date_str: 日期范围。
        eddy_datasets:
          - 兼容旧版: dict，如 {'ACS': acs, 'ACL': acl, 'CS': cs, 'CL': cl}，每个值是“轨迹列表(list of tracks)”；
          - 新版便捷: 字符串列表/元组，如 ['acs','acl','cs','cl']，函数将按配置从 META_tracks 读取对应时间段与区域内的轨迹；
          - None: 并行 worker 模式下默认使用全局 worker_eddy_datasets。
        plot_unrelated_eddies: 是否绘制未与 Argo 交互的涡旋。
        plot_unrelated_argo: 是否额外绘制所有 Argo 剖面位置（空心圆），用于提供基准分布背景。
        save_fig, show_fig: 输出控制。
        skip_save_if_empty: 若为 True 且本图中未绘制任何涡旋（不含底图/Argo点），则跳过保存；默认 False（单次绘图默认不跳过）。
        show_labels: 是否绘制轨迹文本标签（如 ACLXXXX）。
        do_threshold / salinity_threshold / temperature_threshold / depth_interval / depth_merge_tolerance / duplicate_depth_strategy: 传给 calculate_delta_do。
        anomaly_min_depth: (可选) 仅保留异常深度 >= 此值；≤0 不限制；None 表示使用 processing.yml 的 anomaly_min_depth。
        anomaly_color_by: 'delta_do' 或 'do_value'。
        fix_delta_do_colorbar: 若为 True 且按 delta_do 着色，则强制使用 [delta_do_cbar_min, delta_do_cbar_max] 作为色标范围。
        delta_do_cbar_min / delta_do_cbar_max: ΔDO 色标固定范围上下限（仅在 fix_delta_do_colorbar=True 且 anomaly_color_by='delta_do' 时生效）。
        delta_do_cbar_ticks: 自定义 ΔDO 色标刻度列表（None 自动：若只提供上下限则显示两端；若范围>30 添加中点）。
        meta_output_root: 指定 META_tracks 根目录（可覆盖配置默认）。
        save_interacted_eddies (bool): True 时保存本期交互涡旋标签（NPY）；默认 False。
        save_interacting_argo (bool): True 时保存本期交互 Argo 明细（Parquet）；默认 False，输出目录固定为按阈值划分的子目录。

    输出（按阈值子目录 thr{thr}[_depth{d}m] 保存）：
        - All_Tracks_{start}_to_{end}.png
        - Interacted_Eddies_{start}_{end}_thr{thr}[_depth{d}m}.npy（当 save_interacted_eddies=True 时）
        - Interacting_Argo_{start}_{end}_thr{thr}[_depth{d}m}.parquet（当 save_interacting_argo=True 时）
    """
    # --- 0. 确定数据源 ---
    local_eddy_datasets = eddy_datasets
    if circle_enlargement_factor is None:
        circle_enlargement_factor = globals().get('circle_enlargement_factor', 1.2)
    if do_threshold is None:
        do_threshold = _default_delta_do_threshold
    if salinity_threshold is None:
        salinity_threshold = _default_salinity_threshold
    if temperature_threshold is None:
        temperature_threshold = _default_temperature_threshold
    if depth_interval is None:
        depth_interval = _default_depth_interval
    if depth_merge_tolerance is None:
        depth_merge_tolerance = _default_depth_merge_tolerance
    if duplicate_depth_strategy is None:
        duplicate_depth_strategy = _default_duplicate_depth_strategy
    if anomaly_min_depth is None:
        anomaly_min_depth = _cfg_anomaly_min_depth
    # 若未显式提供，则在并行 worker 中从全局共享获得
    if local_eddy_datasets is None:
        try:
            global worker_eddy_datasets
            local_eddy_datasets = worker_eddy_datasets
        except NameError:
            local_eddy_datasets = None

    # 新版便捷用法：若传入的是字符串列表/元组（kinds），先不加载等值线以加速第一轮判定
    lazy_contour_mode = False
    if isinstance(local_eddy_datasets, (list, tuple)):
        kinds = [str(k).lower() for k in local_eddy_datasets]
        local_eddy_datasets = _load_eddy_datasets_for_range(
            kinds=kinds,
            start_date=pd.to_datetime(start_date_str),
            end_date=pd.to_datetime(end_date_str),
            output_root=meta_output_root,
            include_contours=False,
        )
        lazy_contour_mode = True

    start_date = pd.to_datetime(start_date_str)
    end_date = pd.to_datetime(end_date_str)
    
    # --- 1. 动态加载Argo数据 ---
    min_year, max_year = start_date.year, end_date.year
    argo_data_list = []
    for year in range(min_year, max_year + 1):
        try:
            argo_data_list.append(load_argo_data(year, data_dir=argo_path))
        except FileNotFoundError:
            print(f"Warning: Argo data for year {year} not found, skipping.")
    if not argo_data_list:
        print(f"No Argo data found for the period {start_date_str} to {end_date_str}.")
        return

    argo_data = pd.concat(argo_data_list, ignore_index=True)
    
    # --- 2. 预处理Argo数据 ---
    argo_dates = pd.to_datetime(argo_data[['Year', 'Month', 'Day']])
    argo_in_range = argo_data[(argo_dates >= start_date) & (argo_dates <= end_date)]
    needed_argo_data = pd.DataFrame()
    base_argo_positions = pd.DataFrame()
    if not argo_in_range.empty:
        lon_vals = argo_in_range['Longitude'].to_numpy(dtype=float, copy=False)
        lat_vals = argo_in_range['Latitude'].to_numpy(dtype=float, copy=False)
        lon_mask = _region_lon_mask(lon_vals, lonmin, lonmax)
        lat_mask = (lat_vals >= latmin) & (lat_vals <= latmax)
        geo_mask = lon_mask & lat_mask
        argo_in_geo_range = argo_in_range[geo_mask].copy()
        if not argo_in_geo_range.empty:
            # 所有剖面基础位置（避免按深度重复）
            base_argo_positions = (
                argo_in_geo_range.sort_values(['Profile_number','Depth'])
                .groupby('Profile_number', as_index=False)
                .first()[['Profile_number','Longitude','Latitude']]
            )
            anomalies = calculate_delta_do(
                argo_in_geo_range,
                depth_interval=depth_interval,
                do_threshold=do_threshold,
                salinity_threshold=salinity_threshold,
                temperature_threshold=temperature_threshold,
                anomaly_min_depth=anomaly_min_depth,
                depth_merge_tolerance=depth_merge_tolerance,
                duplicate_depth_strategy=duplicate_depth_strategy,
                remove_outliers=True,
                verbose=False
            )
            if not anomalies.empty:
                sort_field = 'delta_do' if (
                    anomaly_color_by == 'delta_do' and 'delta_do' in anomalies.columns
                ) else 'do_value'
                anomalies_sorted = anomalies.sort_values(by=[sort_field], ascending=False)
                anomalies_unique = anomalies_sorted.drop_duplicates(subset='Profile_number', keep='first')
                needed_argo_data = anomalies_unique.rename(columns={'depth': 'Anomaly_depth'})

    argo_by_date = defaultdict(list)
    if not needed_argo_data.empty:
        for _, row in needed_argo_data.iterrows():
            date_key = pd.Timestamp(year=int(row['Year']), month=int(row['Month']), day=int(row['Day']))
            argo_by_date[date_key].append({
                'Profile_number': row.get('Profile_number'),
                'Longitude': float(row.get('Longitude')),
                'Latitude': float(row.get('Latitude')),
                'Year': int(row.get('Year')),
                'Month': int(row.get('Month')),
                'Day': int(row.get('Day')),
                'delta_do': float(row.get('delta_do')) if 'delta_do' in row else np.nan,
                'do_value': float(row.get('do_value')) if 'do_value' in row else (float(row.get('DO')) if 'DO' in row else np.nan),
                'Anomaly_depth': float(row.get('Anomaly_depth')) if 'Anomaly_depth' in row else np.nan,
            })

    # --- 3. 检查所有涡旋轨迹 ---
    # 兼容两种形态：
    #  - 旧形态：ds_data 为 [track_list, ...]
    #  - 新形态：ds_data 为 [(track_list, track_id), ...]
    all_tracks_with_names: list[tuple[list, str, int | None]] = []
    for ds_name, ds_data in local_eddy_datasets.items():
        for item in ds_data:
            if isinstance(item, tuple) and len(item) == 2 and isinstance(item[1], (int, np.integer)):
                track_list, tid = item
            else:
                track_list, tid = item, None
            all_tracks_with_names.append((track_list, ds_name, tid))
    # 第一轮判定时的圆半径：使用原始 circle_enlargement_factor（不再放大）。
    _cef_first_pass = circle_enlargement_factor
    results = [
        check_single_track(
            track,
            argo_by_date,
            start_date,
            end_date,
            ds_name,
            circle_enlargement_factor=_cef_first_pass,
            save_interacting_argo=bool(save_interacting_argo),
        )
        for track, ds_name, tid in all_tracks_with_names
    ]

    # --- 4. 绘图（第一轮：不画等值线） ---
    crosses_dateline = bool(_REGION_CFG.get('crosses_dateline') and (lonmax < lonmin))
    central_lon = 180 if crosses_dateline else 0
    data_crs = ccrs.PlateCarree()
    map_crs = ccrs.PlateCarree(central_longitude=central_lon)

    fig = plt.figure(figsize=(40, 30))
    ax = fig.add_subplot(1, 1, 1, projection=map_crs)
    ax.set_title(
        f"Eddy Tracks and Argo ΔDO Anomalies ({start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')})",
        fontsize=20
    )

    # Basemap features（自然地理背景）
    base_ocean = _BASEMAP_COLORS['ocean']
    base_land = _BASEMAP_COLORS['land']
    coast_color = _BASEMAP_COLORS['coastline']
    grid_color = _BASEMAP_COLORS['grid']

    ax.set_facecolor(base_ocean)
    ax.add_feature(cfeature.OCEAN, facecolor=base_ocean, zorder=0)
    ax.add_feature(cfeature.LAND, facecolor=base_land, edgecolor=coast_color, linewidth=0.5, zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.7, edgecolor=coast_color, zorder=1)

    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color=grid_color, alpha=0.45, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # 地图范围处理（跨日界线时扩展 longitude）
    lon_extent_min = lonmin
    lon_extent_max = lonmax
    if crosses_dateline:
        if lon_extent_max < lon_extent_min:
            lon_extent_max += 360
    ax.set_extent([lon_extent_min, lon_extent_max, latmin, latmax], crs=data_crs)
    
    prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    eddy_colors = {'ACE': prop_colors[1], 'CE': prop_colors[0]}

    any_eddy_drawn = False
    label_texts = []
    label_points_x = []
    label_points_y = []
    for result in filter(None, results):
        has_interaction = result['has_interaction']
        if not has_interaction and not plot_unrelated_eddies:
            continue
        
        is_ace = result['is_ace']
        color = eddy_colors['ACE'] if is_ace else eddy_colors['CE']
        text_color = 'red' if has_interaction else 'black'
        
        ax.plot(
            result['center_lon'],
            result['center_lat'],
            color=color,
            alpha=0.4,
            linestyle='--',
            zorder=4,
            transform=ccrs.Geodetic(),
        )
        any_eddy_drawn = True
        for lon_seg, lat_seg in result['in_range_segments']:
            ax.plot(
                lon_seg,
                lat_seg,
                color=color,
                alpha=0.8,
                linestyle='-',
                zorder=4,
                transform=ccrs.Geodetic(),
            )
        
        if show_labels:
            info = result['text_info']
            text_obj = ax.text(
                info['lon'],
                info['lat'],
                info['text'],
                fontsize=12,
                color=text_color,
                weight='bold',
                zorder=5,
                clip_on=False,
                transform=data_crs,
            )
            label_texts.append(text_obj)
            label_points_x.append(info['lon'])
            label_points_y.append(info['lat'])

        # 第一轮不画等值线；若不是 lazy 模式但已有 contours_to_plot，则保持现状
        if (not lazy_contour_mode) and has_interaction:
            for contour_lon, contour_lat in result['contours_to_plot']:
                ax.plot(
                    contour_lon,
                    contour_lat,
                    color=color,
                    linewidth=1,
                    alpha=0.5,
                    zorder=4,
                    linestyle=':',
                    transform=ccrs.Geodetic(),
                )

    if plot_unrelated_argo and not base_argo_positions.empty:
        ax.scatter(
            base_argo_positions['Longitude'], base_argo_positions['Latitude'],
            facecolors='none', edgecolors='gray', linewidths=0.8, s=36,
            label='All Argo Profiles (baseline)', zorder=2,
            transform=data_crs,
        )

    if not needed_argo_data.empty:
        if anomaly_color_by == 'delta_do' and 'delta_do' in needed_argo_data.columns:
            # 应用固定色标（可选）
            scatter_kwargs = {}
            if fix_delta_do_colorbar:
                scatter_kwargs.update(dict(vmin=delta_do_cbar_min, vmax=delta_do_cbar_max))
            depth_label = (
                f' @ depth ≥ {anomaly_min_depth} m'
                if anomaly_min_depth is not None and anomaly_min_depth > 0
                else ''
            )
            sc = ax.scatter(
                needed_argo_data['Longitude'], needed_argo_data['Latitude'],
                c=needed_argo_data['delta_do'], cmap='Reds', s=70,
                edgecolors='black', linewidths=0.5,
                label=f'ΔDO ≥ {do_threshold} μmol kg⁻¹{depth_label}', zorder=3,
                transform=data_crs,
                **scatter_kwargs
            )
            cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
            cbar.set_label('ΔDO / μmol·kg⁻¹', fontsize=20); cbar.ax.tick_params(labelsize=14)
            # 自动或用户自定义刻度
            if fix_delta_do_colorbar:
                if delta_do_cbar_ticks is not None:
                    cbar.set_ticks(delta_do_cbar_ticks)
                else:
                    # 默认：只显示上下限；若范围较大则加中点
                    rng = delta_do_cbar_max - delta_do_cbar_min
                    if rng > 30:
                        mid = (delta_do_cbar_max + delta_do_cbar_min) / 2
                        cbar.set_ticks([delta_do_cbar_min, mid, delta_do_cbar_max])
                    else:
                        cbar.set_ticks([delta_do_cbar_min, delta_do_cbar_max])
        else:
            color_values = needed_argo_data.get('do_value') if 'do_value' in needed_argo_data.columns else needed_argo_data.get('DO')
            sc = ax.scatter(
                needed_argo_data['Longitude'], needed_argo_data['Latitude'],
                c=color_values, cmap='bwr', s=60, vmin=150, vmax=240,
                edgecolors='black', linewidths=0.5,
                label='Argo DO Anomaly Profiles',
                zorder=3,
                transform=data_crs,
            )
            cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
            cbar.set_label('DO / μmol·kg⁻¹', fontsize=20); cbar.ax.tick_params(labelsize=14)

    legend_elements = [
        Line2D([0], [0], color=eddy_colors['ACE'], lw=2, label='ACE Track'),
        Line2D([0], [0], color=eddy_colors['CE'], lw=2, label='CE Track')
    ]
    handles, labels = ax.get_legend_handles_labels()
    extra_labels = [f"Argo ΔDO ≥ {do_threshold}", 'Argo DO Anomaly Profiles', 'All Argo Profiles (baseline)']
    added = None
    for lab in extra_labels:
        if lab in labels:
            added = handles[labels.index(lab)]
            break
    if added is not None:
        ax.legend(handles=legend_elements + [added], fontsize=18, loc='upper left')
    else:
        ax.legend(handles=legend_elements, fontsize=18, loc='upper left')

    # --- 4.2 第二轮：按需补充绘制等值线（仅对 lazy 模式下有交互的涡旋；使用真实 track_id） ---
    if lazy_contour_mode:
        try:
            need_items = [r for r in results if r and r.get('has_interaction')]
            if need_items:
                need_by_ds: dict[str, dict[int, set[pd.Timestamp]]] = {}
                for r in need_items:
                    ds = r['ds_name']
                    tid = r.get('track_id')
                    if tid is None:
                        continue
                    # 仅针对“半径命中”的候选日期进一步检查等值线
                    dates_req = set(r.get('candidate_dates_for_contour', []) or r.get('dates_in_range', []))
                    if not dates_req:
                        continue
                    need_by_ds.setdefault(ds, {})
                    need_by_ds[ds].setdefault(tid, set()).update(dates_req)

                for ds_name, tid_map in need_by_ds.items():
                    kind_l = ds_name.lower()
                    rich = find_track(
                        kind_l,
                        sorted(tid_map.keys()),
                        region=_current_region_key(),
                        include_contours=True,
                        return_list=False,
                    )
                    # 统一为 {track_id: DataFrame}
                    grouped: dict[int, pd.DataFrame]
                    if isinstance(rich, dict):
                        grouped = {
                            k: pd.DataFrame(v, columns=['track_id','time','center_lon','center_lat','max_lon','max_lat','contour_lon','contour_lat','radius','speed_contour_lon','speed_contour_lat'])
                            for k, v in rich.items()
                        }
                    else:
                        df_rich = pd.DataFrame(rich) if isinstance(rich, list) else rich
                        if isinstance(df_rich, pd.DataFrame) and 'track_id' in df_rich.columns:
                            grouped = {tid: g for tid, g in df_rich.groupby('track_id')}
                        else:
                            only_tid = next(iter(tid_map.keys())) if len(tid_map)==1 else None
                            grouped = {only_tid: df_rich if isinstance(df_rich, pd.DataFrame) else pd.DataFrame(df_rich)}

                    for tid, dates_req in tid_map.items():
                        df_t = grouped.get(tid)
                        if df_t is None or df_t.empty:
                            continue
                        if 'date' not in df_t.columns:
                            if 'time' in df_t.columns:
                                df_t = df_t.copy(); df_t['date'] = convert_date(df_t['time'])
                            else:
                                continue
                        mask = df_t['date'].isin(dates_req)
                        df_needed = df_t[mask]
                        if df_needed.empty:
                            continue
                        color2 = eddy_colors['ACE'] if 'AC' in ds_name else eddy_colors['CE']
                        for _, row in df_needed.iterrows():
                            cl = row.get('contour_lon'); ct = row.get('contour_lat')
                            if isinstance(cl, (list, np.ndarray)) and isinstance(ct, (list, np.ndarray)) and len(cl) >= 3 and len(ct) >= 3:
                                # 二次判定：仅当 Argo 点在等值线内时才绘制该 contour
                                try:
                                    contour_lon_arr = np.asarray(cl, dtype=float)
                                    contour_lat_arr = np.asarray(ct, dtype=float)
                                    center_lon_i = float(row.get('center_lon', contour_lon_arr[0]))
                                    contour_lon_norm = center_lon_i + _minimal_lon_diff_deg(contour_lon_arr, center_lon_i)
                                    contour_lon_norm = np.asarray(contour_lon_norm, dtype=float)
                                    
                                    path = MplPath(list(zip(contour_lon_norm, contour_lat_arr)))

                                    date_norm = row['date'].normalize() if isinstance(row['date'], pd.Timestamp) else pd.Timestamp(row['date']).normalize()
                                    entries = argo_by_date.get(date_norm, [])
                                    draw_it = False
                                    for entry in entries:
                                        try:
                                            pt_lon = float(entry.get('Longitude', entry.get('lon')))
                                            pt_lat = float(entry.get('Latitude', entry.get('lat')))
                                        except Exception:
                                            continue
                                        point_lon_norm = float(center_lon_i + _minimal_lon_diff_deg(pt_lon, center_lon_i))
                                        
                                        if path.contains_point((point_lon_norm, pt_lat)):
                                            draw_it = True
                                            break
                                except Exception:
                                    draw_it = False
                                if draw_it:
                                    ax.plot(
                                        contour_lon_norm,
                                        contour_lat_arr,
                                        color=color2,
                                        linewidth=1,
                                        alpha=0.6,
                                        zorder=4,
                                        linestyle=':',
                                        transform=ccrs.Geodetic(),
                                    )
        except Exception as e:
            print(f"[LazyContour] failed to draw contours lazily: {e}")

    # --- 4.3 标签避让 ---
    orig_xlim = ax.get_xlim()
    orig_ylim = ax.get_ylim()
    region_is_global = _current_region_key().lower() == 'global'
    if region_is_global and adjust_text is not None and label_texts:
        try:
            span_x = max(abs(orig_xlim[1] - orig_xlim[0]), 1e-6)
            span_y = max(abs(orig_ylim[1] - orig_ylim[0]), 1e-6)
            span_ratio_x = span_x / 360.0
            span_ratio_y = span_y / 180.0
            span_ratio = max(span_ratio_x, span_ratio_y, 0.02)
            force_mag = 0.6 * span_ratio
            axis_force = 0.05 * span_ratio
            limit_steps = max(200, int(800 * min(span_ratio, 1.2)))
            adjust_text(
                label_texts,
                x=label_points_x,
                y=label_points_y,
                ax=ax,
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.5, alpha=0.6, shrinkA=6, shrinkB=6),
                expand_text=(1.25, 1.25),
                expand_points=(1.6, 1.6),
                expand_axes=(1.04, 1.04),
                force_text=(force_mag, force_mag),
                force_axes=(axis_force, axis_force),
                lim=limit_steps,
                only_move={'text': 'xy'},
            )
        except Exception as e:
            print(f"[plot_all_tracks_in_range] adjust_text failed: {e}")
    ax.set_xlim(orig_xlim)
    ax.set_ylim(orig_ylim)

    # --- 5. 输出控制 ---
    if save_fig:
        # 使用阈值子目录
        region_slug_for_path = _current_region_key()
        thr_str = f"{do_threshold:g}".replace('.', 'p') if isinstance(do_threshold, (int,float)) else str(do_threshold)
        thr_dir = f"thr{thr_str}"
        if anomaly_min_depth is not None and anomaly_min_depth > 0:
            depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
            thr_dir += f"_depth{depth_str}m"
        output_dir = Path(plots_output_root) / region_slug_for_path / "plot_all_tracks_in_range" / thr_dir
        output_dir.mkdir(exist_ok=True, parents=True)
        base_filename = f"All_Tracks_{start_date_str}_to_{end_date_str}.png"
        save_path = output_dir / base_filename
        # 若开启跳过空图保存且确认为空（没有任何涡旋轨迹绘制），则跳过
        if skip_save_if_empty and not any_eddy_drawn:
            print(f"Skip saving empty figure for {start_date_str} to {end_date_str} (no eddies plotted).")
        else:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")
    
    if show_fig:
        plt.show()
    
    plt.close(fig)

    # 汇总输出
    interacted_eddies = []
    interacting_argo_records = []
    for r in filter(None, results):
        if r['has_interaction']:
            interacted_eddies.append(r['text_info']['text'])
        if save_interacting_argo and r.get('interacting_argo'):
            # 补充轨迹标签信息，便于后续统计
            track_label = r['text_info']['text']
            for rec in r['interacting_argo']:
                rec2 = dict(rec)
                rec2['track_label'] = track_label
                rec2['track_id'] = r.get('track_id')
                interacting_argo_records.append(rec2)

    if save_interacting_argo and interacting_argo_records:
        region_slug_for_path = _current_region_key()
        thr_str = f"{do_threshold:g}".replace('.', 'p') if isinstance(do_threshold, (int,float)) else str(do_threshold)
        thr_dir = f"thr{thr_str}"
        if anomaly_min_depth is not None and anomaly_min_depth > 0:
            depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
            thr_dir += f"_depth{depth_str}m"
        out_dir = Path(plots_output_root) / region_slug_for_path / "plot_all_tracks_in_range" / thr_dir
        out_dir.mkdir(exist_ok=True, parents=True)
        thr_str = f"{do_threshold:g}".replace('.', 'p') if isinstance(do_threshold, (int,float)) else str(do_threshold)
        depth_suffix = ''
        if anomaly_min_depth is not None and anomaly_min_depth > 0:
            depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
            depth_suffix = f"_depth{depth_str}m"
        fname_pq = out_dir / f"Interacting_Argo_{start_date_str}_to_{end_date_str}_thr{thr_str}{depth_suffix}.parquet"
        try:
            df_out = pd.DataFrame(interacting_argo_records)
            df_out.to_parquet(fname_pq, index=False)
            print(f"Interacting Argo saved to: {fname_pq}")
        except Exception as e:
            print(f"[WARN] Failed to save interacting Argo parquet: {e}")

    # 保存交互涡旋标签（每期）为 NPY，并不返回
    try:
        region_slug_for_path = _current_region_key()
        thr_str = f"{do_threshold:g}".replace('.', 'p') if isinstance(do_threshold, (int,float)) else str(do_threshold)
        thr_dir = f"thr{thr_str}"
        if anomaly_min_depth is not None and anomaly_min_depth > 0:
            depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
            thr_dir += f"_depth{depth_str}m"
        out_dir = Path(plots_output_root) / region_slug_for_path / "plot_all_tracks_in_range" / thr_dir
        out_dir.mkdir(exist_ok=True, parents=True)
        thr_str = f"{do_threshold:g}".replace('.', 'p') if isinstance(do_threshold, (int,float)) else str(do_threshold)
        depth_suffix = ''
        if anomaly_min_depth is not None and anomaly_min_depth > 0:
            depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
            depth_suffix = f"_depth{depth_str}m"
        if save_interacted_eddies and interacted_eddies:
            eddies_npy = out_dir / f"Interacted_Eddies_{start_date_str}_to_{end_date_str}_thr{thr_str}{depth_suffix}.npy"
            # 保存为标准 Unicode 字符串数组，避免 object 导致读取需 allow_pickle=True
            labels = sorted(set(interacted_eddies))
            np.save(eddies_npy, np.array(labels, dtype=str))
            print(f"Interacted eddies saved to: {eddies_npy}")
    except Exception as e:
        print(f"[WARN] Failed to save interacted eddies npy: {e}")

    return None

def _load_eddy_datasets_for_range(
    *,
    kinds: list[str],
    start_date: pd.Timestamp,
    end_date: pd.Timestamp,
    region_key: str | None = None,
    output_root: str | Path | None = None,
    include_contours: bool = False,
) -> dict:
    """根据 kinds（如 ['acs','acl']）与时间/区域，从 Parquet+Zarr 动态装载需要的涡旋轨迹列表。

    返回: dict，如 {'ACS': [(track_list, track_id), ...], 'ACL': [...]}。
    其中 track_list 为“旧版轨迹结构”的 list[list]，track_id 为真实的轨迹ID（与 parquet 的 track_id 对齐），
    便于后续懒加载等值线阶段精准定位。
    """
    region_slug = region_key or _current_region_key()
    root = _ensure_meta_tracks_root(output_root)
    region_dir = Path(root) / region_slug
    if not region_dir.exists():
        raise FileNotFoundError(f"区域目录不存在：{region_dir}")

    def _locate_daily_source(kind_l: str) -> tuple[Path, str]:
        daily_file = region_dir / f"{kind_l}_daily.parquet"
        daily_tmp_dir = region_dir / f"{kind_l}_daily_tmp"
        daily_dir = region_dir / f"{kind_l}_daily.parquet"  # 目录形式
        if daily_file.exists() and daily_file.is_file():
            return daily_file, 'file'
        if daily_tmp_dir.exists() and daily_tmp_dir.is_dir():
            return daily_tmp_dir, 'dir'
        if daily_dir.exists() and daily_dir.is_dir():
            return daily_dir, 'dir'
        raise FileNotFoundError(f"未找到 daily 数据：{daily_file} 或 {daily_tmp_dir} / {daily_dir}")

    # dateline 兼容：使用全局 lonmin/lonmax/latmin/latmax
    def _normalize_lon(val: np.ndarray | float) -> np.ndarray | float:
        return (np.asarray(val, dtype=float) + 180.0) % 360.0 - 180.0

    lon_min_cfg = float(lonmin)
    lon_max_cfg = float(lonmax)
    lat_min_cfg = float(latmin)
    lat_max_cfg = float(latmax)
    lon_min_eff = float(_normalize_lon(lon_min_cfg))
    lon_max_eff = float(_normalize_lon(lon_max_cfg))

    raw_span = abs(lon_max_cfg - lon_min_cfg)
    eff_span = (lon_max_eff - lon_min_eff) % 360.0
    is_global_lon = (raw_span >= 359.5) or (eff_span >= 359.5) or np.isclose(eff_span, 0.0, atol=1e-6)

    crosses_cfg = bool(_REGION_CFG.get('crosses_dateline'))
    simple_interval = (not crosses_cfg) and (not is_global_lon) and (lon_min_eff <= lon_max_eff)

    def _build_lon_mask(lon_vals: np.ndarray) -> np.ndarray:
        if lon_vals.size == 0 or is_global_lon:
            return np.ones(lon_vals.size, dtype=bool)
        lon_norm = _normalize_lon(lon_vals)
        if lon_min_eff <= lon_max_eff:
            return (lon_norm >= lon_min_eff) & (lon_norm <= lon_max_eff)
        return (lon_norm >= lon_min_eff) | (lon_norm <= lon_max_eff)

    def _geo_mask(df: pd.DataFrame) -> np.ndarray:
        if df.empty:
            return np.zeros(0, dtype=bool)
        lon_vals = df['center_lon'].to_numpy(dtype=float, copy=False)
        lat_vals = df['center_lat'].to_numpy(dtype=float, copy=False)
        lat_mask = (lat_vals >= latmin) & (lat_vals <= latmax)
        lon_mask = _build_lon_mask(lon_vals)
        return lat_mask & lon_mask

    out: dict[str, list] = {}
    # 预先计算时间边界（闭开区间），避免逐行日期转换：time ∈ [start_ts, end_exclusive)
    start_ts = pd.Timestamp(start_date).normalize()
    end_exclusive = pd.Timestamp(end_date).normalize() + pd.Timedelta(days=1)
    for kind in kinds:
        kind_l = kind.lower()
        if kind_l not in {'acs','acl','cs','cl'}:
            continue
        daily_source, source_type = _locate_daily_source(kind_l)

        # 收集时间+区域内的 track_id
        track_ids: set[int] = set()
        if source_type == 'file':
            pf = pq.ParquetFile(daily_source)
            for rg in range(pf.num_row_groups):
                # 先用统计信息按时间裁剪（使用 'time' 列的 min/max）
                # 同时在不跨日界线时，尝试用经纬度 min/max 做粗滤
                read_rg = True
                try:
                    schema_names = pf.schema.names
                    # time stats
                    if 'time' in schema_names:
                        t_idx = schema_names.index('time')
                        t_stats = pf.metadata.row_group(rg).column(t_idx).statistics
                        if t_stats and t_stats.has_min_max:
                            tmin = pd.to_datetime(t_stats.min)
                            tmax = pd.to_datetime(t_stats.max)
                            # 跳过与 [start_ts, end_exclusive) 无交集的行组
                            if (tmax < start_ts) or (tmin >= end_exclusive):
                                read_rg = False
                    # lat/lon stats（仅当使用简单区间时启用）
                    if read_rg:
                        try:
                            lat_idx = schema_names.index('center_lat')
                            lat_stats = pf.metadata.row_group(rg).column(lat_idx).statistics
                            if lat_stats and lat_stats.has_min_max:
                                lat_min, lat_max_ = float(lat_stats.min), float(lat_stats.max)
                                if (lat_max_ < lat_min_cfg) or (lat_min > lat_max_cfg):
                                    read_rg = False
                        except Exception:
                            pass
                    if read_rg and simple_interval:
                        try:
                            lon_idx = schema_names.index('center_lon')
                            lon_stats = pf.metadata.row_group(rg).column(lon_idx).statistics
                            if lon_stats and lon_stats.has_min_max:
                                lon_min_rg = float(lon_stats.min)
                                lon_max_rg = float(lon_stats.max)
                                lon_norm_min = float(_normalize_lon(lon_min_rg))
                                lon_norm_max = float(_normalize_lon(lon_max_rg))
                                if lon_norm_min <= lon_norm_max:
                                    if (lon_norm_max < lon_min_eff) or (lon_norm_min > lon_max_eff):
                                        read_rg = False
                                else:
                                    # 无法可靠判断跨日期线的行组，保留
                                    pass
                        except Exception:
                            pass
                except Exception:
                    pass
                if not read_rg:
                    continue

                table = pf.read_row_group(rg, columns=['track_id','time','center_lon','center_lat'])
                df = table.to_pandas()
                # 直接用时间戳闭开区间过滤，避免生成 'date'
                mask_time = (df['time'] >= start_ts) & (df['time'] < end_exclusive)
                mask_geo = _geo_mask(df)
                df_sel = df.loc[mask_time.to_numpy() & mask_geo]
                if not df_sel.empty:
                    track_ids.update(df_sel['track_id'].astype(int).unique().tolist())
        else:
            # 目录形式：逐个 part 文件按 row-group 读取与裁剪，避免整文件加载与逐行日期转换
            for part in sorted(daily_source.glob('*.parquet')):
                try:
                    ppf = pq.ParquetFile(part)
                except Exception:
                    continue
                schema_names = ppf.schema.names
                for rg in range(ppf.num_row_groups):
                    read_rg = True
                    try:
                        # time stats
                        if 'time' in schema_names:
                            t_idx = schema_names.index('time')
                            t_stats = ppf.metadata.row_group(rg).column(t_idx).statistics
                            if t_stats and t_stats.has_min_max:
                                tmin = pd.to_datetime(t_stats.min)
                                tmax = pd.to_datetime(t_stats.max)
                                if (tmax < start_ts) or (tmin >= end_exclusive):
                                    read_rg = False
                        # lat/lon stats prune（不跨日界线时）
                        if read_rg:
                            try:
                                lat_idx = schema_names.index('center_lat')
                                lat_stats = ppf.metadata.row_group(rg).column(lat_idx).statistics
                                if lat_stats and lat_stats.has_min_max:
                                    lat_min, lat_max_ = float(lat_stats.min), float(lat_stats.max)
                                    if (lat_max_ < latmin) or (lat_min > latmax):
                                        read_rg = False
                            except Exception:
                                pass
                        if read_rg and simple_interval:
                            try:
                                lon_idx = schema_names.index('center_lon')
                                lon_stats = ppf.metadata.row_group(rg).column(lon_idx).statistics
                                if lon_stats and lon_stats.has_min_max:
                                    lon_min_rg = float(lon_stats.min)
                                    lon_max_rg = float(lon_stats.max)
                                    lon_norm_min = float(_normalize_lon(lon_min_rg))
                                    lon_norm_max = float(_normalize_lon(lon_max_rg))
                                    if lon_norm_min <= lon_norm_max:
                                        if (lon_norm_max < lon_min_eff) or (lon_norm_min > lon_max_eff):
                                            read_rg = False
                                # 若归一化后出现 wrap（lon_norm_min > lon_norm_max），不做截断
                            except Exception:
                                pass
                    except Exception:
                        pass
                    if not read_rg:
                        continue
                    tbl = ppf.read_row_group(rg, columns=['track_id','time','center_lon','center_lat'])
                    df = tbl.to_pandas()
                    mask_time = (df['time'] >= start_ts) & (df['time'] < end_exclusive)
                    mask_geo = _geo_mask(df)
                    df_sel = df.loc[mask_time.to_numpy() & mask_geo]
                    if not df_sel.empty:
                        track_ids.update(df_sel['track_id'].astype(int).unique().tolist())

        # 组装旧版轨迹结构（批量提取以减少重复 I/O；失败时回退逐个提取）
        tracks: list = []
        if track_ids:
            try:
                batch = find_track(
                    kind_l,
                    sorted(track_ids),
                    region=region_key,
                    output_root=output_root,
                    include_contours=include_contours,
                    return_list=True
                )
                # 批量返回为 {track_id: list}
                tracks = [(batch[tid], int(tid)) for tid in sorted(batch.keys())]
            except Exception:
                for tid in sorted(track_ids):
                    try:
                        track_list = find_track(
                            kind_l,
                            tid,
                            region=region_key,
                            output_root=output_root,
                            include_contours=include_contours,
                            return_list=True
                        )
                        tracks.append((track_list, int(tid)))
                    except Exception:
                        continue

        out[kind_l.upper()] = tracks

    return out

def worker_wrapper(args: tuple):
    """multiprocessing worker 包装函数。
    参数:
        args: (start_date_str, end_date_str, plot_unrelated_eddies, skip_save_if_empty, show_labels, save_interacting_argo, save_interacted_eddies, do_threshold, anomaly_min_depth)
    返回: None（结果由子函数写盘保存）。
    """
    start_d, end_d, unrelated_flag, skip_empty, show_labels, save_interacting_argo_flag, save_eddies_flag, do_thr, anom_depth = args
    try:
        plot_all_tracks_in_range(
            start_date_str=start_d,
            end_date_str=end_d,
            plot_unrelated_eddies=unrelated_flag,
            save_fig=True,
            skip_save_if_empty=skip_empty,
            show_labels=show_labels,
            show_fig=False,
            save_interacting_argo=bool(save_interacting_argo_flag),
            save_interacted_eddies=bool(save_eddies_flag),
            do_threshold=do_thr,
            anomaly_min_depth=anom_depth,
        )
        return None
    except Exception as e:
        print(f"!!! ERROR processing period {start_d}: {e}")
        return None

def run_batch_plotting_multiprocessing(
    start_date_str: str,
    end_date_str: str,
    eddy_datasets: dict,
    num_workers: int,
    plot_unrelated_eddies: bool = False,
    skip_save_if_empty: bool = True,
    show_labels: bool | None = None,
    do_threshold: float | None = None,
    anomaly_min_depth: float | None = None,
    save_interacted_eddies: bool = False,
    save_interacting_argo: bool = False,
):
    """
    (批处理控制器) 使用multiprocessing启动一个扁平化的并行绘图作业，并带有进度条。
    
    功能:
        1. 创建一个横跨指定数量核心的进程池。
        2. 使用initializer高效地将涡旋数据共享给所有工作进程。
        3. 按月份切分任务，并使用tqdm实时显示处理进度。

    参数:
        start_date_str (str): 批处理的开始日期 'YYYY-MM-DD'。
        end_date_str (str): 批处理的结束日期 'YYYY-MM-DD'。
        eddy_datasets (dict): 【已加载】的、将被共享给所有进程的涡旋数据集。
        num_workers (int): 需要启动的并行工作进程数（核心数）。
        plot_unrelated_eddies (bool): 是否在批处理中绘制无关涡旋。
        skip_save_if_empty (bool): 批处理默认 True（空图不保存）；当传入 False 时，将把 False 透传至 worker，空图也会保存。
        show_labels (bool | None): 是否绘制轨迹标签；None 表示使用智能判定（全球且 plot_unrelated_eddies=True 时默认 False）。
        do_threshold (float | None): ΔDO 阈值；None 使用配置默认。
        anomaly_min_depth (float | None): 最小深度阈值；None 使用配置默认；≤0 表示不限制。
        save_interacted_eddies (bool): True 时各月份写出 `Interacted_Eddies_*.npy` 并在末尾汇总；默认 False。
        save_interacting_argo (bool): True 时各月份写出 `Interacting_Argo_*.parquet` 并在末尾聚合；默认 False。

        输出:
                - 每月图像写入 `plot_outputs/<region>/plot_all_tracks_in_range/<thr_dir>/`。
                - 当 save_interacted_eddies=True：整期交互涡旋标签汇总保存为带阈值后缀的 NPY：
                    `plot_outputs/<region>/plot_all_tracks_in_range/<thr_dir>/eddy_list_thr{thr}[_depth{d}m].npy`。
                - 当 save_interacting_argo=True 时，额外保存整期交互 Argo 汇总（Parquet）：
                    `plot_outputs/<region>/plot_all_tracks_in_range/<thr_dir>/interacting_argo_all_thr{thr}[_depth{d}m].parquet`。
        清理:
                - 批处理开始前仅清空对应阈值子目录 `<thr_dir>`，不会清空整个 `plot_all_tracks_in_range` 目录。
    """
    print("="*60)
    print("      Multiprocessing Batch Plotting with Progress Bar      ")
    print("="*60)
    print(f"Strategy: Creating a single pool of {num_workers} workers.")
    
    # --- 1. 创建按月切分的任务列表 ---
    month_starts = pd.date_range(start=start_date_str, end=end_date_str, freq='MS')
    region_slug_for_path = _current_region_key()
    # 计算生效阈值并确定阈值子目录名
    eff_do_thr = do_threshold if do_threshold is not None else _default_delta_do_threshold
    eff_anom_depth = anomaly_min_depth if anomaly_min_depth is not None else _cfg_anomaly_min_depth
    thr_str = f"{eff_do_thr:g}".replace('.', 'p') if isinstance(eff_do_thr, (int, float)) else str(eff_do_thr)
    thr_dir = f"thr{thr_str}"
    if eff_anom_depth is not None and eff_anom_depth > 0:
        depth_str = f"{eff_anom_depth:g}".replace('.', 'p')
        thr_dir += f"_depth{depth_str}m"
    effective_show_labels = show_labels
    if show_labels is None:
        if plot_unrelated_eddies and region_slug_for_path.lower() in {'global', 'world', 'global_ocean', 'all'}:
            effective_show_labels = False
        else:
            effective_show_labels = True

    tasks = [
        (
            start_date.strftime('%Y-%m-%d'),
            (start_date + pd.tseries.offsets.MonthEnd(1)).strftime('%Y-%m-%d'),
            plot_unrelated_eddies,
            skip_save_if_empty,
            effective_show_labels,
            save_interacting_argo,
            save_interacted_eddies,
            eff_do_thr,
            eff_anom_depth,
        )
        for start_date in month_starts
    ]
    
    print(f"[*] Created {len(tasks)} monthly plotting tasks to be processed by {num_workers} cores.")

    # 在批量任务开始前清理输出目录
    base_output_dir = Path(plots_output_root) / region_slug_for_path / "plot_all_tracks_in_range"
    output_dir = base_output_dir / thr_dir
    try:
        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"[*] Cleared threshold-specific output directory: {output_dir}")
    except Exception as e:
        print(f"[WARN] Failed to clear output directory {output_dir}: {e}")
    
    # --- 2. 启动进程池并执行任务 ---
    start_time_total = tm.time()
    
    # 使用 initializer 来高效地传递一次大的涡旋数据（子进程写盘，不收集返回）
    # maxtasksperchild=1: 强制每个子进程处理完一个任务后重启，彻底释放 Matplotlib 内存和文件句柄，防止死锁和内存泄漏
    with multiprocessing.Pool(processes=num_workers, initializer=init_worker, initargs=(eddy_datasets,), maxtasksperchild=1) as pool:
        print("[*] Processing tasks...")
        for _ in tqdm(pool.imap_unordered(worker_wrapper, tasks), total=len(tasks)):
            pass
        
    end_time_total = tm.time()
    total_duration_minutes = (end_time_total - start_time_total) / 60
    
    print("\n" + "="*60)
    print("--- All Plotting Tasks Have Finished ---")
    print(f"Total execution time: {total_duration_minutes:.2f} minutes.")
    # 聚合交互涡旋标签（从每月 NPY 汇总）
    interacted_labels: set[str] = set()
    try:
        monthly_eddy_files = sorted(output_dir.glob("Interacted_Eddies_*.npy"))
        for f in monthly_eddy_files:
            try:
                # 优先以非 pickle 方式读取；兼容历史 object 数组文件回退到允许 pickle
                try:
                    arr = np.load(f)
                except Exception:
                    arr = np.load(f, allow_pickle=True)
                interacted_labels.update([str(x) for x in np.asarray(arr).ravel().tolist()])
            except Exception:
                pass
        unique_interacted = sorted(interacted_labels)
        print(f"Total interacted eddies: {len(unique_interacted)}")
        if unique_interacted:
            preview = ", ".join(unique_interacted[:20])
            print(f"Sample (first 20): {preview}{' ...' if len(unique_interacted)>20 else ''}")
        # 不再保存未带阈值后缀的 eddy_list.npy，改为仅保存带阈值后缀版本（见后续输出）
        # 另存带阈值后缀的汇总版本（从文件名解析 thr 与 depth）
        thr_tag = None
        depth_tag = None
        if monthly_eddy_files:
            m = re.search(r"_thr([A-Za-z0-9p]+)(?:_depth([A-Za-z0-9p]+)m)?\\.npy$", monthly_eddy_files[0].name)
            if m:
                thr_tag = m.group(1)
                depth_tag = m.group(2)
        if thr_tag:
            suffix = f"_thr{thr_tag}"
            if depth_tag:
                suffix += f"_depth{depth_tag}m"
            eddy_list_suffixed = output_dir / f"eddy_list{suffix}.npy"
            # 使用标准 Unicode 字符串数组，避免后续读取需要 allow_pickle
            np.save(eddy_list_suffixed, np.array(unique_interacted, dtype=str))
            print(f"Eddy list (with thresholds) saved to: {eddy_list_suffixed}")
    except Exception as e:
        print(f"[WARN] Failed to aggregate eddy labels: {e}")

    # 聚合交互 Argo（Parquet）
    if save_interacting_argo:
        try:
            monthly_argo_files = sorted(output_dir.glob("Interacting_Argo_*.parquet"))
            df_parts = []
            for f in monthly_argo_files:
                try:
                    dfp = pd.read_parquet(f)
                    if isinstance(dfp, pd.DataFrame) and not dfp.empty:
                        df_parts.append(dfp)
                except Exception:
                    pass
            if df_parts:
                df_all = pd.concat(df_parts, ignore_index=True)
                if 'Profile_number' in df_all.columns:
                    if 'date' in df_all.columns:
                        df_all.sort_values(by=['Profile_number', 'date'], inplace=True)
                        df_all = df_all.drop_duplicates(subset=['Profile_number', 'date', 'track_label'], keep='first')
                    elif all(col in df_all.columns for col in ('Year', 'Month', 'Day')):
                        temp_date = pd.to_datetime(df_all[['Year', 'Month', 'Day']])
                        df_all = df_all.assign(_date=temp_date)
                        df_all.sort_values(by=['Profile_number', '_date'], inplace=True)
                        df_all = df_all.drop_duplicates(subset=['Profile_number', '_date', 'track_label'], keep='first')
                        df_all.drop(columns=['_date'], inplace=True)
                thr_str2 = f"{eff_do_thr:g}".replace('.', 'p') if isinstance(eff_do_thr, (int,float)) else str(eff_do_thr)
                depth_suffix2 = ''
                if eff_anom_depth is not None and eff_anom_depth > 0:
                    depth_str2 = f"{eff_anom_depth:g}".replace('.', 'p')
                    depth_suffix2 = f"_depth{depth_str2}m"
                argo_parquet = output_dir / f"interacting_argo_all_thr{thr_str2}{depth_suffix2}.parquet"
                df_all.to_parquet(argo_parquet, index=False)
                print(f"Interacting Argo (all) saved to: {argo_parquet}")
        except Exception as e:
            print(f"[WARN] Failed to aggregate/save interacting Argo parquet: {e}")

def load_combined_eddy_list(
    region: str | None = None,
    thr_dir: str | Path | None = None,
    plots_root: str | Path | None = None,
    include_monthly: bool = True,
    deduplicate: bool = True,
    save_path: str | Path | None = None,
) -> list[str]:
    """读取并合并单个阈值目录下的 `eddy_list*.npy`（及可选月度文件）。

    参数:
        region: 区域 slug；None 时复用当前 `switch_region` 的配置。
        thr_dir: 阈值子目录名称或路径；None 时自动扫描 `thr*` 目录并要求只存在一个候选目录。
        plots_root: 可选自定义 `plot_outputs` 根路径；None 使用配置 `plots_output_root`。
        include_monthly: True 时若目录缺少汇总 `eddy_list*.npy`，会退回加载 `Interacted_Eddies_*.npy`（逐月文件）。
        deduplicate: True 返回去重并排序后的唯一列表；False 按读取顺序返回。
        save_path: 可选输出路径；相对路径时会解析到目标阈值目录中。

    返回:
        list[str]: 聚合后的涡旋标签列表；当没有匹配文件时返回空列表。
    """
    region_slug = region or _current_region_key()
    plots_base = Path(plots_root) if plots_root is not None else Path(plots_output_root)
    base_dir = plots_base / region_slug / "plot_all_tracks_in_range"
    if not base_dir.exists():
        raise FileNotFoundError(f"Plot outputs directory not found: {base_dir}")

    if thr_dir is not None:
        target_dir = Path(thr_dir)
        if not target_dir.is_absolute():
            target_dir = base_dir / target_dir
        if not target_dir.is_dir():
            raise FileNotFoundError(f"阈值目录不存在：{target_dir}")
    else:
        candidate_dirs = sorted(
            p for p in base_dir.iterdir() if p.is_dir() and p.name.lower().startswith('thr')
        )
        if not candidate_dirs:
            print(f"[load_combined_eddy_list] No threshold directories found under {base_dir}.")
            return []
        if len(candidate_dirs) > 1:
            raise ValueError(f"检测到多个阈值目录，请通过 thr_dir 指定其一：{[p.name for p in candidate_dirs]}")
        target_dir = candidate_dirs[0]

    collected_labels: list[str] = []
    npy_files = sorted(target_dir.glob('eddy_list*.npy'))
    if not npy_files and include_monthly:
        npy_files = sorted(target_dir.glob('Interacted_Eddies_*.npy'))
    for npy_path in npy_files:
        try:
            try:
                arr = np.load(npy_path)
            except Exception:
                arr = np.load(npy_path, allow_pickle=True)
            flattened = np.asarray(arr).ravel().tolist()
            labels = [str(item) for item in flattened if str(item)]
            collected_labels.extend(labels)
        except Exception as exc:
            print(f"[load_combined_eddy_list] Failed to read {npy_path}: {exc}")

    if not collected_labels:
        print("[load_combined_eddy_list] No eddy labels were loaded.")
        return []

    combined = sorted(set(collected_labels)) if deduplicate else collected_labels

    if save_path is None:
        save_path = target_dir / "eddy_list_combined.npy"
    else:
        save_path = Path(save_path)
        if not save_path.is_absolute():
            save_path = target_dir / save_path
    save_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(save_path, np.array(combined, dtype=str))
    print(f"[load_combined_eddy_list] Combined eddy list saved to: {save_path}")

    return combined

def plot_argo_hotspots(
    start_year: int,
    end_year: int,
    do_threshold: float | None = None,
    salinity_threshold: float | None = None,
    temperature_threshold: float | None = None,
    depth_interval: float | None = None,
    depth_merge_tolerance: float | None = None,
    duplicate_depth_strategy: str | None = None,
    anomaly_min_depth: float | None = None,
    plot_unrelated_argo: bool = True,
    fix_delta_do_colorbar: bool = True,
    delta_do_cbar_min: float = 10.0,
    delta_do_cbar_max: float = 100.0,
    delta_do_cbar_ticks: list | None = None,
    save_fig: bool = False,
    show_fig: bool = True,
    save_data: bool = True,
    dask_scheduler: str | None = None,
    dask_workers: int | None = None,
    dask_memory_limit: str | None = None,
    use_interacting_argo: bool = False,
    split_plots: bool = False
):
    """以 ΔDO 异常方法绘制多年期 Argo 异常分布。

    流程：
      1. 逐年加载 Argo 年度数据并合并；可利用全局 lonmin/latmin/lonmax/latmax 做空间裁剪；
      2. 用 calculate_delta_do 检测每个剖面潜在 ΔDO 异常；
      3. 若 anomaly_min_depth > 0，则按该阈值过滤；
      4. 每个剖面保留最大 ΔDO 一条记录；
      5. 绘制 ΔDO 异常散点（可选固定色标范围），并可选绘制所有匹配剖面基线位置（空心灰圈）。

    参数:
        start_year / end_year: 年度范围（闭区间）。
        do_threshold / salinity_threshold / temperature_threshold / depth_interval / depth_merge_tolerance / duplicate_depth_strategy:
            传给 calculate_delta_do；当为 None 时从 processing.yml 读取默认值。
        anomaly_min_depth: 仅保留异常深度 >= 该值；≤0 不限制；None 表示使用 processing.yml 的 anomaly_min_depth。
        plot_unrelated_argo: 是否绘制所有匹配剖面基线（被筛掉或无异常的）。
        fix_delta_do_colorbar: 是否固定 ΔDO 色标范围。
        delta_do_cbar_min / delta_do_cbar_max / delta_do_cbar_ticks: 色标范围与刻度配置。
        save_fig / show_fig: 输出控制。
        save_data (bool): True 时保存 anomalies 为 Parquet；False 不保存数据，输出路径固定为 `plot_outputs/<region>/plot_argo_hotspots/`。
        dask_scheduler (str | None): Dask 调度器，'threads'|'processes'|'single'，None 默认 'processes'。
        dask_workers (int | None): Dask worker 数量；None 自动取 min(年度数, CPU)。
        dask_memory_limit (str | None): LocalCluster 模式下的单 worker 内存限制，如 '4GB'；None 不限制。
        use_interacting_argo (bool): 是否读取 run_batch_plotting_multiprocessing 生成的交互 Argo 文件，
                                          并在图中区分交互/非交互 Argo，同时统计交互比例。
                                          默认 False。
        split_plots (bool): 若为 True 且 use_interacting_argo=True，则将交互（菱形）与非交互（圆形）异常点分别绘制在两张图中。
                            默认 False。

    输出:
        - 图像（可选）：`plot_outputs/<region>/plot_argo_hotspots/Argo_DeltaDO_Hotspots_*.png`
        - 异常数据（Parquet，可选）：`plot_outputs/<region>/plot_argo_hotspots/anomalies_{start}_{end}_thr{thr}[_depth{d}m].parquet`
    """
    # 从配置读取默认值
    if do_threshold is None:
        do_threshold = _default_delta_do_threshold
    if salinity_threshold is None:
        salinity_threshold = _default_salinity_threshold
    if temperature_threshold is None:
        temperature_threshold = _default_temperature_threshold
    if depth_interval is None:
        depth_interval = _default_depth_interval
    if depth_merge_tolerance is None:
        depth_merge_tolerance = _default_depth_merge_tolerance
    if duplicate_depth_strategy is None:
        duplicate_depth_strategy = _default_duplicate_depth_strategy
    if anomaly_min_depth is None:
        anomaly_min_depth = _cfg_anomaly_min_depth

    # --- 尝试加载交互 Argo 文件（若启用） ---
    interacting_argo_ids: set[int] = set()
    if use_interacting_argo:
        region_slug_for_path = _current_region_key()
        eff_do_thr = do_threshold
        eff_anom_depth = anomaly_min_depth
        
        thr_str = f"{eff_do_thr:g}".replace('.', 'p') if isinstance(eff_do_thr, (int, float)) else str(eff_do_thr)
        thr_dir_name = f"thr{thr_str}"
        depth_suffix = ""
        if eff_anom_depth is not None and eff_anom_depth > 0:
            d_str = f"{eff_anom_depth:g}".replace('.', 'p')
            thr_dir_name += f"_depth{d_str}m"
            depth_suffix = f"_depth{d_str}m"
            
        interacting_file = Path(plots_output_root) / region_slug_for_path / "plot_all_tracks_in_range" / thr_dir_name / f"interacting_argo_all_thr{thr_str}{depth_suffix}.parquet"
        
        if interacting_file.exists():
            print(f"[*] Loading interacting Argo from: {interacting_file}")
            try:
                df_int = pd.read_parquet(interacting_file)
                if 'Profile_number' in df_int.columns:
                    interacting_argo_ids = set(df_int['Profile_number'].unique())
                print(f"[*] Loaded {len(interacting_argo_ids)} unique interacting profiles.")
            except Exception as e:
                print(f"[WARN] Failed to read interacting Argo file: {e}")
        else:
            print(f"[WARN] Interacting Argo file not found: {interacting_file}")

    print(f"--- Building Argo ΔDO Anomaly Map {start_year}-{end_year} ---")

    # --- 按年份加载策略：串行或并行 ---
    years = list(range(start_year, end_year + 1))
    # Worker 参数打包列表（避免闭包，支持多进程 pickling）
    # 捕获当前区域边界到参数中，避免在 Dask 子进程中重新导入模块后回退默认区域
    current_lon_min = float(lonmin)
    current_lon_max = float(lonmax)
    current_lat_min = float(latmin)
    current_lat_max = float(latmax)
    worker_args_list = [
        (
            y,
            depth_interval,
            do_threshold,
            salinity_threshold,
            temperature_threshold,
            anomaly_min_depth,
            depth_merge_tolerance,
            duplicate_depth_strategy,
            current_lon_min,
            current_lon_max,
            current_lat_min,
            current_lat_max,
        )
        for y in years
    ]

    # Dask-only backend
    # 结果累积容器
    baselines_list: list[pd.DataFrame] = []
    anomalies_list: list[pd.DataFrame] = []

    sched = dask_scheduler or 'processes'
    worker_count = dask_workers or min(len(years), os.cpu_count() or 1)
    print(f"[*] Hotspots Dask mode: scheduler={sched}, workers={worker_count}, years={len(years)}")
    cluster = None
    client = None
    if sched == 'processes':
        try:
            cluster = LocalCluster(n_workers=worker_count, threads_per_worker=1,
                                   memory_limit=dask_memory_limit or None,
                                   silence_logs='CRITICAL')
            client = Client(cluster)
        except Exception as e:
            print(f"[WARN] 创建 LocalCluster 失败，改用 dask.compute scheduler='{sched}': {e}")
            cluster = None
            client = None

    delayed_tasks = [delayed(_hotspot_year_worker)(args) for args in worker_args_list]
    # 显示进度：优先使用 distributed.as_completed + tqdm；否则回退到 ProgressBar
    if client is not None:
        try:
            futures = client.compute(delayed_tasks)
            for fut in tqdm(as_completed(futures), total=len(futures), desc="hotspots(dask)"):
                try:
                    baseline_year, anomalies_year = fut.result()
                except Exception as e:
                    print(f"[hotspots] Dask task failed: {e}")
                    continue
                if not baseline_year.empty:
                    baselines_list.append(baseline_year)
                if not anomalies_year.empty:
                    anomalies_list.append(anomalies_year)
        finally:
            if client:
                client.close()
            if cluster:
                cluster.close()
    else:
        # 非 distributed client 情况：使用 dask.diagnostics.ProgressBar
        try:
            with ProgressBar():
                results = compute(*delayed_tasks, scheduler=sched)
        except Exception:
            results = compute(*delayed_tasks, scheduler=sched)
        for baseline_year, anomalies_year in results:
            if not baseline_year.empty:
                baselines_list.append(baseline_year)
            if not anomalies_year.empty:
                anomalies_list.append(anomalies_year)

    if not baselines_list and not anomalies_list:
        print("No data loaded after filtering; abort.")
        return None

    baseline_profiles = pd.concat(baselines_list, ignore_index=True) if baselines_list else pd.DataFrame()
    anomalies = pd.concat(anomalies_list, ignore_index=True) if anomalies_list else pd.DataFrame()

    # 基线与异常已在加载阶段完成区域过滤，这里只再检查是否为空
    if baseline_profiles.empty:
        print("No baseline profiles after filtering.")
    if anomalies.empty:
        print("No ΔDO anomalies detected.")

    # 绘图
    crosses_dateline = bool(_REGION_CFG.get('crosses_dateline') and (lonmax < lonmin))
    central_lon = 180 if crosses_dateline else 0
    data_crs = ccrs.PlateCarree()
    map_crs = ccrs.PlateCarree(central_longitude=central_lon)

    # 准备数据分组
    anom_interacting = pd.DataFrame()
    anom_others = pd.DataFrame()
    scatter_kwargs = {}
    
    if not anomalies.empty:
        if fix_delta_do_colorbar:
            scatter_kwargs.update(dict(vmin=delta_do_cbar_min, vmax=delta_do_cbar_max))

        anom_others = anomalies.copy()
        if use_interacting_argo and interacting_argo_ids:
            is_interacting = anomalies['Profile_number'].isin(interacting_argo_ids)
            anom_interacting = anomalies[is_interacting].copy()
            anom_others = anomalies[~is_interacting].copy()
            
            # 统计输出
            total_anom = len(anomalies)
            count_int = len(anom_interacting)
            pct = (count_int / total_anom * 100) if total_anom > 0 else 0.0
            print(f"[Plot Info] Anomalies: {total_anom}, Interacting: {count_int} ({pct:.1f}%)")

    # 定义绘图任务
    plots_to_generate = []
    if split_plots and use_interacting_argo:
        plots_to_generate.append({
            'name': 'interacting',
            'title_extra': ' (Interacting)',
            'file_suffix': '_interacting',
            'data_list': [
                {'data': anom_interacting, 'marker': 'D', 'edgecolor': 'blue', 'label': f'Interacting (ΔDO ≥ {do_threshold})', 's': 100, 'zorder': 4}
            ]
        })
        plots_to_generate.append({
            'name': 'non_interacting',
            'title_extra': ' (Non-interacting)',
            'file_suffix': '_non_interacting',
            'data_list': [
                {'data': anom_others, 'marker': 'o', 'edgecolor': 'black', 'label': f'Non-interacting (ΔDO ≥ {do_threshold})', 's': 60, 'zorder': 3}
            ]
        })
    else:
        # 合并模式
        combined_data = []
        if not anom_others.empty:
            label_str = f'ΔDO ≥ {do_threshold} μmol kg⁻¹'
            if use_interacting_argo and interacting_argo_ids:
                label_str = f'Non-interacting (ΔDO ≥ {do_threshold})'
            combined_data.append({'data': anom_others, 'marker': 'o', 'edgecolor': 'black', 'label': label_str, 's': 60, 'zorder': 3})
        
        if not anom_interacting.empty:
            combined_data.append({'data': anom_interacting, 'marker': 'D', 'edgecolor': 'blue', 'label': f'Interacting (ΔDO ≥ {do_threshold})', 's': 100, 'zorder': 4})
            
        plots_to_generate.append({
            'name': 'combined',
            'title_extra': '',
            'file_suffix': '',
            'data_list': combined_data
        })

    for p_cfg in plots_to_generate:
        fig = plt.figure(figsize=(40, 30))
        ax = fig.add_subplot(1, 1, 1, projection=map_crs)
        
        depth_title = (
            f' (depth ≥ {anomaly_min_depth} m)'
            if anomaly_min_depth is not None and anomaly_min_depth > 0 else ''
        )
        thr_title = f' (ΔDO ≥ {do_threshold:g} μmol kg⁻¹)'
        ax.set_title(f'Argo ΔDO Anomalies {start_year}-{end_year}{thr_title}{depth_title}{p_cfg["title_extra"]}', fontsize=20)

        # Basemap features
        base_ocean = _BASEMAP_COLORS['ocean']
        base_land = _BASEMAP_COLORS['land']
        coast_color = _BASEMAP_COLORS['coastline']
        grid_color = _BASEMAP_COLORS['grid']
        ax.set_facecolor(base_ocean)
        ax.add_feature(cfeature.OCEAN, facecolor=base_ocean, zorder=0)
        ax.add_feature(cfeature.LAND, facecolor=base_land, edgecolor=coast_color, linewidth=0.5, zorder=0)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.7, edgecolor=coast_color, zorder=1)
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, color=grid_color, alpha=0.45, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

        # 设定范围（处理跨日界线）
        lon_extent_min = lonmin
        lon_extent_max = lonmax
        if crosses_dateline and lon_extent_max < lon_extent_min:
            lon_extent_max += 360
        ax.set_extent([lon_extent_min, lon_extent_max, latmin, latmax], crs=data_crs)

        if plot_unrelated_argo and not baseline_profiles.empty:
            ax.scatter(
                baseline_profiles['Longitude'], baseline_profiles['Latitude'],
                facecolors='none', edgecolors='gray', linewidths=0.7, s=25,
                label='All Argo Profiles (baseline)', zorder=2, transform=data_crs
            )

        sc = None
        has_anom_plot = False
        for d_cfg in p_cfg['data_list']:
            data = d_cfg['data']
            if data.empty: continue
            has_anom_plot = True
            
            sc_curr = ax.scatter(
                data['Longitude'], data['Latitude'],
                c=data['delta_do'], cmap='Reds', s=d_cfg['s'],
                marker=d_cfg['marker'],
                edgecolors=d_cfg['edgecolor'], linewidths=0.5 if d_cfg['marker']=='o' else 1.0,
                label=d_cfg['label'], zorder=d_cfg['zorder'],
                transform=data_crs,
                **scatter_kwargs
            )
            sc = sc_curr

        if not has_anom_plot and anomalies.empty:
             ax.text(0.5, 0.5, 'No ΔDO anomalies', transform=ax.transAxes, ha='center', va='center', fontsize=24, color='red')

        if sc is not None:
            cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.05)
            cbar.set_label('ΔDO / μmol·kg⁻¹', fontsize=20); cbar.ax.tick_params(labelsize=14)
            if fix_delta_do_colorbar:
                if delta_do_cbar_ticks is not None:
                    cbar.set_ticks(delta_do_cbar_ticks)
                else:
                    rng = delta_do_cbar_max - delta_do_cbar_min
                    if rng > 30:
                        mid = (delta_do_cbar_min + delta_do_cbar_max)/2
                        cbar.set_ticks([delta_do_cbar_min, mid, delta_do_cbar_max])
                    else:
                        cbar.set_ticks([delta_do_cbar_min, delta_do_cbar_max])

        ax.legend(fontsize=18, loc='upper left')

        if save_fig:
            region_slug_for_path = _current_region_key()
            out_dir = Path(plots_output_root) / region_slug_for_path / "plot_argo_hotspots"
            out_dir.mkdir(exist_ok=True, parents=True)
            thr_str = f"{do_threshold:g}".replace('.', 'p')
            depth_suffix = ''
            if anomaly_min_depth is not None and anomaly_min_depth > 0:
                depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
                depth_suffix = f"_depth{depth_str}m"
            fname = out_dir / f"Argo_DeltaDO_Hotspots_{start_year}_{end_year}_thr{thr_str}{depth_suffix}{p_cfg['file_suffix']}.png"
            plt.savefig(fname, dpi=300, bbox_inches='tight')
            print(f"Figure saved: {fname}")
        if show_fig:
            plt.show()
        plt.close(fig)

    # 保存 anomalies 为 Parquet（高效压缩存储）
    if save_data and not anomalies.empty:
        region_slug_for_path = _current_region_key()
        out_dir = Path(plots_output_root) / region_slug_for_path / "plot_argo_hotspots"
        out_dir.mkdir(exist_ok=True, parents=True)
        thr_str = f"{do_threshold:g}".replace('.', 'p')
        depth_suffix = ''
        if anomaly_min_depth is not None and anomaly_min_depth > 0:
            depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
            depth_suffix = f"_depth{depth_str}m"
        pq_path = out_dir / f"anomalies_{start_year}_{end_year}_thr{thr_str}{depth_suffix}.parquet"
        try:
            anomalies.to_parquet(pq_path, index=False)
            print(f"Anomalies saved to: {pq_path}")
        except Exception as e:
            print(f"[WARN] Failed to save anomalies parquet: {e}")

    return None

def _hotspot_year_worker(args: tuple) -> tuple[pd.DataFrame, pd.DataFrame]:
    """模块级 worker，支持 multiprocessing pickling。

    参数 args: (
        year,
        depth_interval,
        do_threshold,
        salinity_threshold,
        temperature_threshold,
        anomaly_min_depth,
        depth_merge_tolerance,
        duplicate_depth_strategy,
        lon_min_bound,
        lon_max_bound,
        lat_min_bound,
        lat_max_bound,
    )
    返回: (baseline_df, anomalies_df)
    baseline_df: 每个剖面第一条记录的基本信息
    anomalies_df: 该年筛选出的 ΔDO 异常（每剖面保留最大 delta_do 一条）
    """
    (
        year,
        depth_interval,
        do_threshold,
        salinity_threshold,
        temperature_threshold,
        anomaly_min_depth,
        depth_merge_tolerance,
        duplicate_depth_strategy,
        lon_min_bound,
        lon_max_bound,
        lat_min_bound,
        lat_max_bound,
    ) = args
    try:
        df_y = load_argo_data(year=year)
    except FileNotFoundError:
        print(f"[hotspots] Missing year {year}, skip.")
        return pd.DataFrame(), pd.DataFrame()
    except Exception as e:
        print(f"[hotspots] Error loading year {year}: {e}")
        return pd.DataFrame(), pd.DataFrame()
    # 地理过滤
    lon_vals = df_y['Longitude'].to_numpy(dtype=float, copy=False)
    lat_vals = df_y['Latitude'].to_numpy(dtype=float, copy=False)
    lon_mask = _region_lon_mask(lon_vals, lon_min_bound, lon_max_bound)
    lat_mask = (lat_vals >= lat_min_bound) & (lat_vals <= lat_max_bound)
    geo_mask = lon_mask & lat_mask
    df_geo = df_y[geo_mask].copy()
    if df_geo.empty:
        return pd.DataFrame(), pd.DataFrame()
    # baseline
    baseline = (
        df_geo.sort_values(['Profile_number','Depth'])
        .groupby('Profile_number', as_index=False)
        .first()[['Profile_number','Longitude','Latitude','Year','Month','Day']]
    )
    anomalies_year = calculate_delta_do(
        df_geo,
        depth_interval=depth_interval,
        do_threshold=do_threshold,
        salinity_threshold=salinity_threshold,
        temperature_threshold=temperature_threshold,
        anomaly_min_depth=anomaly_min_depth,
        depth_merge_tolerance=depth_merge_tolerance,
        duplicate_depth_strategy=duplicate_depth_strategy,
        remove_outliers=True,
        verbose=False
    )
    if anomalies_year.empty:
        return baseline, pd.DataFrame()
    anomalies_year = (
        anomalies_year.sort_values('delta_do', ascending=False)
        .drop_duplicates(subset='Profile_number', keep='first')
    )
    needed_cols = [c for c in ['Profile_number','Longitude','Latitude','Year','Month','Day','depth','delta_do','do_value'] if c in anomalies_year.columns]
    anomalies_year = anomalies_year[needed_cols]
    return baseline, anomalies_year

def calculate_delta_do(
    data: pd.DataFrame,
    depth_col: str = 'Depth',
    do_col: str = 'DO',
    salinity_col: str = 'Salinity',
    temperature_col: str = 'Temperature',
    depth_interval: float | None = None,
    do_threshold: float | None = None,
    salinity_threshold: float | None = None,
    temperature_threshold: float | None = None,
    anomaly_min_depth: float | None = None,
    depth_merge_tolerance: float | None = None,
    duplicate_depth_strategy: str | None = None,
    remove_outliers: bool = True,
    verbose: bool = False
) -> pd.DataFrame:
    """
    参考论文方法计算亚表层异常信号（以 DO 为主）。

    步骤概述：
    1. 按 Profile_number 分组，分别处理每个剖面；
    2. 计算 DO 随深度的一阶导数，定位 DO 的正峰值（代表可能的表层富氧水体俯冲信号）；
    3. 以 DO 峰深度为中心，取窗口 [p-Δp, p+Δp]，用两端点连线构造参考剖面；
    4. 在峰值同一深度计算 ΔDO、ΔSalinity 与 ΔTemperature（原始值减参考线值）；
    5. 以 ΔDO ≥ do_threshold 作为必要条件；如设置了 salinity_threshold 或 temperature_threshold > 0，可附加 |ΔSalinity/ΔTemperature| 过滤；
    6. 若 anomaly_min_depth > 0，仅保留深度不小于该阈值的异常；
    7. 同一剖面内若有相距很近的多个候选深度（常见于峰值上下各一点），按 depth_merge_tolerance（dbar）合并，仅保留 delta_do 较大的记录。

    参数：
        data (pd.DataFrame): 包含多个剖面数据的表；需包含 Profile_number、深度、DO、盐度等列。
        depth_col (str): 深度列名，默认 'Depth'。
        do_col (str): 溶解氧列名，默认 'DO'。
        salinity_col (str): 盐度列名，默认 'Salinity'。
        temperature_col (str): 温度列名，默认 'Temperature'。
        depth_interval (float | None): 深度窗口半宽；None → 全局 `_default_depth_interval`。
        do_threshold (float | None): ΔDO 阈值；None → `_default_delta_do_threshold`。
        salinity_threshold (float | None): 盐度阈值；None → `_default_salinity_threshold`；≤0 不启用过滤。
        temperature_threshold (float | None): 温度阈值；None → `_default_temperature_threshold`；≤0 不启用过滤。
        anomaly_min_depth (float | None): ΔDO 异常最小深度；≤0 表示不做深度过滤；None 表示使用 processing.yml 的 anomaly_min_depth。
        depth_merge_tolerance (float | None): 深度近邻合并阈值；None → `_default_depth_merge_tolerance`；≤0 不合并。
        duplicate_depth_strategy (str | None): 同深度多记录聚合策略；None → `_default_duplicate_depth_strategy`。
        remove_outliers (bool): 基础 QC 与规则过滤，默认 True。
        verbose (bool): 是否打印进度信息（开始处理/已处理/未检测到/总共检测到），默认 False。

    返回：
        pd.DataFrame: 每个满足条件的峰值一行，含
        Profile_number, depth, delta_do, delta_salinity, delta_temperature,
        do_value, salinity_value, temperature_value，
        以及 Year/Month/Day/Longitude/Latitude/Platform_number（若存在）。
        若无满足条件记录，返回空表。
    
    提示:
        可调用 `print_current_processing_defaults()` 查看当前全局默认阈值与处理参数，
        以便核对本函数 None 回退所采用的配置来源。
    """
    
    # 检查必要的列是否存在
    required_cols = [depth_col, do_col, salinity_col, 'Profile_number']
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        print(f"警告：缺少必要的列: {missing_cols}")
        return pd.DataFrame()
    
    # 复制数据以避免修改原始数据
    input_data = data.copy()
    
    # 存储所有剖面的结果
    all_results = []
    
    # 按Profile_number分组处理
    profile_groups = input_data.groupby('Profile_number')
    total_profiles = len(profile_groups)
    if verbose:
        print(f"开始处理 {total_profiles} 个剖面...")
    
    processed_profiles = 0
    
    # 参数回退（放在循环外，避免每个剖面重复判定）
    if depth_interval is None:
        depth_interval = _default_depth_interval
    if do_threshold is None:
        do_threshold = _default_delta_do_threshold
    if salinity_threshold is None:
        salinity_threshold = _default_salinity_threshold
    if temperature_threshold is None:
        temperature_threshold = _default_temperature_threshold
    if depth_merge_tolerance is None:
        depth_merge_tolerance = _default_depth_merge_tolerance
    if duplicate_depth_strategy is None:
        duplicate_depth_strategy = _default_duplicate_depth_strategy
    if anomaly_min_depth is None:
        anomaly_min_depth = _cfg_anomaly_min_depth

    for profile_num, profile_data in profile_groups:
        # 质量控制：移除异常值和质量标记不良的数据
        if remove_outliers:
            # 应用Argo QC标准：仅保留等级为{1,2,5,8}的观测
            for var in [do_col, salinity_col, temperature_col]:
                qc_column_name = f"{var}_Flag"
                if qc_column_name in profile_data.columns:
                    good_qc_flags = ['1', '2', '5', '8', 1, 2, 5, 8]
                    bad_qc_mask = ~profile_data[qc_column_name].isin(good_qc_flags)
                    profile_data.loc[bad_qc_mask, var] = np.nan
            
            # 规则法：移除已知的错误值
            if do_col in profile_data.columns:
                bad_do_mask = profile_data[do_col] <= 1.0
                profile_data.loc[bad_do_mask, do_col] = np.nan
        
        # 移除包含NaN值的行
        drop_subset = [depth_col, do_col, salinity_col, temperature_col]
        profile_data_clean = profile_data.dropna(subset=drop_subset)
        
        if len(profile_data_clean) < 5:  # 需要至少5个数据点来计算导数和峰值
            continue  # 跳过数据点太少的剖面

        # 按深度排序
        profile_data_clean = profile_data_clean.sort_values(by=depth_col).reset_index(drop=True)

        # 处理重复或非严格递增的深度值以避免 np.gradient 内部出现除零 (dx1 或 dx2 = 0)
        # 策略：对同一深度的多条记录按 duplicate_depth_strategy 聚合为单条
        def _best_qc_pick(group: pd.DataFrame) -> pd.Series:
            qccol = f"{do_col}_Flag"
            priority = {1: 0, 2: 1, 5: 2, 8: 3}
            def rank(v):
                try:
                    iv = int(v)
                except Exception:
                    return 999
                return priority.get(iv, 999)
            if qccol in group.columns:
                ranks = group[qccol].apply(rank)
                min_rank = ranks.min()
                picked = group.loc[ranks[ranks == min_rank].index]
                # 并列时优先 DO 非空，随后保留第一条
                picked = picked[pd.notna(picked[do_col])] if do_col in picked.columns else picked
                return picked.iloc[0]
            # 无 QC 列则退化为 first
            return group.iloc[0]

        def _mean_pick(group: pd.DataFrame) -> pd.Series:
            # 对核心变量取均值，其他元数据取首个
            first = group.iloc[0].copy()
            for c in [do_col, salinity_col, temperature_col]:
                if c in group.columns:
                    first[c] = group[c].astype(float).mean()
            return first

        if profile_data_clean[depth_col].duplicated().any():
            strategy = (duplicate_depth_strategy or 'best_qc').lower()
            if strategy not in {'best_qc','first','mean','max','min'}:
                strategy = 'best_qc'
            grouped = list(profile_data_clean.groupby(depth_col, sort=False))
            picked_rows = []
            for depth_val, grp in grouped:
                if len(grp) == 1:
                    picked_rows.append(grp.iloc[0])
                    continue
                if strategy == 'best_qc':
                    picked_rows.append(_best_qc_pick(grp))
                elif strategy == 'first':
                    picked_rows.append(grp.iloc[0])
                elif strategy == 'mean':
                    picked_rows.append(_mean_pick(grp))
                elif strategy == 'max':
                    picked_rows.append(grp.loc[grp[do_col].idxmax()])
                elif strategy == 'min':
                    picked_rows.append(grp.loc[grp[do_col].idxmin()])
            profile_data_clean = pd.DataFrame(picked_rows)
            # 可能破坏原索引，重排并按深度排序
            profile_data_clean = profile_data_clean.sort_values(by=depth_col).reset_index(drop=True)
        depth_series = profile_data_clean[depth_col].values
        # 若仍非严格递增（可能存在逆序或噪声导致深度减小），尝试再次排序并去除 <= 前一值 的点
        # （二次排序保证顺序，随后用 np.diff > 0 过滤）
        if np.any(np.diff(depth_series) <= 0):
            # 先强制排序（已经排序过, 这里是保险）
            profile_data_clean = profile_data_clean.sort_values(by=depth_col).reset_index(drop=True)
            # 过滤掉与前一个深度差 <= 0 的观测
            cleaned_rows = [0]
            last_depth = profile_data_clean.loc[0, depth_col]
            for ridx in range(1, len(profile_data_clean)):
                dval = profile_data_clean.loc[ridx, depth_col]
                if dval > last_depth:  # 严格递增才保留
                    cleaned_rows.append(ridx)
                    last_depth = dval
            profile_data_clean = profile_data_clean.iloc[cleaned_rows].reset_index(drop=True)
            depth_series = profile_data_clean[depth_col].values
        # 若清理后点数不足或仍不严格递增则跳过
        if len(profile_data_clean) < 5 or (len(depth_series) > 1 and np.any(np.diff(depth_series) <= 0)):
            continue

        # 计算深度对DO和Salinity的斜率（一阶导数）
        depth_values = profile_data_clean[depth_col].values  # 已保证严格递增
        do_values = profile_data_clean[do_col].values
        salinity_values = profile_data_clean[salinity_col].values
        temperature_values = profile_data_clean[temperature_col].values

        # 使用中心差分法计算导数
        do_slopes = np.gradient(do_values, depth_values)

        # 重要说明：坐标系差异
        # 图像坐标系（海洋学习惯）：深度向下为正，图像上正斜率表示随深度增加变量增大
        # 数据计算坐标系：np.gradient计算dVar/dDepth，图像正斜率对应数值负斜率(<0)
        # 峰值识别：图像上斜率从正变负为负峰值，从负变正为正峰值
        # 因此代码中的判断条件与图像描述看似相反，但逻辑正确

        # 仅定位 DO 的峰值（保留正峰值作为俯冲信号候选）
        do_peaks = []
        for i in range(1, len(do_slopes) - 1):
            # DO 正峰值：图像上斜率从负变正（数值上从正变负）
            if do_slopes[i-1] > 0 and do_slopes[i+1] < 0:
                do_peaks.append((i, 'positive', depth_values[i]))
            # DO 负峰值：如需扩展也可纳入；当前仅保留正峰值

        if not do_peaks:
            continue

        # 对每个 DO 正峰，按同一深度计算 ΔDO、ΔSalinity、ΔTemperature
        use_salinity_filter = (salinity_threshold is not None and salinity_threshold > 0)
        use_temperature_filter = (temperature_threshold is not None and temperature_threshold > 0)
        profile_results = []

        for do_idx, do_type, target_depth in do_peaks:
            # 定义窗口 [p-Δp, p+Δp]
            depth_lower = max(0, target_depth - depth_interval)
            depth_upper = target_depth + depth_interval

            depth_mask = (depth_values >= depth_lower) & (depth_values <= depth_upper)
            if not np.any(depth_mask):
                continue

            depth_range = depth_values[depth_mask]
            do_range = do_values[depth_mask]
            salinity_range = salinity_values[depth_mask]
            temperature_range = temperature_values[depth_mask]

            if len(depth_range) < 2:
                continue

            # 线性参考剖面（两端点连线）
            do_ref_values = np.interp(
                depth_range, [depth_range[0], depth_range[-1]], [do_range[0], do_range[-1]]
            )
            salinity_ref_values = np.interp(
                depth_range, [depth_range[0], depth_range[-1]], [salinity_range[0], salinity_range[-1]]
            )
            temperature_ref_values = np.interp(
                depth_range, [depth_range[0], depth_range[-1]], [temperature_range[0], temperature_range[-1]]
            )

            # 取最接近目标深度的观测点
            target_idx = np.argmin(np.abs(depth_range - target_depth))

            delta_do = do_range[target_idx] - do_ref_values[target_idx]
            delta_salinity = salinity_range[target_idx] - salinity_ref_values[target_idx]
            delta_temperature = temperature_range[target_idx] - temperature_ref_values[target_idx]

            # 判定：ΔDO 为必需条件；ΔSalinity/ΔTemperature 在各自阈值>0时作为附加过滤（与条件）
            cond = (delta_do >= do_threshold)
            if use_salinity_filter:
                cond = cond and (not np.isnan(delta_salinity)) and (abs(delta_salinity) >= salinity_threshold)
            if use_temperature_filter:
                cond = cond and (not np.isnan(delta_temperature)) and (abs(delta_temperature) >= temperature_threshold)
            if cond:
                result = {
                    'Profile_number': profile_num,
                    'depth': target_depth,
                    'delta_do': delta_do,
                    'delta_salinity': delta_salinity,
                    'delta_temperature': delta_temperature,
                    'do_value': do_values[do_idx],
                    'salinity_value': salinity_range[target_idx],
                    'temperature_value': temperature_range[target_idx]
                }

                # 添加额外的剖面信息（如果存在）
                if 'Year' in profile_data_clean.columns:
                    result['Year'] = profile_data_clean['Year'].iloc[0]
                if 'Month' in profile_data_clean.columns:
                    result['Month'] = profile_data_clean['Month'].iloc[0]
                if 'Day' in profile_data_clean.columns:
                    result['Day'] = profile_data_clean['Day'].iloc[0]
                if 'Longitude' in profile_data_clean.columns:
                    result['Longitude'] = profile_data_clean['Longitude'].iloc[0]
                if 'Latitude' in profile_data_clean.columns:
                    result['Latitude'] = profile_data_clean['Latitude'].iloc[0]
                if 'Platform_number' in profile_data_clean.columns:
                    result['Platform_number'] = profile_data_clean['Platform_number'].iloc[0]

                profile_results.append(result)

        # 剖面内“深度近邻合并”：按 delta_do 降序贪心选取，避免在峰上下方重复取点
        if profile_results:
            if anomaly_min_depth is not None and anomaly_min_depth > 0:
                filtered_results = []
                for rec in profile_results:
                    depth_val = rec.get('depth')
                    if depth_val is None or np.isnan(depth_val):
                        continue
                    if depth_val >= anomaly_min_depth:
                        filtered_results.append(rec)
                profile_results = filtered_results
            if not profile_results:
                continue
            if depth_merge_tolerance is not None and depth_merge_tolerance > 0:
                profile_results.sort(key=lambda r: (np.nan_to_num(r['delta_do'], nan=-np.inf)), reverse=True)
                kept = []
                kept_depths = []
                for rec in profile_results:
                    d = rec['depth']
                    if all(abs(d - kd) >= depth_merge_tolerance for kd in kept_depths):
                        kept.append(rec)
                        kept_depths.append(d)
                # 输出前按深度升序，便于查看
                kept.sort(key=lambda r: r['depth'])
                all_results.extend(kept)
            else:
                all_results.extend(profile_results)
        
        processed_profiles += 1
        if processed_profiles % 100 == 0 and verbose:
            print(f"已处理 {processed_profiles}/{total_profiles} 个剖面...")
    
    if not all_results:
        if verbose:
            print("未检测到满足阈值条件的DO异常信号。")
        return pd.DataFrame()
    
    results_df = pd.DataFrame(all_results)
    if anomaly_min_depth is not None and anomaly_min_depth > 0 and 'depth' in results_df.columns:
        results_df = results_df[results_df['depth'] >= anomaly_min_depth]
    if verbose:
        print(f"总共检测到 {len(results_df)} 个潜在的DO异常信号，来自 {len(results_df['Profile_number'].unique())} 个剖面")

    return results_df

def _export_interacting_argo_worker(args):
    """
    Worker function for export_all_interacting_argo to support multiprocessing.
    Optimized for batch processing by Year with Low Memory Footprint.
    """
    y, kinds, cef, r_key, lmin, lmax, ltmin, ltmax = args

    try:
        # 1. 加载 Argo (按年加载)
        df = load_argo_data(y)
        if df.empty: return [], pd.DataFrame()
        
        # 地理过滤 (粗筛)
        pad = 2.0 
        lon_vals = df['Longitude'].to_numpy(dtype=float)
        lat_vals = df['Latitude'].to_numpy(dtype=float)
        
        if lmax < lmin: # 跨日界线
            lon_mask = (lon_vals >= lmin - pad) | (lon_vals <= lmax + pad)
        else:
            lon_mask = (lon_vals >= lmin - pad) & (lon_vals <= lmax + pad)
            
        lat_mask = (lat_vals >= ltmin - pad) & (lat_vals <= ltmax + pad)
        df_geo = df[lon_mask & lat_mask].copy()
        if df_geo.empty: return [], pd.DataFrame()
        
        # 提取 Baseline (去重)
        baseline = (
            df_geo.sort_values(['Profile_number', 'Depth'])
            .groupby('Profile_number', as_index=False)
            .first()[['Profile_number', 'Longitude', 'Latitude', 'Year', 'Month', 'Day']]
        )
        baseline['date'] = pd.to_datetime(baseline[['Year', 'Month', 'Day']])
        
        # 按天分组 Argo (Dict of DataFrames)
        argo_by_day = {d: grp for d, grp in baseline.groupby('date')}

        # 2. 加载当年涡旋数据 (仅元数据，不含轮廓，节省内存)
        start_d = pd.Timestamp(f"{y}-01-01")
        end_d = pd.Timestamp(f"{y}-12-31")
        
        # Pass 1: Load lightweight metadata
        eddy_data_meta = _load_eddy_datasets_for_range(
            kinds=kinds,
            start_date=start_d,
            end_date=end_d,
            region_key=r_key,
            include_contours=False # 关键：先不加载轮廓
        )
        
        # 3. 构建涡旋按天索引 (仅元数据)
        eddy_by_day = defaultdict(list)
        for ds_name, ds_tracks in eddy_data_meta.items():
            for item in ds_tracks:
                track = item[0] if isinstance(item, tuple) else item
                if not track: continue
                
                # 批量转换时间
                raw_times = [p[1] for p in track]
                converted_times = convert_date(raw_times)
                
                if isinstance(converted_times, pd.Series):
                    ts_index = pd.DatetimeIndex(converted_times)
                else:
                    ts_index = pd.DatetimeIndex([converted_times])
                
                mask = (ts_index >= start_d) & (ts_index <= end_d)
                valid_indices = np.where(mask)[0]
                
                for idx in valid_indices:
                    p = track[idx]
                    d = ts_index[idx].normalize()
                    
                    eddy_by_day[d].append({
                        'track_id': p[0],
                        'lon': p[2],
                        'lat': p[3],
                        'radius': p[8],
                        'ds_name': ds_name,
                        # 'contour_lon': p[6], # Pass 1 中这些是 None 或空
                        # 'contour_lat': p[7]
                    })

        interacting_records = []
        
        # 记录需要进一步检查轮廓的候选者
        # candidates[ds_name][track_id] = set(dates)
        candidates = defaultdict(lambda: defaultdict(set))
        
        # 4. Pass 1: 批量圆筛选
        common_days = sorted(set(argo_by_day.keys()) & set(eddy_by_day.keys()))
        
        # print(f"[Info] Year {y}: Checking {len(common_days)} days with Argo & Eddies...")

        for current_date in common_days:
            day_argo_df = argo_by_day[current_date]
            day_eddies = eddy_by_day[current_date]
            
            argo_lons = day_argo_df['Longitude'].to_numpy(dtype=float)
            argo_lats = day_argo_df['Latitude'].to_numpy(dtype=float)
            
            for eddy in day_eddies:
                e_lon, e_lat, e_rad = eddy['lon'], eddy['lat'], eddy['radius']
                eff_rad = e_rad * cef
                
                # 粗筛
                scale_rough = approximate_degree_length(e_lat)
                # 使用 1.5 倍 buffer，并分别计算经纬度方向的度数阈值
                # 注意：meters_per_degree_lon 随纬度增加而减小，因此同样的米数对应的度数会变大
                # 为安全起见，使用 local scale
                rad_deg_lat = (eff_rad / scale_rough['meters_per_degree_lat']) * 1.5
                # 防止极点附近除零或过大，设置上限（例如 180度）
                if scale_rough['meters_per_degree_lon'] < 1000: # 极靠近极点
                     rad_deg_lon = 180.0
                else:
                     rad_deg_lon = (eff_rad / scale_rough['meters_per_degree_lon']) * 1.5
                
                dlon = np.abs(_minimal_lon_diff_deg(argo_lons, e_lon))
                dlat = np.abs(argo_lats - e_lat)
                
                mask_bb = (dlon < rad_deg_lon) & (dlat < rad_deg_lat)
                if not np.any(mask_bb):
                    continue
                
                # 精细距离
                # 使用 adaptive_distance_m 替代手动平面计算，以解决高纬度畸变问题
                c_lons = argo_lons[mask_bb]
                c_lats = argo_lats[mask_bb]
                
                dists = adaptive_distance_m(c_lons, c_lats, e_lon, e_lat)
                mask_circle = dists <= eff_rad
                
                if np.any(mask_circle):
                    # 命中圆，记录为候选，稍后加载轮廓做精确检查
                    candidates[eddy['ds_name']][eddy['track_id']].add(current_date)

        # 5. Pass 2: 批量加载轮廓并精确匹配 (Batch Loading Optimization)
        # 收集所有需要加载轮廓的 track_ids
        for ds_name, track_map in candidates.items():
            target_ids = list(track_map.keys())
            if not target_ids:
                continue
            
            try:
                # 批量加载：一次性读取该数据集下所有候选涡旋的完整轨迹（含轮廓）
                # find_track 返回 {tid: list_of_rows} (当 len(target_ids) > 1)
                # 或 list_of_rows (当 len(target_ids) == 1)
                # 统一处理为 dict
                batch_res = find_track(
                    ds_name.lower(), 
                    target_ids, 
                    region=r_key, 
                    include_contours=True, 
                    return_list=True
                )
                
                tracks_dict = {}
                if len(target_ids) == 1:
                    # 单个 ID 返回的是 list
                    tracks_dict[target_ids[0]] = batch_res
                else:
                    tracks_dict = batch_res
                
                # 遍历每个涡旋进行多边形检测
                for track_id, full_track_list in tracks_dict.items():
                    dates_set = track_map.get(track_id, set())
                    if not dates_set: continue
                    
                    # 构建时间索引: YYYYMMDD int -> row
                    # full_track_list item: [tid, ymd, clon, clat, mlon, mlat, clon_poly, clat_poly, rad, ...]
                    # Index: 1=ymd, 2=clon, 3=clat, 6=clon_poly, 7=clat_poly, 8=rad
                    
                    # 快速筛选：只处理 dates_set 中的日期
                    # 将 dates_set 转为 YYYYMMDD 整数集合以便快速查找
                    target_ymds = {d.year * 10000 + d.month * 100 + d.day for d in dates_set}
                    
                    for row in full_track_list:
                        ymd = row[1]
                        if ymd not in target_ymds:
                            continue
                            
                        # 提取数据
                        e_lon = row[2]
                        e_lat = row[3]
                        c_lon_poly = row[6]
                        c_lat_poly = row[7]
                        e_rad = row[8]
                        
                        # 还原日期对象用于查找 Argo
                        # ymd is int YYYYMMDD
                        y_ = ymd // 10000
                        m_ = (ymd % 10000) // 100
                        d_ = ymd % 100
                        current_date = pd.Timestamp(year=y_, month=m_, day=d_)
                        
                        # 获取当天的 Argo
                        if current_date not in argo_by_day: continue
                        day_argo_df = argo_by_day[current_date]
                        
                        argo_lons = day_argo_df['Longitude'].to_numpy(dtype=float)
                        argo_lats = day_argo_df['Latitude'].to_numpy(dtype=float)
                        argo_ids = day_argo_df['Profile_number'].to_numpy()
                        
                        # 重复圆筛选 (为了拿到 mask_circle 对应的点)
                        eff_rad = e_rad * cef
                        
                        scale_rough = approximate_degree_length(e_lat)
                        rad_deg_lat = (eff_rad / scale_rough['meters_per_degree_lat']) * 1.5
                        if scale_rough['meters_per_degree_lon'] < 1000:
                             rad_deg_lon = 180.0
                        else:
                             rad_deg_lon = (eff_rad / scale_rough['meters_per_degree_lon']) * 1.5

                        dlon = np.abs(_minimal_lon_diff_deg(argo_lons, e_lon))
                        dlat = np.abs(argo_lats - e_lat)
                        mask_bb = (dlon < rad_deg_lon) & (dlat < rad_deg_lat)
                        
                        if not np.any(mask_bb): continue
                        
                        c_lons = argo_lons[mask_bb]
                        c_lats = argo_lats[mask_bb]
                        c_ids = argo_ids[mask_bb]
                        
                        # 使用 adaptive_distance_m 替代手动平面计算
                        dists = adaptive_distance_m(c_lons, c_lats, e_lon, e_lat)
                        mask_circle = dists <= eff_rad
                        
                        if np.any(mask_circle):
                            circle_hits_idx = np.where(mask_circle)[0]
                            
                            # 多边形检测
                            has_poly = False
                            is_inside = np.zeros(len(circle_hits_idx), dtype=bool)
                            
                            if c_lon_poly is not None and c_lat_poly is not None:
                                c_lon_arr = np.asarray(c_lon_poly)
                                c_lat_arr = np.asarray(c_lat_poly)
                                if c_lon_arr.size >= 3:
                                    has_poly = True
                                    c_lon_norm = e_lon + _minimal_lon_diff_deg(c_lon_arr, e_lon)
                                    verts = np.column_stack((c_lon_norm, c_lat_arr))
                                    path = MplPath(verts)
                                    
                                    p_lons = c_lons[mask_circle]
                                    p_lats = c_lats[mask_circle]
                                    p_lons_norm = e_lon + _minimal_lon_diff_deg(p_lons, e_lon)
                                    points = np.column_stack((p_lons_norm, p_lats))
                                    is_inside = path.contains_points(points)
                            
                            for idx_in_subset, inside in enumerate(is_inside):
                                real_idx = circle_hits_idx[idx_in_subset]
                                if has_poly:
                                    m_str = 'poly' if inside else 'circle'
                                else:
                                    m_str = 'circle'
                                    
                                interacting_records.append({
                                    'Profile_number': c_ids[real_idx],
                                    'Year': y,
                                    'Month': current_date.month,
                                    'Day': current_date.day,
                                    'track_id': track_id,
                                    'ds_name': ds_name,
                                    'method': m_str
                                })
            except Exception as e:
                print(f"[Warn] Batch processing failed for {ds_name} in Year {y}: {e}")
                import traceback
                traceback.print_exc()
                continue

        return interacting_records, baseline
    except Exception as e:
        print(f"[Error] Year {y}: {e}")
        import traceback
        traceback.print_exc()
        return [], pd.DataFrame()

def export_all_interacting_argo(
    start_year: int,
    end_year: int,
    eddy_datasets: dict | list[str] | tuple[str, ...] | None = None,
    circle_enlargement_factor: float | None = None,
    output_path: str | Path | None = None,
    num_workers: int = 1,
):
    """
    计算并导出指定年份范围内所有位于涡旋内部的 Argo 剖面数据。

    功能：
        1. 遍历指定年份的 Argo 数据。
        2. 加载对应的涡旋轨迹数据（META Tracks）。
        3. 判断每个 Argo 剖面是否位于涡旋内部（支持多核并行加速）。
        4. 将所有位于涡旋内的 Argo 剖面信息（包含匹配的涡旋ID、类型等）保存为 Parquet 文件。
        5. 同时保存该区域内所有 Argo 剖面的基础信息（用于后续计算交互率分母）。

    参数:
        start_year (int): 起始年份。
        end_year (int): 结束年份。
        eddy_datasets (list, optional): 指定使用的涡旋数据集列表（如 ['acl', 'acs', 'cyclonic', 'anticyclonic']）。默认为 None，使用所有可用数据集。
        circle_enlargement_factor (float, optional): 涡旋边界放大系数。默认为从 processing.yml 中读取。
        output_path (str | Path, optional): 结果文件保存路径。默认为 `plot_outputs/<region>/statistics/all_interacting_argo_<years>.parquet`。
        num_workers (int): 并行进程数。默认为 1。

    输出:
        生成一个 Parquet 文件，包含所有与涡旋发生交互的 Argo 剖面详细信息。
    """
    if circle_enlargement_factor is None:
        circle_enlargement_factor = globals().get('circle_enlargement_factor', 1.2)
    
    region_slug = _current_region_key()
    if output_path is None:
        out_dir = Path(plots_output_root) / region_slug / "statistics"
        out_dir.mkdir(parents=True, exist_ok=True)
        output_path = out_dir / f"all_interacting_argo_{start_year}_{end_year}.parquet"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"--- Exporting All Interacting Argo (Baseline) {start_year}-{end_year} ---")
    print(f"Output: {output_path}")
    print(f"Workers: {num_workers} (Task granularity: Yearly)")
    
    # 准备任务列表 (Year)
    years = list(range(start_year, end_year + 1))
    
    # 准备参数
    kinds = eddy_datasets if eddy_datasets else ['acs', 'acl', 'cs', 'cl']
    if isinstance(kinds, dict): kinds = list(kinds.keys())
    kinds = [str(k).lower() for k in kinds]
    
    worker_args = [
        (y, kinds, circle_enlargement_factor, region_slug, lonmin, lonmax, latmin, latmax)
        for y in years
    ]
    
    print(f"Total tasks: {len(worker_args)}")
    
    # 执行并行
    all_results = []
    all_profiles = []
    
    if num_workers > 1:
        with multiprocessing.Pool(processes=num_workers, initializer=init_worker, initargs=(None,), maxtasksperchild=1) as pool:
            for res_interactions, res_profiles in tqdm(pool.imap_unordered(_export_interacting_argo_worker, worker_args), total=len(worker_args), desc="Processing Years"):
                all_results.extend(res_interactions)
                if not res_profiles.empty:
                    all_profiles.append(res_profiles)
    else:
        for args in tqdm(worker_args, desc="Processing Years"):
            res_interactions, res_profiles = _export_interacting_argo_worker(args)
            all_results.extend(res_interactions)
            if not res_profiles.empty:
                all_profiles.append(res_profiles)
            
    # 保存交互记录
    if all_results:
        df_out = pd.DataFrame(all_results)
        df_out.to_parquet(output_path, index=False)
        print(f"Saved {len(df_out)} interacting records to: {output_path}")
    else:
        print("No interacting profiles found.")

    # 保存区域内所有 Argo 剖面 (Baseline)
    if all_profiles:
        df_profiles = pd.concat(all_profiles, ignore_index=True)
        # 构造 Baseline 文件名: all_region_argo_{start}_{end}.parquet
        baseline_path = output_path.parent / f"all_region_argo_{start_year}_{end_year}.parquet"
        df_profiles.to_parquet(baseline_path, index=False)
        print(f"Saved {len(df_profiles)} region profiles (Baseline) to: {baseline_path}")
    else:
        print("No region profiles found.")

def calculate_interaction_statistics(
    start_year: int,
    end_year: int,
    eddy_datasets: dict | list[str] | tuple[str, ...] | None = None,
    do_threshold: float | None = None,
    salinity_threshold: float | None = None,
    temperature_threshold: float | None = None,
    depth_interval: float | None = None,
    depth_merge_tolerance: float | None = None,
    duplicate_depth_strategy: str | None = None,
    anomaly_min_depth: float | None = None,
    circle_enlargement_factor: float | None = None,
    save_report: bool = True,
    precomputed_file: str | Path | None = None,
    anomalies_file: str | Path | None = None,
    use_precomputed_anomalies: bool = True,
):
    """
    计算并对比 Argo 剖面与涡旋的交互概率（Baseline vs Anomalies）。
    
    功能：
        1. 加载指定年份的所有 Argo 数据。
        2. 加载预计算的交互记录文件（由 export_all_interacting_argo 生成）。
           注意：必须先运行 export_all_interacting_argo 生成该文件，否则报错。
        3. 计算 Baseline：所有 Argo 剖面中，有多少比例落在涡旋内。
        4. 计算 Anomalies：筛选出 ΔDO 异常剖面，计算其中有多少比例落在涡旋内。
           - 若 use_precomputed_anomalies=True (默认)，会自动尝试查找 plot_argo_hotspots 生成的异常文件。
           - 若找到文件，直接读取；若未找到或读取失败，回退到实时调用 calculate_delta_do 计算。
        5. 输出对比报告。
    """
    # 参数回退
    if do_threshold is None: do_threshold = _default_delta_do_threshold
    if salinity_threshold is None: salinity_threshold = _default_salinity_threshold
    if temperature_threshold is None: temperature_threshold = _default_temperature_threshold
    if depth_interval is None: depth_interval = _default_depth_interval
    if depth_merge_tolerance is None: depth_merge_tolerance = _default_depth_merge_tolerance
    if duplicate_depth_strategy is None: duplicate_depth_strategy = _default_duplicate_depth_strategy
    if anomaly_min_depth is None: anomaly_min_depth = _cfg_anomaly_min_depth
    if circle_enlargement_factor is None: circle_enlargement_factor = globals().get('circle_enlargement_factor', 1.2)

    print(f"--- Calculating Interaction Statistics {start_year}-{end_year} ---")
    
    # --- 0. 尝试加载预计算的交互记录 ---
    interacting_ids = set()
    loaded_precomputed = False
    
    region_slug = _current_region_key()
    if precomputed_file is None:
        # 尝试默认路径
        default_file = Path(plots_output_root) / region_slug / "statistics" / f"all_interacting_argo_{start_year}_{end_year}.parquet"
        if default_file.exists():
            precomputed_file = default_file
            
    if precomputed_file and Path(precomputed_file).exists():
        print(f"[*] Loading precomputed interactions from: {precomputed_file}")
        try:
            df_int = pd.read_parquet(precomputed_file)
            if 'Profile_number' in df_int.columns:
                interacting_ids = set(df_int['Profile_number'].unique())
                loaded_precomputed = True
                print(f"[*] Loaded {len(interacting_ids)} unique interacting profiles.")
        except Exception as e:
            print(f"[WARN] Failed to read precomputed file: {e}")
            
    # 1. 加载 Argo 数据 (用于计算分母和检测异常)
    # 优化：尝试加载预计算的区域 Argo 列表 (Baseline)
    baseline_profiles = pd.DataFrame()
    loaded_baseline_file = False
    
    default_baseline_file = Path(plots_output_root) / region_slug / "statistics" / f"all_region_argo_{start_year}_{end_year}.parquet"
    if default_baseline_file.exists():
        print(f"[*] Loading precomputed region profiles from: {default_baseline_file}")
        try:
            baseline_profiles = pd.read_parquet(default_baseline_file)
            if not baseline_profiles.empty:
                # 确保严格符合当前区域定义 (因为 export 可能包含 padding)
                lons = baseline_profiles['Longitude'].to_numpy()
                lats = baseline_profiles['Latitude'].to_numpy()
                mask_geo = _region_lon_mask(lons, lonmin, lonmax) & (lats >= latmin) & (lats <= latmax)
                baseline_profiles = baseline_profiles[mask_geo].copy()
                
                baseline_profiles['date'] = pd.to_datetime(baseline_profiles[['Year', 'Month', 'Day']])
                loaded_baseline_file = True
                print(f"[*] Loaded {len(baseline_profiles)} region profiles (Strictly inside region).")
        except Exception as e:
            print(f"[WARN] Failed to read region profiles file: {e}")

    argo_geo = pd.DataFrame() # 用于实时计算异常的原始数据
    
    if not loaded_baseline_file:
        print("[*] Loading Argo data (Real-time)...")
        argo_list = []
        for y in range(start_year, end_year + 1):
            try:
                df = load_argo_data(y)
                if not df.empty:
                    argo_list.append(df)
            except Exception:
                pass
        if not argo_list:
            print("No Argo data found.")
            return
        argo_all = pd.concat(argo_list, ignore_index=True)
        
        # 地理过滤
        lon_vals = argo_all['Longitude'].to_numpy(dtype=float)
        lat_vals = argo_all['Latitude'].to_numpy(dtype=float)
        lon_mask = _region_lon_mask(lon_vals, lonmin, lonmax)
        lat_mask = (lat_vals >= latmin) & (lat_vals <= latmax)
        argo_geo = argo_all[lon_mask & lat_mask].copy()
        
        if argo_geo.empty:
            print("No Argo data in region.")
            return

        # 提取 Baseline (所有剖面)
        # 按 Profile_number 去重，保留时间/位置
        baseline_profiles = (
            argo_geo.sort_values(['Profile_number', 'Depth'])
            .groupby('Profile_number', as_index=False)
            .first()[['Profile_number', 'Longitude', 'Latitude', 'Year', 'Month', 'Day']]
        )
        baseline_profiles['date'] = pd.to_datetime(baseline_profiles[['Year', 'Month', 'Day']])
    
    # 提取 Anomalies (异常剖面)
    anomalies = pd.DataFrame()
    
    # 自动查找默认异常文件
    if anomalies_file is None and use_precomputed_anomalies:
        region_slug = _current_region_key()
        out_dir = Path(plots_output_root) / region_slug / "plot_argo_hotspots"
        
        thr_str = f"{do_threshold:g}".replace('.', 'p')
        depth_suffix = ''
        if anomaly_min_depth is not None and anomaly_min_depth > 0:
            depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
            depth_suffix = f"_depth{depth_str}m"
            
        default_anomalies_file = out_dir / f"anomalies_{start_year}_{end_year}_thr{thr_str}{depth_suffix}.parquet"
        if default_anomalies_file.exists():
            anomalies_file = default_anomalies_file
            print(f"[*] Found default precomputed anomalies file: {anomalies_file}")
        else:
            print(f"[*] Default precomputed anomalies file not found: {default_anomalies_file}")
            print("    Will proceed with real-time calculation.")

    if anomalies_file and Path(anomalies_file).exists():
        print(f"[*] Loading precomputed anomalies from: {anomalies_file}")
        try:
            anomalies = pd.read_parquet(anomalies_file)
            # 确保只统计当前区域内的异常（如果文件包含更多区域）
            # 使用 Profile_number 与 baseline_profiles (已过滤区域) 取交集最稳妥
            valid_ids = set(baseline_profiles['Profile_number'])
            anomalies = anomalies[anomalies['Profile_number'].isin(valid_ids)].copy()
            print(f"[*] Loaded {len(anomalies)} anomalies (filtered by region).")
        except Exception as e:
            print(f"[WARN] Failed to read anomalies file: {e}. Falling back to calculation.")
            anomalies = pd.DataFrame() # 触发下方重新计算

    if anomalies.empty:
        print("[*] Detecting anomalies (Real-time calculation)...")
        # 如果没有加载 baseline file，说明 argo_geo 已经准备好了
        # 如果加载了 baseline file，但没有 anomalies file，我们需要重新加载 argo_geo 吗？
        # 是的，因为 calculate_delta_do 需要原始剖面数据 (argo_geo)
        
        if argo_geo.empty:
             print("[*] Reloading Argo data for anomaly calculation...")
             # 这里必须重新加载，因为之前可能跳过了
             argo_list = []
             for y in range(start_year, end_year + 1):
                try:
                    df = load_argo_data(y)
                    if not df.empty:
                        argo_list.append(df)
                except Exception:
                    pass
             if argo_list:
                 argo_all = pd.concat(argo_list, ignore_index=True)
                 lon_vals = argo_all['Longitude'].to_numpy(dtype=float)
                 lat_vals = argo_all['Latitude'].to_numpy(dtype=float)
                 lon_mask = _region_lon_mask(lon_vals, lonmin, lonmax)
                 lat_mask = (lat_vals >= latmin) & (lat_vals <= latmax)
                 argo_geo = argo_all[lon_mask & lat_mask].copy()

        if not argo_geo.empty:
            anomalies = calculate_delta_do(
                argo_geo,
                depth_interval=depth_interval,
                do_threshold=do_threshold,
                salinity_threshold=salinity_threshold,
                temperature_threshold=temperature_threshold,
                anomaly_min_depth=anomaly_min_depth,
                depth_merge_tolerance=depth_merge_tolerance,
                duplicate_depth_strategy=duplicate_depth_strategy,
                remove_outliers=True,
                verbose=False
            )
        else:
            print("[Error] Cannot calculate anomalies: No Argo data available.")
    
    anomaly_ids = set()
    if not anomalies.empty:
        anomaly_ids = set(anomalies['Profile_number'].unique())
    
    # 2. 若未加载预计算文件，则报错返回
    if not loaded_precomputed:
        print(f"[Error] Precomputed interaction file not found or failed to load.")
        print(f"Please run 'export_all_interacting_argo(start_year={start_year}, end_year={end_year}, ...)' first.")
        return

    # 4. 统计
    total_baseline = len(baseline_profiles)
    interacted_baseline = len(interacting_ids)
    pct_baseline = (interacted_baseline / total_baseline * 100) if total_baseline > 0 else 0.0
    
    total_anom = len(anomaly_ids)
    # 异常且交互 = 异常ID集合 与 交互ID集合 的交集
    interacted_anom = len(anomaly_ids & interacting_ids)
    pct_anom = (interacted_anom / total_anom * 100) if total_anom > 0 else 0.0
    
    # 5. 输出报告
    source_str = f"Precomputed File ({Path(precomputed_file).name})" if loaded_precomputed else "Real-time Calculation"
    
    report = (
        f"========================================\n"
        f"Interaction Statistics Report ({start_year}-{end_year})\n"
        f"Region: {_current_region_key()}\n"
        f"Source: {source_str}\n"
        f"Thresholds: ΔDO>={do_threshold}, Depth>={anomaly_min_depth}m\n"
        f"----------------------------------------\n"
        f"[Baseline] All Argo Profiles:\n"
        f"  Total Profiles:       {total_baseline}\n"
        f"  Inside Eddies:        {interacted_baseline}\n"
        f"  Interaction Rate:     {pct_baseline:.2f}%\n"
        f"----------------------------------------\n"
        f"[Subset] ΔDO Anomalies:\n"
        f"  Total Anomalies:      {total_anom}\n"
        f"  Inside Eddies:        {interacted_anom}\n"
        f"  Interaction Rate:     {pct_anom:.2f}%\n"
        f"----------------------------------------\n"
        f"Ratio (Anom Rate / Base Rate): {(pct_anom/pct_baseline if pct_baseline > 0 else 0):.2f}x\n"
        f"========================================"
    )
    
    print(report)
    
    if save_report:
        region_slug = _current_region_key()
        out_dir = Path(plots_output_root) / region_slug / "statistics"
        out_dir.mkdir(parents=True, exist_ok=True)
        
        thr_str = f"{do_threshold:g}".replace('.', 'p')
        depth_suffix = ''
        if anomaly_min_depth is not None and anomaly_min_depth > 0:
            depth_str = f"{anomaly_min_depth:g}".replace('.', 'p')
            depth_suffix = f"_depth{depth_str}m"
            
        fname = out_dir / f"interaction_stats_{start_year}_{end_year}_thr{thr_str}{depth_suffix}.txt"
        with open(fname, 'w') as f:
            f.write(report)
        print(f"Report saved to: {fname}")