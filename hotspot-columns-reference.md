# Hotspot Summary Table 变量参考

`export_hotspot_anomaly_summary_table` 默认输出路径：

```
plot_outputs/<method>/<region>/hotspot_anomaly_summary_table/
    hotspot_anomaly_summary_<start>_<end>_<run_tag>.xlsx
```

例如：`plot_outputs/do/global_ocean/hotspot_anomaly_summary_table/hotspot_anomaly_summary_2002_2023_do50_depth300m.xlsx`

compact_columns=True 后保留的列（约 28 列）。

## 基础信息

| 列 | 类型 | 说明 |
|---|---|---|
| `date` | 日期 | Argo 剖面日期 |
| `Year` | int | 年 |
| `Month` | int | 月 |
| `Day` | int | 日 |
| `Profile_number` | int | Argo 剖面编号 |
| `Platform_number` | int | 浮标平台编号 |
| `Longitude` | float | 剖面经度 |
| `Latitude` | float | 剖面纬度 |

## Argo 异常诊断

| 列 | 类型 | 说明 |
|---|---|---|
| `anomaly_depth_m` | float | 异常出现深度 (m) |
| `primary_value` | float | 主要检测指标值（取决于 `detection_config.method`） |
| `anomaly_score` | float | 异常综合评分。由 `calculate_delta_do` 多指标加权 |
| `delta_do` | float | DO 异常幅度 (μmol/kg)。异常点 DO 相对于背景参考的偏离 |
| `do_value` | float | 异常点的溶解氧值 (μmol/kg) |
| `delta_aou` | float | 表观耗氧量异常幅度 (μmol/kg) |
| `delta_temperature` | float | 温度异常幅度 (°C) |
| `delta_salinity` | float | 盐度异常幅度 (psu) |
| `surface_do_ref` | float | 表层 DO 参考值 (μmol/kg)。剖面 0-100m 内 DO 的中位数 |

## 近岸 DO 形态诊断

这组变量专门捕获近岸剖面 DO 的"先快速降低、后缓慢回升"的 V 形特征。

### 计算逻辑

对每个 Argo 剖面，在异常深度以上（0 ~ anomaly_depth）搜索 DO 最低点：

```
surface_ref   = median(DO[0–100m])
min_do        = min(DO[0 ~ anomaly_depth])
min_depth     = argmin(DO[0 ~ anomaly_depth]) 处的深度
anomaly_do    = DO[anomaly_depth]（若未知则取最近深度）

drop     = surface_ref − min_do           # 表层到最低点的降幅
recovery = anomaly_do − min_do            # 最低点到异常点的恢复量
gap      = anomaly_depth − min_depth      # DO 最低深度与异常深度的间距

nearshore_do_dip = (
    min_do ≤ 50       # DO 最低值很低
    AND drop ≥ 100    # 降幅大（快速下降）
    AND recovery ≥ 100 # 恢复量大（缓慢回升）
    AND gap ≥ 100 m   # 最低与异常有足够垂直距离
)
```

### 变量

| 列 | 类型 | 说明 |
|---|---|---|
| `pre_anomaly_do_min` | float | 异常深度以上 DO 最低值 (μmol/kg) |
| `pre_anomaly_do_min_depth_m` | float | 该最低值所在深度 (m) |
| `surface_to_min_do_drop` | float | 表层 DO 到最低点的降幅 (μmol/kg) |
| `min_to_anomaly_do_recovery` | float | 最低点到异常点的 DO 恢复量 (μmol/kg) |
| `min_to_anomaly_depth_gap_m` | float | DO 最低深度与异常深度的垂直间距 (m) |
| `do_v_shape_score` | float | V 形评分 = min(drop, recovery)。越大越接近近岸型 |

## Heave 诊断

以下变量由 `calculate_glorys_vertical_profile_diagnostics` 计算，
在 GLORYS 垂向剖面上追踪 Argo 异常深度对应的等密面。

### 计算流程

1. 确定 Argo 异常深度处的密度 σ_argo
2. 在 σ_argo − 0.5 到 σ_argo 之间搜索等密线
3. 等密线连通性约束：仅保留在 Argo 点 ±200m 垂向范围内存在过的 σ 面
4. 每条等密线在 ±200 km 窗口内找最浅深度 (z_min)，
   在 Argo 附近 ±50 km 找最深深度 (local_z_max)
5. Heave H = max(local_z_max − z_min)，Heave 峰值 σ = argmax 对应的 σ，
   z_min 取自 σ_peak 同一条等密线

### 变量

| 列 | 类型 | 说明 |
|---|---|---|
| `heave_valid_fraction` | float (0–1) | 局地窗口内 σ 有效数据占比。衡量诊断可靠性 |
| `glorys_heave_sigma_argo` | float | Argo 异常深度处的 σ₀ (kg/m³) |
| `glorys_heave_sigma_peak` | float | Heave 峰值所在 σ₀ (kg/m³)。低值=浅层变形，高值=深层变形 |
| `glorys_heave_zmin` | float | heave 峰值 σ 面在窗口内的最浅深度 (m)。越小越接近海表 |
| `glorys_heave_m` | float | Heave 幅度 (m)：等密线局地最深与窗口最浅的垂直距离 |

## 分类

| 列 | 类型 | 说明 |
|---|---|---|
| `hotspot_type` | int (1/2/3) | 1=通风型（H≥100m 且 z_min<300m），2=深层隔离型，3=近岸 DO 骤降型 |

## 涡旋匹配（可选）

| 列 | 类型 | 说明 |
|---|---|---|
| `meta_inside_eddy` | bool | 该剖面是否位于 META 涡旋内部 |
| `meta_eddy_list` | str | 所在涡旋的编号列表（逗号分隔） |
