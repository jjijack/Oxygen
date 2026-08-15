# OFES 56-event Formation–Organization–Retention Temporal Audit Lock

- 日期:2026-08-16
- 分支:`feat/ofes-mechanism-clinch`(worktree `Oxygen-mechanism`)
- 状态:设计文档;**冻结(以本文件最终提交为准)后不得依据 56 事件结果调整
  规则**。数据全部来自既有完成运行,零新增 OFES 场读取。

## 目标与三个冻结假说

目标不是提高 SCV 命中率,而是区分三个冻结的竞争假说(三个结果都接受,不以
"更像 SCV"为成功标准):

- **H1 front/strain formation → anticyclonic organization/retention**:
  锋面或应变先形成、拉伸和向下输送富氧水团,随后反气旋旋转增强并延长
  异常寿命;成熟 SCV 是这一过程的稳定端元。
- **H2 front/filament only**:事件主要是锋面 intrusion/filament;反气旋只
  改变平流方向,不提高异常持续性或物质相干性。
- **H3 pre-existing rotational carrier**:反气旋载体在异常 onset 时已经
  存在;异常水团进入已有旋转载体,先有载体、后装载水团。

## 输入与冻结定义(全部复用既有运行,不新建分类器)

| 定义 | 来源 |
|---|---|
| 56 个 strict events(background 环带完整)| `event_association.parquet`(grid_scv_v2_51c7a542a042)|
| start / peak / last detected day | population peak 表(start_date/peak_date/end_date)|
| deep-entry 日 = 首个 depth_mean ≥ 500 m 观测日 | lifecycle daily 的 depth_mean |
| 逐日 rotation_dominated 冻结标签 | lifecycle daily 的 rotation_dominated |
| 逐日 strain_dominated = ~rotation_dominated(有效 Ro/strain 行)| 同上(派生)|
| 连续指标 r_share = \|Ro\| / (\|Ro\| + normalized_strain)| lifecycle daily 的 rossby_number / normalized_strain |
| persistent anticyclonic carrier:27/56 | lifecycle event summary 冻结分类 |
| strict SCV-compatible:6/56(只作描述性标签,不参与身份门)| lifecycle 冻结分类 |
| 52-event 3D trajectories | trajectory3d population 完成运行 |
| 56-event ventilation histories | trajectory_ventilation 完成运行 |

**DO 不参与任何 rotation/strain/SCV 身份门,只作响应变量。** 不得新建
"更容易命中"的 SCV 分类。

## 逐日机制表 daily_event_mechanism.parquet(一行一个 event_id × date)

| 类别 | 字段 | 来源 |
|---|---|---|
| 时间 | days_from_start / days_from_peak / normalized_lifecycle_phase / start,deep-entry,peak,end flag | lifecycle daily + population |
| 水团信号 | delta_do_max / delta_do_mean / delta_do_p90;area_km2、equivalent_radius_km、pixel_count | lifecycle daily |
| 水团信号(peak 口径)| water_mass/heave 绝对分数、water_mass_dominated | population peak 表(事件级)|
| 动力 | rossby_number、negative_subsurface_rossby_number、normal/shear/total strain、normalized_strain、r_share、rotation_dominated、strain_dominated | lifecycle daily |
| 锋生 | 仅三例 onset 表(E000002/239/176),不进入总体表,不为本任务自造公式 | onset 完成运行 |
| 形态 | bbox 纵横比(lat/lon 跨度换算 km)、area_km2、equivalent_radius_km | lifecycle daily lon/lat min/max |
| 轨迹/相干 | 事件级 resolved 标记;逐日 particle−core 误差与 active fraction(如 trajectories 表可导)| trajectory3d 完成运行 |
| 通风 | 逐日 anomaly 粒子 MLD/near/outcrop 接触比例 | trajectory_ventilation daily 表 |
| 分裂/合并 | 事件级 daily_object_count;逐日对象 key 数(以可得列为准,如实标注)| population + lifecycle daily |

形态量只用于区分 compact rotating carrier 与 elongated strain filament,
不设置事后二分类门槛,先报告连续分布。

## 主分析 A:阶段转化

- 冻结三阶段:early = 事件最初 ≤3 个 detected days;peak = peak day;
  late = 事件最后 ≤3 个 detected days(不足 3 天用全部)。
- 输出:early→peak→late 类(rotation/strain)转化表;strain→rotation、
  rotation→strain、始终 strain、始终 rotation 的事件数;各类中 persistent
  carrier、McCoy-compatible(peak 口径 9/56)、resolved-downward 的比例;
  start→peak 的连续变化 Δr_share、Δstrain、Δbbox 纵横比、Δdelta_do_max、
  Δwater-mass(peak 口径与 start 日代理)。
- 核心问题:strain 是否通常先于 rotation?还是 rotation 在 start 已存在?
- 不用"DO 首次过阈值"日期推断物理领先(detector accumulation 会制造时滞)。

## 主分析 B:rotation 是否延长异常寿命

- Primary:persistent anticyclonic carriers(27)vs non-persistent(29);
  Secondary:peak rotation(29)vs peak strain(27)。
- 冻结响应量:event lifetime days;peak 后 water-mass 分量归一化衰减斜率
  (无逐日 water-mass → 用 delta_do_max 归一化衰减斜率替代并标注);
  peak 后 ΔDO AUC;降至 peak 50% 所需时间(未降则 right-censored);
  thermohaline fingerprint retained-day fraction(如轨迹 tracer 表可导,
  否则用事件级 45/52 联合指纹保持集合作次级量);particle−core 误差;
  生命周期平均 IoU/bbox 纵横比(以可得列为准)。
- 统计:事件等权,event-level bootstrap(seed 20260729,10000 次),
  报告 median difference + 95% CI;Wilcoxon/Mann–Whitney 只作辅助;
  n=56 不用复杂模型。
- 稳健性回归:retention_metric ~ persistent_rotation + peak_DO_amplitude
  + core_depth + target_sigma0 + event_start_day(OLS,不解释为因果,只查
  rotation–retention 差异是否完全由峰值强度或深度造成)。

## 主分析 C:形成信号与组织信号的先后

- 对象:长度 ≥7 天的事件;连续指标 r_share 与 strain、delta_do 增长率、
  post-peak 衰减率。
- lag 搜索范围冻结:±min(5 days, floor((event_length−1)/3));要求 ≥7 个
  有效 paired days;每事件先得到一个 lag,再在事件层汇总;不得把所有
  event-days pool 在一起。
- 解释规则(冻结):
  - strain 领先 rotation,且 rotation 组寿命更长 → 支持 H1;
  - strain 领先,但 rotation 不改善寿命 → 支持 H2;
  - rotation 在 start 已存在,随后 water-mass 增强 → 支持 H3;
  - 无稳定顺序 → 多机制异质性,不强行归一。

## McCoy 桥接(三阶段,不扫全生命周期)

- 若生产 McCoy 管线成本允许:每事件 start/peak/last 三日各跑 17 点
  footprint,复用现有生产 self-ring 口径(不调任何 McCoy/grid-SCV/
  rotation 阈值)。回答:McCoy-compatible signature 是在 start 已存在还是
  peak/late 才出现;是否更常出现在 strain→rotation 转化后;native
  velocity confirmation 是否同步增强。
- 成本高时先在 19 个既有 peak-McCoy-positive 事件执行;必须标为条件子集,
  不外推到 56 事件,且不阻塞主分析 A/B/C。

## "submesoscale" 措辞纪律

只有同时报告 equivalent radius、aspect ratio/bbox 纵横比、Rossby number、
rotation–strain balance(并与当地 deformation radius 比较,若能可靠取得)
之后才讨论尺度。否则正式措辞只用:
`front/strain-dominated filamentary transport`、
`deforming water-mass intrusion`、
`anticyclonic subsurface rotational carrier`、
`McCoy-compatible profile signature`;
不写"subgrid SCV"或"submesoscale SCV 已证实"。

## 代表案例

- E000002:rotation but non-material counterexample(persistent carrier、
  patch 与 resolved-3D gate 均失败)——"有反气旋旋转 ≠ 刚性物质透镜"。
- E000239:strain-dominated resolved downward。
- E000176:frontogenetic/reversed-surface or upward contrast。
- 总体统计完成后,再从 population 客观选一个"strain→rotation 且 retention
  增强"的候选仅作示意图,不作证明来源。

## 交付物

- `daily_event_mechanism.parquet`、`phase_transition_table.parquet`、
  `retention_comparison.parquet`、`event_lag_diagnostics.parquet`、
  `mechanism_transition_summary.json`、
  `ofes-event-mechanism-transition-report.md`。
- 报告必须明确裁决:是否存在 population-level strain→rotation 顺序;
  rotation 是否提高异常寿命/相干性;McCoy-compatible 信号更像形成初期、
  成熟后期还是无固定阶段;当前最高可用措辞(mature SCV / SCV-like
  rotational carrier / front-filament transport)。

## 报告纪律

正式数字只回填独立结果报告,不回填本 lock;运行前置 = 本 lock 冻结 +
输入表 hash 记录;统计单位是事件,不是 event-day 或 particle。
