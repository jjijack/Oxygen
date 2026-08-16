# OFES grid-SCV v2 Systematic-5-Day Background Sensitivity(冻结协议)

- 日期:2026-08-16
- 分支:`feat/ofes-ke-mechanism`(worktree `Oxygen-ofes`)
- 状态:**独立命名的敏感性协议,冻结于运行前**。不修改
  `ofes-mccoy-grid-scv-v2-analysis-lock.md` 的任何条款;本运行**不声称
  完成 365 日年度 null**——论文主 null 仍是已完成的同日 120–240 km
  事件环带对比。正式数字只回填
  `ofes-grid-scv-v2-results-report.md`(第五节),本文件不回填。

## 一、动机与定位

2026-08-16 试跑结果误读已修正:2003-01-01 的 7 个对象实为 5 个
grid_lens(Tier-1 闭合热盐透镜)、3 个 weak-native、1 个 strong-native、
0 个 well-resolved——检测器在背景海域能够检出闭合透镜乃至原生速度
支持对象,只是不覆盖 56 个 DO 事件核心。因此年度背景抽样值得做。

既有 annual runner 的两个缺陷决定本协议:
1. **occupancy 混合所有对象等级**(grid_lens/underresolved/broad_structure
   合并),不能作为论文的 Tier-1 年度 null;
2. **存储失控**:每 tile-day 的 profile_cache.parquet 占 99.97%,
   1 天 63 tile = 21.7 GiB,全年 ~7.7 TiB。

本协议以系统抽样 + 分 Tier 统计 + 紧凑存储运行年度背景敏感性。

## 二、冻结采样设计

- 日期:2003-01-01 起每 5 天一天,共 **73 天**(doy 1, 6, 11, …, 361)。
  确定性序列,无随机种子;与季节周期(365 天)不对齐(365 不被 5 整除),
  每月覆盖 5–7 天。
- 空间:与既有 annual runner 相同的 63 个 √2·R 规则 tile(相同布局
  函数,相同硬验证)。
- 检测:`detect_ofes_grid_scv_v2_day` 原样调用,**不调任何阈值/规则/
  门链**;DO、事件标签与表层涡目录不参与。
- **不做 Tier-3 逐日寿命追踪**:5 天间隔与逐日关联不兼容。对象表只报
  单日分级统计。
- 统计口径:事件/对象等权;occupancy 分母 = tile-day Voronoi 湿格点
  cell-days(与既有 runner 同口径),分子按 detector tier 分列。

## 三、分 Tier occupancy(冻结定义)

| tier 列 | 定义(对象级 flag) |
|---|---|
| tier1 | `tier1_identity == 'grid_lens'` |
| weak_native | `weak_native_support` |
| strong_native | `strong_native_support` |
| well_resolved | `well_resolved_grid_lens` |

分子 = 该 tier 对象(同日同号中心去重后)的 voxels 按日期+坐标键的
cell-days;四个 tier 共享同一分母。汇总表同时给出四列
`occupancy_fraction`。

## 四、紧凑存储(冻结清单)

每 tile-day 目录**只保留**:`day_summary.json`(含 compact 标记)、
`objects.parquet`、`voxels.parquet`、`node_wet_counts.parquet`、
`seeds.parquet`、`prefilter_scan.parquet`、`profile_only.parquet`、
`node_support.parquet`、`layer_diagnostics.parquet`。
`profile_cache.parquet` 在写完 day_summary 后立即删除。
S5 全量预计 ≤ 20 GiB(对比不紧凑的 1.55 TiB)。

resume 校验与既有 `_ofes_grid_v2_load_day_from_dir` 的区别:不要求
profile_cache 存在,要求 day_summary 携带 `storage_mode == 'compact'`
标记与三哈希一致;仅用于本 S5 runner,不影响既有 v2 目录的 resume 语义。

## 五、并行与基准(冻结流程)

- tile-day 级并行:runner 参数 `tile_day_workers`(每任务 worker_count=1,
  实测 worker 数几乎无增益)。
- 先基准:同一批未完成 tile-day,8/12/16 进程各测一次,报告 tile-day/s
  与墙钟外推;**不预设 24 或 32 一定更快**,以基准选最优进程数。
- 全部 4599 tile-days 完成后一次性聚合(去重 + 分 Tier occupancy)。

## 六、复现门(冻结,运行前置)

2003-01-01(新代码哈希下重算)必须与既有试跑逐值复现:
7 个对象、tier1_identity = 5 grid_lens / 1 underresolved /
1 broad_structure、weak_native_support 3、strong_native_support 1、
well_resolved_grid_lens 0。任何不一致 → 停止,先查原因,不得带病跑
73 天。

## 七、解释规则(冻结)

- 本敏感性只回答"年度背景里各 Tier 对象的多寡与分布(月×纬带×σ0)",
  不回答寿命/路径/闭合动力学;
- 报告句式上限:"S5 抽样下 Tier-1 年率估计 X(CI),对比同日事件环带
  null 与 McCoy 观测层富集"——不写"OFES 全年 SCV 目录完成";
- 与 56 事件 Result C 的关系:S5 背景里 Tier-1 对象多,说明"背景有
  透镜、事件核心不在其中",进一步支持"事件核心不是成熟闭合 SCV"的
  多机制解释;对象少则如实报告零率。
