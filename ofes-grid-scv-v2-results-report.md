# Grid-SCV v2 Results Report(节点 6)

- 日期:2026-08-16
- 分支:`feat/ofes-ke-mechanism`(worktree `Oxygen-ofes`)
- 依据:lock `ofes-mccoy-grid-scv-v2-analysis-lock.md` +
  裁决 `ofes-grid-scv-v2-gate-adjudication.md`(sha `d9911a4d…`,
  冻结于重跑前)。**本文档是最终数字的唯一回填位置**,裁决文档不回填。
- 正式运行:commit `50ed103`,code hash `3b00951a…`,manifest schema v3。

## 一、节点 4 门控(正式,2.90 h)

hard_gates 12/12 全绿,`gate_passed = True`:

| 硬门 | 结果 |
|---|---|
| 19 项解析 + 12 项回归 | PASS |
| 分母完整(141 阳性 / 4480 背景) | 无缺失行 |
| 预筛 recall | 1.0(141/141) |
| stored_pipeline_replay | 1.0(141/141 逐位复现冻结分类) |
| 事件级重放 | 9/9、19/19、6/6、11/11 exact |
| 输入与代码 hash | 完整记录(含裁决文档 sha) |

reported_items(不参与 gate_passed):

| 报告项 | 值 | 含义 |
|---|---|---|
| event_reference_nearest_cell | 0.9504(134/141) | 格内位置/插值敏感性 |
| production_self_ring | 0.6312(89/141) | 生产检测器对冻结阳性的 Tier-0 转移性能 |
| 背景假阳性率(per-position ring) | 0.29%(13/4480) | lock 字面口径 |
| 背景 insufficient 数 | 88 | 生产口径 background_iqr 正常记录 |
| 生产口径事件级 | center 9、any 17、center-vel 6、any-vel 10 | 冻结 9/19/6/11 的网格转移 |
| Tier-1 扩展(141 阳性) | 0/141 | 报告 |
| 支撑列/3×3(7 个翻转) | 5/7、2/7、7/7 | 描述性审计,非门 |

## 二、节点 5 探索性表征(post-gate exploratory)

输出 `grid_scv_v2_51c7a542a042/`,56 事件,56 天日片全复用(11 min)。

- 四级交叉表:none 45、profile_only 11、tier1 0、weak 0、strong 0
  ——**事件核心无一被 Tier-1 透镜包含**
- 网格对象:grid_lens 15、weak_native_support 11、cyclonic_technical 1;
  **well_resolved(物理厚度 ≥100 m)= 0**
- Tier-0 种子(生产口径,56 天全域):profile_only 5623 个
- 环带 null(event-equal、面积加权、同 σ0 节点):
  - tier1:core−ring 均值 −0.0020,bootstrap95 [−0.0053, −0.00005],
    双侧 Wilcoxon p = 0.068
  - weak:均值 −0.00048,p = 0.11
  - strong:全零(对象数 0)
- 解读:lock 预判的 Result C 形态——Tier-0 在网格上大量存在
  (5623 种子、89/141 转移),但冻结的 McCoy 阳性是薄/亚格子结构,
  几乎从不长成事件核心处的 well-resolved 闭合透镜;透镜存在处
  (15 个)也无事件级富集(null)。**不声称 v2 通过 McCoy 等价验证**;
  如实报告:重放口径逐位复现旧管线,生产网格口径转移性能如上。

## 三、证据链完整性

- 失败运行归档:`validation/failed_run_20260815_production_caliber/`
  (生产口径 88/89、18/19、10/11 原始证据,未删除)
- 正式运行:`validation/`(manifest schema v3,两区 gate)
- 长跑期间 track.py 零改动;唯一新提交 `50ed103`

## 四、年度扫描决策(待用户)

lock 的节点 6 剩余项:是否启动 `run_ofes_grid_scv_v2_annual_catalog`
(2003 全年逐日逐 tile 检测,Tier-1 环带 occupancy 启动门)。鉴于
56 事件结果为 Result C 形态,年度扫描的科学价值取决于目标——若为
"生产检测器的全球/全年 SCV 目录统计"则可跑;若为"验证 McCoy 等价"
则不需要(已由节点 4 裁决关闭)。决策待定。

### 2026-08-16 工程试跑完成(63 tile-days = 2003-01-01)

- 运行:`/mnt/w2/scratch/user3/Oxygen-cache/ofes_grid_scv_v2_results/
  annual`(schema 2,manifest complete;`reused_complete_run: False`,
  全量重算)。63/63 tile-day 完成,零错误;**墙钟 5.26 h**。
- **产物逐值核实(2026-08-16 修正,前版误读为"0 Tier-1")**:当日 7 个
  对象中 **5 个 grid_lens(Tier-1 闭合热盐透镜)**、1 个 underresolved、
  1 个 broad_structure;速度支持 **3 个 weak-native、1 个 strong-native、
  0 个 well-resolved**。**正确结论:检测器在背景海域能够检出闭合透镜
  乃至原生速度支持对象,但这些对象没有覆盖 56 个 DO 事件核心。** 这比
  "检测器在 OFES 里基本什么也找不到"重要得多,也使背景抽样更值得做。
  (tier1_identity = 5 grid_lens / 1 underresolved / 1 broad_structure;
  weak_native_support 3、strong_native_support 1、well_resolved 0。)
- **已知缺陷:annual occupancy 混合所有对象等级**。现有 runner 把
  grid_lens、underresolved、broad_structure 的 voxel 合并成一个
  occupancy 无 Tier 维度:01-01 总 occupied cells 25,508(其中 Tier-1
  11,249、weak-native 11,161、strong-native 6,989、well-resolved 0)。
  因此现 annual_occupancy.parquet **不能直接作为论文的 Tier-1 年度
  null**,必须按 detector tier 分开统计后再跑 S5(见第四节 S5 协议)。
- annual_tracks 11(全部 boundary_censored=False,strong/weak/seed-lens
  三类 3/5/30 天档全 0——单日试跑无跨日轨迹,不代表检测器无产出)。
- **全年成本外推**:单 tile-day ≈ 5.01 min(301 s)。全年 = 63 tiles ×
  365 days = 22995 tile-days ≈ **1920 h 串行**(~80 天)。当前 runner
  不并行 tile-day。若需全年目录:
  - 串行不可接受;**必须 tile-day 级并行**(63×365 独立任务,理论
    并行度 = 机器核数;按 16 worker 可压到 ~120 h,~5 天)。
  - **存储是比计算更严重的问题**:1 天 63 tile 占 21.7 GiB,其中
    profile_cache.parquet 占 99.97%(21.712 GiB)。照现状 S5 73 天约
    1.55 TiB、全年约 7.7 TiB。紧凑 runner 只保留 day/tile summary、
    objects、voxels、node wet counts、哈希与完成标记,S5 应在几个 GiB
    以内(见 S5 协议)。
- **决策口径修正**:四个原选项(冬季连续块/改旧 lock/半年/全年)全部
  否决——冬季块季节偏置且不匹配全年事件月份,半年不完整,全年收益不
  匹配成本,改旧 lock 冒充 365 日正式 null 不许可。**采用独立命名的
  systematic 5-day sampled sensitivity**(协议见下)。
