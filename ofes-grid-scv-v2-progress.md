# OFES grid-SCV v2 交接快照(2026-08-15,含外部 AI 审核修复)

## 任务背景

在 `Oxygen-ofes` worktree 的 `feat/ofes-ke-mechanism` 分支上,按用户批准的另一 AI
修复方案重建 McCoy-seeded grid-SCV 检测器 v2。目标不是调高 DO 事件命中率,而是
实现一个真正符合原设计的**分层检测器**,并在看任何 56-event 关联数字前冻结全部
规则(`ofes-mccoy-grid-scv-v2-analysis-lock.md`,458 行,已提交 `141c41e`)。

分层定义: Tier 0 `mccoy_profile_seed`(网格列跑完整 McCoy 门链) → Tier 1
`closed_thermohaline_lens`(种子长成闭合 N2 热盐透镜,含 identity 类) → Tier 2
weak/strong `native_anticyclonic_support`(Nencioli 四规则 + N2 边界 circulation)
→ Tier 3 persistence(描述性追踪,仅门控通过后跑)。

## 提交历史(全部在 feat/ofes-ke-mechanism)

- `d2c356e` Optimize `_ofes_grid_connect_components` overlap search(41/41 等价)
- `141c41e` Lock grid SCV detector repair(v2 锁文件)
- `9e9a779` **Build the McCoy-seeded grid SCV detector**(track.py +3290 行 +
  processing.yml v1.43 `grid_scv_v2:` 段;Tier 0–3 全部代码 + 56-event runner +
  审计/控制 runner + 14 项验证框架)
- `3d1040a` **Repair the grid SCV v2 detector and run gates**(外部 AI 审核
  15 项实现层 bug 修复 + 8 项回归测试 + 年度 runner / Tier-3 tracking)
- **未提交(工作区)**:第二轮复审修复(halo 窗口/逐 voxel 体积/Nencioli
  完整门/节点 4 证据链硬门/年度 runner 5 处/Tier-3 嵌套与 split 规则)

v1 拒绝输出已写 SUPERSEDED 标记:
`Oxygen-cache/ofes_grid_scv_results/{annual_*_localrecert*, event_catalog_final}/SUPERSEDED.txt`

## 外部 AI 审核修复清单(2026-08-15,全部落码)

### 检测器(节点 4 前)

1. **Tier-1 真生长** `_ofes_grid_v2_closed_lens_masks` 重写:种子层锚定后沿相邻
   密度节点上下逐层生长——相邻层必须有自己的闭合轮廓、与上一层掩膜直接或单格
   膨胀重叠、层内同号 spice 非零,任一不满足即停(no-missing-node-bridging)。
   按 node 缓存路径/inside 掩膜。反例修复:三层透镜只放中层种子 → 3 节点对象。
2. **异号永不合并**:同 node 同 path 但异号 seed 各自成层(合并键 = path+sign)。
3. **环带 control 不 clamp**:采样点距最近格点 > 1.2×半对角 → 显式拒绝;
   诊断记录 out_of_window / unique_cells / duplicate_cells。
4. **体积/厚度/质心界面法**:每格点物理厚度 = 最浅节点上界面到最深节点下界面
   (中点界面,边界镜像外推);体积加权质心用真实 voxel 体积;检查 14 用独立
   界面法期望(零厚度倾斜 → 0;200 m 场景 → 200×面积、质心 550),打破自证循环。
5. **Tier-2 身份门**:只有 grid_lens 进入 Nencioli/circulation 门;cyclonic
   technical = 同一套 Nencioli+circulation 门(有限样本 ≥80%)的气旋符号版。

### 运行纪律(节点 4 前)

6. **141/4480 硬门**:audit 函数分母 ≠ 冻结值(141/4480)即 RuntimeError;
   transition 对缺失行/重复 sample_id 报 `gate_failed`;gate 加 count-exact +
   missing==0 硬条件。
7. **resume 门**:抽 `_ofes_grid_v2_load_day_from_dir`——status complete + code/
   protocol/config 三 hash + schema + 窗口边界全一致才复用,按**原始窗口边界**
   重载快照(消除 ±3° 硬编码错位);code hash = **整个 track.py SHA-256**(依赖
   闭包完整);profile cache 序列化 float32→float64(消除 resume 阈值翻转)。
8. **并行**:`_ofes_grid_scv_v2_day_from_snapshot` 把 worker_count 传进
   `_ofes_grid_v2_daily_seeds`(最慢的完整 McCoy 分类阶段)。
9. **配置同步**:processing.yml `version: 1.42` → `"1.43"` + 新键
   `known_background_control_profile_count: 4480`。

### 节点 5 前

10. `run_ofes_grid_scv_v2_event_catalog` 硬前置:节点 4 manifest 必须存在、
    status complete、gate_passed、三 hash 一致,否则拒绝运行。
11. `max_events` 进入运行签名 → 试跑独立目录,不会被正式全量运行误复用。
12. 节点 5 复用节点 4 已校验日片(`_ofes_grid_v2_load_day_from_dir`),只有
    缺失日片才重跑;manifest 记录 reused_validation_day_count。
13. 四级交叉表:event_association 新增 `four_level_class`
    (strong>weak>tier1>profile_only>none),summary 输出 cross-tab。
14. 新增 `run_ofes_grid_scv_v2_annual_catalog`:互斥内域 tile 网格(内域 240 km
    相切、480 km 间隔、无双计湿格点分母)、逐日逐 tile 检测 + resume、跨 tile
    对象去重(同 date 同号中心距离 < 0.5×min(r) 弃小者)、月×1°纬度×σ0
    occupancy null、Tier-3 追踪;启动硬门 = 56-event 目录 complete + 三 hash +
    Tier-1 环带 occupancy 均值非零有限(lock §Run discipline item 7)。
15. 新增 `_ofes_grid_v2_track_objects`(Tier-3):相邻模型日、同 spice 号、中心
    位移 ≤55.66 km、半径/厚度/质心深度相对变化各 ≤0.60;splits/merges 开新
    track;≥3/≥5/≥30 天类统计;三检测家族分别追踪。

### 回归测试(节点 4 gate 新增)

`_ofes_grid_v2_regression_validation` → `regression_validation.json`,8 项:
单 seed 三层生长 / 异号同轮廓不合并 / 越界 control 不 clamp / 零厚度倾斜体积
为 0 / underresolved 不升级 Tier-2 / 141·4480 分母缺样本硬失败 / transition
缺失行 gate_failed / 日片复用在 hash·schema·窗口不一致时拒绝。

## 验证状态(第二轮修复后全量重跑,2026-08-15)

- 解析验证:19/19 PASS
- 回归验证:10/10 PASS
- `py_compile` OK

### 第二轮修复(外部 AI 复审意见,全部落码)

- **事件窗口 halo**:读取域 = 分析域 240 km + 候选点自身 McCoy 背景环带
  230 km + margin 20 = 490 km;`daily_seeds` 只在分析域内列跑完整分类,
  halo 列只供 per-column ring controls(越界拒绝逻辑保留)。合成验证的
  run_day core 移到域中心、含透镜检查用 96 格域(56 格域环带背景几何随
  core 移动会让 N2=0 闭合轮廓不稳定),检查 13 改为异号双透镜(spicy+
  minty 水平相邻,0 轮廓连体也按 lock 不合并)
- **体积逐 voxel**:同节点多个 layer 取 mask 并集;`density_node_count`
  按唯一节点;每个 (node, lat, lon) 独立算上下界面/厚度/体积(非均匀层距
  与局地缺层正确);三维经纬深质心全部 volume 加权;修复 nan*False 清零
- **Nencioli 完整门**:方形周界 8 点(offset b 方形角点,非圆周取整);
  8 点全有效;6/8 切向对齐进入 passed/weak/strong;气旋用对称切向门
  (正对齐 <=2/8)+ cyclonic circulation/core zeta
- **节点 4 硬门**:4480 背景位置逐位 >= 61 个有效 per-position controls;
  9/56、19/56、6/56、11/56 冻结计数逐项 exact 复现(reproduction 硬门);
  manifest 记录 3 个输入 parquet 路径+hash 与 McCoy source archive hash
- **年度 runner 修复**:lat_bands 二维索引维度修复;occupied_frame truth
  判定;√2·R 间隔 tile + 矩形 Voronoi 归属(完整覆盖无双计);occupancy
  分子分母统一 cell-days 口径(cell_key 加 date);start gate 拒绝
  max_events 部分目录
- **Tier-3 修复**:strong 对象嵌套进入 strong/weak/seed-lens 三家族追踪
  (per-family track_id 列);track 延续只在双向 1:1 匹配时发生,split/
  merge 一律全部开新 track
- **回归测试**:新增同节点两 mask 传递 union、非均匀层距逐 voxel 厚度

## 冻结输入(路径均存在且已核验)

- population: `plot_outputs/do/ofes_np30_ke/ofes_delta_do_catalog/20030101_20031231_cf957935d38a/event_diagnostics/ofes_events_21efbe902ab7/event_population/ofes_population_254ae68988a6/population_peak_diagnostics.parquet`(sha b91435f2...)
- 141 阳性/4480 控制: 同根 `event_lifecycle/ofes_lifecycle_f7290df019c2/mccoy_virtual_argo/ofes_mccoy_virtual_argo_a8409b27a056/virtual_profile_diagnostics.parquet`(sha 3a2b4e74...)
- **event_summary(注意:不是 trajectory_3d 的那个 52 行文件!)**: 同 a8409b27a056 目录下 `event_summary.parquet`(56 行)
- 输出根: `/mnt/w2/scratch/user3/Oxygen-cache/ofes_grid_scv_v2_results/`

## 下一步

1. 用户提交第二轮修复 commit(工作区未提交)
2. **节点 4**: `run_ofes_grid_scv_v2_validation(population_path, mccoy_diag_path,
   event_summary_path, output_dir=..., worker_count=16, resume=True)`
   - 14 项解析 + 8 项回归(重跑落盘)→ 4 审计日 → 56 事件日扫描(~5 min/日,
     16 workers ≈ 4.7 hr wall)→ seed_control_audit(141 硬门)→
     matched_background_control(4480 硬门)→ tier_transition_summary
     (attrition 到 weak/strong)
   - 门控: 任一失败即停止
3. **节点 5**(门控通过后): `run_ofes_grid_scv_v2_event_catalog(...)` →
   56-event 关联(primary=Tier-1 core containment;strict=strong-native)+
   四级交叉表 + 同日环带 null + Wilcoxon(单侧预声明 + 透明双侧)
4. **节点 6**: 最终报告 + 年度扫描决策(Tier-3 / annual null 代码已就绪,
   成本 = tile 数×365×单日片时长,由节点 6 决定是否运行)

## 运行纪律(lock 冻结,不可违反)

- 看 56-event 结果前规则已全部冻结(已满足)
- resume 只接受 complete + code/protocol/config 三哈希 + 窗口边界一致;
  worker_count 排除在科学签名外
- DO/event labels 绝不进入 detector gates;绝不以 DO overlap 调门
- 输出原子写入;v1 锁文件与 SUPERSEDED 标记不得改动
