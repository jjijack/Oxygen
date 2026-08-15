# OFES 1/30° 沿等密面速度审计(Isopycnal Velocity Audit)Lock

- 日期:2026-08-16
- 分支:`feat/ofes-mechanism-clinch`(worktree `Oxygen-mechanism`)
- 状态:探索性表征(exploratory characterization),非门控冻结实验;本文件是预注册说明,不设事后通过门槛。

## 背景与定位

GLORYS 1/12° 侧已确认 null(2026-06-26):对 Argo V/I 异常算沿等密面运动学俯冲率
`w_along = u·∂z_σ/∂x + v·∂z_σ/∂y`(z 向下为正;正 = 朝深),单日、±3 天先平均场(欧拉)、
±3 天逐日平均(拉格朗日)三重口径都区分不了 V/I(p=0.65~0.87),且时间平均后 p 更大。
机制复核:V 类比 I 类斜率更陡、流速更快,但流向相对等密面斜率的方向余弦 cos≈0,
把"陡×快"抵消到 w_along≈0;正负占比各 ~50%(硬币)。判读:1/12° 平滑掉了亚中尺度
锋面对准,**实证支持需 OFES 1/30°**(submesoscale-permitting)。

本审计是那条 GLORYS 证据链的 OFES-side 落点:在 1/30° 自由运行场里,56 个 OFES
ΔDO 事件核心处沿等密面俯冲是否可分辨、流向是否对准等密面坡向,以及与动力 regime
的关系。结果用于机制实锤的主路径证据,不是 SCV 等价验证(那已由 grid-SCV v2 裁决
关闭)。无论结论是"1/30° 能分辨对准"还是"1/30° 仍随机对准",都是机制节所需证据:
前者支持亚中尺度俯冲载体,后者把机制责任转给已解析的垂直路径(18 条 resolved
pathways)或未解析过程。

## Estimand

- 主估计量:`w_along`(m/day)在事件核心盘与同日 120–240 km 背景环带的分组分布,
  及因子分解 `w_along = |∇z_σ| × |V_σ| × cos`(斜率 × 等密面流速 × 对准余弦)。
- 次估计量:沿等密面 DO 平流倾向代理 `V_σ·∇_σ(DO)`(向高 DO 方向平流为正)与
  对准余弦 cos_DO。
- 区域口径:core = 事件核心 20 km 盘;ring = 120–240 km 环带(同日同场匹配对照)。

## 输入

- 56-event 关联表:`event_association.parquet`(grid-SCV v2 节点 5 完成运行
  `grid_scv_v2_51c7a542a042/`),列 event_id/peak_date/core_lon/core_lat/
  core_depth_m/target_sigma0;SHA-256 记入 manifest。
- population peak 表:`population_peak_diagnostics.parquet`(59 行),join
  kinematic_regime/contribution_regime。
- 可选 trajectory3d population 表(18 resolved pathways),join 后得到
  resolved_down/up/other 标记;缺失时该项报告列留空。

## 窗口与算法

- 事件日集:peak_date ± 3 天(7 天);核心位置固定(与 GLORYS 口径一致,不做逐日追踪)。
- 场窗口:core ±(240 km 外环 + 25 km 平滑余量),全部交付深度(75 层);
  窗口超出交付域时由 loader 裁切并记录裁切比例,ring 不完整则 ring 指标置 NaN。
- σ0:TEOS-10(`_ofes_sigma0_volume`);目标面 = 事件 target_sigma0;交点选择
  参考 core_depth_m(最近交点,保分支一致);crossing_count 记录;core 盘交点
  可用率 < 0.8 的 event-day 记 NaN 并计数。
- 平滑:NaN-safe 均匀滤波,主口径 smooth_km=10 km(约 3 格点);**映射场落盘**
  (z_σ/u_σ/v_σ/do2_σ 逐 event-day NetCDF),平滑半径与 core/ring 半径的敏感性
  事后从落盘场重算,不重读 OFES。
- 度→米:逐格点 `approximate_degree_length(lat)`;∂z/∂x=(∂z/∂lon)/(米每度经)。
- 时间聚合三重口径(镜像 GLORYS):单日(peak 日)、±3 天先平均场(欧拉)、
  ±3 天逐日平均(拉格朗日)。

## 门(轻,探索性)

- G1 输入:56 行齐全,三输入表 hash 与 manifest 记录一致。
- G2 覆盖:每个 event-day 快照加载出有限 T/S/u/v/do2;无"全部日不可用"事件;
  可用性分项计数报告。
- G3 量级合理性:pooled |w_along| 全有限且 ≤ 500 m/day(物理界);高斜率格点
  |w_along| 量级报告(O(10) m/day 量级为 GLORYS 基准)。
- G4 解析零检验:合成平直等密面 + 均匀流 → w_along = 0 至机器精度(函数内建)。

## 报告项(不设门槛)

- 三重口径 × core/ring 的 w_along、|∇z_σ|、|V_σ|、cos、正占比(硬币检查)、
  V_σ·∇_σ(DO)、cos_DO。
- regime 对照:rotation/strain、resolved down/up/other、water-mass/heave
  (双侧检验 + bootstrap,事件等权)。
- core−ring 配对对照(同场背景)。
- 方向一致性审计(报告,非门):E000239(三维 ensemble 下沉 96.3 m)与 E000176
  (上浮 24.6 m)core w_along 符号。

## 措辞纪律

- 本审计证明"解析场在 1/30° 的可分辨对准"与否,不单独证明俯冲发生(无拉格朗日
  约束);与 52-event 三维轨迹、通风历史互相引用时各自保留原有门控结论。
- 试跑/小样本先验不写入门;正式数字只回填结果报告,不回填本 lock。
