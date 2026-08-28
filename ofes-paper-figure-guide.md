# OFES 机制图库阅读指南

配合 `ofes_figures_preview/` 使用。目录 01–06 为过程与验证产物，07 为投稿主图，
08 为旧版候选图参考，09–10 分别保存核心表格和运行 provenance。全部图的冻结数字
以各 `analysis_summary.json` / 正式报告为准，图中标注与之一致。

grid-SCV S5 背景敏感性已完成；它是每 5 日系统抽样敏感性，不冒充 365 日年度 null。

**通用读图习惯**：事件/异常粒子=红 `#d62728`；对照/背景=蓝 `#4c78a8`；次级对照=青
`#72b7b2`；SCV/ensemble=紫 `#7b61a8`；观测核心=近黑 `#111827`。中位线+分位带
（alpha 0.16）；粒子云细线低透明（alpha ~0.15）。大多数图是**事件等权**口径。

---

## 01_enrichment —— 富集证据

### mccoy_virtual_argo_formal.png（McCoy 门链正式运行）
对 56 个严格 DO50 事件逐个复现 McCoy et al. (2020) 的剖面门链（10-dbar 处理→等密度面
异常→Tukey IQR→spice–N2 共位→高斯形态→弱层结→去第一斜压模动力高度内峰→速度确认）。
- **左上「Separate evidence layers (not a nested funnel)」**：四个独立证据层的计数。
  标题刻意强调"不是漏斗"——各层口径不同，不能把计数当层层筛除。
- **右上「Virtual-Argo intercept sensitivity」**：采样半径/事件等效半径 vs 通过率。
  看曲线随半径怎么掉：通过率对采样位置敏感 = 结构薄/亚格子，是截获敏感性的直接证据。
- **左下**：56 事件按虚拟剖面通过率排序，紫色柱=旧 5-day OFES proxy 子集。看紫柱是否
  集中在高通过率端（两判据的重合度）。
- **右下交叉表**：旧判据 × McCoy center-profile 兼容，「deliberately not equated」——
  两个判据是独立口径，不要互相替代。关键数字：事件 19/56 McCoy-any，背景 20/4480。

### positive_control_known_v2.png（已知 SCV 正对照）
用 McCoy 2020 目录里已知的 SCV 事件当**正对照**：把同一门链放到已知 SCV 上，应当亮；
放到背景上，应当不亮。结果：51 个事件虚拟剖面 4 个兼容 vs 背景 0/240，
Fisher p = 8.5×10⁻⁴。**怎么读**：这是整条 McCoy 线的校准——门链不是"在 OFES 里什么都
点不亮"的空转，它在真 SCV 上亮、在背景上不亮。

---

## 02_water_mass —— 水团分解

### population_phase_space.png（56 事件总体相空间）
四面板：**A** 峰分解（heave DO 贡献的分布）、**B** 次表层动力学（|Ro| × normalized
strain 散点）、**C** 表层衰减、**D** 联合水团指纹（θ/S/spice 对比）。
- **怎么看**：A 面板回答"异常有多少来自密度面抬升(heave)、多少来自水团本身"；
  B 面板看事件在旋转主导区还是应变主导区；D 面板看异常核心与背景的 θ/S/spice 差。
- 叙事角色：给出 56 事件总体的**成分构成**——后续所有机制分析（w_along、MLD、transition）
  都在这张图的坐标系里展开。

---

## 03_formation_transport —— 形成与输送

### trajectory_2d_ensembles.png / trajectory_3d_ensembles.png（三事件定深/三维集成）
E000002、E000239、E000176 三个事件：51 粒子 backward（峰→start）+ forward（start→peak）
释放。紫色细线=backward 粒子云，蓝=forward 粒子云，红线=backward 质心，蓝线=forward
质心，黑线=**观测** ΔDO 核心逐日位置。
- **怎么看**：黑线（观测）与红线/蓝线（模式质心）是否同向同行——这是模式轨迹对观测的
  对齐验证（headline 判据：中位垂直对齐 ≤150 m）。3D 版多一个深度对齐维度。
- 叙事角色：轨迹方法的**可信度前置**——对齐了，后面 56 事件的群体轨迹才敢用。

### trajectory_3d_population.png（56 事件三维群体）
两面板：**Resolved vertical direction**（多少事件的垂直位移被粒子解析）与 **Value added
by resolved w**。
- **怎么看**：第一面板给出可解析垂向路径的计数基础（15 条 downward、3 条 upward）；
  第二面板看解析的 w 相对纯水平平流多解释了多大比例的观测深度变化。
- 叙事角色：F10 单路径图的**总体背景**——单条路径不是孤例。

### trajectory_tracers.png（示踪剂指纹保持）
forward 路径上逐日多示踪剂（θ/S/spice 异常）中位异常随归一化时间的演化 + 总体保持率。
- **怎么看**：曲线是否守住异常符号（保持率门槛 0.75）——守住说明 onset→peak
  过程中热盐指纹与平流携带的水团保持一致。
- 叙事角色：为水团输送提供独立于 DO 幅度的热盐指纹支持。

### trajectory_ventilation.png（通风史回溯）
三组粒子（anomaly 红 / hydrographic control 蓝 / kinematic control 青）在 30/60/90 天
回溯窗口的：直接 MLD 接触率、近 MLD(25 m)率、等密度面出露机会。
- **怎么看**：红线是否系统高于蓝/青线（分位带不重叠=强信号）。anomaly 粒子有更强
  的近期通风史 → 异常水团在 onset 前几个月接触过冬季混合层。
- 叙事角色：**潜沉的上游证据**——与 F5（核心在冬季 MLD 下）、F10（单路径）合成完整
  的"通风→脱钩→潜沉"链条。

### w_validation.png（垂向速度验证）
OFES w 场的逐日连续方程残差与反号对照。**先于一切 w_along 分析的门槛**：w 数据可信，
w_along（沿等密面下沉）的结果才可解释。不是科学结果图，是数据质量证明。

### onset_regime_evolution.png（onset 标量演化）
三个客观机制代表事件 pre-start→peak 的逐日标量：frontogenesis 率、密度梯度、归一化
kinematics。
- **怎么看**：onset 之前 frontogenesis/密度梯度是否先增强——若先于 ΔDO 增强，则
  单事件上支持“锋生组织在先、异常在后”的候选解释。群体 lag 审计未发现统一的
  strain→rotation 日级时滞，因此本图不用于推广固定两阶段顺序。

### E000002 / E000054 / E000176_onset_maps.png（单事件 onset 地图）
单个事件的 onset 日地图：等密度面 DO 场 + 水平 frontogenesis 场 + 流场箭头。
- **怎么看**：异常核心出现在锋面/流轴的什么位置，frontogenesis 大值区与异常核是否
  空间重合。三个事件 = 三种机制形态（旋转型/应变型）的典型代表。

---

## 04_organization_retention —— 组织与滞留

### event_lifecycle.png（群体生命周期）
三面板：**carrier 层级**（persistent anticyclonic carrier / SCV-compatible /
surface-obscured SCV 计数柱，虚线=28 的半数参考）、**旋转日分数散点**（x=旋转主导日
比例，y=旋转日中的弱/反转表层比例，紫=SCV 兼容）、**群体合成**（次表层 |Ro| vs 表层
核心加权 |Ro| 随归一化时间）。
- **怎么看**：左柱看 carrier 有多少（正式生命周期口径 27/56）；中图看 SCV 兼容事件是否集中在高旋转日
  比例区；右图两条 |Ro| 线随时间的分离/汇合 = 深/浅动力耦合的演化。
- 叙事角色：把 56 事件压缩成"carrier 身份 + 生命周期"两维。

---

## 05_exemplar_events —— 单事件案例

### E000002/034/122/222/276_diagnostic_evolution.png（单事件逐日诊断）
单个事件的逐日时间序列：peak 深度、归一化 kinematics（|Ro|、strain 等）。
- **怎么看**：一个事件从 start 到 peak 的深度与动力学轨迹；五事件并排对比可看
  rotation 型 vs strain 型事件的演化差异。

### E000002/034/122/222/276_isopycnal_context.png（单事件等密度面上下文）
单个事件在等密度面坐标下的五类 context 地图（等密度面深度−参考等）。
- **怎么看**：事件核心在等密度面地形（涡/锋面位势）里的位置——核心坐在透镜的哪一侧、
  上游是什么结构。

---

## 06_auxiliary —— 辅助验证

### annual_catalog_overview.png（年度目录总览）
四面板：(a) 全年 ΔDO 发生（空间热图）、(b) DO50 事件峰位置、(c) delivered-depth 分布、
(d) 锁定正式 top-5。
- **怎么看**：DO50 事件的季节/空间分布、深度分布是否集中于 500–900 m 目标层。

### hosoda_benchmark.png（Hosoda 观测基准）
用 Hosoda et al. (2021) 文献观测日期（2003 年 4–5 月，KE 区）在 OFES 场上做高 DO/低盐
双等密度面集成检验。
- **怎么看**：同一 2003 OFES 配置是否能再现 Hosoda 报道的水团/filament 图景。这是
  文献一致性与管线正对照，不是独立模拟或独立事件的外部验证。

### negative_control_comparison.png（负对照）
负对照事件的三面板：精确分解、动力学与层结、匹配分数。
- **怎么看**：负对照事件**没有**同样的分解结构/动力信号——证明检测器的阳性不是
  方法学伪影。与 01 的 positive_control 构成正负闭环。

---

## 07_manuscript_figures —— 论文主图

最终主图由 plot_ofes_manuscript_figures.py 从已完成的 parquet/JSON 产物重组；
脚本只负责读取、统计已有结果和绘图，不重新运行 OFES。

### Figure1_global_scv_ke.png（全球 SCV 与 KE 集中）
三面板保留互斥的全球地图、log-x OR forest 和 KE/OFES 分析区。244 个对象
统一称为 DO-evaluable McCoy SCVs，四类分别为 DO50 17、DO35–50 22、
DO20–35 28、below DO20 177；类别互斥且总和为 244，DO50/DO35+/DO20+
分别为 17/39/67。OR 直接读取正式 sweep（5.86、13.01、20.22；META
约 1.1–1.2）。KE 面板同时标出 KE frame（140–180°E, 25–45°N）与
OFES domain（140–170°E, 25–45°N），并标注 16 of 17 DO50 carriers
occurred in the Kuroshio Extension；这表示空间集中性，不是区域背景概率检验。

### Figure2_ofes_water_mass_mccoy.png（OFES 水团与 McCoy 证据）
四面板依次显示 161 quality events 的 water-mass/heave 分解、绝对水团贡献
分数、联合 θ–S 状态和 event-core McCoy enrichment。56 strict events 高亮；
水团贡献分数中位数为 86.2%（54/56 dominated）和 78.2%（138/161 dominated）。
McCoy 面板使用 19/56、事件等权差 0.144、bootstrap 95% CI [0.076, 0.221]，
主图标注配对单侧 p；双侧审计值留在 caption/结果报告，不再挤入面板。θ 与 S
是一个联合热盐状态的两个坐标，不是三条独立证据；总标题同时概括水团主导与
McCoy-compatible profile enrichment。

### Figure3_ventilation_downward.png（通风史与向下路径）
左侧是正式 30-day trajectory-complete paired subset 的三项 event-equal
ventilation contrasts（direct、within 25 m、isopycnal outcrop），分别读取
hydrographic n=28 与 kinematic n=27 的 CI/p。每个事件先作 anomaly−control，
再事件等权；三个指标相关，六个点不是六次独立验证。右侧三个并列诊断框各自
保留分母：displacement-classifiable 24/56（21 down、3 up），strict
resolved pathways 18（15 down、3 up），以及 w_along n=19、daily mean
主要值 +4.8 m d⁻¹、13/19 downward、nominal p=0.040。它们是 related but non-independent
diagnostics with distinct eligibility criteria，不是漏斗或可相加计数。

### Figure4_E000073_case.png（E000073 配对路径案例）
四面板展示 anomaly 与 hydrographic-control ensembles 的中位线/IQR 以及代表
粒子。轨迹在事件峰值条件下初始化并向后积分，但按自然时间顺序显示：
Trajectories were initialized at the event peak and integrated backward;
reconstructed histories are displayed in chronological order。最早端称
earliest reconstructed position/depth within the integration window，不称
release/start/source。E000073 初始连续 direct contact 为 13 d、累计为 20 d；
Jan 26 与 Feb 1–6 是 intermittent re-encounters，Feb 7 以后才是 final
detachment and major descent。hydrographic control 的正式含义是 0 of 51
trajectory days with direct MLD contact，不是 51 个 control 粒子；图中蓝色点线为
control trajectory 自身路径上的 MLD 中位数，不能用 anomaly-path MLD 替代核验。

### Figure5_rotational_organization.png（反气旋组织与一致 retention tendency）
三面板使用 nested hierarchy（56 strict → 27 persistent anticyclonic
rotational carriers → 6 SCV-compatible → 1 surface-obscured SCV-compatible）、
per-stage McCoy expression（start/peak/last = 10/56、19/56、9/56）和三种
retention reference frames。6/6 与 1/1 是定义内包含关系，不是独立证据；
peak 的 19/56 只能表述为 McCoy-compatible profile expression is most
concentrated at the anomaly peak。三种 carrier−non-carrier decay-slope
差异方向一致，但主文口径为 consistent but statistically limited tendency
toward slower decay，不宣称显著性或因果。

Figure 6 的 surface-eddy 说明与正式数值位于 ofes-surface-eddy-results-report.md；
它使用同一 16-event 分母呈现 15/16 same-sign surface–deep Ro correspondence
对 0/16 closed-SSH containment。

## 08_legacy_paper_candidates —— 旧版参考图（不作为主图）

### F1_enrichment_main.png（富集主图）
左：事件核心 vs 同天背景的 McCoy 兼容率均值柱 + 逐事件配对灰线；右：事件等权均值差
的 bootstrap 分布（10,000 次），橙带=95% CI。
- **怎么看**：左图柱高差 0.144 且配对线大多向右上走；右图分布整体在 0 右侧、CI
  [0.076, 0.221] 不含 0；预设方向的单侧配对 p=9.8×10⁻⁵，透明双侧值
  p=1.97×10⁻⁴；19/56 事件至少包含一条
  McCoy-compatible 虚拟剖面。**这是全论文富集主张的核心图**。

### F3_global_ofes_bridge.png（全球↔OFES 桥接图）
三面板：**(a)** 全球 carrier 率分组柱（log y）——All Argo / META / McCoy SCV 三组在
DO20/35/50 的 ΔDO carrier 频率，柱顶标 n/分母；(b) OR 森林图（log x）——SCV vs
all-Argo 的 OR 5.9/13.0/20.2（红），META 对照方块（青，OR≈1.1）；(c) OFES 56 事件
核心 vs 同天背景的 McCoy 兼容率。
- **怎么看**：左/中回答"SCV 里 DO 极端多不多"（4.5–19 倍载率、OR 5.9–20.2，
  且随阈值收紧单调爬升）；
  右回答对偶问题"DO 极端核里 SCV 签名多不多"（Δ=0.144）。两问互为对偶，suptitle 点明
  "From global statistics to mechanism"。META 行是尺度对照：表面涡尺度几乎不富集，
  富集是 SCV（次表层小尺度）专属。
- 叙事角色：论文"全球统计→OFES 机制"两部分的**衔接图**，数字全部来自冻结 parquet。

### F2_E000193_exemplar.png（代表性全链案例）
E000193 是同时呈现 McCoy-compatible、persistent carrier 与 resolved-downward
信息的代表事件之一。四面板：(a) ΔDO 时间线+旋转状态（灰带=峰日）；(b) 峰日 17 点虚拟 Argo 足迹
（红=McCoy 兼容，4/17）；(c) stage bridge（start 0 / peak 0.235 / last 0）；
(d) 逐日核心轨迹（颜色=旋转状态）。
- **怎么看**：四面板从上到下看——异常怎么长、门链在哪一天亮、核心怎么走。一个事件
  串起整条证据链。

### F4_walong_three_calibers.png（沿等密面下沉三口径）
三个口径（单日/场平均±3d/逐日 Lagrangian）下 core vs ring 的 w_along 连线图，红=resolved
-downward 事件。**看 0 线**：全体事件的配对 Wilcoxon p=0.916/0.574/0.401，说明
沿等密面下沉不是均匀铺在所有事件核心上的背景运动；结合红色 resolved 子集与 F10–F13，
可看出强下沉通道集中在可解析事件中。

### F5_winter_mld_vs_core.png（核心深度 vs 冬季 MLD）
56 点散点，对角线=核心深度=冬季最大 MLD。全部点在对角线**下方**（核心更深），
中位超出 +212 m，范围 108–437 m。
- **怎么看**：所有事件核心都在"当地冬季混合层直接够不到"的深度之下，说明异常到达
  核心深度必须经历脱离表层连通层的输送过程。它本身不是通风时间证明；与
  matched-control trajectory ventilation 联用后才支持近期通风—潜沉解释。

### F6_transition_and_lag.png（相态转化与滞后）
左：early→late 相态转化四类计数，柱上标注 carrier 数；右：strain→r_share 与
strain→ΔDO-growth 的最优滞后直方图。
- **怎么看**：19 个 always-rotation 事件全部是 carrier，27 个 carrier 中 23 个在
  start-day 已为 rotation；右图中位滞后 0 天（符号检验 p=1.0/0.65）。结果支持
  旋转载体常在事件早期存在，同时显示群体中没有统一的 strain→rotation 先后顺序。

### F7_retention_forest.png（滞留三口径森林图）
三行误差棒：ΔDO proxy(fixed-site)、water-mass(fixed-site)、water-mass(moving-core)，
carrier−non-carrier 归一化衰减斜率**中位差**（event-equal bootstrap 95% CI）。
图上同时标出 carrier/non-carrier 样本量与双侧 Mann–Whitney p：17/18,p=0.050；
23/25,p=0.127；23/25,p=0.122。bootstrap CI 与秩检验对应不同 estimand，均透明报告。
- **怎么看**：三个口径的 carrier−non-carrier 差都为正，moving-core 水团项沿载体
  接近守恒；组间检验仍属 suggestive（moving-core MW p=0.122）。一致方向支持
  反气旋旋转对水团信号保持的组织作用。

### F8_mccoy_stage_bridge.png（阶段桥接）
per_stage（各阶段自己位置采足迹）vs fixed_peak（峰位固定采）双面板：start/peak/last
三柱=兼容率均值±SEM，柱顶标 n/56。
- **怎么看**：per_stage 口径下信号在 start 已存在（10/56），peak 增至 19/56，last
  为 9/56。start→peak 事后配对 p=0.042，start→last p=0.624：部分载体早期已有，
  McCoy-compatible 形态在异常峰值期最集中。

### F9_global_scv_regional_structure.png（全球 SCV carrier 区域结构）
按冻结 basin 字段拆成 KE、KE 外太平洋、大西洋、印度洋和南大洋，分别显示
DO20/35/50 在 DO-evaluable McCoy SCV 内的 carrier fraction，柱顶为载体数/分母。
- **怎么看**：DO50 的 17 个载体中 16 个位于 KE；DO20/35 在 KE 外太平洋、南大洋和
  大西洋仍有较弱分布。这张图展示**载体的区域集中性**，不是区域背景校正 OR；它把
  全球统计自然连接到 KE 的 OFES process analogue。

### F10_case_pathway_E000073.png（单条潜沉路径）
E000073 粒子 z+0_r1_a270 的四面板：(a) 水平轨迹（紫云=全部 anomaly 粒子，红段=MLD
接触期）；(b) 深度-时间（粒子 vs MLD，黄带=混合层内，01-14 MLD 骤减脱钩）；
(c) σ0 沿程（26.395→26.495，近等密面）；(d) 密度-深度相空间。
- **怎么看**：b 面板是主线——13 天混合层接触 → MLD 164→23 m 骤减把粒子"留在"层结内
  → 37 天下沉 433 m 到 579 m 峰核。c/d 面板显示沿程密度指纹基本保持（Δσ0=+0.100）。
  ensemble 灰线跟随=代表性路径。**这是把统计图落成具体过程的一张图**。

### F11_pathway_control_comparison.png（路径对照并排）
E000073 异常粒子与 hydrographic control 四面板双线对照；对照在同日、同深度释放，
来自 120–240 km 环带并匹配同一目标 σ0（约 26.5）。(a) 水平轨迹（红=异常、蓝=对照，
紫/蓝云=各自 51 粒子 ensemble）；
(b) 深度-时间；(c) σ0 沿程；(d) 密度-深度相空间。
- **怎么看**：两条粒子密度几乎相同（Δσ0 都 ≤0.10）——差异**只在通风史**：异常粒子
  泡混合层 20 天（51 条 control 全部 0 天）并下沉 433 m；对照粒子在层结内几乎不动
  （68 m）。对照组把"通风→潜沉"从单条趣闻变成配对实验。
- 叙事角色：F10 的对照组，trajectory_ventilation（群体口径）的个体版。

### F12_pathway_event_array.png（多事件路径阵列）
三条 resolved-downward 事件（E000133 rotation / E000164 strain / E000172 strain）
的 2×3 阵列：上排=水平轨迹（紫云=ensemble、红线=代表粒子、红粗段=接触期），下排=
深度-时间（红=代表粒子、青=MLD、灰=ensemble 中位、黄带=接触期）。
- 展示粒子按固定规则选择：每个异常 ensemble 中 direct-MLD-contact 天数最多者；若全为
  0，则取有效记录最长者。该规则写入图题，图不代表随机粒子抽样。
- **怎么看**：三个事件形态互补——E000164 泡混合层 65 天、13→579 m（极端沉降）；
  E000172 接触 40 天、149→603 m；E000133 零接触、层结内 529→603 m 缓沉。
  表明 resolved-downward 路径可跨 rotation/strain 两种动力语境复现，F10 不是孤例。
- 叙事角色：F10/F11 的多事件推广。

### F13_pathway_3d.png（3D 潜沉路径）
E000073 的双视角 3D 轨迹（lon-lat-depth）：左=斜俯视（看路径水平形态），右=低角度
侧视（看深度演变）。红=异常粒子（加粗段=混合层内 20 天）、蓝=对照粒子、青虚线=沿
程 MLD、紫云=anomaly ensemble、每 10 天一个日期标记。
- **怎么看**：侧视图最直观——红色粒子从 ~150 m（混合层内）一路向下滑到 579 m 峰核，
  蓝对照在层结内几乎平走；沿程 MLD 线在 01-14 骤减后与粒子分道扬镳。z 轴向下为正
  （海面向下），所以"下沉"在图上表现为 z 增大。
- 注意：lon/lat 用原始度坐标（35°N 处 1°lon≈91 km、1°lat≈111 km），x 方向有
  ~18% 视觉压缩；保存为 png 后视角固定，如需交互旋转可出 jupyter 版。

### F14_surface_containment.png（PET 包含率）
两面板：(a) 四组柱（strict-56 effective/speed、quality-161 effective/speed），
蓝=PET 可分析分母、红=被 PET 涡包含数（全 0）；strict 为 31/56，quality
敏感性为 58/161，其余分别是 25、103 个不可判定。(b) 仅对 31 个可分析严格事件
绘制最近 PET 涡中心距离/有效半径直方图（log x）。
- **怎么看**：strict 与 quality 两个口径中包含数分别为 0/31、0/58；严格事件到
  最近目录涡的中位距离为 5.9×有效半径，说明事件核心与闭合 SSH 目录存在
  系统脱钩。结合 F16 的表层 Ro 同号结果，这不是“没有表层动力”，而是相关
  次表层载体没有组织成传统目录要求的闭合 SSH 涡。

### F15_surface_null_two_sided.png（双侧 null 森林图）
左：两口径森林图（annual 月×1° 与同日 120–240 km 环带）的 core−null 均值与
bootstrap 95% CI——**两个 CI 全在 0 左**（annual p ≈ 1.7×10⁻⁶；同日环带
p ≈ 8.9×10⁻⁵，均为双侧，paired n=31）；右：逐事件 core−null 直方图
（annual 30/31 负、ring 20/31 负，其余为零，无正值）。
- **怎么看**：PET 可分析事件核心优先位于闭合 SSH 目录涡之外，占有率显著低于
  两种背景。这是“系统脱钩”而非少数案例未匹配。图中使用 2026-08-14
  审计修正后的双侧正式口径（原单侧 greater p≈1.0 掩盖了对侧信号）。

### F16_surface_rotation_pet.png（rotation×PET 交叉）
左：29 个峰日旋转主导事件中 16 个 PET 可分析且均未包含目录涡（0/16），
13 个不可判定；右：表层核心加权 Ro 与深层极性同号 26/29。
- **怎么看**：这是“表层有旋转印记、闭合 SSH 目录却脱钩”的关键桥梁。相关
  次表层旋转载体能够产生同号表层 Ro 响应，却未组织成 PET 所要求的闭合 SSH
  等高线涡，从而为全球 SCV 强关联与 META 弱关联的反差提供模式侧解释。

### F17_glorys_reproduction_stratification.png（GLORYS 重现分层）
三个 ΔDO 阈值下比较 GLORYS reproduced(89)与 missed(71)的观测 SCV 载体率；
另 85 个无可比较重现判定的剖面不进入分层。
- **怎么看**：DO20/35/50 在 missed 组均更高（47.9/29.6/18.3% vs
  19.1/11.2/3.4%）。全球 SCV 富集并非只由模式里容易重现的大尺度对象驱动；
  携带强 DO 异常的 SCV 反而更容易落在模式分辨率盲区。该图为探索性分层。

### F18_analysis_funnels.png（分析分母漏斗）
左：McCoy 4084 目录 → 4066 时间窗 → 263 BGC/DO 联合分析交集 → 245 当前 DO 可算
→ 244 冻结主 OR 分母；右：OFES 161 quality 普查 → 56 strict DO50 → 52 个 3D
population → 34 条轨迹完整 → 28 个 30 天配对通风事件。
- **怎么看**：全球 244/245 的唯一差异是 P32111；把它作为非载体加入只把 DO50 OR
  从 20.22 改为 20.13。OFES 后段缩小来自时间域和对照可用性，而非事后机制筛选。

### F19_ventilation_stratification.png（通风史分层森林图）
三指标（direct/near/outcrop）× 三子集（all/strain/rotation，vs hydro 对照）
的事件等权差森林图，标注 p 值和分母：全体 n=28、strain n=14、rotation n=14。
- **怎么看**：全事件高质量轨迹的 anomaly−control 差为正；strain 行在该子集中
  更高（direct +0.150,p=0.015）。regime 分层在全部可算事件中会减弱,因此主结论是
  异常粒子的近期通风史增强,分组差异作为过程异质性的描述。

### F20_formation_retention_division.png（形成-滞留散点）
54 事件散点：x=anomaly−control 30 天 MLD 接触率差、y=滞留衰减斜率
（moving-core，越正=保持越好）；紫=rotation、橙=strain、空心=非
trajectory-complete diagnostic subset、黑圈=carrier。
- **怎么看**：carrier（黑圈）整体位于更慢的 moving-core 衰减一侧,而近期通风史
  跨组异质。图把水团进入深层与后续旋转组织放在两个互补坐标上,同时保留轨迹
  完整性差异。

### F21_walong_resolved_subset.png（w_along resolved 子集）
左：population flag 的 resolved-down 子集（n=19）三口径 core/ring w_along 柱
（+7.3/+3.9/+4.8，
Lagrangian p=0.040 nominal）；右：观测垂直方向评估（21 降 / 3 升 / 32 低于
25 m 门槛）。
- **怎么看**：可分类事件以下沉为主,强 core−ring 对准集中在 resolved 子集；总体
  事件并没有共同的运动学对准背景。这里的 n=19 与轨迹严格门下 15 条
  high-confidence resolved-downward 使用不同判定链，不能混成同一分母；p=0.040
  为探索性子集结果。

### F22_grid_scv_s5_background.png（grid-SCV S5 背景可达性）
左：四个 detector tier 的逐月日均对象数；中：按 σ0 的背景占据率，灰带为 56 个
DO50 事件的 5–95% 密度范围；右：按事件月份、纬带和 σ0 匹配后，56 个事件核心的
预期命中数与实际 0 命中。
- **怎么看**：背景有 549 个 Tier-1 对象，但面积占据率仅 0.227%；事件匹配后的
  Tier-1 预期命中数只有 0.237，出现 0 个的参考概率约 79%。因此 0/56 说明成熟闭合
  grid-SCV 不是事件核心的常见解析形态，却不足以否定 SCV-like/亚格子载体。
- 叙事角色：解释 grid detector 零重合为什么不推翻全球 McCoy SCV 富集，并把主证据
  放回剖面富集、原生速度支持子集与 Lagrangian 过程链。
