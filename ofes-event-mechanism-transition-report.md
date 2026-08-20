# OFES 56-event Formation–Organization–Retention 时间序审计结果报告

- 日期:2026-08-16
- 分支:`feat/ofes-mechanism-clinch`(worktree `Oxygen-mechanism`)
- 依据:lock `ofes-event-mechanism-transition-analysis-lock.md`(冻结设计)。
  本文档是正式数字的唯一回填位置,lock 不回填。
- 完成运行:`/mnt/w2/scratch/user3/Oxygen-cache/ofes_mechanism_transition/
  ofes_mechanism_transition_be1b4e1863c5/`(manifest schema 1,complete)。
  数据全部复用既有完成运行(lifecycle/population/trajectory/tracers/
  ventilation/McCoy),零新增 OFES 读取。

## 一、主分析 A:阶段转化(early→late 类,56 事件)

| 转化 | n | persistent carrier | SCV-compatible | McCoy-any | resolved-down |
|---|---|---|---|---|---|
| strain→strain | 19 | 0 | 0 | 6 | 6 |
| rotation→rotation | 19 | 19 | 4 | 8 | 6 |
| strain→rotation | 8 | 2 | 1 | 3 | 4 |
| rotation→strain | 10 | 6 | 1 | 2 | 3 |

- **rotation 在 early 阶段(前 ≤3 个 detected 日)已存在的事件 =
  29/56**(19 始终 + 10 转出);strain 起始 = 27/56(19 始终 + 8 转入)。
  两类起始形态几乎各半。
- strain→rotation(8)只有 2/8 成为 persistent carrier——**后期转入旋转
  并不带来滞留身份**;persistent carrier 几乎全部(19/19)在 always-
  rotation(early 与 late 都为 rotation)。**口径注意(2026-08-16 修正)**:
  early 三日分类 ≠ 真正 start-day——按逐日 rotation_dominated 的
  start-day 口径,persistent carrier 中 start-day 即为 rotation 的是
  23/27,early 三日口径为 25/27;"19/19"只是 always-rotation 组内的
  数,不应读成"carrier 从第一天起就旋转"。SCV-compatible(4/6)与
  McCoy-any(8/19)集中于 always-rotation。
- start→peak 连续变化:Δr_share 中位 0.0(正占比 0.48,无系统增长);
  Δnormalized_strain 中位 +0.011;Δbbox 纵横比中位 +0.18(对象向 peak
  变得更细长);ΔDO 中位 +7.17(构造性增长)。

## 二、主分析 C:形成-组织 lag(事件级)

| 对 | n(≥7 天事件) | 中位 lag(天) | 正/负/零 | 符号检验 p |
|---|---|---|---|---|
| strain → r_share | 26 | 0.0 | 11/11/4 | 1.0 |
| strain → DO 增长 | 23 | 0.0 | 8/11/4 | 0.648 |

无 population-level 的"strain 先、rotation 后"顺序;lag 分布对称且中位
为 0。限制:事件中位长度 ~10 天、日分辨率、k_max ≤ 5,只能排除 ≥1 天的
稳定时滞,不能排除天内或未解析过程。

### 30 天通风史精确数字(2026-08-16 回填,数据源 ventilation 运行
`bd4e029be65e` 的 group_comparison,事件等权差 vs σ0-matched
hydrographic control,bootstrap 95% CI + 配对 Wilcoxon):

| 子集(30 天窗口) | direct MLD contact | near MLD | outcrop |
|---|---|---|---|
| 全事件(n=28) | **+0.109 [0.041,0.190] p=0.0049** | +0.129 p=0.0010 | +0.116 p=0.0017 |
| strain(n=14) | **+0.150 [0.045,0.277] p=0.015** | +0.177 p=0.009 | +0.165 p=0.006 |
| rotation(n=14) | +0.068 p=0.158 | +0.082 p=0.046 | +0.067 p=0.158 |

- 60/90 天窗口 n≤11,统计力崩溃;**群体级"近期通风"的定量证据止于
  30 天**;90 天只有案例级(路径图 F11)。
- **口径注意(2026-08-16)**:上表来自 ventilation 运行中 30 天控制匹配
  完整的事件对(n=28,strain/rotation 各 14)。扩展到全部 54 个可算事件后,全事件
  通风差仍约 +0.11(稳健),但 strain/rotation 分层减弱(0.121 vs
  0.103)——分层只在高质量轨迹子集成立,未通过子集的 rotation 事件
  通风史反而最高(0.185)。**论文只主张"异常粒子整体有近期通风史",
  不主张群体级 strain/rotation 分层**。
- carrier 的近期通风差整体偏弱(0.093 vs non-carrier 0.131),提示这类
  旋转载体更多整合较早通风或沿等密面远距离输入的水团;rotation 事件内部
  仍保留明显异质性。

## 三、主分析 B:retention 对比(事件级,event-level bootstrap)

| 响应量 | carrier(27)vs 其余(29)中位 | diff | bootstrap 95% CI | MW p |
|---|---|---|---|---|
| lifetime(days) | 11 vs 6 | +5 | [−3, +7] | 0.74 |
| post-peak 衰减斜率(归一/day) | −0.0139 vs −0.0261 | +0.012 | [−0.003, +0.051] | **0.050** |
| post-peak AUC | 8.75 vs 4.19 | +4.57 | [−0.94, +6.42] | 0.117 |
| DO 指纹保持日占比 | 1.0 vs 1.0 | 0 | [0, 0] | 0.31 |

Secondary(peak rotation 29 vs 27):lifetime +5(p=0.75);衰减斜率 +0.011
(p=0.125);post-peak AUC +5.5,CI [1.05, +7.38],p=0.015。
time_to_half 53/56 右删失(检测期内未降到 peak 50%),不参与比较。

**解释纪律(2026-08-16 修正)**:AUC 受观测持续时间直接影响,不得写成
"衰减速度显著更慢";直接衰减证据只有 slope——代理版 p=0.050 仅名义
边缘(CI 跨 0),真实水团项版 p=0.127,retention 结论保持"有限"。

**真实水团项版,两种几何口径(2026-08-16 按审核修正)**:

- **fixed-site persistence(Eulerian,运行
  `ofes_daily_water_mass_afa024fe7a6e`,peak 日与 population 逐位一致
  56/56 零差)**:固定 peak 核心位置 + 事件 core_depth_m 参考深度。
  carrier 组归一化水团项衰减斜率 −0.118/day vs 非 carrier −0.254/day
  (慢约 2.2 倍,方向与 ΔDO 代理一致,median diff +0.136,bootstrap CI
  [−0.001, +0.377],MW p=0.127);post-peak AUC 0.94 vs 0.72(p=0.18);
  wm_peak 75.6 vs 65.6 μmol kg⁻¹。测的是同一固定站点上的信号持续,
  **不是移动载体的 retention**。
- **moving-core(运行 `ofes_daily_water_mass_256fb2a16c4e`,56 任务 0
  错误,peak 日目标 σ0 与冻结值 56/56 零差)**:水平位置逐日跟随对象峰,
  目标 σ0 冻结为事件 peak 值,参考深度沿前一日等密面交点连续选择。
  carrier 组水团项衰减斜率 −0.0025/day vs 非 carrier −0.0885/day
  (median diff +0.086,bootstrap CI [0.004, 0.190],MW p=0.122);
  wm_peak 52.5 vs 42.0 μmol kg⁻¹(p=0.094;峰值值与 fixed-site 不同是
  参考深度口径所致——moving-core 用逐日连续交点深度,fixed-site 用
  事件 core_depth_m,不矛盾)。**随载体移动时水团信号几乎不衰减**
  (斜率比 fixed-site 小一个量级以上)——这是移动参考系下 retention 的
  直接证据,但组间差异仍未达显著(与 FDR 纪律一致)。

三个口径(ΔDO 代理、fixed-site、moving-core)方向全部一致:carrier 组
衰减更慢。moving-core 的斜率差 CI 不含 0 但 MW p=0.122——如实记录,
retention 结论保持"有限";AUC 类响应受观测持续时间直接影响,不作
独立证据。另外全 650 观测日 w_along(`ofes_walong_6f4cf2337fd1`,
0 错误、56/56 全天数)已落盘,可做逐日俯冲强度时序。

稳健性回归(lifetime ~ persistent_carrier + peak_ΔDO + core_depth +
sigma0 + start_doy,n=56):carrier β=+3.1,p=0.23(ns);**peak_ΔDO
β=+0.43,p=3e-8**;core_depth β=+0.089,p=1e-3;其余 ns。事件寿命主要由
异常幅度与核心深度决定,**rotation 的寿命优势在控制后被解释掉**;
carrier 组衰减斜率名义上减半(p=0.050,CI 跨 0)——AUC 类响应受观测
持续时间直接影响,不作"衰减更慢"的独立证据(见解释纪律)。

## 四、明确裁决(lock 要求)

1. **是否存在 population-level strain→rotation 顺序?不存在。** lag 中位
   0、符号检验 p=1.0/0.65;转化双向对称(8 vs 10);Δr_share 无系统增长。
2. **rotation 是否提高异常寿命/相干性?部分、且有限。** 寿命无独立贡献
   (幅度/深度主导);直接衰减证据(斜率)仅名义边缘(p=0.050,CI 跨 0;
   真实水团项版 fixed-site p=0.127、moving-core p=0.122),AUC 类响应
   受观测持续时间影响不作独立证据;指纹保持两组无差别。
3. **McCoy-compatible 信号更像哪个阶段?** 三阶段桥接两版完成(全部 56
   事件 × start/peak/last,168 请求,peak 日与完成运行一致性检查 56/56
   零错配;先前的 19-event 条件子集运行 `14d0247c0a06` 保留为溯源,不再
   引用):

   | 版本 | start | peak | last | start–last 配对 p |
   |---|---|---|---|---|
   | fixed_peak(Eulerian sensitivity,`adcffa512438`) | 4/56(0.023) | 19/56(0.148) | 6/56(0.087) | 0.058 |
   | **per_stage(正式,`99a6ed364260`)** | **10/56(0.072)** | 19/56(0.148) | **9/56(0.106)** | **0.624** |

   fixed_peak 版三个阶段统一用 peak 日位置的 17 点足迹——"4→19→6"只
   说明信号经过固定 peak 位置时强。per_stage 版改用 start/peak/last
   各自逐日位置后,start 已有 10/56 兼容,peak 增至 19/56,last 为 9/56。
   **McCoy 型剖面形态可在最早检测阶段出现,并在异常峰值期最集中**;
   start→peak 的事后配对比例检验 p=0.042(any-event exact p=0.049),
   start→last 则无差异(p=0.624)。结果支持"部分载体早期已建立并在
   peak 增强",不支持把所有事件概括为后期才成熟或统一 carrier-first。
   域边事件(E000222 等)230 km 控制环落出交付域的
   控制点显式剔除并计数(生产路径 peak 日不受影响,一致性零错配)。
   E000171(唯一 surface-obscured SCV-compatible)peak 日 17/17 兼容 +
   17/17 原生速度确认——表层盲区 SCV 端元的典范案例。
4. **论文级机制综合:** `front/strain-dominated filamentary transport`
   与 `anticyclonic subsurface rotational carrier` 是 OFES KE 中并存的
   两种常见形态;前者突出水团形成与输送,后者突出组织与保持。strict
   SCV-like 是 always-rotation 组中最成熟的端元(4/19),事件窗内两类作用
   可以同时发生。

**H1/H2/H3 群体图景**:半数事件保持 strain 语境,而 anticyclonic
organization 是最大的连贯载体类;27 个 carrier 中 23 个在 start-day
已为 rotation,per_stage McCoy 桥接则从 start 的 10/56 增至 peak 的
19/56。最简洁的机制解释是:**通风水团经锋面/应变输送进入次表层后,
反气旋旋转在相当一部分事件中承担组织与保持;部分载体很早已经存在,
峰值期的 SCV 型热盐结构最集中。**

## 五、对论文主线的收束

**SCV 主线的 OFES 侧富集证据(已存在于完成运行,本文档升为主叙事数字):**

- 全球观测:P(DO50 | SCV)= 6.97% vs 普通 Argo 0.369%,OR = 20.22。
- OFES 同构:事件核心 McCoy 型虚拟剖面通过率 vs 同日背景 controls =
  14.9% vs 0.45%(背景 20/4480,与全球 0.369% 同量级),pooled OR =
  38.77。论文主估计量为**事件等权通过率差 0.144、bootstrap 95% CI
  [0.076, 0.221]、paired Wilcoxon p = 9.84e-5**;pooled OR 用于直观展示
  富集量级(952 个足迹与 4480 个控制具有事件内相关性)。全球
  OR 20.22 与 OFES OR 38.77 的估计口径、条件概率方向与抽样设计不同,
  不可比较"孰强孰弱",只陈述"同方向富集"。
- 载体证据:27/56 persistent anticyclonic carrier;6/56 严格
  SCV-compatible,其中 6/6 都是 carrier——**注意这是定义内包含关系
  (SCV-compatible 定义含持久旋转),不是独立证据**;独立信息是 1/6
  surface-obscured(E000171,表层盲区)与 11/56 原生速度确认(反气旋且
  rotation>strain)。表层弱/反极性在 carrier 内 4/27 persistent——
  **表层涡目录(META)盲区的直接解释**。
- 输送证据:15 条 resolved-downward 中 9 条由旋转载体承担,其轨迹积分与
  核心 w_along 在同一流场上方向收敛(G1),构成可解析下沉子集的互补证据。
- 全链条 exemplar(SCV-compatible ∩ persistent carrier ∩ resolved
  downward):E000276、E000073、E000267;E000276 为 blind top-5 第 3。

**结构层级的综合**:grid-SCV Tier-1 在事件核心闭合为 0,而 McCoy 剖面
签名与反气旋旋转载体仍显著存在。56 个事件由此勾勒出从锋面水团异常、
变形中的热盐透镜到成熟 SCV 的连续体;全球 SCV 富集抓住了这个连续体中
最稳定、最容易长期保存异常的端元,OFES 则补出了它的形成与组织过程。
S5 背景抽样进一步校准成熟透镜在全年背景中的可达性。

**代表案例**:E000002 = rotation but non-material counterexample;
E000239 = strain-dominated resolved downward;E000176 = frontogenetic/
upward contrast;全链条三例作示意图候选(E000276 优先,blind 案例天然
成立)。

## 六、待执行

- 年度 grid-SCV 背景线:S5 系统抽样敏感性(独立协议)已冻结,紧凑
  runner 已提交,基准与 73 天运行进行中(见 ofes worktree 结果报告
  第四节)。
- McCoy 桥接已完成(全 56 事件,fixed_peak 与 per_stage 两版,见裁决
  3);19-event 子集运行保留为溯源。
- 全 650 观测日 w_along 与逐日 water-mass/heave 分解均已完成并回填
  (见主分析 B;fixed-site 与 moving-core 两口径)。
- 161 事件 quality 普查已完成(见第七节)。

## 七、161 事件 quality 普查(2026-08-16 补跑完成)

- 运行:`event_diagnostics/ofes_events_4447b51c2486`
  (`diagnose_ofes_ranked_events`,candidate_count=161 = 全部
  quality-eligible),0.91 h,1290 event-days,manifest complete,
  selected_events 161 行。
- **一致性**:与旧 population 运行(59 事件)的公共事件在 8 个关键峰值列
  (target_sigma0、fixed-depth/water-mass/heave 分解、θ/AOU 等密面对比、
  Ro、normalized_strain)上**逐位零差(59/59)**——新旧 code hash 不同但
  诊断机器输出不变,transition 报告全部数字立得住。
- **诊断通过**:153/161(all_sampled_days_passed);56 严格事件 56/56 全
  通过,其余 105 事件 97/105。
- **普查口径峰值动力学**:

| 子集 | n | 反气旋(负 Ro) | rotation-dominant | 双条件 | wm fraction 中位 |
|---|---|---|---|---|---|
| 56 严格 | 56 | 52/56 | 29/56 | 27/56 | 71.7 |
| 其余 105 | 105 | 63/105 | 39/105 | 29/105 | 71.5 |
| 合计 | 161 | 115/161 | 68/161 | 56/161 | 71.6 |

- **水团异常分量在 161 普查上依旧主导**(wm fraction 中位 71.6,p25/p75
  60.1/89.2;heave 中位 11.1)——56 事件的分解结论在全部 quality 事件上
  复现,不是 deep 子集的偶然。
- 严格子集的反气旋占比(52/56)高于非严格(63/105)——deep-sensitivity
  选择口径天然富集反气旋事件,这是选择效应,不是新证据;普查口径
  (115/161)是更公允的总体估计。
- 旧 5 事件运行(ofes_events_21efbe902ab7 顶层 summary 5 行)是最初
  blind top-5 溯源,保留不动。
