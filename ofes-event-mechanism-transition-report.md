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

- **rotation 在 onset 已存在的事件 = 29/56**(19 始终 + 10 转出);strain
  起始 = 27/56(19 始终 + 8 转入)。两类起始形态几乎各半。
- strain→rotation(8)只有 2/8 成为 persistent carrier——**后期转入旋转
  并不带来滞留身份**;persistent carrier 几乎全部(19/19)在 onset 就是
  rotation。SCV-compatible(4/6)与 McCoy-any(8/19)集中于 always-rotation。
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

## 三、主分析 B:retention 对比(事件级,event-level bootstrap)

| 响应量 | carrier(27)vs 其余(29)中位 | diff | bootstrap 95% CI | MW p |
|---|---|---|---|---|
| lifetime(days) | 11 vs 6 | +5 | [−3, +7] | 0.74 |
| post-peak 衰减斜率(归一/day) | −0.0139 vs −0.0261 | +0.012 | [−0.003, +0.051] | **0.050** |
| post-peak AUC | 8.75 vs 4.19 | +4.57 | [−0.94, +6.42] | 0.117 |
| DO 指纹保持日占比 | 1.0 vs 1.0 | 0 | [0, 0] | 0.31 |

Secondary(peak rotation 29 vs 27):lifetime +5(p=0.75);衰减斜率 +0.011
(p=0.125);**post-peak AUC +5.5,CI [1.05, +7.38],p=0.015**。
time_to_half 53/56 右删失(检测期内未降到 peak 50%),不参与比较。

**真实水团项版(650 观测日逐日 water-mass/heave 分解,运行
`ofes_daily_water_mass_afa024fe7a6e`,peak 日与 population 逐位一致
56/56 零差)**:carrier 组归一化水团项衰减斜率 −0.118/day vs 非 carrier
−0.254/day(慢约 2.2 倍,方向与 ΔDO 代理一致,median diff +0.136,
bootstrap CI [−0.001, +0.377],MW p=0.127);post-peak AUC 0.94 vs 0.72
(p=0.18);wm_peak 75.6 vs 65.6 μmol kg⁻¹。真实分量上的滞留优势方向一致
但弱于代理版、未达显著——如实记录,retention 结论保持"有限"。
另外全 650 观测日 w_along(`ofes_walong_6f4cf2337fd1`,0 错误、56/56
全天数)已落盘,可做逐日俯冲强度时序。

稳健性回归(lifetime ~ persistent_carrier + peak_ΔDO + core_depth +
sigma0 + start_doy,n=56):carrier β=+3.1,p=0.23(ns);**peak_ΔDO
β=+0.43,p=3e-8**;core_depth β=+0.089,p=1e-3;其余 ns。事件寿命主要由
异常幅度与核心深度决定,**rotation 的寿命优势在控制后被解释掉**;但
peak-rotation 组的 post-peak AUC 仍显著更大(衰减更慢),carrier 组衰减
斜率减半(p=0.050,CI 跨 0)。

## 四、明确裁决(lock 要求)

1. **是否存在 population-level strain→rotation 顺序?不存在。** lag 中位
   0、符号检验 p=1.0/0.65;转化双向对称(8 vs 10);Δr_share 无系统增长。
2. **rotation 是否提高异常寿命/相干性?部分、且有限。** 寿命无独立贡献
   (幅度/深度主导);衰减速度确实更慢(AUC p=0.015 显著、斜率 p=0.05
   边缘);指纹保持两组无差别。
3. **McCoy-compatible 信号更像哪个阶段?** 三阶段桥接已完成(全部 56
   事件 × start/peak/last,168 请求,运行 `ofes_mccoy_stage_bridge_adcffa512438`,
   peak 日与完成运行一致性检查 56/56 零错配;先前的 19-event 条件子集运行
   `14d0247c0a06` 保留为溯源,不再引用):start 4/56、peak 19/56(构造)、
   last 6/56;mean 兼容占比 0.023 / 0.148 / 0.087;start–last 配对
   Wilcoxon p=0.058。**信号是 peak-成熟期特征,不是 onset 特征**;持续到
   late 的 6 个事件中 5 个是 rotation_day_fraction=1.0 的全旋转 persistent
   carrier(另 1 个是 E000176,rotation 占比 0.5 的锋生对照)——**SCV 型
   形态的滞留与反气旋载体重合**。E000171(唯一 surface-obscured
   SCV-compatible)peak 日 17/17 兼容 + 17/17 原生速度确认——表层盲区
   SCV 端元的典范案例。
4. **当前最高可用措辞:** `front/strain-dominated filamentary transport`
   与 `anticyclonic subsurface rotational carrier` 是 OFES KE 中并存的
   两种常见形态(起始各 ~50%);strict SCV-like 端元是 always-rotation
   组的少数子集(4/19);无两阶段时序。

**H1/H2/H3 裁决(全部为报告,不设门槛)**:H1(锋面先、反气旋后)不被支持
(无时序);H2(纯锋面)部分支持(半数事件始终 strain),但 rotation 的衰减
优势反对"反气旋完全无关";H3(载体先行)形态上常见(29/56 onset 即
rotation),但无"后装载水团"的时序证据。→ **多机制异质性,不强行归一**
(lock 冻结规则)。

## 五、对论文主线的收束

**SCV 主线的 OFES 侧富集证据(已存在于完成运行,本文档升为主叙事数字):**

- 全球观测:P(DO50 | SCV)= 6.97% vs 普通 Argo 0.369%,OR = 20.22。
- OFES 同构:事件核心 McCoy 型虚拟剖面通过率 vs 同日背景 controls =
  14.9% vs 0.45%(背景 20/4480,与全球 0.369% 同量级),pooled OR =
  **38.77**(Fisher p = 1.55e-87),事件等权通过率差 0.144,bootstrap95
  [0.076, 0.221],paired Wilcoxon **p = 9.84e-5**。全球富集在模式里被
  同方向、更强地复现。
- 载体证据:27/56 persistent anticyclonic carrier;6/56 严格
  SCV-compatible(**6/6 都是 carrier**,1/6 surface-obscured);11/56
  原生速度确认(反气旋且 rotation>strain);表层弱/反极性在 carrier 内
  4/27 persistent——**表层涡目录(META)盲区的直接解释**。
- 输送证据:15 条 resolved-downward 中 9 条由旋转载体承担,且其核心
  运动学 w_along 显著为正(G1)。
- 全链条 exemplar(SCV-compatible ∩ persistent carrier ∩ resolved
  downward):E000276、E000073、E000267;E000276 为 blind top-5 第 3。

**限制的重新定位**:grid-SCV Tier-1 闭合 0 不否定 SCV 重要性,而是
"1/30° 网格分辨不了闭合透镜动力学"——观测层富集最强与模式层闭合最难是
同一枚硬币(SCV 是形成连续体的稳定端元;McCoy 定义本身含锋面潜沉、
重层化隔离与亚中尺度卷起)。本审计的"无两阶段时序"只否定"总体逐日
时序规律",降为"形成路径异质"的次级陈述,不动摇富集主线。

**代表案例**:E000002 = rotation but non-material counterexample;
E000239 = strain-dominated resolved downward;E000176 = frontogenetic/
upward contrast;全链条三例作示意图候选(E000276 优先,blind 案例天然
成立)。

## 六、待执行

- 年度 grid-SCV 背景线:工程试跑完成后按成本外推决策(见 ofes worktree
  结果报告第四节)。
- McCoy 桥接已完成(全 56 事件,见裁决 3);19-event 子集运行保留为溯源。
- 全 650 观测日 w_along 与逐日 water-mass/heave 分解均已完成并回填
  (见主分析 B)。
- 可选大件:161 个 quality-eligible 事件中非 deep 的 ~100 个补跑 population
  诊断(入口 `diagnose_ofes_ranked_events`),事件选择口径待定。
