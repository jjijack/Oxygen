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
3. **McCoy-compatible 信号更像哪个阶段?** 三阶段桥接未跑(lock 允许的
   条件子集,待执行);现有线索:SCV-compatible 与 McCoy-any 集中于
   always-rotation 组,与载体形态共现,不构成阶段证据。
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

- 全球观测:成熟次表层相干反气旋(SCV)是深层高氧异常最强、最稳健的单一
  富集结构(P(DO50|SCV)=6.97% vs 背景 0.369%,OR=20.22)。
- OFES:异常水团由通风水团经锋面/应变与可分辨垂向输送进入温跃层,承载
  形态在锋面 filament 与反气旋旋转载体之间大致各半并存;SCV-like 端元
  是旋转组内少数子集;无 population-level 形成顺序。
- 网格 detector 负结果:解释"为什么不能把所有异常直接等同于成熟 SCV"
  (形成阶段、形变过程与模式分辨率共同限制)。
- 代表案例:E000002 = rotation but non-material counterexample;
  E000239 = strain-dominated resolved downward;E000176 = frontogenetic/
  upward contrast。"strain→rotation 且 retention 增强"的示意图候选可
  从 strain→rotation(8)∩ carrier(2)中选,只作示意图不作证明。

## 六、待执行

- McCoy 三阶段桥接(start/peak/last 三日 footprint,19-event 条件子集
  优先),不阻塞上述结论。
- 年度 grid-SCV 背景线:工程试跑完成后按成本外推决策(见 ofes worktree
  结果报告第四节)。
