# OFES 表面涡(PET)检测、追踪与 56-event 关联结果报告

- 日期:2026-08-24
- 分支:`feat/ofes-surface-eddies`(worktree `Oxygen-ofes-surface`)
- 依据:lock `ofes-surface-eddy-analysis-lock.md`(检测)与
  `ofes-surface-eddy-association-lock.md`(关联)。本文档是正式数字的唯一
  回填位置,lock 不回填。
- 正式运行:`ofes-np30-2003-v4`(检测 manifest complete);原关联后处理 commit
  `59ef7b1`,本次可分析分母审计修正记录于 schema-v2 association manifest。
- 输入全部 SHA-256 记录于 `surface_eddy_summary.json`;检测环境
  `ofes-pet`(Python 3.10.20,pyeddy_tracker 3.6.1)。

## 一、检测与追踪(全年 365 日)

- 365/365 日 manifest complete;逐日对象 **614 个**(反气旋 294 /
  气旋 320),中位振幅 0.090 m(p05–p95:0.008–0.313 m),
  中位有效半径 40.6 km(p05–p95:10.7–64.2 km)。
- 轨迹 **64 条**:persistent(≥30 天)5 条、long(10–29 天)12 条、
  short(2–9 天)27 条、untracked 20 条;虚拟观测占比 4.5%。
- filter-valid 域 = 全网格的 18.2%(700-km Bessel 高通核足迹约束,
  边界裁剪核语义按 lock),检测对象零 boundary-censored。

## 二、与 56 个严格 DO50 事件的关联

- 严格分母 56;峰核 filter-valid 且 ocean-valid 的分析合格数 **31**。
- **分析合格峰核被 PET 涡包含:0/31**(有效轮廓与速度轮廓双口径、
  任何极性);其余 25 个峰核在 PET 滤波有效域外,记为不可判定而非阴性。
  可计算的 peak-footprint 与 PET 轮廓重叠中位 **0.0**(n=32;该足迹诊断的
  可用性门与峰核 analysis-eligible 门不同)。
- 161 个 quality-eligible 敏感性人群中 58 个峰核在 PET 可分析域内:
  包含数 **0/58**(双口径),其余 103 个不可判定。
- 31 个分析合格严格事件到最近 PET 涡中心的距离中位为
  **5.87×有效半径**;全体严格事件的 9.22×不作为主口径,因为混入了
  滤波域外峰核。
- 29 个峰日旋转主导事件中 16 个 PET 可分析:包含数 **0/16**;
  另 13 个不可判定。

结论:在 PET 可分析域内,这些深层 DO50 峰核**未被识别为闭合 SSH 表面涡**。
PET 全年检出 614 对象、64 轨迹,说明检测器并非空转;该结果限定于冻结的
PET/META-family 闭合 SSH 目录,不等同于“无任何表层动力表达”。

## 三、双侧 null 检验(2026-08-14 审计修正后的正式口径)

事件级 core−null(事件等权均值、10,000 次种子 bootstrap 95% CI、
配对 Wilcoxon):

| null | mean | bootstrap 95% CI | 双侧 Wilcoxon p | 个体负/事件 |
|---|---|---|---|---|
| 同天 120–240 km 环带(面积加权) | −0.0190 | [−0.0291, −0.0105] | 8.86×10⁻⁵ | 20/31 |
| 年尺度 月×1° 纬度分层 | −0.0099 | [−0.0140, −0.0066] | 1.69×10⁻⁶ | 30/31 |

- **两个口径的 95% CI 都完全落在 0 左侧**:31 个 PET 可分析事件核心的
  闭合 SSH 涡占有率显著低于背景;这是目录占有率的空间负关联,不直接判定
  其物理成因。
- **口径修正说明**:原 summary 只报告单侧 greater p(annual 0.999999、
  ring 0.999956),表面读作"无正关联"。2026-08-14 独立审计指出该单侧
  口径掩盖了对侧;本文档以双侧为正式口径,less 单侧 p 分别
  ≈8.5×10⁻⁷ / 4.4×10⁻⁵,双侧 ≈1.7×10⁻⁶ / 8.9×10⁻⁵。
- 个体层面 20/31(环带)、30/31(年度)事件为负,其余分别为 11/31、
  1/31 零差,无一可分析事件为正偏(max = 0.0)。

## 四、旋转与 McCoy 交叉表

- rotation_29 与表面 Ro:深极性反气旋 27 / 气旋 2;
  **表层核心加权 Ro 同极性 26/29**。深层旋转信号在表层 Ro 上有影子,
  但 16 个 PET 可分析峰核均未形成闭合 SSH 等高线涡(0/16),另 13 个
  因滤波有效域限制不可判定。
- McCoy 事件级字段与 lock 冻结计数逐值一致:
  center 9/56、center+velocity 6/56、any 19/56、any+velocity 11/56。
  PET 包含数与 McCoy 兼容数之间的交叉无变异可用(包含数全 0),
  交叉表如实保留零行。

## 五、解读与叙事位置

- **"闭合 SSH 目录脱钩"的模式侧证据**:McCoy 2020 SCV 的定义特征之一
  是表层表达弱或被遮挡。本线显示 PET 可分析的 OFES 深层 DO50 峰核未被
  闭合 SSH 涡捕获(0/31),且距最近目录涡中位 5.87R;它支持“表层目录可能
  漏掉关键次表层结构”,但不单独证明卫星看不见的具体物理原因。
- 与 mechanism 线共同构成次表层组织证据:事件侧 McCoy-compatible
  信号相对背景富集(19/56),全球观测中 SCV 对 DO50 的关联最强
  (OR 20.2),而本线显示同类 OFES 深层事件在 PET 可分析域内未进入闭合
  SSH 涡目录。
  三条证据合起来支持"表层目录难以捕获关键次表层载体"这一主叙事。
- **适用范围**:本线是 META 家族的 SSH 涡实现,直接回答表面涡目录能否
  捕获这些深层事件;形成与输送的因果解释由通风史、三维轨迹和生命周期线
  提供。两部分组合后,表层盲区与次表层过程各自有独立证据支撑。

## 六、图件

- 历史运行图:`plots/event_containment_and_contour_comparison.png`、
  `rotation_pet_expression.png`、三案例图与四工程日轮廓图。
- 论文候选图由受版本控制的 `plot_ofes_surface_eddy_paper_figures.py` 从上述
  审计输出生成;权威输出位于本 worktree 的
  `plot_outputs/do/ofes_np30_ke/surface_eddy/paper_figures/`,并同步到整合图包:
  - `F14_surface_containment.png`:可分析分母包含率 + 最近涡距离分布
  - `F15_surface_null_two_sided.png`:双侧 null 森林图(两口径)
  - `F16_surface_rotation_pet.png`:rotation×PET 三态交叉(0/16 eligible、
    13 unassessable、surface Ro 同极性 26/29)
