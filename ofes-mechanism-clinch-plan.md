# 溶解氧异常机制证据计划(Mechanism Evidence Plan)

- 日期:2026-08-16
- 分支:`feat/ofes-mechanism-clinch`(worktree `Oxygen-mechanism`)
- 状态:历史执行路线图；A--D 与 S5 均已完成。写作数字以 `paper-outline.md`
  和各 `*-results-report.md` 为准，不以本计划中的运行中描述为准。
- 定位:回答"导致 ΔDO 异常的主要机制是什么"的证据路线图。与 `feat/ofes-ke-mechanism`
  的 grid-SCV v2 线共享 OFES 数据与事件表,但代码/运行隔离(年度跑期间那条线冻结
  track.py)。

## 一、现状:证据分级盘点(截至 2026-08-16)

**高可信度证据(可直接写入论文)**

1. 异常是真实的:δDO 经 WOA 背景验证(P29812 +108 ≈ 110.6;多浮标证 ~2 周相干斑块)。
2. 水团项主导:OFES 56 事件 54/56 water-mass dominated,水团绝对分数中位数 0.862。
3. 自上补充为主:峰 vs OMZ 核位置判别,90% 自上,N.Pac 98.7%;深层补充只是
   I-BG 大西洋/南大洋深层小子群;AOU 是判决工具。
4. 通风接触机会富集:30 天轨迹 MLD 接触/近接触/outcrop 相对两组 matched controls
   增强；完整配对子集为 hydro n=28、kinematic n=27(p=0.001–0.008),
   56/56 指完成运行覆盖,不是统计检验分母；60/90 天受左截断限制。
5. OFES 18 条 resolved pathways(15 下沉、3 上浮,51/52 数值门);E000002 是
   关键反例(patch gate 与 resolved-3D gate 均失败)。
6. GLORYS 1/12° 沿等密面俯冲 null:三重时间口径都区分不了 V/I(p=0.65–0.87),
   机制 = 流向相对等密面随机(cos≈0)。
7. McCoy 观测层富集(DO20/35/50 载体 67/39/17 of 244,OR 5.86–20.22)但动力学层
   未闭合(grid-SCV v2:生产 Tier-0 转移 89/141、事件核心 Tier-1 闭合 0、
   well-resolved 0;重放口径 141/141 已由裁决关闭)。

**缺口(证据目标)**

- G1 **1/30° 沿等密面对准**:已完成。总体 core−ring 为 null,可解析下沉子集
  在两个互补运动学口径上收敛。
- G2 **背景统计**:已用独立 S5 协议完成 73 个系统日期抽样。Tier-1 背景
  occupancy 0.227%;事件匹配后的 56 核心预期命中数 0.237,P(0)≈0.79。
- G3 **V/I 轴 OFES 类比**:观测侧本地冬季 MLD 阈值已落 code;OFES 侧可做
  heave-可达性(事件核深 vs 本地冬季 MLD)类比。
- G4 **E000002 来源**:沿等密面输送/形变/混合/core switching 四选,不能再拟合
  刚性下沉 parcel。

**明确不可能/已排除(不浪费算力)**

- tracer tendency 闭合:交付只有 9 变量、无平流诊断输出,做不了。
- PV 作年龄/通风轴:已证非年龄轴(纬度-季节记忆),排除。
- w_cross 散度统计:GLORYS 侧噪声太大;OFES 侧最多描述性报告。

## 二、主路径假说(分级证据写法,不是定论)

> ΔDO 异常以水团信号为主,由"自上通风 + 沿流/沿等密面输送"到达,承载结构异质
> (反气旋旋转俘获 27/56 persistent carrier,strain/deformation 27/56);观测层的
> McCoy 型 SCV 富集强,但 1/30° 动力学层闭合失败,支持"亚格子 SCV 型载体是重要
> 观测类而非可分辨动力学实体";区域机制分家(KE 两机制、SO SAMW/AAIW、NA SPMW、
> 印度洋锋面),不做全球单一机制。

## 三、执行顺序(昂贵步骤先行,都"不管结论如何都需要")

| 步 | 内容 | 状态 | 成本估计 |
|---|---|---|---|
| A | w_along 审计(G1)| **完成** | 0.13 h;结果见 ofes-walong-results-report.md |
| B | grid-SCV 背景统计(G2)| **S5 完成** | 73 天 × 63 tile,4599/4599 |
| C | OFES 冬季 MLD + heave 可达性(G3)| **完成** | 0.16 h;结果见 ofes-winter-mld-results-report.md |
| D | Formation–organization–retention 时间序审计 | **完成** | 结果见 ofes-event-mechanism-transition-report.md |
| E | E000002 沿等密面 DO 通量(G4)| 降级为次要 | 总体时序优先(另一个 AI 裁决) |
| F | 论文机制节收束 | **可启动** | 当前工作 |

- A 核心结论:总体 core−ring null(三重口径 p=0.40–0.86,正占比≈硬币,与
  GLORYS 1/12° null 同形态);**resolved-downward 15 事件子集核心 w_along
  显著为正(欧拉/拉格朗日 p=0.032/0.018)**——同一 OFES 流场的两个互补诊断
  (轨迹位移 vs 运动学对准)在可分辨俯冲子集上收敛;同源,非独立证据。
  俯冲是子集机制,不是总体机制。(2026-08-16 平滑口径修正后重跑,数字以
  结果报告为准。)
- A/B 并行(A 在机制 worktree,B 在 ofes worktree,不冲突)。
- C 复用 A 的 T/S 读取模式(每 event-day 最近单元剖面)。
- D 完全复用 A 落盘的映射场(逐 event-day nc),零新增 I/O。

## 五、两阶段框架(2026-08-16 采纳的科学方向)

> 原工作假说是锋面/应变过程先**形成并向下输送**水团异常,随后反气旋旋转/
> SCV-like 结构可能参与**组织与保存**——把"亚中尺度过程 vs SCV"从二选一
> 改成可检验的两阶段框架；第二阶段是否成立由衰减与对照诊断决定，不在此预设。
> 网格 SCV 检测失败不推翻全球 SCV 主线,反而提示 OFES 更可能解析了形成/变形
> 阶段,没有完整解析最终闭合、长寿命的 SCV 透镜。

完成后的裁决是:群体层面没有稳定的 strain→rotation 日级时滞,形成、输送
与旋转组织可重叠发生；因此正文使用“互补作用/连续过程”,不再把固定两阶段
先后顺序当作已支持结论。

现有证据映射:

- 阶段 1(形成—输送):18 resolved pathways(15 下)+ G1 运动学收敛 +
  onset 锋生(E000176 锋生型正例、E000239 strain-非锋生反例)+ 水团项在
  deep-entry 增强(E000002)。
- 阶段 2(组织—滞留):27/56 persistent anticyclonic carrier、6/56
  SCV-compatible、McCoy 观测层富集(OR 5.86–20.22);grid-SCV v2 剖面形态
  大量存在(Tier-0 5623 种子、89/141 转移)但闭合透镜 0——"形态学信号在、
  动力学闭合不在"。

T1–T3 已并入正式 transition 审计(lock
`ofes-event-mechanism-transition-analysis-lock.md`,已完成,见其结果报告):

- T1 阶段时序:15 条 resolved-down 事件中,rotation onset 相对 deep-entry
  的先后(lifecycle 逐日 rotation_dominated/depth_mean,650 观测日)。
  lag > 0 = 两阶段时序支持。
- T2 滞留:persistent carrier vs 非 carrier 的事件寿命/DO 指纹持续差异。
- T3 锋生-下沉关联:三例逐日锋生率与下沉速率的个例级并列表(56 事件无
  逐日锋生量,只做三例)。


## 四、与论文结构对应

- 观测统计(Argo/META/McCoy 富集)→ 已定稿部分。
- OFES process analogue:水团分解 → 18 pathways → 通风历史 → **A(w_along)** →
  **B(年度背景)** → G3/G4。
- 机制节最后按"主路径假说 + 已排除项 + 剩余缺口"三层写,与用户此前认可的
  结构一致。
