# 论文整体逻辑梳理存档

- 日期:2026-08-16
- 状态:工作文档,随分析推进更新
- 定位:物理海洋学(涡旋对 DO 的物理输运,DO 为示踪剂)

## 一、核心叙事线

> 全球 Argo 样本显示次表层 DO 极端**选择性富集于 SCV**:载体率随极端强度单调增强
> (OR 5.9→20.2),而表面涡关联始终很弱(OR≈1.1);最强载体高度集中在 KE。
> 1/30° OFES 因而针对这一观测热点提供
> 过程侧解释:56 个深层事件以水团异常为主,McCoy-compatible 剖面信号相对
> 背景显著富集(19/56 vs 20/4480),同时在 SSH 表面涡场合格事件中零包含并显著回避
> 表面涡。高质量回溯显示异常粒子近期接触混合层的比例高于匹配对照,可解析
> 下沉子集的轨迹与沿等密面运动学方向收敛;反气旋次表层旋转载体则在三个
> retention 口径中都表现出更慢的水团信号衰减。形成、输送和组织在事件窗内
> 可以重叠发生,共同构成从通风水团到表层不可见次表层载体的连续过程。

一句话:**"全球识别 SCV 为最强的组织载体,OFES 揭示富氧水团如何从通风层
进入深层并被表层不可见的次表层结构组织和保持。"**

## 二、证据分层

**S 层(承重墙,强统计+稳健性)**:
- 全球样本 SCV 载体率 OR 阶梯 5.86→13.01→20.22(阈值安全+EKE 重匹配+置换+区域);
  DO50 载体 16/17 位于 KE，区域集中性作为全球统计→OFES 的桥梁
- 全球 META 对照 OR≈1.1–1.2
- OFES 对偶富集 Δ=0.144 [0.076,0.220],p=9.8×10⁻⁵ vs 背景 20/4480
- OFES 表面目录脱钩(PET 0/31 analysis-eligible strict events + 双侧负 null
  两口径 CI 全在 0 左;不否认表层 Ro 的动力表达)
- 冬季 MLD 几何 56/56(中位 +212 m,108–437 m)

**A 层(过程证据,群体统计+高置信子集+案例)**:
- 通风史回溯(见第三节三个新分析)
- 15 条高置信 resolved-downward pathways + F10–F13 案例
- 载体早期存在:27 个 persistent carriers 中 23 个在 start-day 已为 rotation;
  McCoy-compatible 信号 start/peak/last 为 10/19/9 个事件,峰值期最集中
- 滞留方向一致:ΔDO proxy、fixed-site persistence、moving-core 三口径均显示
  carrier 衰减更慢,moving-core 水团信号接近守恒(组间差异为 suggestive)
- 轨迹可信链(2D/3D 对齐、指纹保持、w_validation、Hosoda、正负对照)

**B 层(范围与分辨率边界)**:
- w_along 总体无 population-wide 对准,但 resolved-downward 子集方向收敛
- 机制线单变量检验经 BH-FDR 后不保留显著项,主统计资产集中在富集和通风线
- GLORYS 对观测 McCoy SCV 的重现不完整,突出模式分辨率边界
- Zhu 与 grid-SCV 不覆盖事件核心,说明成熟闭合透镜不是 56 事件的普遍形态;
  S5 背景抽样用于校准 detector 的背景可达性
- OFES 单年单区定位为 KE process analogue

## 三、过程线升级分析(2026-08-16 完成)

### 1. w_along resolved 子集 + 小位移事件
- 总体 null 复现:p=0.916/0.574/0.401
- resolved 子集(n=19):single +7.3 (p=0.210)、field_mean +3.9 (p=0.080)、
  **daily_mean +4.8 (13/19 向下,p=0.040 nominal)**
- 失败模式:观测方向可评估 24/56,其中 **21 向下 / 3 向上**;未评估 32 个是
  位移 <25 m(太小,非方向相反)
- **建议表述**: "In the 19 resolved events Lagrangian core w_along is
  +4.8 m/day (p=0.040 nominal), signs consistent across calibers; the
  classifiable vertical displacements are predominantly downward
  (21 of 24), while 32 events remain below the 25-m displacement threshold."
  (p=0.040 属探索性子集结果)

### 2. 通风接触率(30 天窗口统计力最高)
事件等权差(anomaly vs 对照),bootstrap CI + 配对 Wilcoxon:

| 子集 | direct | near MLD | outcrop |
|---|---|---|---|
| 全事件 vs hydro (n=28) | +0.109 p=0.005 | +0.129 p=0.001 | +0.116 p=0.002 |
| 全事件 vs kin (n=27) | +0.101 p=0.008 | +0.120 p=0.002 | +0.107 p=0.003 |
| strain 子集 (n=14) | +0.150 p=0.015 | +0.177 p=0.009 | +0.165 p=0.006 |
| rotation 子集 (n=14) | +0.068 p=0.158 | +0.082 p=0.046 | +0.067 p=0.158 |

**群体级近期通风的主估计量来自 trajectory-complete paired subset 的 30 天窗口;
60/90 天用于案例和时间范围补充。**

## 四、形成与滞留的互补角色(工作综合)

异常粒子相对匹配对照具有稳定的近期通风优势;在 30 天完整控制配对子集中,
strain 组的接触率提升更强,而滞留三口径都指向 anticyclonic carrier。
两条线合成一个形成—组织分工图景:

> **锋面/应变过程为富氧水进入深层提供形成与输送路径,反气旋旋转则为其中
> 一部分事件提供组织与保持。** 两种作用在群体中并存,事件窗内不要求固定
> 的先后顺序。

文献对应:CE 冷淡水跨锋输送(Chen2021,83% subduction patches)是应变/锋生
过程;反气旋旋转与 SCV-like 组织和三个口径的一致较慢衰减相联系。

**口径修正(2026-08-16,F20 绘制时发现)**:strain/rotation 的通风分层
(direct contact +0.150 vs +0.068,strain p=0.015)来自 30 天
group_comparison 中控制匹配完整的 **28 个事件**(strain 14 / rotation 14)。
扩展到全部 54 个可算事件后分层减弱(strain 0.121 vs rotation 0.103);未通过子集里
rotation 反而最高(0.185)。carrier 的通风史整体偏弱(0.093 vs non-carrier
0.131)。**写作时分层结论必须限定在 trajectory-complete paired subset**,F20 图上以
实心/空心同时呈现两个子集。

## 五、结构提案(8 节)

1. **Introduction** — SCV 生态学(McCoy 2020)+ 涡-生地化 + 次表层 DO 极端;
   问题:谁承载?如何形成?落点:表面涡目录与次表层结构的观测空隙
2. **Data & Methods** — 全球 Argo+ΔDO、McCoy 目录+重现能力审计、OFES NP30+
   56 事件、虚拟 Argo、轨迹验证链、预注册纪律(BH-FDR、lock、双侧口径);
   **prediction table**(机制→预测→检验→结果,把三个 null 写成预测命中);
   明示全球 4084→4066→263→245/244 与 OFES 161→56→52→34→28 漏斗
3. **全球统计** — OR 阶梯 vs META≈1.1;KE/Pacific-outside-KE/Atlantic/
   Indian/SO 区域载率结构;EKE 重匹配、置换与区域稳健性;
   V/I 分类轴(观测 heave+冬季 MLD)作第 5 节的观测前奏;模式重现能力讨论
4. **OFES 事件群体** — 对偶富集(P(SCV|DO) 呼应全球 P(DO|SCV));表面目录脱钩
   (PET 0/31 eligible+双侧负 null);剖面签名与成熟闭合透镜的层级差异;水团分解
5. **形成:通风-脱钩-潜沉** — 三层递进:几何(F5)→ 群体通风史(分析 3)
   → 案例(F10–F13);机制落点 CE 冷淡水跨锋(σ0≈26.4,Chen2021/Nagano2016)
6. **组织与滞留** — start-day rotation、peak-concentrated McCoy 信号与
   滞留三口径(F7);与第 5 节合成形成-组织连续过程
7. **尺度连接与过程分辨率** — w_along 的群体背景与 resolved 子集、
   GLORYS/OFES 的 SCV 重现能力、grid-SCV 与 S5 背景;说明轨迹为何能从
   网格平均场中提取出集中在可解析子集的下沉通道
8. **Discussion** — Keutgen2026 差异化(他们 V/I=未来 reemergence/碳、EKE
   热点;我们=过去通风史/分类、SCV 结构);分类学贡献;OMZ 外推

## 六、待决策项

1. META/PET 双对照线:分散在 3、4 节作对照(推荐),不单列
2. V/I 分类轴:第 3 节末、以"承载者的通风史"引出(推荐),作 3→5 伏笔
3. FDR 0 存活:第 7 节正文一句话 + 附录表
4. "56/56 在冬季 MLD 下"的同义反复风险:靠通风史回溯(分析 3)扛,几何事实
   只作背景;措辞"ventilated before onset and within 30 days for a
   significant excess of particles"
5. 全球 244 SCV 与 OFES 56 事件非同一样本:桥接靠对偶富集,不声称再现样本

## 七、解释范围

- 15 条高置信 resolved-downward 路径代表可解析子集,其余事件保留为异质路径
- w_along 子集 p=0.040 标为 exploratory;主结果是轨迹与运动学方向收敛
- rotation 在多数 carrier 的最早检测日已存在,但不据此声称严格因果先后
- 30 天通风史作为群体主估计量,60/90 天承载案例与更长时间范围

## 八、正文 evidence table（审稿口径）

| 主张 | 主估计量 | 统计读数 | 分析地位 | 层级 |
|---|---|---|---|---|
| SCV 对强 DO 异常的关联远强于普通 Argo | DO50 SCV vs All-Argo OR | 20.22,p=1.4×10⁻¹⁶ | 同分母阈值扫描;EKE/空间/置换另作稳健性 | S |
| META 表层涡关联弱 | DO20/35/50 META OR | 1.13/1.22/1.14 | 与 SCV 同阈值对照 | S |
| OFES 事件核心富集 McCoy-compatible 剖面 | 事件等权 core−background | +0.144 [0.076,0.220],双侧 p=9.8×10⁻⁵ | 锁定事件/背景口径 | S |
| OFES 事件与闭合 SSH 涡目录脱钩 | PET local/annual null | 双侧 p=8.9×10⁻⁵/1.7×10⁻⁶ | 冻结检测器;双侧方向由审计补报 | S(观测事实) |
| 异常粒子近期通风史增强 | 30 天 direct MLD contact | +0.109,p=0.0049,n=28 | trajectory-complete 配对子集 | A |
| 可解析下沉子集运动学收敛 | Lagrangian w_along | +4.8 m d⁻¹,p=0.040,n=19 | exploratory subset | A |
| carrier 的水团信号衰减更慢 | 三个 reference frames | p=0.050/0.127/0.122,方向一致 | exploratory;不声称显著延寿 | A |
| 成熟闭合网格 SCV 并非普遍可解析 | Zhu/grid-SCV containment | 0 event core containment | 分辨率/方法边界 | B |
