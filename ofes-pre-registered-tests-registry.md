# 预注册检验注册表(机制收口三 lock)+ BH-FDR 校正

- 日期:2026-08-16
- 分支:`feat/ofes-mechanism-clinch`(worktree `Oxygen-mechanism`)
- 目的:在论文写作前锁定"哪些 p 值可以叫显著"。这是纪律性文档,
  不是新计算——所有数字来自既有完成运行与三份 lock。
- 规则:主检验集 = 三 lock 预注册的科学主张检验(claims);解析验证
  (G4 6.6e-12)、稳健性回归的控制项、描述性判定(56/56 在 MLD 之下,
  无 p 值)不纳入 FDR 分母。

## 一、主检验集(BH-FDR 结果)

| # | 检验 | 来源 lock | p | 单测名义 | FDR 后 |
|---|---|---|---|---|---|
| 1 | w_along core−ring 单日 | G1 lock | 0.916 | ns | ns |
| 2 | w_along core−ring 欧拉场 | G1 lock | 0.574 | ns | ns |
| 3 | w_along core−ring 逐日平均 | G1 lock | 0.401 | ns | ns |
| 4 | lag strain→r_share 符号检验 | transition lock | 1.00 | ns | ns |
| 5 | lag strain→DO 符号检验 | transition lock | 0.648 | ns | ns |
| 6 | retention lifetime(carrier vs 其余) | transition lock | 0.74 | ns | ns |
| 7 | retention 衰减斜率 | transition lock | 0.050 | 边缘 | ns |
| 8 | retention post-peak AUC(primary) | transition lock | 0.117 | ns | ns |
| 9 | McCoy start–last 配对 | transition lock | 0.058 | 边缘 | ns |

(2026-08-16 修正:w_along 三口径 p 值取自最终 3 格平滑运行
`ofes_walong_16e0d8158d21` 的正式 tests 块——单日 0.916、欧拉场
0.574、逐日平均 0.401;初版表误抄了旧 5 格敏感性数字,已更正。)

**BH-FDR(q=0.05 与 q=0.10):0/9 存活。** 本表登记的三份机制收口 lock
中的九个检验在多重比较校正下无一显著；这不包括另行预注册且显著的
McCoy 富集与 matched-control 通风史检验。这与既有措辞纪律一致
(transition 报告已用"有限/
无两阶段时序"),但把约束显式化:**论文中机制线不得出现"显著"字样,
一律用描述性/一致性措辞**(如 p=0.050 的衰减斜率只能写"名义上边缘,
未通过多重比较校正";AUC 受观测持续时间直接影响,不得写成"衰减速度
显著更慢"——直接衰减证据是 slope,代理版仅名义边缘、真实水团项版
p=0.127)。

## 二、探索性检验(不进 FDR 分母,论文须标注 exploratory)

| 检验 | p | 备注 |
|---|---|---|
| resolved-down 15 事件核心 w_along 为正 | 0.032/0.018 | 子集事后选择,探索性 |
| peak-rotation post-peak AUC(secondary) | 0.015 | transition 报告已标 Secondary |
| 真实水团项 retention 衰减 | 0.127 | 事后扩展运行(用户要求"昂贵的事都跑") |
| 水团项 post-peak AUC | 0.18 | 同上 |
| 稳健性回归 carrier β | 0.23 | 控制项,非主张 |

## 三、论文的显著性资产(任何 FDR 下存活,主叙事)

| 检验 | p | 来源 |
|---|---|---|
| OFES 事件核心 McCoy 通过率 vs 同日背景,OR=38.77 | Fisher 1.55e-87 | grid-SCV 线完成运行 |
| 事件等权通过率差 0.144 | paired 单侧 9.84e-5；双侧 1.97e-4 | 同上 |
| 全球 P(DO50\|SCV) OR=20.22 | 极强 | 检测线既有结论 |

**含义:论文的统计显著主张集中在富集线;机制线的角色是"富集同构的
机制解释",全部以描述性/一致性措辞呈现。** 这与"SCV 主线 + 机制补充"
的叙事结构一致,写作时无需再为机制线寻找显著性。

## 四、写作期说明(不需要新计算)

- 模式代表性限定:OFES 冬季 MLD 中位 358 m vs 观测(HT ~278 m)——写入
  limitation;单年 2003、自由运行无同化,机制陈述加"在 OFES 配置内"。
- 循环论证声明一句:carrier 身份门只用动力学量(rotation_day_fraction),
  响应量只用示踪量(ΔDO/水团项),两者构造上独立;DO 不参与任何身份门。
- grid-SCV S5 已完成:背景 Tier-1 occupancy 0.227%,事件匹配后的 56 核心
  预期命中数 0.237、P(0)≈0.79。该项是背景可达性/功效校准,不进入上述
  九个机制检验的 FDR 分母。
