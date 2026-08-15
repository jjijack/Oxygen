# OFES 冬季 MLD 审计(Winter MLD Audit)Lock

- 日期:2026-08-16
- 分支:`feat/ofes-mechanism-clinch`(worktree `Oxygen-mechanism`)
- 状态:探索性表征;本文件是预注册说明,不设事后通过门槛。

## 目的

观测侧 V/I 轴用"逐点本地冬季 MLD 阈值"(heave + 冬季混合层可达)分类。OFES 侧
需要定量同构的量:56 个事件核心深度相对其核心处本地冬季混合层的距离。这是
"自上补充"论点的 OFES-side 落点——若所有核心都在本地冬季 MLD 之下数百米,
则异常水团必然经潜沉/沿等密面输运从通风层到达(观测侧 90% 自上的定量同构)。

## Estimand

- 每事件:core_depth_m − winter_max_MLD_core(核心处 2003-01-01..03-31 逐日
  MLD 的冬季最大值),以及 winter_median_MLD、有效日数、封顶日数。
- 与 30 天轨迹通风可达性(event_horizon_summary 的
  median_minimum_depth_minus_mld_m,anomaly 组,30 天)交叉:两套独立可达性
  (快照冬季混合 vs 轨迹过程)的对照表。

## 算法

- 每 (event, winter day):核心最近网格单元取 T/S 剖面(加载窗口 ±0.05°,
  最近单元按大圆距离)→ TEOS-10 σ0 → `_mld_from_sigma`(参考深度 10 m,
  密度阈值 0.03 kg m⁻³,与通风历史同口径)→ 单日 MLD;阈值在剖面内从未达到
  的日记封顶(capped)。
- 冬季窗口固定 2003-01-01..03-31:交付只有 2003 一个完整年,所有事件共用
  同一冬季是快照启发式的诚实选择;2002 冬不可用。

## 门(轻)

- G1 输入 56 行、hash 记录。
- G2 每事件冬季有效剖面天数(不封顶)≥ 45/90 的事件数报告。
- G3 MLD 量级合理性:冬季最大 MLD 分布报告(KE 文献量级 ~100–400 m)。

## 报告项(不设门槛)

- core_depth_m − winter_max_MLD 的分布(全部应为正值且数百米量级;若出现
  负值 = 核心在冬季混合层内,同样如实报告)。
- 与 rotation/strain、resolved down/up、water-mass/heave 的交叉表。
- 快照 winter_max_MLD 与 30 天轨迹 median_minimum_depth_minus_mld_m 的
  一致性表(两口径测不同对象:冬季混合极值 vs 30 天过程接触;不一致不互否)。

## 措辞纪律

- 快照 MLD 启发式只给"核心相对通风层的深度尺",不构成通风证明;通风证据
  保留在 56-event 通风历史(MLD/outcrop matched-control)的原门控结论内。
