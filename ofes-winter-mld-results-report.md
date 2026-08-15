# OFES 冬季 MLD 审计结果报告(机制证据 G3)

- 日期:2026-08-16
- 分支:`feat/ofes-mechanism-clinch`(worktree `Oxygen-mechanism`)
- 依据:lock `ofes-winter-mld-audit-lock.md`(冻结于运行前)。本文档是正式
  数字的唯一回填位置,lock 不回填。
- 完成运行:`/mnt/w2/scratch/user3/Oxygen-cache/ofes_winter_mld_results/
  ofes_winter_mld_b63817508eb8/`(manifest schema 1,status complete)。

## 一、门与运行

| 门 | 结果 |
|---|---|
| G1 输入(56 行 + hash) | PASS |
| G2 覆盖(5040 event-days,0 错误,56/56 有效天数 ≥45) | PASS |
| G3 量级(冬季最大 MLD 分布,见下) | 报告 |

耗时 0.16 h(8 workers);MLD 口径与 56-event 通风历史一致(参考 10 m,
密度阈值 0.03 kg m⁻³),最近单元剖面(格距 0 m,无插值)。

## 二、主结果:核心相对本地通风层的深度尺

| 量 | 值 |
|---|---|
| 核心处冬季最大 MLD | 中位 358 m(p10 247,p90 446) |
| core_depth_m − 冬季最大 MLD | 中位 +212 m,范围 +108 ~ +437 m,**0 个为负** |
| 30 天轨迹可达性(median min depth−MLD ≤ 0) | 5/54 事件 |

**全部 56 个事件核心都在本地冬季混合层之下 100–450 m**:快照口径只确立深度尺前提:**当地冬季混合不能在最终核心深度直接通风**;
结合 matched-control ventilation histories(30 天 MLD/近 MLD/outcrop 接触
相对两组对照显著富集,p=0.001–0.005),支持**此前或远端通风后再经向下
和/或侧向输送**。快照 MLD 单独不构成俯冲证明。

OFES 冬季 MLD 中位 358 m 比 GLORYS(~240 m)/HT(~278 m)深,符合自由运行
OFES 的 KE 深混合特征;本报告只用作深度尺,不作模式间 MLD 评判。

## 三、regime 交叉(报告项)

core − winter max MLD 中位(m):rotation 206 / strain 223;resolved_down
231 / rest 210;resolved_up 264 / rest 213;water_mass 212 / heave 407
(n=2,不足判)。各组都在同量级,无强结构——"位于通风层之下 100–450 m"
是全体事件的共同前提,不是某一 regime 的特征。

## 四、与 w_along(G1)的衔接

- 全体核心在通风层之下(G3)+ 总体无沿等密面对准富集(G1)→ 向下输送
  由解析 w 或未解析过程承担。
- resolved-downward 15 事件同时具备:拉格朗日下沉 + 核心运动学 w_along
  为正(G1)+ 核心位于通风层下 ~230 m(G3)→ 这批事件是"冬季通风层水团经沿等密面坡向下输运"的**最强 internally
  consistent subset**;注意 G1 与轨迹使用同一套 OFES 流场,属同源互补
  诊断,不是独立证据(见 w_along 结果报告)。
