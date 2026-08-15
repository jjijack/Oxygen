# Grid-SCV v2 Gate Adjudication

- 日期:2026-08-16
- 分支:`feat/ofes-ke-mechanism`(worktree `Oxygen-ofes`)
- 地位:对 `ofes-mccoy-grid-scv-v2-analysis-lock.md`(下文简称 lock)节点 4
  门控的运行后裁决。**本文档不改 lock**;lock 的解析/回归/分母/预筛/
  重放要求维持原判,本文档只对"mapped cells 141/141"一项作出拆分裁决。
- 冻结方式:本文档在正式重跑前提交冻结;重跑后的最终数字写入独立结果
  报告或 `ofes-grid-scv-v2-progress.md`,**不回填本文档**。

## 一、原假设与失败记录

lock 的 seed audit 要求:"Tier-0 must reproduce 141/141 at the exact
positions and at the mapped cells; any failure is reported with its full
failure stage and is a stop-and-report gate."

2026-08-15 两次正式运行均未通过该项。失败 manifest(归档于
`validation/failed_run_20260815_production_caliber/`,代码 hash
`28b2dc434d05879a3539c3c4c434a525795fc7c57f4b5c78d12ee58b5cd86277`,
即 commit `1a71768` 的 track.py)记录的生产口径数字:

| 指标 | 失败 manifest 值 |
|---|---|
| exact-position(样本位置双线性 + 自身局地环带,混合线) | 88/141(0.6241) |
| nearest-cell(最近网格单元 + 自身局地环带,生产口径) | 89/141(0.6312) |
| 事件级(失败 manifest 记录的是 exact-position 混合线) | center 9/9、any-event **18/19**、center-vel 6/6、any-vel **10/11** |
| 事件级(生产口径,自同一运行 nearest-cell 线可算) | center 9/9、any-event **17/19**、center-vel 6/6、any-vel **10/11** |
| 背景假阳性率(per-position ring) | 0.29%(13/4480) |
| 背景不足 61 控制数 | 88 |

## 二、裁决:两套口径的混装

原 lock 把"**旧虚拟管线重放**"和"**生产网格转移**"混装在同一 141/141
要求中,本次裁决将其拆开。

- **旧虚拟管线口径**(55e24d4 虚拟 Argo 管线):同事件所有剖面共享**事件
  核心环带** reference,控制剖面在采样点**双线性插值**,剖面自身也在
  采样点双线性插值。141 个冻结阳性、事件级 9/19/6/11、背景 20 个阳性
  都是这个口径下的数字。该口径经 diag5 逐位复现验证
  (pycnocline 25.760484172932138 与存储值逐位一致)。
- **生产网格口径**(lock Tier-0 定义:"same-day, **same-position**
  120--240 km ring";track.py `_ofes_grid_v2_ring_controls` /
  `_ofes_grid_v2_tier0_column`):每个网格列围绕**自身位置**建局地环带,
  控制剖面取**最近格点**。生产检测器(节点 5 起)走这个口径。

两者在背景定义(事件核心 vs 列自身)与剖面来源(双线性插值 vs 网格列/
最近格点)上都不相同,因此**任何正确实现都无法让生产口径在 141 个
位置逐位复现冻结分类**。2026-08-15 的诊断进一步确认:

- 同一事件环带下,最近网格单元剖面 vs 采样点双线性剖面:7 个翻转
  (134/141),全部为**格内位置与插值敏感性**——样本离单元中心
  0.75--1.56 km,单元中心剖面与"双线性@单元中心"剖面逐位相同
  (证明剖面构建正确),而样本处通过、中心处失败;
- 7 个翻转的失败 stage:gaussian(lens_extent ×2、missing_lens_limits、
  no_matching_spice_n2_peak)、dynamic_height、profile_offset_qc ×2;
- 这是双线性插值的固有性质(四节点插值,非模式解析出的亚格子物理),
  不构成实现缺陷。

因此 lock 的"at the mapped cells 141/141"在**生产口径**下不可达;在
**重放口径**下(exact positions)141/141 已达且逐位一致。

## 三、三证据线

节点 4 的 seed audit 同时记录三条线:

1. **stored_pipeline_replay**(硬门):事件核心环带 + 样本位置双线性剖面。
   必须 141/141,事件级必须复现 9/19/6/11。这是旧虚拟管线结果的重放
   验证。
2. **event_reference_nearest_cell**(报告):同一事件环带 + 最近网格单元
   剖面。度量格内位置/插值敏感性。裁决依据值为 134/141(2026-08-15
   预演,正式值以结果报告为准)。
3. **production_self_ring**(报告):网格列自身局地环带 + 单元剖面,即
   生产检测器的实际行为。度量生产网格转移性能。失败 manifest 值为
   89/141(正式值以结果报告为准)。

背景控制(4480)按 lock 字面恢复 per-position ring 口径,假阳性率与
不足控制数均为报告项(见下)。

## 四、硬门集合(hard_gates)

以下各项必须全部通过,`gate_passed` 才为 True:

1. 19 项解析验证全过;
2. 12 项回归验证全过;
3. 分母完整:已知阳性恰 141、背景控制恰 4480,行缺失为 0;
4. 预筛 recall = 1.0(141/141);
5. 旧管线重放 exact Tier-0 = 1.0(141/141);
6. 事件级 9/19/6/11 重放一致(重放口径);
7. 输入文件与代码 hash 完整记录。

## 五、报告项集合(reported_items)

以下各项写入 manifest 与审计输出,但**不参与 gate_passed**:

- event_reference_nearest_cell 复现率(裁决依据 134/141);
- production_self_ring 复现率(失败 manifest 89/141);
- 背景假阳性率(生产口径 per-position ring;失败 manifest 0.29%);
- 背景 insufficient 控制数(失败 manifest 88)。**说明**:4480 行必须
  完整、缺失行必须为零(硬门);但某些控制位置自身局地环带不足 61 条
  QC 通过控制,是生产检测器应记录的 `background_iqr` 失败,数量只报告,
  不要求为零;
- 四支撑列 / 3×3 邻域召回(描述性审计,复核后记录,不注册为成功门);
- 生产口径事件级计数(自失败运行 nearest-cell 线:center 9、any 17、
  center-vel 6、any-vel 10);
- Tier-1 / Tier-2 升级率(按 lock 原义报告,不设要求)。

## 六、节点 5/6 定位

- 节点 5(event catalog)在 hard_gates 全过后,以 **post-gate exploratory
  characterization** 运行;
- **不得**把生产检出是否达到 19/9 设为新门槛;生产数字与冻结数字的差
  异如实报告;
- 在正式记录中**不得**声称"v2 已通过 McCoy 等价验证";可陈述:v2 完整
  复现旧虚拟管线(重放口径 141/141),生产网格口径转移性能为报告值。

## 七、运行纪律

- 本次裁决冻结后,正式重跑在**全新的 validation 目录**从零运行;
- 旧 validation 目录整体归档为
  `validation/failed_run_20260815_production_caliber/`,不删除;
- 长跑开始后不再修改 track.py;
- 日片生产路径的真实变更(检测器行为改动)必须重跑;统计汇总、绘图、
  报告类更新不触发重跑。
