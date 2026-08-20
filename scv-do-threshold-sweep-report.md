# McCoy SCV ΔDO 阈值安全富集结果报告(全球)

- 日期:2026-08-16
- 分支:`main`(worktree `Oxygen`)
- 数据源:冻结 parquet
  `plot_outputs/do/global_ocean/scv_do_threshold_sweep/
  scv_do_threshold_sweep_2002_2023_depth300m.parquet`(OR 与 Wilson CI)、
  `mccoy_scv_delta_do_thresholds_2002_2023_depth300m.parquet`(阈值安全表)、
  `screen_mccoy_scvs_against_glorys/mccoy_glorys_miss.parquet`(检测器
  fidelity)。
- 本文档是正式数字的回填位置;sweep 输出的 csv/parquet/txt 不回填。

## 一、主结果:OR 阶梯

统一 `P(DO_t | group)`,三组:All DO-evaluable Argo(n=230117)、
META-matched DO-evaluable Argo(n=93161)、DO-evaluable McCoy SCVs
(冻结主分析 n=244;一剖面敏感性见第四节)。

| 阈值 | All Argo 载率 | META 载率 | SCV 载率(载/分母) | SCV OR [95% CI] | Fisher p(nominal) |
|---|---|---|---|---|---|
| DO20 | 0.0606 | 0.0680 | 0.2746 (67/244) | **5.86 [4.42, 7.77]** | 6.4e-26 |
| DO35 | 0.0144 | 0.0175 | 0.1598 (39/244) | **13.01 [9.22, 18.35]** | 2.7e-28 |
| DO50 | 0.0037 | 0.0042 | 0.0697 (17/244) | **20.22 [12.30, 33.26]** | 1.4e-16 |

- OR 随极端强度单调爬升(5.9→13.0→20.2);META 表面涡对照全程 ≈1.1–1.2。
- Fisher p 是 nominal 描述,不作空间独立性或因果证据。

## 二、区域结构:最强载体集中在 KE

按冻结 source table 的 `basin` 字段分盆地,并从 Pacific 中单列
25–45°N、140–180°E 的 KE 区。下表是**McCoy SCV 内的 carrier fraction**,
不是各区相对背景的区域 OR:

| 区域(n,各阈值相同) | DO20 | DO35 | DO50 |
|---|---:|---:|---:|
| KE(n=54) | 35/54 (64.8%) | 29/54 (53.7%) | **16/54 (29.6%)** |
| Pacific outside KE(n=50) | 13/50 (26.0%) | 3/50 (6.0%) | 0/50 |
| Atlantic(n=64) | 7/64 (10.9%) | 3/64 (4.7%) | 1/64 (1.6%) |
| Southern Ocean(n=51) | 12/51 (23.5%) | 4/51 (7.8%) | 0/51 |
| Indian(n=26) | 0/26 | 0/26 | 0/26 |

DO50 的 17 个 carrier 中 16 个位于 KE,说明全球样本中的最强关联具有
显著区域集中性。该结构加强了选择 KE OFES 作为 process analogue 的依据;
它不支持把 raw OR 表述为各海盆强度均匀的普遍效应。区域背景校正仍由既有
空间标准化、分层置换和 EKE 重匹配分析承担。

## 三、模式重现能力与 DO 富集(探索性分层)

245 个 DO-evaluable McCoy 剖面中,160 个具有可直接比较的 GLORYS
重现判定(89 reproduced、71 missed;其余 85 不进入该分层)。观测 DO
异常在 GLORYS-missed 组反而更常见:

| 阈值 | GLORYS reproduced | GLORYS missed | Fisher p(nominal) |
|---|---|---|---|
| DO20 | 17/89 (19.1%) | 34/71 (47.9%) | 1.5e-4 |
| DO35 | 10/89 (11.2%) | 21/71 (29.6%) | 0.0046 |
| DO50 | 3/89 (3.4%) | 13/71 (18.3%) | 0.0026 |

这不是 McCoy 目录完整度估计,也不改写观测 DO 标签;它表明全球富集并非
只由模式中容易重现的大尺度 SCV 驱动。相反,携带强 DO 异常的观测 SCV
更常落在 GLORYS 无法重现的组中,与 OFES/GLORYS 对细薄次表层结构的
分辨率边界相呼应。该分层未控制区域、密度和季节,作为探索性分辨率诊断。

## 四、分母漏斗与 244/245 敏感性

- McCoy 原始目录 4084 条,其中 4066 条位于 2002–2023 时间窗;
- 263 条进入 BGC/DO 联合分析交集;
- 当前阈值管线判定其中 245 条 profile-level DO-evaluable;
- 冻结主 OR 沿用既有联合表的 244 条 cohort。

唯一差异是 `2017:32111`。它在当前阈值管线中可计算且为 non-carrier,
但旧联合管线记录为 `argo_no_signal`,因此未进入冻结主分母。把它作为
non-carrier 加回 n=245 后,SCV OR 为 DO20 5.83、DO35 12.94、DO50
20.13(主结果分别为 5.86、13.01、20.22),结论逐项不变。正文以 n=244
作为预先冻结的主估计量,n=245 明确标为单剖面保守敏感性,不再混用。

## 五、口径说明

- "全球"指全球目录上的分析范围,不表示各盆地载率均匀;最强 DO50 结果
  由 KE hotspot 主导。
- 相关图件:F3_global_ofes_bridge.png(OR 阶梯+桥接)、
  F9_global_scv_regional_structure.png(区域结构)、
  F17_glorys_reproduction_stratification.png(模式重现分层)、
  F18_analysis_funnels.png(分母漏斗),输出在 Oxygen-mechanism
  `plot_outputs/do/ofes_np30_ke/ofes_paper_figures/`。
