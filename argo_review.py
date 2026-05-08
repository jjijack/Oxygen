#!/usr/bin/env python
"""交互式 Argo-GLORYS 剖面复核工具（独立脚本）。

运行方式::

    # 本地有显示器
    python argo_review.py

    # SSH 转发（需要 X11）
    ssh -X user@host
    python argo_review.py

    # headless 环境（图片保存到文件，无法交互点击）
    python argo_review.py

修改下方 CONFIG 区域的参数来调整时间范围、区域、检测方法等。
"""

import os
import sys

# 自动选择最佳 backend
def _pick_backend():
    # SSH X11 转发或本地 DISPLAY
    if os.environ.get('DISPLAY'):
        for backend in ('TkAgg', 'Qt5Agg', 'GTK3Agg'):
            try:
                import matplotlib
                matplotlib.use(backend)
                import matplotlib.pyplot as plt
                plt.figure()  # 验证能否真正创建窗口
                plt.close()
                return backend
            except Exception:
                continue
    # headless fallback
    import matplotlib
    matplotlib.use('Agg')
    return 'Agg'

_backend = _pick_backend()
import matplotlib.pyplot as plt
import track

# ──────────────────────── CONFIG ────────────────────────
REGION = 'kuroshio_extension'
START_TIME = '2014-05-01'
END_TIME = '2014-05-31'
DETECTION_METHOD = 'do'
DETECTION_KWARGS = dict(
    do_threshold=50.0,
    anomaly_min_depth=300.0,
    cbar_min=10.0,
    cbar_max=100.0,
)
REVIEW_MODE = 'all'           # 'filtered' | 'selected' | 'all'
HORIZONTAL_VARIABLE = 'vorticity'
NEEDED_DEPTH = 0
XMIN, XMAX = -400.0, 400.0
YMIN, YMAX = 0.0, 1000.0
SAVE_CLICKED_FIGURES = False
VERBOSE = True
# ──────────────────────── END CONFIG ────────────────────


def main():
    track.switch_region(REGION)
    det_cfg = track.make_detection_config(DETECTION_METHOD, **DETECTION_KWARGS)

    if _backend == 'Agg':
        print("[ArgoReview] headless 模式：批量生成并保存所有 filtered profile 的 detail 图。")
        print(f"  输出目录: plot_outputs/{DETECTION_METHOD}/{REGION}/plot_interactive_argo_glorys_review/\n")

        reviewer = track.plot_interactive_argo_glorys_review(
            start_time=START_TIME,
            end_time=END_TIME,
            region=REGION,
            detection_config=det_cfg,
            review_mode=REVIEW_MODE,
            horizontal_variable=HORIZONTAL_VARIABLE,
            needed_depth=NEEDED_DEPTH,
            xmin=XMIN, xmax=XMAX,
            ymin=YMIN, ymax=YMAX,
            save_clicked_figures=True,
            show_fig=False,
            verbose=VERBOSE,
        )

        # 遍历所有可点击的 profile，逐个生成 detail 图
        targets = reviewer.profiles.loc[reviewer.profiles['review_status'] == 'filtered']
        if REVIEW_MODE == 'all':
            targets = reviewer.profiles
        elif REVIEW_MODE == 'selected':
            targets = reviewer.profiles.loc[reviewer.profiles['is_anomaly']]

        total = len(targets)
        print(f"共 {total} 个 profile 待处理。\n")
        for i, row_idx in enumerate(targets.index, 1):
            try:
                reviewer.open_profile_by_row_index(row_idx)
                plt.close('all')
                if i % 50 == 0 or i == total:
                    print(f"  进度: {i}/{total}")
            except Exception as exc:
                print(f"  [{i}/{total}] 跳过 row {row_idx}: {exc}")
                continue

        print(f"\n完成。图片已保存到上述目录。")
        return

    # GUI 模式：弹窗交互
    reviewer = track.plot_interactive_argo_glorys_review(
        start_time=START_TIME,
        end_time=END_TIME,
        region=REGION,
        detection_config=det_cfg,
        review_mode=REVIEW_MODE,
        horizontal_variable=HORIZONTAL_VARIABLE,
        needed_depth=NEEDED_DEPTH,
        xmin=XMIN, xmax=XMAX,
        ymin=YMIN, ymax=YMAX,
        save_clicked_figures=SAVE_CLICKED_FIGURES,
        verbose=VERBOSE,
        show_fig=True,
    )
    print(f"\n共 {len(reviewer.profiles)} 个 profile 可复核。")
    print("点击地图上的点查看详情，关闭主窗口退出。\n")
    plt.show()


if __name__ == '__main__':
    main()
