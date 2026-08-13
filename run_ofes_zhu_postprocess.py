#!/usr/bin/env python
"""Finish the locked Zhu catalog with event association and McCoy bridge."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

import track


def _wait_for_catalog(catalog_root: Path, poll_seconds: int) -> Path:
    while True:
        manifests = sorted(catalog_root.glob('*/manifest.json'))
        for manifest_path in manifests:
            manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
            if manifest.get('status') == 'complete':
                return manifest_path.parent
            if manifest.get('status') == 'failed':
                raise RuntimeError(f'Zhu catalog failed: {manifest_path}: {manifest.get("error")}')
        time.sleep(max(10, int(poll_seconds)))


def _find_event_population(root: Path) -> Path:
    candidates = []
    for path in root.glob('do/ofes_np30_ke/ofes_delta_do_catalog/*/event_diagnostics/*/event_population/*/population_peak_diagnostics.parquet'):
        manifest_path = path.parent / 'manifest.json'
        if not manifest_path.exists():
            continue
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        if manifest.get('status') == 'complete' and len(pd.read_parquet(path)) == 59:
            candidates.append(path)
    if not candidates:
        raise FileNotFoundError('No completed 59-event OFES population diagnostics found.')
    return max(candidates, key=lambda path: path.stat().st_mtime_ns)


def _find_mccoy_summary(root: Path) -> Path | None:
    candidates = []
    patterns = (
        'do/ofes_np30_ke/ofes_delta_do_catalog/*/event_diagnostics/*/event_population/*/event_lifecycle/*/mccoy_virtual_argo/*/event_summary.parquet',
        'do/ofes_np30_ke/ofes_delta_do_catalog/*/event_diagnostics/*/event_population/*/event_onset/*/mccoy_virtual_argo/*/event_summary.parquet',
    )
    for pattern in patterns:
        for path in root.glob(pattern):
            manifest_path = path.parent / 'manifest.json'
            if not manifest_path.exists():
                continue
            manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
            if manifest.get('status') == 'complete':
                candidates.append(path)
    return max(candidates, key=lambda path: path.stat().st_mtime_ns) if candidates else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog-root', type=Path, default=Path('plot_outputs/do/ofes_np30_ke/zhu_scv'))
    parser.add_argument('--project-root', type=Path, default=Path('.'))
    parser.add_argument('--poll-seconds', type=int, default=60)
    parser.add_argument('--skip-bridge', action='store_true')
    args = parser.parse_args()
    zhu_run = _wait_for_catalog(args.catalog_root, args.poll_seconds)
    population = _find_event_population(args.project_root)
    association = track.run_ofes_zhu_event_association(zhu_run, population)
    print(f'completed Zhu event association: {association["run_dir"]}', flush=True)
    if args.skip_bridge:
        return
    bridge = track.run_ofes_zhu_mccoy_bridge(zhu_run)
    print(f'completed Zhu McCoy bridge: {bridge["run_dir"]}', flush=True)
    mccoy_summary = _find_mccoy_summary(args.project_root)
    if mccoy_summary is not None:
        cross = track.build_ofes_zhu_mccoy_event_crosstab(
            association['run_dir'] / 'event_association.parquet',
            mccoy_summary,
        )
        print(f'completed Zhu McCoy event crosstab: {cross["run_dir"]}', flush=True)


if __name__ == '__main__':
    main()
