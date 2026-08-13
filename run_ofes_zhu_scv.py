#!/usr/bin/env python
"""Run the locked Zhu OFES grid-SCV catalog."""

from __future__ import annotations

import argparse
from pathlib import Path

import track


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start-date', default='2003-01-01')
    parser.add_argument('--end-date', default='2003-12-31')
    parser.add_argument('--output-dir', type=Path, default=None)
    parser.add_argument('--no-resume', action='store_true')
    parser.add_argument('--skip-preflight', action='store_true')
    args = parser.parse_args()
    result = track.run_ofes_zhu_scv_catalog(
        args.start_date,
        args.end_date,
        output_dir=args.output_dir,
        resume=not args.no_resume,
        validate_preflight=not args.skip_preflight,
    )
    manifest = result['manifest']
    print(
        f"completed {manifest.get('completed_days', 0)}/"
        f"{manifest.get('total_days', 0)} days: {result['run_dir']}",
        flush=True,
    )


if __name__ == '__main__':
    main()
