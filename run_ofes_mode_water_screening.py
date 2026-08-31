"""Thin command-line runner for the frozen OFES mode-water source screen."""
from __future__ import annotations

import argparse
import contextlib
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
            stream.flush()
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=None)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--memory-limit', default=None)
    parser.add_argument('--log-file', type=Path, default=None)
    parser.add_argument('--preflight-only', action='store_true')
    parser.add_argument('--validate-only', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--cache-only', action='store_true')
    args = parser.parse_args()
    if args.log_file:
        args.log_file.parent.mkdir(parents=True, exist_ok=True)
    stream = open(args.log_file, 'a', buffering=1) if args.log_file else None
    try:
        if stream:
            print(f'[runner] parent={ROOT}', file=stream, flush=True)
        with contextlib.redirect_stdout(_Tee(sys.stdout, stream) if stream else sys.stdout), contextlib.redirect_stderr(_Tee(sys.stderr, stream) if stream else sys.stderr):
            import track
            track.switch_region('kuroshio_extension', verbose=False)
            result = track.run_ofes_mode_water_source_screening(
                args.output_dir,
                workers=args.workers,
                memory_limit=args.memory_limit,
                overwrite=args.overwrite,
                preflight_only=args.preflight_only,
                validate_only=args.validate_only,
                cache_only=args.cache_only,
            )
            if not args.validate_only and not args.preflight_only:
                print('MODE_WATER_SCREENING_COMPLETE', result.get('output_dir', result.get('preflight')))
    except Exception:
        traceback.print_exc(file=stream or sys.stderr)
        return 1
    finally:
        if stream is not None:
            stream.close()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
