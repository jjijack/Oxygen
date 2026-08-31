"""按选定日期运行 OFES SCV reverse enrichment producer 或 cache reducer。"""

import argparse
import hashlib
import json
import resource
import time
from pathlib import Path

import track


BASELINE = 'c0b49ca5dcb16e2b45694fc72f4a63be3e10e7fd'


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--date', action='append', required=True, help='OFES date; repeat for multiple dates.')
    parser.add_argument('--scope', choices=('primary_300_1000', 'deep_500_1000'), action='append', default=None)
    parser.add_argument('--output-dir', type=Path, default=None)
    parser.add_argument('--cache-only', action='store_true')
    parser.add_argument('--skip-object-layer', action='store_true')
    parser.add_argument('--approved-code-sha256', action='append', default=None, help='Explicitly approved producer code SHA-256; cache-only only.')
    args = parser.parse_args()
    if args.approved_code_sha256 and not args.cache_only:
        parser.error('--approved-code-sha256 is only valid with --cache-only.')
    scopes = tuple(args.scope or ('primary_300_1000', 'deep_500_1000'))
    output = args.output_dir or track._ofes_scv_reverse_paths()['output']
    if args.cache_only:
        result = track.reduce_ofes_scv_reverse_enrichment_cache(output, scopes=scopes, dates=args.date, write_outputs=True, approved_code_sha256=args.approved_code_sha256)
        print(json.dumps({'status': 'cache_only_complete', 'output_dir': str(output), 'scopes': list(scopes), 'science_signature': result.get('science_signature')}, default=str))
        return
    progress_path = output / 'production_progress.jsonl'
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    for date in args.date:
        started = time.time()
        start_record = {'event': 'start', 'date': str(date), 'scopes': list(scopes), 'epoch_seconds': started}
        with progress_path.open('a', encoding='utf-8') as handle:
            handle.write(json.dumps(start_record) + '\n')
        print(json.dumps(start_record), flush=True)
        try:
            track.produce_ofes_scv_reverse_enrichment_day(date, scopes=scopes, output_dir=output, write_object_layer=not args.skip_object_layer)
        except Exception as error:
            failure = {'event': 'failure', 'date': str(date), 'elapsed_seconds': time.time() - started, 'peak_rss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'error': repr(error)}
            with progress_path.open('a', encoding='utf-8') as handle:
                handle.write(json.dumps(failure) + '\n')
            print(json.dumps(failure), flush=True)
            raise
        finished = {'event': 'finish', 'date': str(date), 'elapsed_seconds': time.time() - started, 'peak_rss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        with progress_path.open('a', encoding='utf-8') as handle:
            handle.write(json.dumps(finished) + '\n')
        print(json.dumps(finished), flush=True)
    s5_root = track._ofes_scv_reverse_paths()['s5']
    manifest_path = (s5_root / 'manifest.json') if (s5_root / 'manifest.json').exists() else (s5_root / 'days' / 'manifest.json')
    source_manifest = json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path.exists() else {}
    fixed_dates = source_manifest.get('dates', args.date)
    run_status = 'complete' if set(args.date) == set(fixed_dates) else 'preflight_complete'
    identity = track._ofes_scv_reverse_input_identity()
    identity['runner_code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    track._ofes_scv_reverse_atomic_json({'status': run_status, 'baseline_commit': BASELINE, 'requested_dates': args.date, 'fixed_dates': fixed_dates, 'scopes': list(scopes), 's5_expected': {'dates': 73, 'tile_days': 4599, 'object_days': 949, 'tier1_object_days': 549, 'weak_native_object_days': 309, 'strong_native_object_days': 153, 'well_resolved_object_days': 98}, 'input_identity': identity}, Path(output) / 'producer_manifest.json')
    print(json.dumps({'status': 'producer_complete', 'output_dir': str(output), 'dates': args.date, 'scopes': list(scopes)}, default=str))


if __name__ == '__main__':
    main()
