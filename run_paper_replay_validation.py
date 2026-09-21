#!/usr/bin/env python
"""Replay saved small-table bootstrap results and test validation-report failures.

No detector, matching, raw fields, or scientific producers are invoked. Functions
are AST-loaded to avoid importing track.py's unrelated scientific dependencies.
"""
from __future__ import annotations
import argparse
import ast
import contextlib
import hashlib
import io
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd


def load_functions(path, names, namespace):
    tree = ast.parse(path.read_text())
    nodes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name in names]
    assert {n.name for n in nodes} == set(names)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), namespace)


def main(source: Path, package: Path, formal: Path, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    ns = dict(np=np, pd=pd, hashlib=hashlib)
    track_names = ['_matched_set_cells', '_mantel_haenszel_odds_ratio', '_matched_dependency_components',
                   '_bootstrap_matched_or_by_cluster', '_stable_analysis_seed', '_scv_matched_scope_mask',
                   '_format_detection_value']
    load_functions(source/'track.py', track_names, ns)
    ns.update(track=SimpleNamespace(**{k: ns[k] for k in track_names}), THRESHOLDS=(20.,35.,50.),
              SCOPES=('Global','Global excluding KE'), BOOTSTRAP_ITERATIONS=2000, RANDOM_SEED=42)
    load_functions(source/'run_paper_ready_delivery.py', ['restore_bootstrap_input_order', 'bootstrap_facts',
                   'bootstrap_status', 'summarize_queue', 'pair_queue'], ns)
    queue = pd.read_csv(package/'tables/qualified_analysis_queue.csv')
    # Original runner's left merges and row filters preserve this cached cohort order.
    # Recover by identities only; never search permutations against interval endpoints.
    keys = ['match_set_id','profile_number','is_scv']
    if formal is None:
        linked = queue.copy()
        ns["restore_bootstrap_input_order"](linked)
    else:
        original = pd.read_parquet(formal, columns=keys)
        assert not original.duplicated(keys).any() and not queue.duplicated(keys).any()
        original['_source_order'] = np.arange(len(original))
        linked = queue.drop(columns=['bootstrap_input_order'], errors='ignore').merge(
            original, on=keys, how='left', validate='one_to_one', sort=False)
        assert linked['_source_order'].notna().all() and len(linked)==504
        linked['bootstrap_input_order'] = linked['_source_order'].rank(method='first').astype(int)-1
        linked = linked.drop(columns='_source_order')
        if 'bootstrap_input_order' in queue:
            np.testing.assert_array_equal(queue.bootstrap_input_order, linked.bootstrap_input_order)
    expected = pd.read_csv(package/'tables/scv_matched_effects.csv')
    recovered = ns['restore_bootstrap_input_order'](linked)
    results = []
    baseline = None
    for name, frame in [('original',recovered), ('export_order',linked),
                        ('shuffled',linked.sample(frac=1,random_state=20260922))]:
        actual = ns['summarize_queue'](frame, 'primary_metadata_supported')
        if baseline is None:
            baseline=actual
            cols = [c for c in expected if pd.api.types.is_numeric_dtype(expected[c])]
            np.testing.assert_allclose(actual[cols].to_numpy(float), expected[cols].to_numpy(float), rtol=1e-13, atol=1e-13)
        else:
            pd.testing.assert_frame_equal(actual, baseline, check_exact=True)
        results.append(dict(test=name,passed=True,rows=len(actual)))
        print(name+' replay passed', flush=True)
    # Missing and duplicate ordinals must not silently create new sampling orders.
    for bad in [linked.drop(columns='bootstrap_input_order'), linked.assign(bootstrap_input_order=0)]:
        try: ns['restore_bootstrap_input_order'](bad)
        except ValueError: pass
        else: raise AssertionError('Invalid order accepted')
    fresh=ns['pair_queue'](linked.drop(columns='bootstrap_input_order'),linked.qualification_pass.astype(bool))
    np.testing.assert_array_equal(fresh.bootstrap_input_order,np.arange(len(fresh)))
    details=dict(catalog_rows_tested=0,status_counts={},confirmed_mismatch_checks=0,confirmed_mismatch_failures=0,
        confirmed_mismatch_source_status_counts={},metadata_supported_checks=0,metadata_supported_failures=0,
        special_rows_648_649_1110=[])
    for fail in (False,True):
        test_ns=dict(pd=pd,Path=Path,json=json,track=SimpleNamespace(__file__=str(source/'track.py')),
            synthetic_validation=lambda: [dict(test=f'synthetic_{i}',passed=not(fail and i==2)) for i in range(6)],
            real_263_validation=lambda: ([],details),screen_compatibility_validation=lambda: dict(passed=True))
        load_functions(source/'run_matcher_source_fix_validation.py',['main'],test_ns)
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            code=0
            try: test_ns['main'](Path(tmp))
            except SystemExit as exc: code=exc.code
            payload=json.loads((Path(tmp)/'matcher_source_fix_validation.json').read_text())
            note=(Path(tmp)/'matcher_source_fix_note.md').read_text()
            assert code==int(fail) and payload['synthetic_passed']==6-int(fail)
            assert f"Synthetic cases: {6-int(fail)}/6 passed" in note
        results.append(dict(test='matcher_report_failure' if fail else 'matcher_report_success',passed=True,exit_code=code))
    # Write only the order-augmented queue and validation artifacts, never effects.
    linked.to_csv(output/'qualified_analysis_queue.csv', index=False)
    baseline.to_csv(output/'replayed_effects.csv', index=False)
    (output/'replay_checks.json').write_text(json.dumps(dict(checks=results, formal_order_source=str(formal) if formal else 'recorded bootstrap_input_order',
        recovery='Original cohort identities' if formal else 'Recorded exported ordinals',
        invalid_order_checks=2, new_queue_order_check=True),indent=2)+'\n')
    print('All bounded replay and report checks passed',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--package',type=Path,required=True)
    p.add_argument('--formal',type=Path,help='Original cohort for one-time recovery; omit when queue includes order')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); main(a.source,a.package,a.formal,a.output)
