"""第 4 步（条件保持）共同支持审计：复现整核匹配与 4B 局地块检查，并写入持久归档。

整核：`build_ofes_grid_lens_conditional_retention` 的第一阶段（氧核种子与平移候选）与匹配规则，另记逐透镜覆盖与种子几何。
局地块：内外用同一个圆盘模板（按 km 半径取水平柱）；内侧为透镜氧核内整块都在的块，外侧沿用整核的平移候选，
块内全部种子在同一 σ0 上都留在透镜外 DO20 核内才算完整块；块级性质按块重新汇总，沿用冻结的 5 个性质卡尺。
两种设计都没有形成可比样本，因此没有正式前向比较（见 OFES.ipynb 第 4 步收口小节）。

用法：python audit_ofes_lens_retention_support.py [reverse_cache_dir] [--workers 8]
"""
import argparse
import json
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

import track

REVERSE_CACHE = '/mnt/w2/scratch/user3/Oxygen-cache/ofes_scv_reverse_enrichment_v1_20260831'
OUTPUT = Path('plot_outputs/do/ofes_np30_ke/grid_lens_conditional_retention/support_audit')
FORWARD_DAYS = 30
CONTROLS_PER_LENS = 3
CALIPER_SD = 1.0
PROPERTIES = ['depth_m', 'oyashio_ratio', 'do2', 'delta_do', 'core_thickness_m']
CENTRAL_RADII_KM = (4.0, 7.5, 11.0, 15.0)
TILING_RADII_KM = (7.5, 11.0)


def column_components(columns) -> tuple[int, int]:
    """水平柱集合的 8 邻接连通块数与最大连通块柱数。"""
    columns = np.asarray(sorted(set(columns)), dtype=int)
    if not len(columns):
        return 0, 0
    image = np.zeros(tuple(columns.max(axis=0) - columns.min(axis=0) + 1), dtype=bool)
    image[columns[:, 0] - columns[:, 0].min(), columns[:, 1] - columns[:, 1].min()] = True
    labels, count = ndimage.label(image, structure=np.ones((3, 3)))
    return int(count), int(np.bincount(labels.ravel())[1:].max())


def template(radius_km: float, lat: float, grid: dict, lon_step: float) -> list[tuple[int, int]]:
    """以中心柱为原点、水平距离不超过 radius_km 的网格偏移。"""
    scale = track.approximate_degree_length(lat)
    dy = grid['step'] * float(scale['meters_per_degree_lat']) / 1000.0
    dx = lon_step * float(scale['meters_per_degree_lon']) / 1000.0
    reach = (int(radius_km // dy), int(radius_km // dx))
    return [(i, j) for i in range(-reach[0], reach[0] + 1) for j in range(-reach[1], reach[1] + 1)
            if (i * dy) ** 2 + (j * dx) ** 2 <= radius_km ** 2]


def inside_blocks(lens_seeds: pd.DataFrame, offsets: list, grid: dict, lon_step: float, tile: bool) -> list[set]:
    """透镜氧核内整块都在的局地块；tile=False 只取离氧核质心最近的一块，True 时由近及远铺满互不重叠的块。"""
    columns = set(zip(lens_seeds['source_lat_index'].astype(int), lens_seeds['source_lon_index'].astype(int)))
    full = [c for c in columns if all((c[0] + i, c[1] + j) in columns for i, j in offsets)]
    scale = track.approximate_degree_length(float(lens_seeds['lat'].mean()))
    dy = grid['step'] * float(scale['meters_per_degree_lat']) / 1000.0
    dx = lon_step * float(scale['meters_per_degree_lon']) / 1000.0
    centre = (lens_seeds['source_lat_index'].mean(), lens_seeds['source_lon_index'].mean())
    full.sort(key=lambda c: ((c[0] - centre[0]) * dy) ** 2 + ((c[1] - centre[1]) * dx) ** 2)
    used, blocks = set(), []
    for c in full:
        block = {(c[0] + i, c[1] + j) for i, j in offsets}
        if block & used:
            continue
        used |= block
        blocks.append(block)
        if not tile:
            break
    return blocks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('reverse_cache_dir', nargs='?', default=REVERSE_CACHE)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    settings = track._ofes_grid_lens_retention_settings()
    scope_root = Path(args.reverse_cache_dir).expanduser().resolve() / 'primary_300_1000'
    source = track.load_ofes_grid_lens_source_comparison()
    core = track.load_ofes_grid_lens_core_oxygen_contrast()
    units = track.load_ofes_grid_lens_units()
    objects = units['objects'].set_index('object_id')
    lenses = source['lens_sample'].copy()
    lenses['forward_available_days'] = (pd.Timestamp(settings['data_end_date']) - lenses['date']).dt.days.astype(int)
    eligible = lenses.loc[lenses['forward_available_days'].ge(FORWARD_DAYS)]
    tasks = [(date, part, objects.loc[part['object_id']], scope_root, core['end_members'], True, settings)
             for date, part in eligible.groupby('date', sort=True)]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        stage = list(executor.map(track._ofes_grid_lens_retention_candidate_day, tasks))
    patches = pd.concat([result['patches'] for result in stage], ignore_index=True)
    seeds = pd.concat([result['seeds'] for result in stage], ignore_index=True)
    for name in ('source_lat_index', 'source_lon_index', 'seed_index'):
        seeds[name] = seeds[name].astype('Int64')

    whole = OUTPUT / 'whole_core'
    whole.mkdir(parents=True, exist_ok=True)
    selected, scales = track._ofes_grid_lens_retention_match(patches, seeds, CONTROLS_PER_LENS, CALIPER_SD)
    track._ofes_grid_lens_retention_balance(patches, selected, scales).to_csv(whole / 'initial_balance.csv', index=False)
    patches.to_parquet(whole / 'candidate_patches.parquet')
    selected.to_csv(whole / 'control_selection.csv', index=False)
    scales.rename('lens_sd').to_csv(whole / 'matching_scales.csv')
    lens_patches = patches.loc[patches['material'].eq('lens')].set_index('lens_track_id')
    lenses['oxygen_core_seeds'] = lenses['track_id'].map(lens_patches['seeds']).fillna(0).astype(int)
    lenses['controls_selected'] = lenses['track_id'].map(selected.groupby('lens_track_id').size()).fillna(0).astype(int)
    lenses['status'] = np.select(
        [lenses['forward_available_days'].lt(FORWARD_DAYS), lenses['oxygen_core_seeds'].eq(0), lenses['controls_selected'].eq(0)],
        ['insufficient_forward_window', 'no_oxygen_core_seed', 'insufficient_whole_core_common_support'], default='matched',
    )
    lenses.to_csv(whole / 'lens_sample.csv', index=False)

    match = list(track._OFES_GRID_LENS_RETENTION_MATCH)
    core_lenses = lens_patches.loc[lens_patches['seeds'].gt(0)]
    candidates = patches.loc[patches['material'].eq('control') & patches['seeds'].gt(0) & patches['jet_side'].eq('south')]
    lens_seeds = seeds.loc[seeds['material'].eq('lens') & seeds['particle_group'].eq('core')]
    candidate_seeds = seeds.loc[seeds['unit_id'].isin(candidates['unit_id'])]
    coverage = []
    for lens_id in core_lenses.index:
        part = candidates.loc[candidates['lens_track_id'].eq(lens_id)]
        z = (part[match].astype(float) - core_lenses.loc[lens_id, match].astype(float)) / scales
        ratio = part['seeds'] / core_lenses.loc[lens_id, 'seeds']
        passed = z[PROPERTIES].abs().le(CALIPER_SD).all(axis=1)
        own = lens_seeds.loc[lens_seeds['lens_track_id'].eq(lens_id)]
        own_columns = list(zip(own['source_lat_index'], own['source_lon_index']))
        kept = [column_components(zip(frame['source_lat_index'], frame['source_lon_index']))
                for _, frame in candidate_seeds.loc[candidate_seeds['unit_id'].isin(part.loc[passed, 'unit_id'])].groupby('unit_id')]
        coverage.append({
            'lens_track_id': lens_id, 'lens_seeds': int(core_lenses.loc[lens_id, 'seeds']),
            'lens_columns': len(set(own_columns)), 'lens_components': column_components(own_columns)[0],
            'lens_largest_component_columns': column_components(own_columns)[1],
            'lens_q90_km': float(core_lenses.loc[lens_id, 'horizontal_q90_km']),
            'south_candidates': len(part), 'pass_properties': int(passed.sum()),
            'pass_all': int(z.abs().le(CALIPER_SD).all(axis=1).sum()),
            **{f'pass_{name}': float(z[name].abs().le(CALIPER_SD).mean()) if len(part) else np.nan for name in match},
            'kept_ratio_median': float(ratio.median()) if len(part) else np.nan,
            'property_pass_kept_ratio_median': float(ratio[passed].median()) if passed.any() else np.nan,
            'property_pass_kept_largest_component_median': float(np.median([value[1] for value in kept])) if kept else np.nan,
            'property_pass_kept_components_median': float(np.median([value[0] for value in kept])) if kept else np.nan,
        })
    pd.DataFrame(coverage).to_csv(whole / 'coverage_by_lens.csv', index=False)

    local = OUTPUT / 'local_blocks'
    local.mkdir(parents=True, exist_ok=True)
    grid = track._ofes_grid_lens_grid(pd.Timestamp(lenses['date'].iloc[0]))
    lon_step = float(np.median(np.diff(grid['lon'])))
    seed_groups = {unit: frame for unit, frame in candidate_seeds.groupby('unit_id')}
    rows, pairs = [], []
    for tile, radii in ((False, CENTRAL_RADII_KM), (True, TILING_RADII_KM)):
        for radius_km in radii:
            for lens_id, own in lens_seeds.groupby('lens_track_id'):
                offsets = template(radius_km, float(own['lat'].mean()), grid, lon_step)
                blocks = inside_blocks(own, offsets, grid, lon_step, tile)
                outside = matched = 0
                for block in blocks:
                    members = own.loc[[c in block for c in zip(own['source_lat_index'].astype(int), own['source_lon_index'].astype(int))]]
                    index = set(members['seed_index'])
                    for unit in candidates.loc[candidates['lens_track_id'].eq(lens_id), 'unit_id']:
                        hit = seed_groups[unit].loc[seed_groups[unit]['seed_index'].isin(index)]
                        if len(hit) < len(index):
                            continue
                        outside += 1
                        z = (hit[PROPERTIES].mean() - members[PROPERTIES].mean()) / scales[PROPERTIES]
                        matched += int(z.abs().le(CALIPER_SD).all())
                        if tile:
                            pairs.append({'radius_km': radius_km, 'lens_track_id': lens_id, 'candidate_id': unit,
                                          'block_seeds': len(index), **{f'z_{name}': value for name, value in z.items()}})
                rows.append({'blocks': 'tiled' if tile else 'central', 'radius_km': radius_km, 'lens_track_id': lens_id,
                             'template_columns': len(offsets), 'inside_blocks': len(blocks),
                             'outside_full_blocks': outside, 'matched_pairs': matched})
    pd.DataFrame(rows).to_csv(local / 'block_support.csv', index=False)
    pd.DataFrame(pairs).to_csv(local / 'full_block_residuals.csv', index=False)

    manifest = {
        'analysis': 'ofes_grid_lens_conditional_retention_support_audit',
        'status': 'insufficient_common_support',
        'git_commit': subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip(),
        'script': Path(__file__).name,
        'inputs': {
            'reverse_cache_dir': str(scope_root.parent),
            **{name: {'output_dir': str(Path(result['output_dir']).resolve()), 'updated_at_utc': result['manifest'].get('updated_at_utc')}
               for name, result in (('source_comparison', source), ('core_oxygen_contrast', core), ('lens_units', units))},
        },
        'parameters': {
            'forward_days': FORWARD_DAYS, 'controls_per_lens': CONTROLS_PER_LENS, 'caliper_sd': CALIPER_SD,
            'whole_core_matching_variables': match, 'local_block_properties': PROPERTIES,
            'central_radii_km': list(CENTRAL_RADII_KM), 'tiling_radii_km': list(TILING_RADII_KM),
            'control_max_radius_km': settings['control_max_radius_km'], 'control_spacing_km': settings['control_spacing_km'],
            'matching_scales': {name: float(value) for name, value in scales.items()},
        },
        'counts': {
            'lenses': int(len(lenses)), 'forward_eligible': int(len(eligible)), 'with_oxygen_core': int(len(core_lenses)),
            'whole_core_matched': int(selected['lens_track_id'].nunique()) if len(selected) else 0,
            'candidates_evaluated': int(patches['material'].eq('control').sum()),
        },
        'updated_at_utc': pd.Timestamp.now(tz='UTC').isoformat(),
    }
    (OUTPUT / 'manifest.json').write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
