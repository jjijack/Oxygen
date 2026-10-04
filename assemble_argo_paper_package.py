#!/usr/bin/env python
"""Assemble the Argo Letter package from existing plot outputs.

The script only copies: every figure and table is produced by a track.py
function into its own plot_outputs directory, and the mapping below records
which output becomes which figure or table of the paper.  Figures go to
``argo_paper_package/figures``; their source-data CSVs and the tables go to
``argo_paper_package/tables``; the Markdown table summary becomes
``argo_paper_package/tables.md``.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
import track

REGION = 'global_ocean'
# Paper figure -> (producing function's output directory, figure file stem).
FIGURES = {
    'Fig1': ('plot_argo_do_occurrence_overview', 'argo_do_occurrence_overview'),
    'Fig2': ('plot_argo_do_occurrence_by_eddy_setting', 'argo_do_occurrence_by_eddy_setting'),
    'Fig3': ('plot_scv_matched_anchors', 'scv_matched_anchors'),
    'Fig4': ('plot_scv_case_local_context/ke_minty_spicy', 'scv_case_local_context'),
    'FigS1': ('plot_argo_do_occurrence_maps', 'argo_do_occurrence_maps'),
    'FigS2': ('do_density_coordinate_sensitivity', 'do_density_coordinate_sensitivity_do50_depth300m'),
    'FigS3': ('plot_scv_matching_diagnostics', 'scv_matching_diagnostics'),
    'FigS4': ('plot_scv_glorys_representation', 'scv_glorys_representation'),
    'FigS5': ('plot_argo_ke_sampling_by_year', 'argo_ke_sampling_by_year'),
    'FigS6': ('plot_argo_lens_case/pmove_2014', 'argo_lens_case'),
    'FigS7': ('plot_scv_case_profiles/ke_minty_spicy', 'scv_case_profiles'),
    'FigS8': ('plot_scv_case_profiles/atlantic_spicy', 'scv_case_profiles'),
}
# Source data written by producers other than the figure's plotter.
EXTRA_SOURCES = {
    'figS2_density_coordinate_details.csv': (
        'do_density_coordinate_sensitivity',
        'do_density_coordinate_sensitivity_2002_2023_do50_depth300m_tol0p0001_details.csv',
    ),
}
TABLES_DIR = 'export_scv_matched_effect_tables'
TABLES = {
    'tableS1_matched_results.csv': 'matched_results.csv',
    'tableS2_ke_only_matched_results.csv': 'ke_matched_results.csv',
    'tableS3_core_alignment_permutation_summary.csv': 'core_alignment_permutation_summary.csv',
    'tableS4_common_grid_matched_results.csv': 'common_grid_matched_results.csv',
    'tableS5_core_oxygen_contrast.csv': 'core_oxygen_contrast.csv',
    'tableS6_positive_anchors.csv': 'positive_anchors.csv',
}
TABLES_MARKDOWN = 'matched_effect_tables.md'


def source_name(figure: str, name: str) -> str:
    """Prefix a plotter's source-data file with the figure number; panel-lettered names attach directly."""
    prefix = figure[0].lower() + figure[1:]
    return f'{prefix}{name}' if re.match(r'^[a-h]+_', name) else f'{prefix}_{name}'


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--package', type=Path, default=REPO / 'argo_paper_package',
                        help='package directory (default: argo_paper_package next to this script)')
    args = parser.parse_args()
    outputs = track.make_detection_config('do').output_dir('', REGION)
    figures_dir, tables_dir = args.package / 'figures', args.package / 'tables'
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    copies = []
    for figure, (directory, stem) in FIGURES.items():
        source_dir = outputs / directory
        copies.append((source_dir / f'{stem}.png', figures_dir / f'{figure}.png'))
        copies += [(csv, tables_dir / source_name(figure, csv.name)) for csv in sorted(source_dir.glob('*.csv'))
                   if figure != 'FigS2']
    copies += [(outputs / directory / name, tables_dir / target) for target, (directory, name) in EXTRA_SOURCES.items()]
    copies += [(outputs / TABLES_DIR / name, tables_dir / target) for target, name in TABLES.items()]
    copies.append((outputs / TABLES_DIR / TABLES_MARKDOWN, args.package / 'tables.md'))

    missing = [str(source) for source, _ in copies if not source.exists()]
    if missing:
        raise FileNotFoundError('Missing plot outputs:\n  ' + '\n  '.join(missing))
    for source, target in copies:
        shutil.copy2(source, target)
        print(f'{source.relative_to(outputs)} -> {target.relative_to(args.package)}')


if __name__ == '__main__':
    main()
