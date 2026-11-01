"""Real-data count concordance, without rerunning mapping or counting.

Default: render the saved bin counts. --refresh recomputes them from the retained
normalized Parquet matrices. All nonzero-union entries are included, unsampled.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = HERE/'raw'
MATRICES = ROOT/'benchdata/large-scale-2026-11-01/matrices'
AXES = ROOT/'experiments/2026-11-01-competitor-stages/raw'
CACHE = RAW/'concordance-figure-data.json'
SOURCES = [RAW/'matrix-comparison.json', AXES/'canonical-genes.json', AXES/'canonical-cells.json']
MODES = [('count', 'Conservative + UMI correction (default)'), ('usa', 'Optional USA + exact UMIs')]
STATES = [('S', 'Spliced'), ('U', 'Unspliced'), ('A', 'Ambiguous')]


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def source_hashes():
    return {str(p.relative_to(ROOT)): digest(p) for p in SOURCES}


def pearson(n, sx, sy, xx, yy, xy):
    return float((xy-sx*sy/n)/math.sqrt((xx-sx*sx/n)*(yy-sy*sy/n)))


def prepare():
    import pyarrow.parquet as pq
    expected = json.loads(SOURCES[0].read_text())
    ng, nc = [len(json.loads(p.read_text())) for p in SOURCES[1:]]
    bins, extent = 96, 4.25  # common axes: log10(1 + count), through 17,782 counts
    paths = [MATRICES/f'500m-{tool}-{state}.parquet'
             for tool in ('piscem-count', 'ours-count', 'ours-usa') for state, _ in STATES]
    inputs = {str(p.relative_to(ROOT)): digest(p) for p in paths}

    def load(tool, state):
        table = pq.read_table(MATRICES/f'500m-{tool}-{state}.parquet')
        g = table['gene_id'].combine_chunks().to_numpy().astype(np.uint64)
        c = table['cell_id'].combine_chunks().to_numpy().astype(np.uint64)
        v = table['count'].combine_chunks().to_numpy().astype(np.uint64)
        assert (g < ng).all() and (c < nc).all() and (v > 0).all()
        keys = g*nc+c
        assert (keys[1:] > keys[:-1]).all()
        assert np.log10(1+int(v.max())) < extent
        logs = np.log1p(v.astype(np.float64))
        positions = np.floor(logs/math.log(10)/extent*bins).astype(np.int64)
        return keys, v, logs, positions

    report = {}
    for state, _ in STATES:
        kx, x, lx, bx = load('piscem-count', state)
        sx, xx = int(x.sum()), int(np.dot(x, x))
        for mode, _ in MODES:
            ky, y, ly, by = load(f'ours-{mode}', state)
            sy, yy = int(y.sum()), int(np.dot(y, y))
            assert max(sx, sy) < 2**32  # bounds uint64 count products/sums above
            common, ix, iy = np.intersect1d(kx, ky, assume_unique=True, return_indices=True)
            n = len(kx)+len(ky)-len(common)
            del common
            xy = int(np.dot(x[ix], y[iy]))
            minimum = int(np.minimum(x[ix], y[iy]).sum())
            overlap = minimum/(sx+sy-minimum)
            full_r = pearson(ng*nc, sx, sy, xx, yy, xy)
            ref = expected[f'500m-ours-{mode}_vs_500m-piscem-count'][state]
            assert sx == ref['b_total'] and sy == ref['a_total'] and n == ref['nonzero_union']
            assert math.isclose(full_r, ref['full_matrix_pearson'], abs_tol=1e-12)
            assert math.isclose(overlap, ref['weighted_jaccard'], abs_tol=1e-12)
            raw_r = pearson(n, sx, sy, xx, yy, xy)
            log_r = pearson(n, float(lx.sum()), float(ly.sum()), float(np.dot(lx, lx)),
                            float(np.dot(ly, ly)), float(np.dot(lx[ix], ly[iy])))
            hist = np.bincount(bx[ix]*bins+by[iy], minlength=bins*bins).reshape(bins, bins)
            only_x = np.ones(len(x), dtype=bool); only_x[ix] = False
            only_y = np.ones(len(y), dtype=bool); only_y[iy] = False
            hist[:, 0] += np.bincount(bx[only_x], minlength=bins)
            hist[0, :] += np.bincount(by[only_y], minlength=bins)
            assert int(hist.sum()) == n and hist[0, 0] == 0
            report[f'{mode}-{state}'] = dict(
                histogram=hist.tolist(), nonzero_union=n, shared_zeros_excluded=ng*nc-n,
                nonzero_union_raw_pearson=raw_r, nonzero_union_log1p_pearson=log_r,
                full_matrix_raw_pearson=full_r, weighted_jaccard=overlap,
                alevin_fry_molecules=sx, constellation_molecules=sy,
            )
            print(mode, state, 'r(raw)', round(raw_r, 4), 'r(log1p)', round(log_r, 4), flush=True)
            del ky, y, ly, by, ix, iy, only_x, only_y
    CACHE.write_text(json.dumps(dict(
        source_sha256=source_hashes(), parquet_sha256=inputs, genes=ng, cells=nc,
        pairs=500_000_000, bins=bins, log10p1_extent=extent,
        x='piscem + alevin-fry cr-like WTA / exact UMI', y='Constellation',
        population='all gene-cell entries nonzero in either compared state matrix; no sampling',
        correlations='Pearson on raw or log1p counts over the nonzero union; shared zeros excluded',
        histogram_axes='x then y; uniformly binned in log10(1 + molecule count)',
        rows=report,
    ), indent=2)+'\n')


def render():
    data = json.loads(CACHE.read_text())
    assert data['source_sha256'] == source_hashes(), 'Source summaries changed; rerun with --refresh'
    rows = data['rows']; extent = data['log10p1_extent']
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'text.color': '#17283B', 'axes.labelcolor': '#435069',
        'xtick.color': '#435069', 'ytick.color': '#435069',
        'axes.edgecolor': '#CBD5DF', 'axes.spines.top': False, 'axes.spines.right': False,
        'svg.fonttype': 'none', 'svg.hashsalt': 'constellation-real-count-concordance',
    })
    fig = plt.figure(figsize=(12, 10), facecolor='white')
    fig.text(.055, .965, 'CONSTELLATION  /  REAL-DATA CONCORDANCE', fontsize=11, weight='bold', color='#087F8C')
    fig.text(.055, .92, 'Count agreement with piscem + alevin-fry', fontsize=24, weight='bold')
    fig.text(.055, .881, '500M read pairs · PBMC 10k v3', fontsize=13, color='#526175')
    vmax = max(max(max(row) for row in values['histogram']) for values in rows.values())
    ticks = [0, 1, 10, 100, 1000, 10000]
    tick_positions = np.log10(1+np.array(ticks))
    labels = ['0', '1', '10', '100', '1k', '10k']
    cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('white')
    for i, (mode, title) in enumerate(MODES):
        bottom = .515 if i == 0 else .135
        fig.text(.10, bottom+.305, title, fontsize=13, weight='bold', color='#087F8C')
        for j, (state, name) in enumerate(STATES):
            record = rows[f'{mode}-{state}']
            hist = np.array(record['histogram'], dtype=np.int64)
            assert hist.sum() == record['nonzero_union']
            ax = fig.add_axes([.10+.28*j, bottom, .22, .264])
            im = ax.imshow(np.ma.masked_equal(hist.T, 0), origin='lower', extent=[0, extent, 0, extent],
                           cmap=cmap, norm=LogNorm(vmin=1, vmax=vmax), interpolation='nearest', aspect='equal')
            ax.plot([0, extent], [0, extent], '--', linewidth=.8, color='#657180', alpha=.7)
            ax.set_xticks(tick_positions, labels)
            ax.set_yticks(tick_positions, labels)
            ax.tick_params(labelsize=8.5, length=3, pad=4)
            if i == 0: ax.set_title(name, fontsize=13, weight='bold', pad=10)
            stats = (f"Pearson r = {record['nonzero_union_raw_pearson']:.4f}\n"
                     f"Overlap = {100*record['weighted_jaccard']:.1f}%")
            ax.text(.035, .96, stats, transform=ax.transAxes, va='top', fontsize=9.3,
                    linespacing=1.4, bbox=dict(facecolor='white', edgecolor='none', alpha=.90, pad=3))
    cax = fig.add_axes([.925, .135, .012, .644])
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.set_label('Gene–cell entries per bin', fontsize=10, labelpad=8)
    colorbar.ax.tick_params(labelsize=8.5)
    fig.text(.028, .49, 'Constellation molecule counts', rotation=90, va='center', fontsize=12)
    fig.text(.49, .075, 'piscem + alevin-fry molecule counts', ha='center', fontsize=12)
    fig.text(.055, .025, 'Dashed line = equal counts', fontsize=9.5, color='#526175')
    fig.text(.945, .025, 'github.com/gwcmowry/constellation', ha='right', fontsize=9.5, color='#526175')
    for extension in ('png', 'svg'):
        fig.savefig(HERE/f'figures/real_data_concordance.{extension}', dpi=200, facecolor='white',
                    metadata={'Date': None} if extension == 'svg' else {})
    plt.close(fig)
    print('Rendered 2400×2000 real-data concordance PNG and SVG.')


def render_simple():
    """One agreement measure, with a full 0–100% axis and both count policies."""
    data = json.loads(CACHE.read_text())
    assert data['source_sha256'] == source_hashes(), 'Source summaries changed; rerun with --refresh'
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 12,
        'text.color': '#17283B', 'axes.labelcolor': '#435069',
        'xtick.color': '#435069', 'ytick.color': '#17283B',
        'axes.edgecolor': '#CBD5DF', 'axes.spines.top': False,
        'axes.spines.right': False, 'axes.spines.left': False, 'axes.axisbelow': True,
        'svg.fonttype': 'none', 'svg.hashsalt': 'constellation-simple-count-agreement',
    })
    fig = plt.figure(figsize=(10, 7.5), facecolor='white')
    fig.text(.065, .955, 'CONSTELLATION  /  REAL-DATA AGREEMENT', fontsize=11, weight='bold', color='#087F8C')
    fig.text(.065, .90, 'How closely do the counts match?', fontsize=25, weight='bold')
    fig.text(.065, .852, 'Constellation vs piscem + alevin-fry · 500M read pairs', fontsize=13.5, color='#526175')
    ax = fig.add_axes([.20, .22, .75, .51])
    positions = [0, 1.5, 3]
    settings = [('count', 'Conservative (default)', '#087F8C', -.19),
                ('usa', 'USA (optional)', '#91CBCE', .19)]
    values = {}
    for mode, label, color, offset in settings:
        values[mode] = {}
        for y, (state, _) in zip(positions, STATES, strict=True):
            value = 100*data['rows'][f'{mode}-{state}']['weighted_jaccard']
            assert 0 <= value <= 100
            values[mode][state] = value
            ax.barh(y+offset, 100, height=.31, color='#EEF2F5', zorder=2)
            ax.barh(y+offset, value, height=.31, color=color, zorder=3)
            ax.text(value-1.8, y+offset, f'{value:.1f}%', ha='right', va='center',
                    color='white' if mode == 'count' else '#153C42', weight='bold', fontsize=12.5)
    ax.set_xlim(0, 100)
    ax.set_ylim(3.7, -.7)
    ax.set_yticks(positions, [name for _, name in STATES])
    ax.tick_params(axis='y', length=0, pad=12, labelsize=13)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.tick_params(axis='x', length=4, pad=7, labelsize=11)
    ax.set_xlabel('Count overlap', fontsize=12, labelpad=12)
    ax.grid(axis='x', color='#E6EBF0', linewidth=.8)
    fig.legend(handles=[Patch(facecolor=color, label=label) for _, label, color, _ in settings],
               loc='lower left', bbox_to_anchor=(.195, .765), ncol=2, frameon=False,
               fontsize=12, handlelength=1.6, columnspacing=2)
    fig.text(.065, .08, '100% = identical counts',
             fontsize=14, weight='bold', color='#087F8C')
    fig.text(.065, .026, 'github.com/gwcmowry/constellation', fontsize=9.5, color='#526175')
    for extension in ('png', 'svg'):
        fig.savefig(HERE/f'figures/count_agreement_simple.{extension}', dpi=200, facecolor='white',
                    metadata={'Date': None} if extension == 'svg' else {})
    plt.close(fig)
    (RAW/'simple-concordance-figure-data.json').write_text(json.dumps({
        'source': 'raw/concordance-figure-data.json', 'source_sha256': digest(CACHE),
        'metric': 'weighted Jaccard of matching gene-cell counts, expressed as a percentage',
        'definition': '100 * sum(min(Constellation, alevin-fry)) / sum(max(Constellation, alevin-fry))',
        'percentages': values, 'pairs': data['pairs'], 'genes': data['genes'], 'cells': data['cells'],
        'reference': data['x'], 'output_pixels': [2000, 1500],
        'matplotlib_version': matplotlib.__version__,
    }, indent=2)+'\n')
    print('Rendered simplified 2000×1500 count-overlap PNG and SVG.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh', action='store_true', help='Recompute binned data from local Parquet matrices')
    args = parser.parse_args()
    if args.refresh or not CACHE.exists(): prepare()
    render()
    render_simple()
