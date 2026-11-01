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
    fig.text(.055, .881, '500M pairs · PBMC 10k v3 · 58,395 genes × 11,769 cells', fontsize=13, color='#526175')
    vmax = max(max(max(row) for row in values['histogram']) for values in rows.values())
    ticks = [0, 1, 10, 100, 1000, 10000]
    tick_positions = np.log10(1+np.array(ticks))
    labels = ['0', '1', '10', '100', '1k', '10k']
    cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('white')
    for i, (mode, title) in enumerate(MODES):
        bottom = .55 if i == 0 else .18
        fig.text(.11, bottom+.285, title, fontsize=13, weight='bold', color='#087F8C')
        for j, (state, name) in enumerate(STATES):
            record = rows[f'{mode}-{state}']
            hist = np.array(record['histogram'], dtype=np.int64)
            assert hist.sum() == record['nonzero_union']
            ax = fig.add_axes([.11+.28*j, bottom, .20, .24])
            im = ax.imshow(np.ma.masked_equal(hist.T, 0), origin='lower', extent=[0, extent, 0, extent],
                           cmap=cmap, norm=LogNorm(vmin=1, vmax=vmax), interpolation='nearest', aspect='equal')
            ax.plot([0, extent], [0, extent], '--', linewidth=.8, color='#657180', alpha=.7)
            ax.set_xticks(tick_positions, labels)
            ax.set_yticks(tick_positions, labels)
            ax.tick_params(labelsize=8.5, length=3, pad=4)
            if i == 0: ax.set_title(name, fontsize=13, weight='bold', pad=10)
            stats = (f"r(raw) = {record['nonzero_union_raw_pearson']:.4f}\n"
                     f"r(log1p) = {record['nonzero_union_log1p_pearson']:.4f}\n"
                     f"Count overlap = {record['weighted_jaccard']:.3f}")
            ax.text(.035, .96, stats, transform=ax.transAxes, va='top', fontsize=9.3,
                    linespacing=1.4, bbox=dict(facecolor='white', edgecolor='none', alpha=.90, pad=3))
    cax = fig.add_axes([.925, .18, .012, .61])
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.set_label('Gene–cell entries per bin (log color scale)', fontsize=10, labelpad=8)
    colorbar.ax.tick_params(labelsize=8.5)
    fig.text(.028, .49, 'Constellation molecules per gene–cell entry', rotation=90, va='center', fontsize=12)
    fig.text(.49, .123, 'piscem + alevin-fry molecules per gene–cell entry', ha='center', fontsize=12)
    fig.text(.055, .083, 'All nonzero-union entries shown; shared zeros excluded. Axes use log10(1 + count); dashed line is equality.', fontsize=9.3, color='#526175')
    fig.text(.055, .057, 'r = Pearson over that union. Count overlap = Σmin / Σmax. Alevin-fry uses exact UMIs; counting policies differ.', fontsize=9.3, color='#526175')
    fig.text(.055, .031, 'Agreement is not biological ground truth. Data & protocol: github.com/gwcmowry/constellation', fontsize=10, color='#526175')
    for extension in ('png', 'svg'):
        fig.savefig(HERE/f'figures/real_data_concordance.{extension}', dpi=200, facecolor='white',
                    metadata={'Date': None} if extension == 'svg' else {})
    plt.close(fig)
    print('Rendered 2400×2000 real-data concordance PNG and SVG.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh', action='store_true', help='Recompute binned data from local Parquet matrices')
    args = parser.parse_args()
    if args.refresh or not CACHE.exists(): prepare()
    render()
