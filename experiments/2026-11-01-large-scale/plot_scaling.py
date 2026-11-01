"""Shareable runtime-versus-read-count figure from saved native workflow timings."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SOURCE = HERE/'raw/summary.json'
DATA = json.loads(SOURCE.read_text())
OUT = HERE/'figures'
OUT.mkdir(exist_ok=True)
METHODS = [
    ('ours', 'Constellation', '#087F8C', 'o'),
    ('piscem', 'piscem + alevin-fry', '#3D64AE', 's'),
    ('kallisto', 'kallisto + bustools', '#B56A29', 'D'),
]


def duration(seconds):
    rounded = int(seconds + .5)
    minutes, seconds = divmod(rounded, 60)
    return f'{minutes}m {seconds:02d}s' if minutes else f'{seconds}s'


plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 12,
    'text.color': '#17283B', 'axes.labelcolor': '#435069',
    'xtick.color': '#435069', 'ytick.color': '#435069',
    'axes.edgecolor': '#CBD5DF', 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.axisbelow': True,
    'svg.fonttype': 'none', 'svg.hashsalt': 'constellation-read-scaling-share',
})
fig = plt.figure(figsize=(10, 7.5), facecolor='white')
fig.text(.065, .955, 'CONSTELLATION  /  BENCHMARK', fontsize=11, weight='bold', color='#087F8C')
fig.text(.065, .90, 'Read mapping + counting runtime', fontsize=25, weight='bold')
fig.text(.065, .854, '100M → 500M read pairs · PBMC 10k v3 · lower is faster', fontsize=13.5, color='#526175')
ax = fig.add_axes([.12, .29, .595, .495])
ax.set_xlim(50, 520)
ax.set_ylim(0, 28)
ax.set_xticks([100, 200, 300, 400, 500])
ax.set_yticks([0, 5, 10, 15, 20, 25])
ax.set_xlabel('Read pairs (millions)', fontsize=13, labelpad=12)
ax.set_ylabel('Runtime (minutes)', fontsize=13, labelpad=12)
ax.grid(axis='y', color='#E6EBF0', linewidth=.8)
ax.tick_params(length=4, pad=7)
plotted = {}
for key, name, color, marker in METHODS:
    seconds = [DATA['cases'][f'{n}m-{key}']['pipeline_seconds'] for n in (100, 500)]
    minutes = [s/60 for s in seconds]
    plotted[key] = {'pairs': [100_000_000, 500_000_000], 'seconds': seconds}
    ax.plot([100, 500], minutes, color=color, marker=marker, markersize=8,
            linewidth=3, markeredgecolor='white', markeredgewidth=1.2,
            zorder=4 if key == 'ours' else 3)
    ax.text(1.045, minutes[1]/28, f'{name}\n{duration(seconds[1])}',
            transform=ax.transAxes, color=color, fontsize=13,
            weight='bold', va='center', linespacing=1.5, clip_on=False)

speedup = DATA['speedups']['500']['piscem']['map_and_count']
fig.text(.065, .151, f'{speedup:.2f}× faster than piscem + alevin-fry at 500M pairs',
         fontsize=14.5, color='#087F8C', weight='bold')
fig.text(.065, .108, 'Ryzen 9 7900X · 12 physical cores · prebuilt splici indexes · one run per configuration',
         fontsize=9.8, color='#526175')
fig.text(.065, .080, 'Constellation: conservative + UMI correction. Competitor exact-UMI policies differ.',
         fontsize=9.8, color='#526175')
fig.text(.065, .052, 'Only 100M and 500M were measured here; lines connect observations. Index building excluded.',
         fontsize=9.8, color='#526175')
fig.text(.065, .023, 'Data & full protocol: github.com/gwcmowry/constellation',
         fontsize=9.8, color='#526175')

for extension in ('png', 'svg'):
    path = OUT/f'runtime_vs_read_pairs.{extension}'
    fig.savefig(path, dpi=200, facecolor='white', metadata={'Date': None} if extension == 'svg' else {})
plt.close(fig)
(HERE/'raw/scaling-figure-data.json').write_text(json.dumps({
    'source': 'raw/summary.json', 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'runtime': 'native mapping + complete native counting; prebuilt indexes',
    'constellation_policy': 'conservative intersection + directional UMI correction',
    'series': plotted, 'replicates': 1, 'intermediate_sizes_measured': False,
    'output_pixels': [2000, 1500], 'matplotlib_version': matplotlib.__version__,
}, indent=2)+'\n')
print('Wrote 2000×1500 PNG and vector SVG runtime-versus-read-pairs figures.')
