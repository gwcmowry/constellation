"""Shareable precision/recall figure from controlled source-origin holdouts.

No new benchmarks or model fitting. Read-level gene mapping is distinct from
downstream USA/UMI resolution and from agreement on biological PBMC data.
"""
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TRUTH = ROOT/'experiments/2026-10-30-consensus/raw'
BUILD = ROOT/'experiments/2026-10-31-build-count/raw'
FILES = [TRUTH/'heldout-accuracy.json', TRUTH/'heldout150-accuracy.json',
         BUILD/'validation.json', BUILD/'quality.json']
holdouts = {length: json.loads(path.read_text()) for length, path in zip((91, 150), FILES[:2], strict=True)}
regression = json.loads(FILES[2].read_text())
state_check = json.loads(FILES[3].read_text())
METHODS = [
    ('constellation-consensus', 'Constellation', '#087F8C', 'o'),
    ('piscem-forward', 'piscem', '#3D64AE', 's'),
    ('kallisto-forward', 'kallisto', '#B56A29', 'D'),
]

assert regression['mapping_bytes_identical'] and regression['pipeline_mapping_bytes_identical']
matrices = {k: v for k, v in regression['pipeline_count_output_bytes_identical'].items() if k.endswith('.mtx')}
assert len(matrices) == 7 and all(matrices.values())
rows = []
for group, length in enumerate((91, 150)):
    for index, (key, name, color, marker) in enumerate(METHODS):
        source = holdouts[length][key]
        p = source['positive']
        assert p['n'] == 360_000 and p['correct_unique'] + p['wrong_unique'] == p['unique']
        assert math.isclose(p['unique_precision'], p['correct_unique']/p['unique'], abs_tol=1e-14)
        assert math.isclose(p['unique_recall'], p['correct_unique']/p['n'], abs_tol=1e-14)
        assert source['random/exact']['n'] == 30_000 and source['random/exact']['mapped'] == 0
        if key == 'constellation-consensus':
            # Later state-preserving output independently reproduced these gene calls.
            check = state_check['heldout' if length == 91 else 'heldout150']['counts']
            assert check['intersection_assigned'] == p['unique']
            assert check['intersection_correct'] == p['correct_unique']
        rows.append(dict(length=length, method=key, name=name, color=color, marker=marker,
                         y=group*4+index, precision=100*p['unique_precision'],
                         recall=100*p['unique_recall'], positives=p['n'], unique=p['unique'],
                         correct_unique=p['correct_unique'], wrong_unique=p['wrong_unique'],
                         candidate_recall=100*p['candidate_set_recall'],
                         intergenic_assigned=100*source['intergenic/exact']['assignment_rate'],
                         random_assigned=source['random/exact']['mapped']))
ours = [r for r in rows if r['method'] == 'constellation-consensus']
peers = [r for r in rows if r['method'] == 'piscem-forward']
recall_difference = [a['recall']-b['recall'] for a, b in zip(ours, peers, strict=True)]

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11,
    'text.color': '#17283B', 'axes.labelcolor': '#435069',
    'xtick.color': '#435069', 'ytick.color': '#435069',
    'axes.edgecolor': '#CBD5DF', 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.spines.left': False, 'axes.axisbelow': True,
    'svg.fonttype': 'none', 'svg.hashsalt': 'constellation-mapping-quality-share',
})
fig = plt.figure(figsize=(10, 7.5), facecolor='white')
fig.text(.065, .955, 'CONSTELLATION  /  VALIDATION', fontsize=11, weight='bold', color='#087F8C')
fig.text(.065, .90, 'Mapping accuracy on controlled holdouts', fontsize=23, weight='bold')
fig.text(.065, .854, 'Known-source simulations · 360,000 positive reads per length',
         fontsize=12.5, color='#526175')

for panel, (metric, left, title, definition, limits, ticks, decimals) in enumerate([
    ('precision', .25, 'Unique-gene precision', 'Correct / unique calls', (99.90, 100.00), [99.90, 99.95, 100.00], 3),
    ('recall', .66, 'Correct unique recall', 'Correct / all positive reads', (80, 92), [80, 85, 90], 2),
]):
    fig.text(left, .80, title, fontsize=14, weight='bold')
    fig.text(left, .772, definition, fontsize=9, color='#526175')
    ax = fig.add_axes([left, .22, .245, .525])
    ax.set_xlim(*limits)
    ax.set_ylim(6.7, -.7)
    ax.set_xticks(ticks)
    ax.set_yticks([r['y'] for r in rows], [f"{r['name']} · {r['length']} bp" for r in rows] if panel == 0 else [])
    if metric == 'precision': ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlabel(f'{metric.capitalize()} (%)', fontsize=11, labelpad=10)
    ax.grid(axis='x', color='#E6EBF0', linewidth=.8)
    ax.axhline(3, color='#E6EBF0', linewidth=.8)
    ax.tick_params(axis='y', length=0, pad=10)
    ax.tick_params(axis='x', length=4, pad=7, labelsize=10)
    for r in rows:
        ax.scatter(r[metric], r['y'], color=r['color'], marker=r['marker'], s=72,
                   edgecolor='white', linewidth=1, zorder=3)
        ax.annotate(f"{r[metric]:.{decimals}f}%", (r[metric], r['y']),
                    xytext=(8, 0), textcoords='offset points', ha='left', va='center',
                    fontsize=11.5, color=r['color'],
                    weight='bold' if r['name'] == 'Constellation' else 'normal',
                    annotation_clip=False)
    for tick, r in zip(ax.get_yticklabels(), rows):
        tick.set_color(r['color'])
        if r['name'] == 'Constellation': tick.set_weight('bold')

fig.text(.065, .085, f"Constellation: {ours[0]['precision']:.3f}–{ours[1]['precision']:.3f}% unique-gene precision",
         fontsize=14, color='#087F8C', weight='bold')
fig.text(.065, .025, 'github.com/gwcmowry/constellation', fontsize=9.5, color='#526175')

for extension in ('png', 'svg'):
    fig.savefig(HERE/f'figures/mapping_quality.{extension}', dpi=200, facecolor='white',
                metadata={'Date': None} if extension == 'svg' else {})
plt.close(fig)
(HERE/'raw/quality-figure-data.json').write_text(json.dumps({
    'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in FILES},
    'rows': rows, 'constellation_policy': 'read-level consensus, as enabled by map --splicing-state',
    'precision_definition': 'correct unique source-gene calls / unique calls on positive reads',
    'recall_definition': 'correct unique source-gene calls / all positive reads',
    'recall_difference_vs_piscem_percentage_points': recall_difference,
    'optimization_regression': {'pairs': 10_000_000, 'mapping_bytes_identical': True, 'seven_matrices_bytes_identical': True},
    'limitations': ['Controlled simulated source origins, not biological PBMC error rates.',
                    'Read-level mapping, not downstream molecule/UMI accuracy.',
                    'Error variants share source windows; no independent-replicate error bars are claimed.',
                    'Shared sequence can make the sampled source origin unidentifiable.'],
    'output_pixels': [2000, 1500], 'matplotlib_version': matplotlib.__version__,
}, indent=2)+'\n')
print('Wrote controlled mapping-quality PNG and SVG from saved counts.')
