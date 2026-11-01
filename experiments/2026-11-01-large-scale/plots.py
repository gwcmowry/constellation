"""Render README figures/tables from saved results; never run native benchmarks.

Use the UV-managed benchmark Python with matplotlib==3.10.7 installed.
PNG files are embedded in Markdown; SVG files retain vector text for export.
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = HERE / 'raw'
FIGURES = HERE / 'figures'
BUILD_RAW = ROOT / 'experiments/2026-11-01-competitor-stages/raw'
WORKFLOWS = [
    ('ours', 'count', 'Constellation · conservative', '#087F8C'),
    ('ours', 'usa', 'Constellation · USA', '#087F8C'),
    ('piscem', 'count', 'piscem + alevin-fry', '#3D64AE'),
    ('kallisto', 'count', 'kallisto + bustools', '#B56A29'),
]
BUILD_NAMES = {'ours': 'Constellation', 'piscem': 'piscem', 'kallisto': 'kallisto'}


def load_data():
    summary = json.loads((RAW / 'summary.json').read_text())
    build = json.loads((BUILD_RAW / 'summary.json').read_text())
    rows = {}
    for size in (100, 500):
        rows[size] = []
        for tool, mode, label, color in WORKFLOWS:
            case = summary['cases'][f'{size}m-{tool}']
            count = case['stages']['usa']['seconds'] if mode == 'usa' else case['count_seconds']
            total = case['map_seconds'] + count
            peak = max(case['stages']['map']['rss_gib'], case['stages']['usa']['rss_gib']) if mode == 'usa' else case['pipeline_peak_rss_gib']
            rows[size].append(dict(tool=tool, mode=mode, label=label, color=color,
                                  mapping=case['map_seconds'], counting=count, total=total,
                                  rss=peak, throughput=size / total,
                                  build_plus_pipeline=build['build'][tool]['seconds'] + total))
    return summary, build, rows


def runtime_table(rows):
    lines = ['| Pairs | Workflow | Map (s) | Count (s) | Total (s) | Peak RSS (GiB) |',
             '|---|---|---:|---:|---:|---:|']
    for size in (100, 500):
        for r in rows[size]:
            lines.append(f"| {size}M | {r['label']} | {r['mapping']:.2f} | {r['counting']:.2f} | **{r['total']:.2f}** | {r['rss']:.2f} |")
    return '\n'.join(lines)


def build_tables(build, rows):
    lines = ['| Builder | Build wall time (s) | Peak RSS (GiB) | Index size (GB) |',
             '|---|---:|---:|---:|']
    for tool, label in BUILD_NAMES.items():
        b = build['build'][tool]
        lines.append(f"| {label} | **{b['seconds']:.2f}** | {b['rss_gib']:.2f} | {build['index_bytes'][tool]/1e9:.3f} |")
    lines += ['', '**Calculated build + mapping + counting totals:** the saved build time plus',
              'the measured sample workflow, not a new combined timing run.', '',
              '| Workflow | 100M pairs (s) | 500M pairs (s) |', '|---|---:|---:|']
    for a, b in zip(rows[100], rows[500], strict=True):
        lines.append(f"| {a['label']} | {a['build_plus_pipeline']:.2f} | {b['build_plus_pipeline']:.2f} |")
    return '\n'.join(lines)


def render_build_section(figure_prefix='figures/', source_prefix='../2026-11-01-competitor-stages/'):
    _, build, rows = load_data()
    p, k = [build['build'][t]['seconds'] / build['build']['ours']['seconds'] for t in ('piscem', 'kallisto')]
    return f'''### Index construction and first-sample cost

![Index build wall time, peak memory and index size for the three native builders]({figure_prefix}index_build.png)

[Vector SVG]({figure_prefix}index_build.svg). Prior measurements from this same
Ryzen 9 7900X host, twelve physical cores, one fresh build per tool from the same
prepared Ensembl 93 splici FASTA (k31). No builds or large read tests were rerun
to make these figures. Reference preparation from genome/GTF is excluded.

{build_tables(build, rows)}

Constellation's measured build is {p:.2f}× faster than piscem's and {k:.2f}×
faster than kallisto's. Build cost is independent of input read volume **for a
fixed reference and configuration** and can be amortized across samples.
The indexes retain different information: Constellation uses sparse k31/s20
gene/state keys; the competitors build denser target/path graphs. These are
native workflow costs, not equal-information index-construction workloads.
GB denotes decimal disk size; GiB denotes binary memory size.

[Build protocol and original measurements]({source_prefix}README.md) ·
[Raw build results]({source_prefix}raw/summary.json).
'''


def render_quality_section(figure_prefix='figures/', report_prefix='../'):
    truth_root = ROOT/'experiments/2026-10-30-consensus/raw'
    sources = [json.loads((truth_root/name).read_text()) for name in ('heldout-accuracy.json', 'heldout150-accuracy.json')]
    ours = [d['constellation-consensus']['positive'] for d in sources]
    gaps = [100*(d['constellation-consensus']['positive']['unique_recall']-d['piscem-forward']['positive']['unique_recall']) for d in sources]
    return f'''### Mapping quality on known-source holdouts

![Controlled holdout gene-assignment precision and correct unique recall at 91 bp and 150 bp]({figure_prefix}mapping_quality.png)

[Sharing PNG]({figure_prefix}mapping_quality.png) · [Vector SVG]({figure_prefix}mapping_quality.svg).
The same consensus gene-decision rule used by `map --splicing-state` reaches
**{100*ours[0]['unique_precision']:.3f}% precision at 91 bp** and
**{100*ours[1]['unique_precision']:.3f}% at 150 bp**. Correct unique recall is
{abs(gaps[0]):.2f} percentage points below piscem at 91 bp and {gaps[1]:.2f} points
above it at 150 bp, so these results do not establish zero accuracy/recall tradeoff.

Each length has 360,000 positives: 90,000 source windows generated from genome/GTF
without mapper code, each with four error conditions. Precision is correct unique source-gene
calls divided by unique calls on positives; recall divides by all positives.
These are read-level controlled scores before UMI resolution, not biological
PBMC error rates. Optional USA resolution is not shown. Intergenic and random-DNA
controls and per-error breakdowns are in the [full holdout report]({report_prefix}2026-10-30-consensus/README.md).

Separately, the [build/count optimization audit]({report_prefix}2026-10-31-build-count/raw/validation.json)
found byte-identical assignments and all seven count matrices on 10M real pairs
before versus after those implementation changes. This checks output preservation
within our conservative pipeline, not identity with competitor outputs.
'''


def render_concordance_section(figure_prefix='figures/', data_prefix='raw/', include_table=True):
    data = json.loads((RAW/'concordance-figure-data.json').read_text())
    section = f'''### Real-data count agreement

![Count overlap with piscem and alevin-fry for spliced, unspliced and ambiguous counts, comparing default and optional USA modes]({figure_prefix}count_agreement_simple.png)

[Sharing PNG]({figure_prefix}count_agreement_simple.png) · [SVG]({figure_prefix}count_agreement_simple.svg).
Each bar answers: **how closely do the count values match for the same genes
in the same cells?** Higher is closer; **100% means identical counts everywhere**.
For example, counts of 10 versus 8 at one gene/cell give 80% overlap. The overall
score sums the smaller counts and divides by the sum of the larger counts.

Both Constellation modes are compared with the same piscem/alevin-fry result
from **500M PBMC read pairs**. Optional USA gives higher overlap in all three
layers here. The default also corrects UMIs; USA and alevin-fry use exact UMIs,
so the counting rules differ. This measures **agreement, not mapping accuracy**.

[Detailed scatterplots and correlations (PNG)]({figure_prefix}real_data_concordance.png) ·
[Detailed SVG]({figure_prefix}real_data_concordance.svg).

'''
    if include_table:
        section += '| Constellation mode | Spliced overlap | Unspliced overlap | Ambiguous overlap |\n|---|---:|---:|---:|\n'
        for mode, label in [('count', 'Conservative default'), ('usa', 'Optional USA')]:
            values = [100*data['rows'][f'{mode}-{state}']['weighted_jaccard'] for state in 'SUA']
            section += f'| {label} | {values[0]:.1f}% | {values[1]:.1f}% | {values[2]:.1f}% |\n'
        section += f'''
<details>
<summary>Technical scatterplots and correlation definitions</summary>

![Detailed gene-by-cell count concordance]({figure_prefix}real_data_concordance.png)

These density plots use the common {data['genes']:,}-gene × {data['cells']:,}-cell
axes. Every entry nonzero in either compared state matrix is included, without
sampling or summing genes across cells. Shared zeros are excluded. The dashed
line is equal counts; color shows the number of entries in each log-scaled bin.

Pearson correlation is calculated on raw counts and `log1p` counts over that
nonzero union. Log counts reduce the influence of highly abundant entries and
expose low-count discrepancies. These have a different population from the
full-matrix correlations that include shared zeros. Count overlap is the weighted
Jaccard score, `sum(min) / sum(max)`; shared zeros contribute nothing to either sum.
Both mapping evidence and counting policy contribute to differences. This is
concordance rather than independent biological truth.

'''
        section += '| Constellation mode | Layer | Raw-count r | Log1p-count r | Count overlap |\n|---|---|---:|---:|---:|\n'
        for mode, label in [('count', 'Conservative + UMI correction'), ('usa', 'USA, exact UMIs')]:
            for state in 'SUA':
                row = data['rows'][f'{mode}-{state}']
                section += f"| {label} | {state} | {row['nonzero_union_raw_pearson']:.4f} | {row['nonzero_union_log1p_pearson']:.4f} | {row['weighted_jaccard']:.4f} |\n"
        section += '\n</details>\n\n'
    section += f'[Binned data, matrix hashes and correlation definitions]({data_prefix}concordance-figure-data.json).\n'
    return section


def update_root_readme():
    summary, _, rows = load_data()
    ours_500 = summary['cases']['500m-ours']
    ours_100 = summary['cases']['100m-ours']
    agreement = json.loads((RAW/'matrix-comparison.json').read_text())['500m-ours-usa_vs_500m-piscem-count']
    mass = sum(v['a_total'] + v['b_total'] for v in agreement.values())
    l1 = sum(v['l1'] for v in agreement.values())
    overlap = (mass - l1) / (mass + l1)
    prefix = 'experiments/2026-11-01-large-scale/'
    block = f'''<!-- BENCHMARK RESULTS START -->
### 100M and 500M paired reads

Real PBMC 10k v3 reads, Ryzen 9 7900X / 64 GB RAM, twelve physical cores.
**One timed pass per workflow and size**, run serially with prebuilt indexes.
All methods receive identical compressed FASTQs and the same 11,769-cell list.
Spliced/unspliced/ambiguous output is enabled throughout. Timings include native
index loading, mapping output and the complete correction/grouping/counting
pipeline through matrix writing. Index builds and input preparation are excluded
from these first tables; build costs are reported separately below.

![Mapping plus counting wall times at 100M and 500M pairs, with conservative and USA Constellation results separated]({prefix}figures/pipeline_runtime.png)

[Vector SVG]({prefix}figures/pipeline_runtime.svg). Hatched segments are counting;
each panel has its own linear time axis. No error bars are implied by one run.

For sharing: [runtime versus read count (PNG)]({prefix}figures/runtime_vs_read_pairs.png)
or [vector SVG]({prefix}figures/runtime_vs_read_pairs.svg).

{runtime_table(rows)}

At 500M pairs, the conservative Constellation workflow is
**{summary['speedups']['500']['piscem']['map_and_count']:.2f}× faster than piscem/alevin-fry** and
**{summary['speedups']['500']['kallisto']['map_and_count']:.2f}× faster than kallisto/bustools** in these runs.
The conservative row includes directional UMI correction. The optional USA row
uses `--resolution usa --no-correct-umis` on the same mappings. Competitor rows
use exact UMIs with their native resolution rules; their statistical computations
and mapping evidence are not identical.
Competitor versions: piscem 0.24.0, alevin-fry 0.19.0, kallisto 0.52.0,
bustools 0.45.1.

{render_quality_section(prefix+'figures/', 'experiments/')}

![Mapping-and-counting throughput and peak resident memory at both input sizes]({prefix}figures/throughput_memory.png)

[Vector SVG]({prefix}figures/throughput_memory.svg). Throughput uses **all input
pairs divided by mapping + counting wall time**. Peak memory is the largest
native stage RSS, not the sum of stages. Constellation's mapping stays near
{ours_500['stages']['map']['rss_gib']:.2f} GiB; its default counting grows from
{ours_100['stages']['count']['rss_gib']:.2f} to {ours_500['stages']['count']['rss_gib']:.2f} GiB.
At 500M, {ours_500['count_metrics']['phase_seconds']['ingest']:.2f} of
{ours_500['count_seconds']:.2f} counting seconds go into serial ingestion/partitioning.

The optional USA output has S/U/A Pearson correlations of
**{agreement['S']['full_matrix_pearson']:.4f} / {agreement['U']['full_matrix_pearson']:.4f} /
{agreement['A']['full_matrix_pearson']:.4f}** with piscem/alevin-fry at 500M,
computed on the full gene–cell matrices including shared zeros, with weighted
count overlap **{overlap:.4f}**.
All six mapping runs processed the exact requested input counts and all 24
state matrices passed validation. These measurements establish a lead among
the tested configurations on this dataset; agreement is not biological ground
truth, and broader SOTA claims require more methods and independent datasets.

{render_concordance_section(prefix+'figures/', prefix+'raw/', include_table=False)}

{render_build_section(prefix+'figures/', 'experiments/2026-11-01-competitor-stages/')}

[Full scaling protocol, outputs and bottlenecks]({prefix}README.md) ·
[Raw timings]({prefix}raw/summary.json) · [Matrix comparisons]({prefix}raw/matrix-comparison.json) ·
[Reproducible figure script]({prefix}plots.py).
<!-- BENCHMARK RESULTS END -->
'''
    path = ROOT / 'README.md'
    body = path.read_text()
    if '<!-- BENCHMARK RESULTS START -->' in body:
        start = body.index('<!-- BENCHMARK RESULTS START -->')
        end = body.index('<!-- BENCHMARK RESULTS END -->', start) + len('<!-- BENCHMARK RESULTS END -->')
    else:
        start = body.index('## Current Performance\n') + len('## Current Performance\n\n')
        end = body.index('The [build-to-matrix experiment]', start)
        block += '\n### Earlier benchmarks and implementation history\n\n'
    path.write_text(body[:start] + block.rstrip() + '\n\n' + body[end:].lstrip('\n'))


def make_figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    _, build, rows = load_data()
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                         'text.color': '#17283B', 'axes.labelcolor': '#435069',
                         'xtick.color': '#435069', 'ytick.color': '#17283B',
                         'axes.titleweight': 'bold', 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.spines.left': False,
                         'axes.edgecolor': '#CBD5DF', 'axes.axisbelow': True,
                         'svg.fonttype': 'none', 'svg.hashsalt': 'constellation-benchmarks-2026-11-01'})
    saved = []

    def pale(color):
        return tuple(0.42 * value + 0.58 for value in to_rgb(color))

    def style(ax):
        ax.grid(axis='x', color='#E6EBF0', linewidth=.8)
        ax.tick_params(axis='y', length=0, pad=10)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.set_ylim(3.65, -.65)

    def heading(fig, title, subtitle):
        fig.text(.04, .965, title, fontsize=21, weight='bold', va='top')
        fig.text(.04, .905, subtitle, fontsize=11, color='#526175', va='top')

    def save(fig, stem):
        for extension in ('png', 'svg'):
            path = FIGURES / f'{stem}.{extension}'
            fig.savefig(path, dpi=180, facecolor='white', metadata={'Date': None} if extension == 'svg' else {})
            saved.append(str(path.relative_to(ROOT)))
        plt.close(fig)

    labels = [item[2] for item in WORKFLOWS]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.1))
    fig.subplots_adjust(left=.235, right=.955, bottom=.24, top=.75, wspace=.36)
    heading(fig, 'Mapping to splicing-state matrices',
            'PBMC 10k v3 · 12 physical cores · prebuilt indexes · same compressed inputs')
    for panel, (ax, size) in enumerate(zip(axes, (100, 500), strict=True)):
        data = rows[size]
        maximum = max(r['total'] for r in data)
        for y, r in enumerate(data):
            ax.barh(y, r['mapping'], color=r['color'], height=.58)
            ax.barh(y, r['counting'], left=r['mapping'], color=pale(r['color']),
                    edgecolor=r['color'], linewidth=.5, hatch='////', height=.58)
            ax.text(r['total'] + maximum * .025, y, f"{r['total']:,.1f} s", va='center', fontsize=11, weight='bold')
        ax.set_yticks(range(4), labels if panel == 0 else [])
        ax.set_title(f'{size}M pairs', loc='left', pad=14, fontsize=15)
        ax.set_xlabel('Wall time (seconds)', labelpad=10)
        ax.set_xlim(0, maximum * 1.30)
        style(ax)
    fig.legend(handles=[Patch(facecolor='#526175', label='Mapping'),
                        Patch(facecolor='#E1E6EC', edgecolor='#526175', hatch='////', label='Complete counting pipeline')],
               loc='lower left', bbox_to_anchor=(.23, .80), ncol=2, frameon=False, fontsize=11)
    fig.text(.04, .055, 'Independent linear time axes. One measured run per row; no index build or input preparation included.', fontsize=10, color='#526175')
    fig.text(.04, .018, 'Conservative: directional UMI correction. USA: exact UMIs. Competitor resolution policies differ.', fontsize=10, color='#526175')
    save(fig, 'pipeline_runtime')

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.1))
    fig.subplots_adjust(left=.235, right=.955, bottom=.24, top=.75, wspace=.36)
    heading(fig, 'Throughput and memory as read volume grows',
            'Complete native mapping + counting · same 12-core host · one measured run per size')
    for panel, (ax, metric, title, unit) in enumerate(zip(axes,
            ('throughput', 'rss'), ('Throughput · higher is better', 'Peak RSS · lower is better'),
            ('Million input pairs / second', 'Peak resident memory (GiB)'), strict=True)):
        maximum = max(r[metric] for data in rows.values() for r in data)
        for size, offset in ((100, -.17), (500, .17)):
            for y, r in enumerate(rows[size]):
                ax.barh(y+offset, r[metric], height=.29,
                        color=pale(r['color']) if size == 100 else r['color'],
                        edgecolor=r['color'], linewidth=.5, hatch='///' if size == 100 else None)
                ax.text(r[metric] + maximum*.025, y+offset, f"{r[metric]:.2f}", va='center', fontsize=10)
        ax.set_yticks(range(4), labels if panel == 0 else [])
        ax.set_title(title, loc='left', pad=14, fontsize=14)
        ax.set_xlabel(unit, labelpad=10)
        ax.set_xlim(0, maximum*1.22)
        style(ax)
    fig.legend(handles=[Patch(facecolor='#D4DCE5', edgecolor='#526175', hatch='///', label='100M pairs'),
                        Patch(facecolor='#526175', label='500M pairs')],
               loc='lower left', bbox_to_anchor=(.23, .80), ncol=2, frameon=False, fontsize=11)
    fig.text(.04, .055, 'Throughput = all input pairs / (mapping + counting wall time). Memory = maximum RSS of any native stage.', fontsize=10, color='#526175')
    fig.text(.04, .018, 'Prebuilt indexes; construction is excluded. USA and conservative rows reuse the same mappings.', fontsize=10, color='#526175')
    save(fig, 'throughput_memory')

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.1))
    fig.subplots_adjust(left=.15, right=.96, bottom=.31, top=.73, wspace=.40)
    heading(fig, 'A reusable index: build once per reference',
            'Prior measurements on the same Ryzen 9 7900X · 12 physical cores · prepared splici FASTA · k31')
    colors = [WORKFLOWS[i][3] for i in (0, 2, 3)]
    for panel, (ax, title, unit) in enumerate(zip(axes,
            ('Build wall time', 'Peak build memory', 'Final index size'),
            ('Seconds', 'GiB', 'GB (decimal)'), strict=True)):
        values = [(build['build'][tool]['seconds'], build['build'][tool]['rss_gib'], build['index_bytes'][tool]/1e9)[panel] for tool in BUILD_NAMES]
        for y, (value, color) in enumerate(zip(values, colors, strict=True)):
            ax.barh(y, value, height=.55, color=color)
            ax.text(value + max(values)*.025, y, f'{value:,.2f}' if panel < 2 else f'{value:.3f}', va='center', fontsize=11, weight='bold')
        ax.set_yticks(range(3), list(BUILD_NAMES.values()) if panel == 0 else [])
        ax.set_title(title, loc='left', pad=14, fontsize=14)
        ax.set_xlabel(unit, labelpad=10)
        ax.set_xlim(0, max(values)*1.35)
        style(ax)
        ax.set_ylim(2.65, -.65)
    fig.text(.04, .13, 'One saved fresh build per tool. No benchmarks rerun. Build cost is independent of read volume for a fixed reference/configuration.', fontsize=10, color='#526175')
    fig.text(.04, .077, 'Index content differs: sparse gene/state keys (Constellation) versus denser target/path graphs (competitors).', fontsize=10, color='#526175')
    fig.text(.04, .025, 'All builders start from the same prepared FASTA. Genome/GTF reference preparation is excluded.', fontsize=10, color='#526175')
    save(fig, 'index_build')

    sources = [RAW/'summary.json', RAW/'matrix-comparison.json', BUILD_RAW/'summary.json']
    manifest = {'renderer': 'matplotlib', 'matplotlib_version': matplotlib.__version__,
                'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                'figures': saved, 'rows': rows, 'builds': build['build'], 'index_bytes': build['index_bytes'],
                'native_benchmarks_rerun': False,
                'build_plus_pipeline_is_sum_of_separate_measurements': True}
    (RAW/'figure-data.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Rendered {len(saved)} PNG/SVG artifacts from saved measurements.')


if __name__ == '__main__':
    make_figures()
    update_root_readme()
