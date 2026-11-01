"""Write the final measured report after all six native workflows and comparisons."""
import ast,json,subprocess
from pathlib import Path
from plots import render_build_section, update_root_readme
ROOT=Path(__file__).resolve().parents[2];RAW=Path(__file__).parent/'raw'
subprocess.run(['python3',str(RAW.parent/'summarize.py')],stdout=(RAW/'summary.log').open('w'),check=True)
s=json.loads((RAW/'summary.json').read_text());assert len(s['cases'])==6
comparisons=json.loads((RAW/'matrix-comparison.json').read_text());assert len(comparisons)==6
ec_bytes={size:next(f['bytes'] for f in json.loads((RAW/f'{size}m-kallisto.artifacts.json').read_text())['files'] if f['path']=='map/matrix.ec') for size in (100,500)}
text='''## Measured results

**One timed pass per tool and size**, same genuine paired reads, twelve cores,
prebuilt indexes. Combined times are mapping plus the complete native counting
pipeline. They exclude input preparation, index construction, validation and
optional additional counting policies.

![Mapping plus counting wall times, with conservative and USA results separated](figures/pipeline_runtime.png)

[Vector SVG](figures/pipeline_runtime.svg). Counting is hatched; panels use
independent linear time axes. The USA row reuses the same Constellation mapping.

For sharing: [runtime versus read count (PNG)](figures/runtime_vs_read_pairs.png)
or [vector SVG](figures/runtime_vs_read_pairs.svg). This line chart uses the
conservative Constellation workflow and only the two measured input sizes.

| Pairs | Workflow | Mapping | Counting | Combined | Peak pipeline RSS |
|---|---|---:|---:|---:|---:|
'''
for size in (100,500):
 for tool,title in [('ours','Constellation, conservative'),('piscem','piscem + alevin-fry'),('kallisto','kallisto + bustools')]:
  c=s['cases'][f'{size}m-{tool}'];text+=f"| {size}M | {title} | {c['map_seconds']:.2f}s | {c['count_seconds']:.2f}s | {c['pipeline_seconds']:.2f}s | {c['pipeline_peak_rss_gib']:.2f} GiB |\n"
text+='\n'
for size in (100,500):
 ratios=s['speedups'][str(size)];text+=f"At {size}M, Constellation's complete native workflow is **{ratios['piscem']['map_and_count']:.2f}× faster than piscem/alevin-fry** and **{ratios['kallisto']['map_and_count']:.2f}× faster than kallisto/bustools** in these runs.\n\n"
text+='''The default Constellation count includes directional UMI edit correction;
the competitor settings here use exact UMIs and different resolution rules.
Those are native workflow comparisons, not identical statistical computations.
The extra Constellation USA count keeps exact UMIs and native-style spliced
preference, while mapping evidence still differs from the competitors.

| Pairs | Our USA count | Mapping + USA count | USA count peak RSS |
|---|---:|---:|---:|
'''
for size in (100,500):
 c=s['cases'][f'{size}m-ours'];q=c['stages']['usa'];text+=f"| {size}M | {q['seconds']:.2f}s | {c['usa_pipeline_seconds']:.2f}s | {q['rss_gib']:.2f} GiB |\n"
text+='\n'+render_build_section()+'\n'
text+='''
## Scaling and remaining bottlenecks

![Mapping-and-counting throughput and peak memory at 100M and 500M pairs](figures/throughput_memory.png)

[Vector SVG](figures/throughput_memory.svg). Throughput uses all input pairs;
memory is the largest native stage RSS. Construction is excluded from both.

| Workflow | Mapping-time multiplier, 500M / 100M | Counting-time multiplier | Combined multiplier |
|---|---:|---:|---:|
'''
for tool,label in [('ours','Constellation'),('piscem','piscem + alevin-fry'),('kallisto','kallisto + bustools')]:
 r=s['scaling'][tool];text+=f"| {label} | {r['map_seconds']:.2f}× | {r['count_seconds']:.2f}× | {r['pipeline_seconds']:.2f}× |\n"
a,b=[s['cases'][f'{n}m-ours'] for n in (100,500)]
text+=f'''
Constellation mapping throughput is {a['million_pairs_per_second']:.3f}M pairs/s
at 100M and {b['million_pairs_per_second']:.3f}M pairs/s at 500M. Mapping RSS stays
near {b['stages']['map']['rss_gib']:.2f} GiB. Counting RSS grows from
{a['stages']['count']['rss_gib']:.2f} to {b['stages']['count']['rss_gib']:.2f} GiB;
the counter retains evidence across the dataset and does not have constant
memory with input size. No process reached the 46 GiB guard.

The 500M default count spends {b['count_metrics']['phase_seconds']['ingest']:.2f}s
ingesting/partitioning evidence, {b['count_metrics']['phase_seconds']['resolve_and_collapse']:.2f}s
resolving/correcting molecules, {b['count_metrics']['phase_seconds']['gather']:.2f}s
gathering/sorting counts and {b['count_metrics']['phase_seconds']['write']:.2f}s
writing. Input ingestion is now the main counter bottleneck. The current
[`count-splicing` implementation](../../crates/constellation-cli/src/cmd_count/splicing.rs)
parses, corrects barcodes and partitions evidence in one serial loop before
parallel cell resolution. Parallel ingestion or a fused mapping/counting path
that avoids the CSTRAD2 write/read round trip are the clearest next targets.
Controlling retained evidence memory is a separate requirement; fusion alone
does not make that memory bounded. These changes are not implemented here.

Our mapping CPU cost per pair is nearly unchanged between these two runs, while
wall throughput is lower at 500M. Input bytes, filesystem I/O, cache conditions
and the selected read distribution differ. The logs retain CPU and I/O counters;
these single runs do not isolate a unique cause for the wall-time difference.

Alevin-fry's collation grows from 1.12s at 100M to 47.89s at 500M. Its own
phase log localizes most of the increase to gather: 0.47s versus 41.33s;
scatter takes 0.49s versus 6.07s. Both runs use the optimized collator, eleven
spool files and ten gather workers, so this is not evidence of a switch from
in-memory to external collation. The 500M stage records substantial filesystem
I/O and less than one core of average CPU utilization. Cache, spool access and
output waiting are plausible contributors; this experiment does not separate
them. The default 2 GiB collation budget is retained. Inspection of the pinned
source shows both record-count planning and buffer tuning capped at 2 GiB when
using `--memory-limit`, so simply raising that option is not a demonstrated
remedy. See [phase/source diagnosis](raw/alevin-collation-diagnostic.json).

Bustools' state-count metadata is not a universal constant: its `matrix.ec`
grew from 18.9 MB in the prior 10M run to {ec_bytes[100]/1e6:.1f} MB at 100M and
{ec_bytes[500]/1e6:.1f} MB at 500M. The dictionary changes
as new compatibility sets are observed. Thus the old 17s empty-input setup
measurement cannot simply be subtracted at every input size. Its 100M barcode
correction also used only 4.26 CPU-seconds during 12.01 wall seconds, showing
substantial waiting in that stage; wall figures include output I/O/system effects.

## Molecule outputs and agreement

Forward-compatible mapping assignments, before whitelist filtering and molecule
resolution (denominator: all input pairs):

| Pairs | Constellation | Piscem | Kallisto |
|---|---:|---:|---:|
'''
for size in (100,500):
 rates=[s['cases'][f'{size}m-{tool}']['forward_assignment_fraction']*100 for tool in ('ours','piscem','kallisto')]
 text+=f"| {size}M | {rates[0]:.3f}% | {rates[1]:.3f}% | {rates[2]:.3f}% |\n"
text+='''
Piscem's orientation-matched figure comes from permit-list correction statistics;
its raw mapper total also includes reverse-orientation mappings. Mapping
assignment is not the same as a retained, correctly assigned molecule.

| Pairs | Evidence / counter | Spliced | Unspliced | Ambiguous | Total molecules |
|---|---|---:|---:|---:|---:|
'''
for size in (100,500):
 for tool,mode,title in [('ours','count','Ours, conservative + UMI correction'),('ours','usa','Ours, exact-UMI USA'),('piscem','count','piscem + alevin-fry'),('kallisto','count','kallisto + bustools')]:
  m=s['cases'][f'{size}m-{tool}']['matrix_totals'][mode]
  text+=f"| {size}M | {title} | {m['S']['total']:,} | {m['U']['total']:,} | {m['A']['total']:,} | {m['total']:,} |\n"
text+='''
These totals are not an accuracy ranking. The exact-UMI intersection baseline
is available in our `molecules_before_umi_correction` metric, which separates
the USA resolution change from disabling UMI correction:

| Pairs | Intersection before UMI correction | USA exact UMIs | Policy-only net change |
|---|---:|---:|---:|
'''
for size in (100,500):
 c=s['cases'][f'{size}m-ours'];before=c['count_metrics']['molecules_before_umi_correction'];after=c['usa_metrics']['molecules_after_umi_correction']
 text+=f"| {size}M | {before:,} | {after:,} | {after-before:+,} ({(after/before-1)*100:+.2f}%) |\n"
text+='''
Matrix agreement with piscem/alevin-fry, including implicit zeros on the common
58,395-gene × 11,769-cell axes. S/U/A columns are Pearson correlations. Weighted
overlap is the sum of coordinate-wise minimum counts divided by the sum of
maximum counts across all three states (weighted Jaccard):

| Pairs | Comparison | S | U | A | Weighted overlap |
|---|---|---:|---:|---:|---:|
'''
for size in (100,500):
 for name,title in [('ours-count','Ours conservative'),('ours-usa','Ours USA'),('kallisto-count','Kallisto/bustools')]:
  m=comparisons[f'{size}m-{name}_vs_{size}m-piscem-count']
  total=sum(v['a_total']+v['b_total'] for v in m.values());l1=sum(v['l1'] for v in m.values());overlap=(total-l1)/(total+l1)
  text+=f"| {size}M | {title} | {m['S']['full_matrix_pearson']:.4f} | {m['U']['full_matrix_pearson']:.4f} | {m['A']['full_matrix_pearson']:.4f} | {overlap:.4f} |\n"
text+='''
L1 discrepancies, weighted overlaps and gene-total correlations are retained
alongside Pearson values. Agreement does not establish biological truth, and
these real-read tests do not replace the previous independent source-origin
accuracy tests.

[All timing/resource results](raw/summary.json),
[matrix comparisons](raw/matrix-comparison.json),
[input and binary provenance](raw/provenance.json), and
[paired-input validation](raw/input-validation.json).

'''
p=RAW.parent/'README.md';body=p.read_text()
if '<!-- RESULTS -->' in body:
 body=body.replace('<!-- RESULTS -->',text)
else:
 start=body.index('## Measured results');end=body.index('\n## Protocol',start)
 body=body[:start]+text+body[end:]
p.write_text(body)
update_root_readme()
for p in RAW.parent.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
(RAW/'completion.json').write_text(json.dumps({'all_six_native_workflows_complete':True,'matrix_comparison_sets':len(comparisons),'exact_input_pairs':[100000000,500000000],'production_code_changed':False},indent=2)+'\n')
print('Final report written.')
