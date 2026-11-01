"""Check completion, exact input sizes and retained artifacts after the suite."""
import ast
import hashlib
import json
from pathlib import Path

import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
RAW = Path(__file__).parent / 'raw'
DATA = ROOT / 'benchdata/large-scale-2026-11-01'
summary = json.loads((RAW / 'summary.json').read_text())
assert len(summary['cases']) == 6
runs = [p for p in RAW.glob('*.run.json') if p.name.startswith(('100m-', '500m-'))]
assert len(runs) == 22
for path in runs:
    result = json.loads(path.read_text())
    assert result['exit_code'] == 0 and not result['stopped_by_memory_guard'], path
for name, case in summary['cases'].items():
    metrics = case['mapping_metrics']
    assert metrics.get('num_reads', metrics.get('n_processed')) == case['pairs'], name
    assert max(stage['rss_gib'] for stage in case['stages'].values()) < 46

# The same directory also retains nine earlier normalizer-validation fixtures.
files = [p for p in (DATA / 'matrices').glob('*.parquet') if p.name.startswith(('100m-', '500m-'))]
assert len(files) == 24
for path in files:
    parquet = pq.ParquetFile(path)
    assert parquet.schema_arrow.names == ['gene_id', 'cell_id', 'count']
    assert str(parquet.schema_arrow.field('count').type) == 'uint64'
    label, state = path.stem.rsplit('-', 1)
    expected = json.loads((RAW / f'{label}.matrix-summary.json').read_text())[state]['nonzero']
    assert parquet.metadata.num_rows == expected, path

provenance = json.loads((RAW / 'provenance.json').read_text())
digest = hashlib.sha256((DATA / 'constellation.frozen').read_bytes()).hexdigest()
assert digest == provenance['binaries']['ours']['sha256']
comparisons = json.loads((RAW / 'matrix-comparison.json').read_text())
assert len(comparisons) == 6
for comparison in comparisons.values():
    assert set(comparison) == set('SUA')
    for values in comparison.values():
        assert -1 <= values['full_matrix_pearson'] <= 1
        assert 0 <= values['weighted_jaccard'] <= 1
for path in RAW.parent.glob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
assert not any((DATA / 'work').iterdir())
audit = {
    'native_stage_exit_codes_all_zero': True,
    'native_stages': len(runs),
    'exact_mapping_input_counts_verified': 6,
    'validated_parquet_state_matrices': len(files),
    'binary_snapshot_checksum_matches': True,
    'all_peaks_below_guard': True,
    'no_pipeline_scratch_remaining': True,
    'matrix_comparison_sets': len(comparisons),
}
(RAW / 'final-audit.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps(audit, indent=2))
