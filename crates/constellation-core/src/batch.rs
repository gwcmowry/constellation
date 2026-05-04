use crate::sketch::SketchRecord;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SketchBucket {
    pub start: usize,
    pub end: usize,
}

pub fn make_sketch_buckets(records: &[SketchRecord], target_size: usize) -> Vec<SketchBucket> {
    let target_size = target_size.max(1);
    (0..records.len())
        .step_by(target_size)
        .map(|start| SketchBucket {
            start,
            end: (start + target_size).min(records.len()),
        })
        .collect()
}
