//! Compatibility with spliced and intronic reference classes, not isoform identity.

/// Stable on-disk codes. Unknown is uninformative; Conflict is sticky.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum SplicingState {
    Spliced = 1,
    Unspliced = 2,
    Both = 3,
    Conflict = 4,
    Unknown = 8,
}

impl SplicingState {
    pub fn from_code(code: u8) -> Option<Self> {
        Some(match code {
            1 => Self::Spliced,
            2 => Self::Unspliced,
            3 => Self::Both,
            4 => Self::Conflict,
            8 => Self::Unknown,
            _ => return None,
        })
    }
    pub(super) fn from_class_bits(bits: u8) -> Self {
        Self::from_code(bits & 3).unwrap_or(Self::Unknown)
    }
    /// Intersect constraints from seeds/reads belonging to the same gene/molecule.
    pub fn intersect(self, other: Self) -> Self {
        use SplicingState::*;
        match (self, other) {
            (Conflict, _) | (_, Conflict) => Conflict,
            (Unknown, x) | (x, Unknown) => x,
            (a, b) => Self::from_code((a as u8) & (b as u8)).unwrap_or(Conflict),
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::Spliced => "S",
            Self::Unspliced => "U",
            Self::Both => "A",
            Self::Conflict => "C",
            Self::Unknown => "X",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::SplicingState::*;
    #[test]
    fn compatibility_is_order_independent_and_distinguishes_conflict() {
        let states = [Spliced, Unspliced, Both, Conflict, Unknown];
        for a in states {
            for b in states {
                for c in states {
                    assert_eq!(a.intersect(b), b.intersect(a));
                    assert_eq!(a.intersect(b).intersect(c), a.intersect(b.intersect(c)));
                }
            }
        }
        assert_eq!(Both.intersect(Spliced), Spliced);
        assert_eq!(Both.intersect(Unspliced), Unspliced);
        assert_eq!(Spliced.intersect(Unspliced), Conflict);
        assert_eq!(Conflict.intersect(Both), Conflict);
        assert_eq!(Unknown.intersect(Unspliced), Unspliced);
    }
}
