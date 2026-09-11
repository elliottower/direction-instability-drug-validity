# Superseded: outputs from the defective matched-contrast code

**Do not use for scientific reporting.**

`matched_well_contrast.csv` (2026-09-08) came from a version of
`experiments/08c_matched_well_contrast.py` whose plate-to-source mapping used
`np.searchsorted` without mapping back through the sorter. On unsorted plate
codes that returns positions in the sorted view rather than real rows, which
invalidated the permutation null. The defect was found by an adversarial test
on deliberately unordered codes; a toy test had passed only because its codes
happened to be pre-sorted.

`within_vs_across_source.csv` (2026-09-08) came from `08b`, which compared
within-source raw wells against across-source consensuses. Those are different
estimands at unequal aggregation depth, so the comparison supports a direction
but not a magnitude. It was superseded by the matched raw-well design.

Neither has been regenerated with corrected code. Retained for provenance only.
