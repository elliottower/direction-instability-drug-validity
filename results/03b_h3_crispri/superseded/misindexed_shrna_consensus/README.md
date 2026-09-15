# Superseded: mis-indexed shRNA consensus directions

**Do not use for scientific analysis.** Retained for provenance.

`h3_crispri_results.json` in this directory was computed with target consensus directions indexed
by metadata row position rather than by signature id. The two source files carry
the same 154,993 signature ids in different order and agree at 115 positions, so
each target's direction was the mean of metadata-position-selected unrelated
shRNA signatures. Deviation 9 in `DEVIATION_LOG.md` records the mechanism, the
evidence and the scope.

| | sha256 |
|---|---|
| superseded, this file | `1e96ef13ad1dd1629b2b4f9ebfd90087bde3c2f444f184c9d6071e5be2ea68bd` |
| corrected, one directory up | not yet produced |

The corrected version requires the Replogle K562 h5ad, which is not held locally or on any Modal volume checked; the file at the canonical path is still the superseded one until that rerun happens.

Raw direction instability never uses a target direction and is unaffected in
both versions.
