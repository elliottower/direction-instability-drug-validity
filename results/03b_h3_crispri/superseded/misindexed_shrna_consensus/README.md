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
| corrected, one directory up | `f00f8428071178bd2317139c5a1f534e458931684fb6da65c750f9f1e5456eb3` |

The corrected version was produced in commit `f822fb1` from the scPerturb
`ReplogleWeissman2022_K562_essential.h5ad` (Zenodo 10044268, sha256
`412fd0df8c4ccea9f4db91cd88033c49200838b29d40945e48574be588b48789`). Only the
shRNA fields differ between the two versions; the CRISPRi directions never used
the defective path.

Raw direction instability never uses a target direction and is unaffected in
both versions.
