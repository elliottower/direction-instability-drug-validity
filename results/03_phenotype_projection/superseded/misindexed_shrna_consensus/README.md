# Superseded: mis-indexed shRNA consensus directions

**Do not use for scientific analysis.** Retained for provenance.

`phenotype_projection_results.json` in this directory was computed with target consensus directions indexed
by metadata row position rather than by signature id. The two source files carry
the same 154,993 signature ids in different order and agree at 115 positions, so
each target's direction was the mean of metadata-position-selected unrelated
shRNA signatures. Deviation 9 in `DEVIATION_LOG.md` records the mechanism, the
evidence and the scope.

| | sha256 |
|---|---|
| superseded, this file | `65e5d10e272037987384f89e6208de478fa17f6b8fb946add893e5a24c2d80c4` |
| corrected, one directory up | `fd69e26fc9a3917323065b631688baeab8b283f735c8bf5b16210ba67bd21425` |

Raw direction instability never uses a target direction and is unaffected in
both versions.
