# Superseded: mis-indexed shRNA consensus directions

**Do not use for scientific analysis.** Retained for provenance.

`localization_results.json` in this directory was computed with target consensus directions indexed
by metadata row position rather than by signature id. The two source files carry
the same 154,993 signature ids in different order and agree at 115 positions, so
each target's direction was the mean of metadata-position-selected unrelated
shRNA signatures. Deviation 9 in `DEVIATION_LOG.md` records the mechanism, the
evidence and the scope.

| | sha256 |
|---|---|
| superseded, this file | `36808886a93ab4b801225c26e1699f693686a2318ed80c5cd185989f5cd3acd4` |
| corrected, one directory up | `9c7ada30ad633486c92ee844c48203c7ae2c133a4a6419be5e399b8b91b90c04` |

Raw direction instability never uses a target direction and is unaffected in
both versions.
