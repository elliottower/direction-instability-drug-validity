# Superseded implementation

**Do not run. Retained because the implementation manifest pins it.**

`build_h3_bundle.py` built the H3 cohort bundle from `lincs_shrna.npz`, chosen over
the GCTX rebuild on the grounds that the rebuild "reproduced their values exactly
but not their coordinate order". Deviation 11 establishes that the ordering was the
other way round: the extraction's columns do not follow the identifiers it
declares, and the rebuild's do. Amendment 2, frozen at `fd1ae8d`, retires both the
extraction and this script in A9.

Two further defects, neither of which affected a registered artifact because no
registered command ran this file:

- Line 48 asserts that the two extractions share a gene axis and passes on the
  defective file, because both declare the same identifiers. Amendment 2's A8
  registers a second layer for exactly this reason.
- The eligibility floor counted metadata rows rather than distinct signature ids,
  so a duplicated row weighted one signature twice. The three scripts that do run
  count distinct ids.

The bundle the analyses use comes from `experiments/modal_h3_rebuild.py`.
