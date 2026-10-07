# Does the library-version difference between the two environments account for the residual?

**Date:** 2026-10-07
**Status:** DRAFT
**Commit SHA:** PENDING
**Kind:** forensic diagnostic. No registered statistic is computed and no
confirmatory claim rests on it. It asks how the deposited numbers were produced.
**Relates to:** Deviation 14, which records the residual;
`experiments/PREREG_H3_S1S3_AMENDMENT_5.md` (frozen `04fba74`), whose gate rule
this cannot change; `experiments/PREREG_H3_S1S3_DIAGNOSTIC_D9_INDEXING.md` (frozen
`e684d16`), the candidate this follows.

## The question

**Does running the legacy route under the library versions that produced the
deposited records remove the residual difference from those records?**

Categorical, over all 795 drugs and all three invariant quantities.

## Foreknowledge

The residual is known in full. The authoritative route differs from the deposited
`proj_shrna` by at most 1.139e-06 and the retired route by at most 9.755e-07, both
below a 1e-7 relative floor that three candidates have failed to explain: the
float32 mean accumulator, a float32 projection path, and the pre-fix shRNA consensus
indexing of Deviation 9.

The version difference was found before this was written and is why it is written.
`f822fb1` of 2026-09-15 produced `results/03b_h3_crispri/h3_crispri_results.json` by
running `experiments/03b_h3_crispri_ground_truth.py`, and no Modal wrapper runs that
script: it ran locally, under the `uv.lock` committed at that revision, which pins
numpy 2.4.6, pandas 2.3.3, scipy 1.18.0 and anndata 0.12.19. The reconstruction ran
in the container, whose image pins numpy 2.1.3, pandas 2.2.3, scipy 1.14.1 and
anndata 0.11.4, and the audit records those two versions in its own provenance.
`uv.lock` is byte-identical at `f822fb1` and at `04fba74`, so the local resolution has
not moved since.

numpy's pairwise summation blocks its accumulation by a size threshold, and its dot
products dispatch to whichever BLAS the wheel carries. Both changed between 2.1 and
2.4. A difference of this order is what that would produce.

## Hypotheses

**H1.** Under the library versions `uv.lock` pins at `f822fb1`, the legacy route
reproduces the deposited records more closely than it does under the versions the
audit was measured with.

**H2.** It does not, and the residual is not a library-version effect.

H1 carries the diagnostic. Neither outcome changes which artifact is authoritative.

## Implementation, fixed before it runs

**One image, one run, and no search.** A single Modal image pins numpy 2.4.6, pandas
2.3.3, scipy 1.18.0 and anndata 0.12.19, the versions `uv.lock` resolves at
`f822fb1`. The legacy aggregation and the legacy reference are the same code the
sealed audit ran, on the same pinned inputs, by the digests the audit records.

**No other version combination is tried under this plan.** If H1 does not hold, that
is the answer; a second image with different pins would be fitting an environment to
an observed difference, and this plan does not authorize one. Any further version
work needs its own frozen plan saying what it is for.

The interpreter version is not recoverable: `uv.lock` records `requires-python
>= 3.10` and does not pin an interpreter, and nothing records which Python ran on
2026-09-15. The image keeps the Python the audit used, so the interpreter is held
constant and only the libraries move. A residual that survives under matched
libraries is therefore not shown to be free of an interpreter effect.

`P_shrna`, `E_shrna` and `D` are computed for all 795 drugs and compared against the
deposited records. Every per-drug error is written, not only the maxima.

**Where the stage lives.** `stage_environment_diagnostic` is added to
`experiments/modal_h3_execute.py`, which Amendment 5 pins at
`3bf0050c09f76cf77722164bf58015f49fde315f57455905e47e9f2bb68032e6`. That digest
describes the file as it stood at the freeze and not afterwards, and the file now
hashes to something else. The alternative was a second copy of the retired
aggregation in a new file, and two copies of the code under test can drift, which is
worse than a wrapper whose change is written down. Amendment 5's own evidence does not
rest on this file: the sealed audit is pinned outside itself in
`registry/frozen/reproduction_legacy_pin.json` and the gate enforces that pin on every
run.

## Inference criteria

| outcome | holds when |
|---|---|
| H1 holds | for all three quantities, the maximum absolute error against the deposited records falls below 1e-08 |
| H1 holds in part | the maximum absolute error falls by at least a factor of ten for at least one quantity, without reaching 1e-08 |
| H2 holds | no quantity's maximum absolute error falls by as much as a factor of two |

**The residual is described as a library-version effect only under H1.** A fall at
one quantity, or a maximum that merely looks smaller, decides nothing.

## What this cannot establish

Matching the recorded library versions does not reconstruct the machine: the BLAS a
wheel carries depends on the platform it was built for, and the deposited records were
produced on an Intel Mac rather than in a Debian container. H1 would show that the
libraries account for the difference; H2 would not rule out an environment effect of
another kind, and would leave the residual where Amendment 5 leaves it.

## Disposition

The result is preserved either way and recorded in Deviation 14. Amendment 5 is
frozen and does not change on the outcome: the deposited records remain the
comparison target, the flat tolerance remains on the legacy comparison, and the
elementwise rule remains on the production comparison. Nothing here regenerates the
deposited records.

## Maximum claim under this diagnostic

Confirmatorily, nothing. The paper may report, as a statement about provenance rather
than about H3, whether the library versions recorded for the two environments account
for the difference between the deposited values and the reconstructable routes, at the
precision stated, over all 795 drugs and all three invariant quantities.
