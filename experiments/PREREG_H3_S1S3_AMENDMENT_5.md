# Amendment 5: the two comparisons one tolerance was carrying

**Date:** 2026-10-07
**Status:** FROZEN. The R0-R7 execution of 2026-10-07 refused at this gate and is
void and preserved. No R1-R7 association, resampling result or criterion outcome
has been computed under any construction. 279 tests pass under `python` and under
`python -O`; the logs are pinned below.
**Commit SHA:** PENDING
**Amends:** `experiments/PREREG_H3_REFERENCE_DISCORDANCE.md` (frozen `7f57136`),
whose reproduction tolerance this splits, and
`experiments/PREREG_H3_S1S3_AMENDMENT_3.md` (frozen `03175d0`), whose B9 it
narrows. `experiments/PREREG_H3_S1S3_AMENDMENT_4.md` (frozen `8edf227`) is
untouched.
**Discovery record:** Deviation 14 in `DEVIATION_LOG.md`.

## What refused

The registered driver ran once on 2026-10-07 and refused at its own reproduction
gate before computing any association:

```
P_shrna reproduces to 1.14e-06, tolerance 1e-06
```

That execution is void under B9 and is preserved as a failed execution. It is not
reclassified, and nothing from it is reported.

One drug of 795 crossed the flat tolerance: levofloxacin, whose deposited
`proj_shrna` is 8.6574612530 and whose error is 1.1385e-06, 2.21 float32 ULPs at
that magnitude.

## What the diagnosis established, before anything was changed

Two explanations were proposed and both were refuted by measurement. The retired
loader takes the mean of a float32 array, which numpy accumulates in float32; the
float32 and float64 accumulations of the same signatures differ from each other by
about 1e-9 and both sit about 1e-6 from the deposited values. A float32 projection
path was the second; no dtype route reproduces the deposited values below a 1e-7
relative floor.

What holds instead: context membership is identical between routes for every drug
examined, and among the nine drugs the localization diagnostic examined the
difference is spread rather than concentrated in one context pair, with 6 of
levofloxacin's 21 pairs above half the maximum per-pair difference of 5.658e-07.
**The diagnostics found no localized context-pair discrepancy among the examined
focus drugs and rejected the two tested dtype explanations; the remaining residual
has not been explained.**

A third candidate was named and tested under a plan frozen before it ran,
`experiments/PREREG_H3_S1S3_DIAGNOSTIC_D9_INDEXING.md`: the pre-fix shRNA consensus
indexing of Deviation 9. H2 holds. Over all 795 drugs the pre-fix reference differs
from the deposited records by up to 7.818 in `P_shrna` and 0.448 in `E_shrna`, against
9.755e-07 and 9.843e-08 for the corrected reference, and no drug of the cohort is
closer under it.

That outcome was fixed by provenance rather than discovered by the measurement, and
the registration says so plainly: the deposited records this gate compares against
are the corrected artifact, recomputed with the corrected function at `f822fb1` on
2026-09-15, so a defect whose fix produced them cannot account for a residual measured
against them. The 7.818 and 0.448 are the same two figures Deviation 9 reports for the
corrected reconstruction against the superseded artifact, measured in the other
direction. The diagnostic is recorded for what it does establish, that the deposited
reference directions are the corrected consensus, and no further candidate for the
residual is named.

**Route to route, over the whole cohort.** The two reconstructable routes are closer
to each other than either is to the deposited values in the cohort's maxima and
medians, by a factor of 2.5 for `P_shrna`, 4.5 to 5.4 for `E_shrna` and 31 to 38 for
`D`. Per drug it does not hold universally: the routes are mutually closer for 568 of
795 drugs on `P_shrna`, 611 on `E_shrna` and 778 on `D`. The nine-drug localization
diagnostic had suggested five to ten times for every drug examined, and the cohort
does not carry either part of that.

## The two comparisons

The registration's flat absolute tolerance describes two float64 computations from
one extraction. Amendment 3's B1 then made a different artifact authoritative, and
the gate went on applying one rule to two comparisons that are not the same
comparison.

**D1. Legacy reproduction.** The retired extraction, read on its own declared axis,
and the retired aggregation against the deposited records, under the frozen flat
absolute `RECON_TOL` of 1e-06. This is the comparison the registration describes and
the tolerance is unchanged for it.

**D2. Production equivalence.** The authoritative GCTX-derived route against the same
deposited records, under Amendment 2's elementwise numerical-equivalence rule,
`|a-b| <= atol + rtol|b|` with `atol = 1e-8` and `rtol = 1e-6`, for every drug and
every invariant quantity. The compared statistics are float64; the rule supplies the
error model for values derived from float32 source data, and is not a claim that
either side is a float32 number.

Both are measured in
`results/03d_h3_reference_discordance/two_layer_gate_measurement.json`, over all 795
drugs:

| quantity | legacy max error | legacy residual | production max error | production residual | production over flat 1e-06 |
|---|---|---|---|---|---|
| P_shrna | 9.755e-07 | 0.122 | 1.139e-06 | 0.195 | 1 |
| E_shrna | 9.843e-08 | 0.446 | 9.365e-08 | 0.421 | 0 |
| D | 7.298e-08 | 0.128 | 7.496e-08 | 0.129 | 0 |

The legacy route reproduces the deposited analysis to the frozen flat tolerance on
every drug and every quantity, nothing over. The authoritative route satisfies the
elementwise rule on every drug at a fifth of its allowance for P, and fails the flat
rule for one. **The frozen tolerance is therefore kept, not loosened**: it governs
the comparison it was registered for, which passes.

**D3. Where the legacy layer is established, and how it is trusted.** Once, outside
production, in a sealed audit this amendment pins. It is not re-run inside a
production invocation: reading the retired extraction there would reintroduce the
known-bad axis declaration that Amendment 3's B2 put a check in the loader to keep
out.

A record consumed on its own word is not evidence, so the audit is a trust boundary
and is treated as one. Its bytes are hashed and compared with the expected digest in
`registry/frozen/reproduction_legacy_pin.json` before they are parsed, by the same
rule Amendment 4 applies to the analysis bases: an artifact carrying its own expected
digest attests to itself. The audit must then declare the schema and the rule the pin
names rather than any of its own, must have been measured on the inputs the pin names,
and must carry the digest and row count of the complete per-drug comparison behind its
summaries. Its predicate, its count of drugs over the tolerance, and its list of those
drugs are three statements of one fact, and a record that contradicts itself is refused
whichever of them is read. If the retired route ever fails its own frozen tolerance,
that is a finding and is not amended away; the gate says so in its own refusal message.

**D3a. What the production gate checks about the cohort.** The deposited records must
carry no duplicated drug identifier, because keying by drug silently keeps the last of
a pair and compares against it. Every drug of the arm must have a deposited record, as
a stated refusal rather than an incidental key error. The arm must hold exactly the
cohort registered for that comparison, so a deposited record that went uncompared is a
refusal rather than a silence. The two CRISPRi fields are both present or both absent,
because a measured difference must not disappear because one field of a pair was
serialized.

**D4. Every element, with no way round it.** No aggregate satisfies the rule for an
element, no drug is excepted, no drug-specific tolerance exists, and no drug leaves
the cohort. Levofloxacin in particular receives no special treatment and is not
removed: it is neither the largest `P_shrna` in the cohort, 8.66 against 16.95, nor
the largest error in ULPs, 2.21 against 3.29 for pefloxacin. It is the one drug whose
magnitude and numerical path put it across a flat threshold.

**D5. What the `D` check cannot establish.** `D` is one minus the mean pairwise
cosine, and a cosine does not change when either argument is rescaled. Multiplying
any signature by any positive scalar therefore leaves `D` exactly unchanged, so the
reproduction gate on `D` cannot detect a scaling error of any size in any signature.
This is a property of the measure, not of this data, and the registration did not
state it.

No other derived cosine statistic is given that job. Amplitude is established
upstream, where it can be: the pinned source checksums, the elementwise comparison of
every reconstructed raw signature against the pinned GCTX, which gate v2 holds at a
maximum absolute difference of 0.0 over 42,498 signatures, the exact membership and
context counts, and the per-signature finite and nonzero-norm checks. The reproduction
gate then establishes that verified inputs flow through the intended aggregation and
projection, which is what it can establish and all of it.

**D6. What the gate records.** Per drug and per quantity: the computed value, the
deposited value, the absolute and relative error, the allowance, and the normalized
residual `|a-b| / (atol + rtol|b|)`. Per quantity: the maxima of each, the worst drug
by residual and by absolute error, the counts over each rule, and the complete list
of drugs over each. A non-finite computed or deposited value is refused by name
before any comparison, so it can neither dominate a maximum silently nor be skipped.
Every drug of this cohort is registered as carrying both shRNA fields, so a record
missing one is refused as a cohort defect rather than passed over, which is what the
previous implementation did.

**D7. Freeze discipline.** `tests/test_reproduction_gate_layers.py` requires that a
difference inside the elementwise allowance passes and one over it refuses; that 39
drugs reproducing exactly cannot carry a fortieth that does not; that a failing
legacy layer refuses before the production layer is reached; that a record missing a
shRNA field is a cohort defect; that a non-finite value never reaches the comparison;
that every drug and quantity is reported with its allowance; and that the measured
legacy layer holds on the real 795-drug cohort while the production layer is the one
the flat rule rejects, for levofloxacin alone.

On the trust boundary D3 describes, it requires refusal when the audit is absent,
when no pin is supplied, when the audit's hash is not the pinned one, when one byte
of it is altered, when a replacement satisfying every predicate is substituted, when
its drug count is not the pinned one, when its predicate contradicts its count or its
list of failing drugs, when its reported maximum does not clear the tolerance it
claims to meet, when it was measured on other inputs, when the complete per-drug table
behind it is not the pinned one, and when it declares a rule of its own. On the cohort,
it requires refusal for a duplicated deposited drug, for a deposited record the arm
never compares, for an arm drug with no deposited record, and for either CRISPRi field
present without its pair. The flat rule is strict at exactly the tolerance, in the
per-drug row and in the audit alike. The suite runs under `python` and `python -O`.

## Pins

The sealed legacy audit is pinned from outside itself, in
`registry/frozen/reproduction_legacy_pin.json`, which the gate reads to authenticate
it. Its digest is restated here for the record and is not read from this file.

| artifact | sha256 |
|---|---|
| `results/03d_h3_reference_discordance/legacy_reproduction_audit.json` | `0810be5b88bbe429ca313c4d3d2b0e0072e05feaf5f90ef15a8654fd9562d693` |
| `results/03d_h3_reference_discordance/legacy_reproduction_per_drug.json` | `4194f88c10291cc81500196b09bd2462922c31d6feaf14a7962f9afcd6529d0a` |
| `results/03d_h3_reference_discordance/route_to_route_all_drugs.json` | `b3b97cc7ed9decb8e6bf9ec6b4073913edb504313690c3031226c3416526c226` |
| `results/03d_h3_reference_discordance/two_layer_gate_measurement.json` | `27608a93e848f30c324f89d247ea4451c9b02eb4e71bb5b01207757124306bb5` |
| `results/03d_h3_reference_discordance/reproduction_error_localized.json` | `c46878d529084526d2607a64519adcedddb205949a704be438e65651a46a5663` |
| `results/03d_h3_reference_discordance/reproduction_tolerance_diagnosis.json` | `d0c4615e015e8ed5663b50821455d498dfd1a5a50ff3d8a27dbad166ff832105` |
| `results/03d_h3_reference_discordance/reproduction_tolerance_power.json` | `e4867903ee43de6df2fce4e592837f14a8474c1d0f2a4ec948feff2c3b0cd883` |
| `results/03d_h3_reference_discordance/d9_indexing_diagnostic.json` | `6e81efac52a00bc7edabf792ae797bc5fb9d28e6da9df0aee962e8eaaad513a1` |
| `registry/frozen/reproduction_legacy_pin.json` | `d5a443fe76b8d1e03c00879a8358d2e2fe75885fad1d3f185f4b2820f7b30111` |
| `experiments/03d_h3_reference_discordance.py` | `caf3a1da1e945d7a87dd6c1b0388ea149d7578e481d22f944536ab362bf9e7db` |
| `experiments/modal_h3_execute.py` | `3bf0050c09f76cf77722164bf58015f49fde315f57455905e47e9f2bb68032e6` |
| `tests/test_reproduction_gate_layers.py` | `9f722a31f20a15eab9ee231382940d0040954daddab9923098659f933494a8c8` |
| `results/03d_h3_reference_discordance/test_log_at_amendment_5.txt` | `956357f53307375d16205bba343565eb7074995e719e10e7272f19a5e659f276` |
| `results/03d_h3_reference_discordance/test_log_at_amendment_5_optimized.txt` | `37e68e721105fd1186d83f115f59dcdde7bb6862b110499f7fb1fbdffbd83750` |

The audit records the identity of the code that measured it, and the pin carries that
record: `analysis_code_sha256` is this version of `03d` and `measurement_code_sha256`
the stage that ran. The gate compares the audit against the pin, not against whatever
is on disk, so later work on the analysis code does not silently invalidate the
quarantined evidence or quietly re-authorize it.

## Foreknowledge at the time of this amendment

The refusal and the full numerical diagnosis were known before this was written. In
particular: which drug crossed the flat tolerance and by how much, the per-drug error
distribution for all 795 drugs on all three invariant quantities, that the elementwise
rule holds on all of them, and the normalized residuals under it.

`atol = 1e-8` and `rtol = 1e-6` are Amendment 2's, frozen at `fd1ae8d` before this
refusal, and are not fitted to place the observed error inside a boundary. **Their
application to this gate was chosen after the flat rule refused.** That is the
foreknowledge this amendment carries and it is not reduced by the constants
predating it.

Three injected-error ladders are recorded with the diagnosis. They are diagnostic
illustrations of what each rule refuses and were not used to choose `atol` or `rtol`.
Two of the three cannot move `D` at all, for the reason D5 gives, and the first
cannot move `E`; that limit is stated in the artifact rather than left for a reader
to find.

R0.4 and R0.6 have been computed, under Amendment 4 and on its frozen bases. No
R0.5, R1 to R7 association, bootstrap, permutation or criterion outcome exists.

The D9 diagnostic and the cohort route-to-route comparison were both run before this
text was finalized, and both outcomes are reported above whether or not they
supported what the draft said. The route-to-route claim in the earlier draft did not
survive them.

## Maximum claim under this registration

Confirmatorily, the paper may say: in a sealed and externally pinned legacy audit, the
deposited shRNA quantities reproduce from the retired extraction and the retired
aggregation with absolute error below 1e-06 for every drug and every invariant
quantity; and under the authoritative source-derived route, all 795 drug-level
comparisons for `P_shrna`, `E_shrna` and `D` satisfy Amendment 2's elementwise rule
`|a-b| <= 1e-8 + 1e-6|b|`, at maximum normalized residuals of 0.1953, 0.4211 and
0.1288 respectively.

On the two routes: the paper may say that over the registered cohort the two
reconstructable routes agree with each other more closely than either agrees with the
deposited values, in the maximum and the median of every invariant quantity, by a
factor of 2.5 for `P_shrna`, 4.5 to 5.4 for `E_shrna` and 31 to 38 for `D`. It may not
say this holds for every drug, because it holds for 568, 611 and 778 of 795.

It may not say that the authoritative route reproduces the deposited values to 1e-06
absolute, because for one drug it does not. It may not say the residual difference
between the deposited values and either modern route has been explained; it has been
bounded, shown among the nine drugs examined not to concentrate in one context pair,
shown not to arise from the aggregation dtype, and shown not to arise from the
pre-fix consensus indexing, and no further account of it is offered. It may not say that the
`D` reproduction check establishes anything about signature amplitudes.
