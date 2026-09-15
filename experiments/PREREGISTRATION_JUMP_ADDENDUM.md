# Does cross-source rank preservation depend on how many sources a compound was profiled in?

**Registered:** (commit SHA recorded on freeze; pushed before any outcome is computed)
**Author:** Elliot Tower
**Amends:** `PREREGISTRATION.md` and `PREREG_ROBUSTNESS_CHECKS.md` (R6), commit `249abaf`
**Scope:** JUMP-CP Cell Painting only. Adds one confirmatory hypothesis, freezes a
cohort-exclusion rule, and reclassifies analyses whose outcomes are already known.

## Foreknowledge of data or evidence

Known before writing this: on a cohort that wrongly retained eight plate-level
positive controls, a 252-feature morphological ablation gave Spearman
$\rho = 0.9977$ over 25,000 compounds; its size-matched random null placed the
observed value below all 200 draws ($p_{MC} = 1/201$); leave-one-source-out
correlations ran 0.9975–0.9977; the median absolute percentile-rank shift was
1.04%. A matched raw-well contrast produced a positive median difference under
code later shown to contain an indexing defect that invalidated its permutation
null.

Not computed, and not inspected, by anyone: any correlation stratified by source
count; any result on the control-excluded cohort; any result using the
194-feature frozen rule.

Cohort sizes and source-count distributions below were read from well metadata
before this registration was frozen. No direction instability was computed to
obtain them. That is design information and is reported in
`results/qc_control_audit/control_audit.json`.

## Control-exclusion rule

The assembled profile table carries no control-role column, so controls are
identified by plate coverage. **Exclude DMSO (`JCP2022_033924`) and every
compound present on at least 50% of eligible COMPOUND plates.** Applied to the
release this names exactly eight compounds, each on 93.6–100.0% of the 1,557
eligible plates:

    JCP2022_012818  JCP2022_025848  JCP2022_035095  JCP2022_037716
    JCP2022_046054  JCP2022_050797  JCP2022_064022  JCP2022_085227

The most frequent retained compound sits at 6.9% of plates, so no compound lies
near the threshold. No other frequency-based exclusion is permitted. The eight
identifiers are frozen in `results/qc_control_audit/excluded_control_ids.txt`.

## Research questions or hypotheses

The eligible cohort is dominated by minimally covered compounds: of 24,992
compounds observed in at least five sources, 23,348 are observed in exactly
five and 1,644 in six or more. The registered R6 bins (5–6, 7–8, 9–10) cannot
support a contrast, because the 9–10 bin holds six compounds. The bins below
are set from those counts alone.

If cross-source consistency is estimated from few sources, its sampling
variance is larger, and a ranking built from it should be more sensitive to
perturbation of the feature space. That predicts an interaction between source
coverage and ablation robustness.

**H1.** The raw-versus-ablated Spearman correlation is higher among compounds
observed in six or more sources than among compounds observed in exactly five.

**H2.** Within the five-source stratum, the correlation still meets the
registered preservation threshold of $\rho \geq 0.83$, so any coverage effect
is a matter of degree rather than a failure confined to sparse compounds.

H1 carries the design. H2 exists so that the registration can return a result
that strengthens the score rather than only one that embarrasses it: if H2
holds, rank preservation is not an artifact of well-covered compounds.

## Sample size

23,348 compounds at five sources; 1,644 at six or more. Both strata support a
Spearman correlation and a compound-level bootstrap. The design cannot resolve
coverage beyond seven sources: 296 compounds have seven or more and 6 have nine
or more, so no claim is made about high-coverage behavior.

## Inference criteria

| hypothesis | holds when |
|---|---|
| H1 | $\Delta_\rho = \rho_{\geq 6} - \rho_{5} > 0$ and its 95% bootstrap CI excludes zero |
| H2 | $\rho_{5} \geq 0.83$ |

$\Delta_\rho$ is estimated by 10,000 compound-level bootstrap resamples drawn
independently within each stratum, seed 20260910. Both strata are reported with
their correlations and interval regardless of direction. The three registered
R6 bins are additionally reported descriptively, with the 9–10 bin's size stated
alongside it.

**All hypotheses are void if the control-exclusion rule names a number of
compounds other than eight.**

**H1 is void if either stratum holds fewer than 500 compounds after exclusion.**

**H2 is void if the ablation feature list fails to resolve against the release
schema.**

## Analyses whose status this registration fixes

| analysis | status |
|---|---|
| 194-feature ablation, frozen rule at `249abaf` | post-foreknowledge specification check; not confirmatory |
| 252-feature ablation | unregistered sensitivity analysis |
| R6 three-bin stratification | registered descriptive analysis, executed on the corrected cohort |
| H1 coverage contrast | confirmatory, registered here before computation |
| matched raw-well contrast | repaired sensitivity analysis; its direction was seen under defective code |

No analysis in this table other than H1 and H2 is reported as confirmatory.

## Design and procedure

Ablation removes a frozen, hashed feature list from each compound's per-source
consensus profile; direction instability is recomputed and the two rankings are
compared by Spearman correlation. Consensus profiles are per-source medians of
replicate wells on COMPOUND plates of the pinned release
(`profiles_var_mad_int.parquet`, 803,853 rows, 3,180 features). Caches are
cleared before the run.

**Existing data:** yes; the release is public and pinned. **Data collection:**
N/A — no new data are collected. **Blinding:** N/A — no treatment is assigned.
**Missing data:** wells with non-finite feature values are excluded at load and
their count is reported. **Outliers:** none are removed beyond the control rule
above. **Exploratory analyses:** any analysis not named in this file is
exploratory and is labelled as such.
