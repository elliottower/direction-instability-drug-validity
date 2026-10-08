# Pre-Registration Deviation Log

## Entry 1: Gene-set method and parameters for H3/H4/H5

**Date:** 2026-07-07
**Status:** Pre-specified before running experiments 03, 04, 05
**Applies to:** H3 (phenotype projection), H4 (localization), H5 (transport stability)

### Context

The original pre-registration (commit 8d74bd0) froze hypotheses, decision
thresholds, and scorer functions (geometry/bracket_norm.py) but did not specify
how to construct the pathway gene sets needed by experiments 03 and 04. The
experiment scripts were listed in the frozen-scorer section but were never
written. This entry closes that gap by declaring the gene-set method and all
analytic parameters before any experiment is run.

### Declared analytic choices

**Primary gene-set method (H3, H4): genetic perturbation connectivity.**
For each drug with an annotated target gene (from the Drug Repurposing Hub via
frozen_drug_labels.json), the "on-target" signal is defined by the matched
shRNA knockdown signature of that target gene in LINCS L1000
(lincs_shrna.npz). This avoids researcher degrees of freedom in gene-set
construction (no GO-BP union, no arbitrary pathway database selection).

- **H3 phenotype_direction:** The consensus shRNA knockdown signature
  (mean across all cell lines and hairpins) of the drug's annotated target
  gene, used as the phenotype_direction argument to
  phenotype_projected_bracket().

- **H4 region_mask:** The top N=100 genes by absolute z-score in the
  consensus shRNA knockdown signature of the drug's annotated target gene.
  N=100 is fixed at ~10% of the 978 landmark genes, within the
  benchmark-optimal range for L1000 connectivity methods (Xie et al. 2019,
  Briefings in Bioinformatics 21(6):2194).

- **H5 cross-validation:** Leave-one-cell-line-out. For each fold, compute
  transport_stable_bracket and raw direction_instability on the training
  cell lines. Held-out outcome is cosine similarity between the held-out
  cell line's signature and the training-set consensus (mean) direction.
  This matches the direction-based estimand of DI. Spearman rank
  correlation between bracket score and held-out cosine is the test
  statistic.

**Sensitivity analysis (H3, H4): Drug Repurposing Hub MOA grouping.**
As a secondary check, repeat H3 and H4 using MOA class membership as the
pathway definition: for each MOA class with >=10 drugs, compute the mean
drug signature across the class as the phenotype_direction (H3) or use the
top-100 genes of that mean signature as the region_mask (H4). If both
primary and sensitivity analyses agree on pass/fail, the result is robust.
If they disagree, report both and flag the discrepancy.

**Inclusion criteria:**
- Drugs must have >= 5 cell lines in the compound LINCS data
- Target gene must have >= 3 shRNA hairpins in lincs_shrna.npz
- Only drugs with non-null target annotation are eligible for H3/H4

### Justification

Using genetic perturbation connectivity rather than pathway databases
eliminates the largest unregistered knob (gene-set construction method and
size). The shRNA data is already in the LINCS ecosystem, uses the same 978
landmark genes, and provides a direct empirical definition of "on-target"
without requiring curator judgment about which GO terms are relevant.

### What this does NOT change

- Hypotheses H3, H4, H5 and their decision thresholds remain exactly as
  specified in PREREGISTRATION.md (commit 8d74bd0)
- Scorer functions in geometry/bracket_norm.py remain unchanged
- All other pre-registered decisions (H1, H2, H6) are unaffected


## Entry 2: H5 fold threshold adjustment for fewer than 5 valid folds

**Date:** 2026-07-07
**Status:** Pre-specified before running experiment 05
**Applies to:** H5 (transport stability)

The pre-registered criterion is "transport-stable bracket outpredicts raw in
at least 4 of 5 folds." If fewer than 5 cell lines produce valid folds
(i.e., >=20 drugs with data in both training and held-out), the threshold
adjusts to ceil(n_valid_folds * 0.8). This maintains the 80% win-rate
requirement while avoiding discarding the experiment when a rare cell line
has too few drugs.

In practice, the LINCS data has 71 cell lines and the first 5 alphabetically
all produced valid folds, so this fallback was not triggered.


## Entry 3: H3 and H4 shared dependency on shRNA consensus

**Date:** 2026-07-07
**Status:** Documented before interpreting results
**Applies to:** H3 (phenotype projection), H4 (localization)

H3 uses the shRNA consensus signature as the phenotype_direction vector.
H4 uses the top-100 genes of that same shRNA consensus as the region_mask.
These are not independent tests — both derive from the same genetic
perturbation object. The paper must state this shared dependency and not
present H3 and H4 as two independent confirmations of the validity ladder.
H5 (transport stability) is fully independent of H3/H4 since it uses no
gene sets or shRNA data


## Entry 4: H5-full replaces original 5-fold H5

**Date:** 2026-07-08
**Status:** Pre-registered extension replaces original 5-fold analysis
**Applies to:** H5 (transport stability)

### Context

The original H5 used 5 alphabetically-selected cell lines. The full
leave-one-out (H5-full) was pre-registered in PREREGISTRATION_EXTENDED.md
(f3c88ed) with two co-primary criteria: (1) TS wins >= 80% of valid folds,
(2) Wilcoxon signed-rank p < 0.05 on per-fold delta-rho.

### Result

Both co-primary criteria satisfied: 66/66 folds (100%), Wilcoxon
p = 8.2e-13, mean delta-rho = +1.10 [1.00, 1.19]. Frequency-improvement
correlation rho = -0.32 (p = 0.009) confirms the rare/common regime
distinction generalizes. The paper reports H5-full as the primary result
with the original 5-fold table retained as representative examples.

### What this changes

- H5 is now reported as 66/66 folds (not 5/5)
- Wilcoxon signed-rank is a co-primary criterion alongside fold-win fraction
- The "5 alphabetical cell lines" limitation is removed from the paper

## Entry 5: H3-CRISPRi convergent-validity outcome (c)

**Date:** 2026-07-08
**Status:** Pre-registered outcome realized; interpretation follows frozen decision tree
**Applies to:** H3-CRISPRi (convergent-validity check, PREREGISTRATION_EXTENDED.md)

### Pre-registration reference

Committed at f3c88ed (tagged prereg-extended-psb) on 2026-07-08 11:56:58,
before experiments ran (~12:00-12:15). The three-outcome decision logic was
frozen as:
- (a) rho_CRISPRi >= 0.376: convergent validity confirmed
- (b) 0 < rho_CRISPRi < 0.376: robustness confirmed with K562 confound explanation
- (c) rho_CRISPRi <= 0: interpret via K562-only confound first

### Result

Outcome (c) realized: rho_CRISPRi = -0.12, p = 0.16 (n = 131).

Secondary criterion fired as pre-specified: raw bracket correlates
positively with CRISPRi enrichment (rho = +0.27, p = 0.002) while
projected bracket does not — the reverse of the shRNA pattern. This
indicates the CRISPRi ground truth itself is confounded by K562
cell-type-specific effects, consistent with the pre-registered
explanation.

### Interpretation (following frozen decision tree)

The shRNA-based H3 result (rho = 0.38) is not refuted but not reinforced.
CRISPRi does not provide independent convergent validity, attributable
to the pre-registered single-cell-line limitation.

Post-hoc power analysis (added after initial commit): n=131 has >99%
power to detect rho=0.376 (alpha=0.05, one-sided). The CRISPRi 95% CI
[-0.29, 0.05] excludes the shRNA effect size entirely; the two
correlations are significantly different (z=5.47, p<10^-7). The CRISPRi
result is therefore a genuine divergence from the shRNA finding, not an
underpowered null. The K562-only confound remains the pre-registered
explanation for this divergence.

### What this does NOT change

- H3 shRNA result and its Holm-Bonferroni survival are unaffected
- The pre-registered interpretation was followed exactly as frozen
- No post-hoc re-analysis or threshold adjustment was performed

**Note added 2026-09-22.** The shRNA values in this entry (rho = 0.376 and 0.38,
the z-test against it, and the matched comparison) were computed under the
indexing defect recorded in Deviation 9; the entry is left as written. Corrected
values are in Deviation 9. Two readings in this entry are not supported by what
was tested: that the raw correlation "indicates the CRISPRi ground truth itself is
confounded by K562 cell-type-specific effects", and that the divergence is
"attributable to" the single-cell-line limitation. No analysis of cell type was
run. The K562 account is one candidate among several, registered for test in
`experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`.

---

## Deviation 6: H5 win test compared oppositely-oriented statistics

**Pre-registered:** Transport-stable instability outpredicts raw direction
instability in >= 80% of valid leave-one-cell-line-out folds, evaluated as
`ts_rho > raw_rho` on Spearman correlations against held-out cosine.

**Actual:** H5 is reported as not confirmed.

**Why:** Direction instability measures inconsistency and held-out cosine
measures agreement, so a working score must correlate negatively with the
outcome. Transport-stable instability subtracts a Frechet-variance term that
also grows with inconsistency and grows faster, reversing the sign. The
registered comparison of signed correlations is therefore satisfiable by
orientation rather than by prediction, and returns 66/66 for any sufficiently
large penalty weight. The registration omitted an orientation convention; that
omission is the deviation.

Compared on predictive strength, |ts_rho| > |raw_rho| in 10 of 66 folds, and
the transport-stable variant loses in the two largest (A375: 0.155 against
0.480, n = 8,328; A549: 0.201 against 0.494, n = 8,888).

**When:** Post-analysis.

### What this does NOT change

- The 66-fold computation itself is unchanged; only the comparison is corrected
- Raw direction instability predicts held-out cell-line agreement, mean Spearman
  rho = -0.602, correctly signed in all 66 folds
- The permutation null (p_perm < 10^-3) still establishes that the prediction is
  drug-specific rather than algebraic
- No other hypothesis depends on H5

Reproduced by `experiments/audit_h5_orientation_and_hdac.py`; output in
`results/audit_h5_orientation_and_hdac/audit.json`.

---

## Deviation 7: Figure 5 rebuilt from committed source data

**Pre-registered:** N/A -- figure construction was not specified.

**Actual:** `fig5_hdac()` previously hardcoded six direction-instability values
(0.14-0.45) and six cell-line counts that reproduce no stored artifact. It now
reads `results/fig5_hdac_source/fig5_hdac_source.json`, built by
`experiments/build_fig5_source.py` from this repository's own toxicity-correction
results, and asserts at build time that every plotted value equals its source row.

Selectivity is recorded as three ordered categories rather than six untied ranks:
the pan-HDAC inhibitors differ in potency rather than isoform breadth, and the
relative selectivity of tubacin and PCI-34051 is assay-dependent. All 11
cross-category comparisons are concordant, Kendall tau_b = 0.856. No significance
test is reported; six compounds do not support one.

**When:** Post-analysis.

---

## Deviation 8: the registered rebuild names a loader the deposited artifact did not use

**Pre-registered:** rebuild the LINCS matrices "with the existing loader",
pinned as `drug-perturbation-geometry@1dc20a2`, `data/lincs_loader.py`
(`PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md`, frozen `f288507`).

**Actual:** the loader is used to read signature metadata. Per-cell-line
aggregation follows `experiments/03_phenotype_projection.py`, which produced the
deposited artifact.

**Why:** the loader's `get_consensus_signatures` selects one signature per cell
line, by highest `distil_ss` where that column exists and otherwise the first
row. The deposited H3 artifact was not built that way: it averages every
signature for each drug-cell pair. Reproducing the deposited `D`, `P` and `E` to
the registered tolerance of 1e-6 is only possible with the mean. The registered
phrase did not distinguish the loader's metadata functions from its consensus
function, and the two were not the same procedure.

Separately, `GSE92742_Broad_LINCS_sig_info.txt.gz` carries no `distil_ss`
column, so the loader's primary selection rule is unavailable on this release
regardless.

**When:** before any statistic was computed. The first extraction attempt failed
on the missing column rather than producing a cohort.

### What this does NOT change

- The registered hypotheses, statistics, criteria and seeds are untouched
- The per-drug reproduction requirement is unchanged and remains the check that
  the rebuilt cohort is the deposited one
- The analysis is void if reproduction fails, as registered

---

## Deviation 9: shRNA consensus directions were indexed by metadata row, not signature id

**Pre-registered:** H3 tests whether phenotype-projected instability correlates
with on-target enrichment, where the on-target direction is the consensus shRNA
knockdown signature of the drug's annotated target gene (`PREREGISTRATION.md`,
commit `249abaf`).

**Actual:** the deposited H3 artifact was computed against consensus directions
that are not those genes' hairpin means. `build_shrna_consensus` mapped each
signature id to its row position in `lincs_shrna_siginfo.csv.gz` and used that
number to index the signature matrix in `lincs_shrna.npz`. The two files carry
the same 154,993 signature ids in different order, agreeing at 115 of 154,993
positions. Each target's consensus was therefore the mean of
metadata-position-selected unrelated shRNA signatures. The compound side of the same scripts indexes by npz position
and is unaffected.

**How it was identified:** the rebuild-validation gate of
`PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md` (frozen `f288507`) failed. That
registration required rebuilt values to reproduce the deposited artifact within
1e-6, and correct reconstruction does not. S1, S2 and S3 were therefore void
before any of them was computed, and none was. The rebuild reproduced raw
direction instability to 7.5e-08 and the per-drug
cell-line counts exactly, while projected instability and enrichment differed by
up to 7.82 and 0.448. Replicating the metadata-row indexing reproduced the
deposited values to 1.9e-07 and 1.3e-07, which identifies how they were
produced.

**Scope:** `03_phenotype_projection.py` (H3), `04_localization.py` (H4) and
`03b_h3_crispri_ground_truth.py` (the CRISPRi convergent-validity check) carry
the same function. A repository-wide audit found no other instance: every other
position map is built from the signature matrix's own ids. Raw direction
instability never uses a target direction, so H1, H2, H5 and the held-out
cell-line prediction do not depend on it.

**Corrected results, H4 and CRISPRi:** both were recomputed with the corrected
function (commits `fde2d9b` and `f822fb1`, 2026-09-15). H4's macro-averaged AUROC
gap is +0.0005, 95% CI [−0.0256, +0.0282], against +0.018 deposited; H4 was not
confirmed before or after, and the per-class gaps behind the "domain-dependent"
reading also move. The CRISPRi convergent-validity check gives rho = −0.1247
(p = 0.156, n = 131) for projected instability against CRISPRi enrichment,
against −0.12 deposited; its CRISPRi directions never used the defective path.
On the 114 drugs with both references, the shRNA projected rho is 0.3398 against
0.42 deposited.

**Provenance:** each superseded artifact is kept byte-for-byte under
`superseded/misindexed_shrna_consensus/` beside its corrected replacement, with a
README recording both hashes.

| artifact | superseded sha256 | corrected sha256 |
|---|---|---|
| `results/03_phenotype_projection/phenotype_projection_results.json` | `65e5d10e2720…` | `fd69e26fc9a3…` |
| `results/04_localization/localization_results.json` | `36808886a93a…` | `9c7ada30ad63…` |
| `results/03b_h3_crispri/h3_crispri_results.json` | `1e96ef13ad1d…` | `f00f84280711…` |

The corrected H3 and H4 artifacts were produced from `lincs_subset.npz` (sha256
`2ad0f5d30ab8…`) and `lincs_shrna.npz` (`4a990e5072a4…`) on the Modal volume
`drug-perturbation-vol`; the CRISPRi artifact additionally from the scPerturb
`ReplogleWeissman2022_K562_essential.h5ad` (Zenodo 10044268, sha256
`412fd0df8c4c…`).

**Corrected result, H3:** with consensus directions built by signature id,
Spearman rho between projected instability and on-target enrichment is 0.3172
(p = 4.8e-20, n = 795), against 0.3756 deposited. Raw direction instability
against enrichment is -0.0580, against -0.0433 deposited. Both pre-registered
criteria are still met: |rho_projected| > 0.3 and |rho_raw| < 0.15.

**When:** post-analysis, before the registered sensitivity analyses were
computed.

### What this does NOT change

- H1, H2, H5 and the 66-fold held-out prediction result, none of which use a
  target direction
- The registered criteria themselves, which H3 still meets after correction

Meeting H3's registered numerical criterion is not the same as confirming the
construct it was meant to test. Convergence with the independent CRISPRi reference
is absent: under CRISPRi, projected instability is not positively associated with
alignment and raw instability is (rho = +0.2723, p = 0.00165, n = 131). H3's
interpretation as separating on-target from off-target consistency is unresolved
and is the subject of `experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`.

### What it sharpens

Projected instability and enrichment are constructed from the same target
direction, so an association between them is available whether or not that
direction is the annotated target. The erroneous construction shows that a
shared direction can generate a substantial association without correct
drug-target matching: it produced a *higher* correlation (0.3756) than the
correct directions (0.3172). It is not itself a draw from the registered
target-permutation null, which preserves the multiset of valid target
directions, and unrelated mixtures need not share their norm or covariance
structure. The corrected association remains subject to the same shared-axis
concern, and S2 has not been run.

---

## Deviation 10: the registered description of the PRISM release does not match the release

**Pre-registered:** R5's toxicity measure is "minus the median, across cell lines, of
the PRISM Repurposing 19Q4 primary-screen log-fold-change (2.5 µM,
replicate-collapsed)", and a drug is mapped "by exact name match after
normalization, against the `Drug.Name` and `Synonyms` fields of the pinned compound
list" (`experiments/PREREG_H3_REFERENCE_DISCORDANCE.md`, frozen `7f57136`).

**Actual:** the release was retrieved from figshare 9393293 on 2026-09-30 and
differs from that description in two ways.

- Its treatment table carries `name` and no synonym field, so exact-name matching
  uses `name` alone. The registered rule is otherwise unchanged: identifiers first,
  then exact names normalized by case folding and whitespace, and no source outside
  the pinned files.
- Its doses are not uniformly 2.5 µM. The primary screen holds one treatment per
  compound, 4,147 in the HTS screen at 0.03–5.6 µM and 539 in MTS004 at 2.29–3.12
  µM. Each compound is therefore screened at one dose, which is recorded per drug
  rather than filtered on.

The registered field names and dose were written from the 24Q2 public release,
which is the one that was on disk when the registration was drafted, and which
covers 5 of the 131 CRISPRi-arm drugs.

**When:** at retrieval, before any log-fold-change value was summarized and before
any R5 statistic was computed.

### What this does NOT change

- The mapping rule's substance: identifiers first, then exact names, sourced only
  from the pinned release, with the accepted and rejected tables written before any
  value is summarized
- R5's gate, its quantities, its intervals or its readings
- Coverage, counted at mapping and before any value was summarized: 121 of the 131
  CRISPRi-arm drugs map over 38 targets, against a registered minimum of 60 drugs
  over 20 targets

Recorded in `results/03d_h3_reference_discordance/r5_prism_mapping.json` at run
time, with the dose of every mapped drug.

---

## Deviation 11: the shRNA extraction stores correct values under a wrong gene axis

**Pre-registered:** the reconstruction gate of
`experiments/PREREG_H3_S1S3_CORRECTED_BASELINE.md` (frozen `7f57136`) requires an
independent rebuild from the GEO GCTX to reproduce the pinned extraction, joined on
identifiers, for both the compound matrices and the shRNA target consensuses.

**Actual:** the gate ran on 2026-10-02 and failed on its shRNA half. The forensics
that followed establish what the disagreement is.

- Every one of the 14,656 shRNA signatures differs from the rebuild as stored, with
  a median maximum absolute difference of 6.63.
- The two versions have identical norms, to sixteen decimal places, and a median
  cosine of 0.004. Equal norms with no alignment rules out a rescaling and points
  at a permutation; the two bullets below are what establish one.
- Each signature holds the same multiset of values, and no rebuilt vector matches a
  different extraction row, so the permutation is of the gene axis rather than of
  the signature labels.
- One bijective permutation maps the extraction onto the rebuild for **all 14,656
  signatures with a maximum absolute difference of 0.0**.
- That permutation is the only one: all 978 gene columns carry distinct value
  profiles across the 14,656 signatures, and the closest pair of columns, genes
  1788 and 2523, differs by 5.63 elementwise against a gate tolerance of 1.73e-05.
  Two identical columns would have admitted a second reconciling map; none exist.
- The order that permutation implies matches none of the orders the pipeline uses:
  not the gene-info file's row order, not the landmark ids sorted as strings, not
  the gene symbols alphabetically.
- `lincs_subset.npz` and `lincs_shrna.npz` declare the *same* gene axis, and the
  compound half of the gate passed against that declared axis before the shRNA half
  ran. The compound matrix follows the declared labels; the shRNA matrix does not.

**What this means:** `lincs_shrna.npz` holds the right numbers under the wrong
labels. Any quantity pairing a drug signature with a target direction — H3's
signless target-axis alignment and projected dispersion, and H4's localization —
was computed between vectors in different gene orders. The Deviation 9 replacement
values for H3 and H4 are invalidated rather than merely shifted: they are not
measurements of their stated quantities, because that deviation fixed which
signatures were averaged into each consensus and never which gene each column
carried.

**Scope:** H3 and H4, and the shRNA fields of the CRISPRi records
(`proj_shrna`, `enrich_shrna`). Raw direction instability, H1, H2, H5 and the
66-fold held-out prediction never pair the two matrices and are unaffected. The
CRISPRi directions are built from the Perturb-seq file into the compound matrix's
own coordinates and are unaffected.

**How it was identified:** the registered reconstruction gate, which exists because
Deviation 9 was an identifier-versus-position error and a gate comparing values by
position could not have seen it. The gate's shRNA half was added after a review
observed that a compound-only comparison cannot test the object Deviation 9
corrupted.

**Why the checks already in the code could not see it.** A `Reference` carries
`positions`, an index into the 978-landmark order, so every matrix is assumed to
follow one order and a positional axis cannot disagree with itself.
`experiments/build_h3_bundle.py:48` does compare the two files' declared gene axes
and asserts they match — the assertion passes, because both declare the same
identifiers. The same script reproduces the deposited `D`, `P` and `E` to within
1e-6 and asserts on it; that check cannot detect the defect either, because it
reads the mislabeled matrix on both sides of the comparison. Only a comparison
against an independent parse of the pinned GCTX, which is what the gate performs,
reaches the values rather than the labels.

**When:** before any S1-S3 statistic, and before any R0-R7 statistic, was computed.

**Artifacts:** `results/03c_h3_sensitivity/shrna_axis_recovery.json` carries the
recovered permutation and its verification;
`gate_shrna_diagnostic.json`, `gate_permutation_test.json` and `gate_axis_test.json`
carry the forensics above; `shrna_axis_uniqueness.json` carries the column-profile
comparison behind the uniqueness of the recovered permutation.

### What this does NOT change

- H1, H2, H5 and the 66-fold held-out prediction, none of which pair the two
  matrices
- The CRISPRi directions, or any quantity built from them alone
- The R5 mapping, which is frozen and uses neither matrix

**Superseded in part by Deviation 12.** The evidence above stands. The reading
drawn from it does not: the compound matrix carries the same permutation, so the
quantities said here to pair vectors in different gene orders did not.

## Deviation 12: both extractions carry the same wrong gene axis, and the paired quantities hold

**Pre-registered:** Amendment 2 to the corrected-baseline registration, frozen
`fd1ae8d`, specifies a repair on the reading that the shRNA extraction's gene axis
is mislabeled while the compound extraction's is not. Its A7 requires the
reconstruction gate to pass on both halves before any statistic runs.

**Actual:** the gate ran on 2026-10-03 against the preserved canonical artifacts.
Two of its three comparisons passed and one failed.

| comparison | result |
|---|---|
| canonical shRNA signatures and consensuses against the source | passes: 14,656 signatures at maximum absolute difference 0.0, 258 targets, stored directions at 1.1e-16 |
| compound rebuild's drug-by-cell matrices against the source | passes: 3.1e-06 across the 795-drug cohort |
| retained `lincs_subset.npz` against the source | fails: every one of the 795 drugs at exactly 20.0, the distance between the clipping bounds |

The forensics that followed establish what the compound disagreement is.

- It is present at the level of single signatures, so no grouping or aggregation is
  responsible. Coverage is complete at 41,643 of 41,643 cohort signatures, the two
  copies of the compound metadata are identical, and negation does not explain it.
- Each signature holds the same multiset of values as the source signature its label
  names, tested by sorting: 40 of 40, and no row is monotonic.
- Matching each declared column against the source columns by its profile across
  400 signatures recovers a bijection: **978 of 978 columns, two of them in place.**
- That order is the same order the shRNA extraction carries. Compared against the
  permutation of Deviation 11 across all 978 columns, **no position differs.**

**What this means:** both extractions declare one gene axis and hold another, and it
is the same other. A cosine is invariant under a permutation applied to both
operands, so every quantity computed from the two files together sits in one
coordinate system and is unaffected by the mislabeling. Raw direction instability
recomputed from the extraction and from the pinned source both reproduce the
deposited values, to 7.3e-09 and 3.0e-08 across twenty drugs, which is what that
invariance predicts.

Deviation 11's evidence stands and its reading does not. H3 and H4 were not
computed between vectors in different gene orders. The reading rested on the
compound half of the gate appearing to pass, which came from a gate version with a
known intersection defect that wrote no record.

**Amendment 2's premise is therefore refuted.** The amendment is frozen and is not
edited. Its repair — rebuilding the target directions from the GEO rebuild — would
change no reported quantity, because the quantity it was meant to repair was never
computed across mismatched axes. Its A8 layers, its supersession discipline and the
source comparison it requires of the gate stand on their own and are kept.

**What is affected.** Any quantity pairing an extraction with a reference labeled
outside it, and any operation that selects a column by gene name.
`experiments/03b_h3_crispri_ground_truth.py` maps symbols through the declared gene
ids and places the Perturb-seq reference in that order while the drug signatures sit
in the permuted order, and `geometry/references.py`'s `without_gene` drops a
coordinate by symbol. At discovery neither had been measured. The pooled C0 arm has
since been recomputed, below; the registered C1-K562, C1-RPE1 and C1-GW arms and
R7f have not run.

**How it was identified:** the registered reconstruction gate, on the comparison
Amendment 2 added — the retained compound extraction against an independent read of
the pinned source. The earlier gate compared that file only against the rebuild.

**When:** before any S1-S3 statistic and before any R0-R7 statistic was computed.

**Artifacts:** `results/03c_h3_sensitivity/gate_run_2026-10-03/reconstruction_gate.json`
is the gate's own report, the first a failing gate has written.
`compound_discrepancy_diagnostic.json`, `compound_row_test.json`,
`compound_column_match.json`, `compound_defect_scope.json` and
`shared_axis_confirmed.json` carry the forensics above, and
`input_pin_check.json` records that the five inputs matched the registration's pins
before the run.

**A count that does not match the record.** `lincs_subset.npz` holds 167,266
signatures. Deviation 9 describes the source files as carrying 154,993. The two
statements are about the run each was written from and have not been reconciled.

### The CRISPRi arm, measured

The registered convergent-validity check pairs a drug signature with a CRISPRi
reference placed by gene symbol, so it is the one quantity here whose operands sit
in different coordinate systems. Four routes on one fixed 131-drug, 41-target
cohort, with the permutation recovered from the pinned source before any CRISPRi
quantity was computed:

| route | projected rho | raw rho |
|---|---|---|
| the deposited configuration | -0.1247 (p = 0.16) | +0.2723 |
| corrected, extraction relabeled | **+0.5410** (p = 2.6e-11) | -0.1920 |
| corrected, through the verified rebuild | +0.5410 | -0.1920 |
| control: one permutation applied to both operands | +0.5410 | -0.1920 |

The deposited configuration reproduces the deposited values to sixteen decimals,
the two corrected paths agree, and the control returns the corrected value, which
it must, because no cosine changes when one permutation is applied to both
operands.

The deposited C0 reversal did not survive coordinate correction. The registered
K562-only confound explanation addressed an artifact of the reference coordinates
and is not applied. The registered H3 criterion is projected `|rho| > 0.3` and raw
`|rho| < 0.15`: corrected, the CRISPRi arm meets the first and not the second,
0.541 and 0.192, so it satisfies one of the two conditions. The shRNA arm
satisfies both.

This covers the pooled C0 construction on the deposited cohort. The registered
R0-R7 arms built on C1-K562, C1-RPE1 and C1-GW have not run.

### What this does NOT change

No deposited value has been replaced. Quantities computed solely between the
compound and shRNA extractions are invariant to their shared permutation and
reproduce the deposited results in the checks reported above: H1, H2, H3, H4, H5,
raw direction instability and the 66-fold held-out prediction. Quantities involving
an externally labeled reference or a coordinate selected by gene name: corrected C0
is measured and reported below, and the rest await computation and are not to be
interpreted until recomputed. No stored deposited value has been overwritten; one
reported scientific quantity now has a corrected replacement estimate.

- The canonical shRNA artifacts, which the gate verified against the source
- The R5 mapping, which is frozen and uses neither matrix

### The affected surface, named

`experiments/03d_h3_reference_discordance.py:604` builds the symbol list from the
extraction's declared `gene_ids`, and that list places five external references:
C0, C1-K562, C1-RPE1, C1-GW and C1-GW-phenotype-positive. Every externally mapped
CRISPRi arm in R0-R7 therefore inherits the mismatch, not only the pooled
construction of `03b`. The same list feeds R7f's target-gene sensitivity at
`03d:950`, where `Reference.without_gene` removes the coordinate *named* for a
target rather than the one holding that target's values.

The deposited manuscript reports CRISPRi quantities, and they are the exception to
the paragraph above. The deposited C0 arm has now been recomputed; the remaining
R0-R7 external-reference quantities and the symbol-selected sensitivity of R7f
await computation.

`03d` calls no load-time axis check, so the script most directly concerned with the
CRISPRi comparison could consume the known-bad declaration. The check is added
there.

### Scope of the generalization

The shared permutation was recovered from 400 compound signatures and matches the
shRNA permutation, itself verified on all 14,656 shRNA signatures, at every one of
978 positions. Applied to **all 41,643 cohort compound signatures**, it reconciles
every one of them with the pinned source at a maximum absolute difference of
**0.0**, with no failures. The invariance claim therefore rests on the cohort
rather than on the block the permutation was recovered from. The multiset test
covered 40 signatures and the direction-instability reproduction 20 drugs; both are
now corollaries of the exact all-cohort agreement.

## Deviation 13: released Replogle pseudobulk sources carry undefined matched-landmark coordinates, and R0.4 and R0.6 used the wrong gene space

**Date:** 2026-10-06
**Status:** found while running R0.4; R0.4 and R0.6 have run, R0.5 and R1-R7 have not
**Applies to:** R0.2, R0.4, R0.5, R0.6, R0.8, R4, R6a, R7c, and every quantity built
from C1-RPE1, C1-GW or C1-GW-phenotype-positive. Seven matched landmark coordinates
are affected across two source files: one in the RPE1 release and six in the
genome-wide release. The RPE1 mechanism is established below; the genome-wide
coordinates are recorded as deposited `+inf` values with no mechanism attributed.

### What was found in the RPE1 release

`rpe1_normalized_bulk_01.h5ad` and `RPE1_essential_single_cell.h5ad` hold `+inf`
in two gene columns and in no others: ATF3 (`ENSG00000162772`) and CCL2
(`ENSG00000108691`), 253 values over 240 of 2,679 bulk rows and 7,913 values over
7,913 cells. The K562-essential release is clean and holds neither gene on its
axis.

The mechanism is in the released normalization. Replogle et al. z-normalize each
gene within each GEM group against that group's non-targeting control cells. In
group 46 the 158 control cells all hold 0.0 for ATF3, and in group 9 the 180
control cells all hold 0.0 for CCL2, so the divisor is zero. The raw counts are
finite for every gene in every group; no gene has constant raw values. The
reconstruction of the published procedure from the raw files does not reproduce
the deposited matrix — 42 genes sampled, 10,404,475 values compared, maximum
absolute difference 1824.3968881650242 on the finite positions, zero of 42 genes
agreeing within 0.1 — while the non-finite positions agree exactly. No deposited
value is therefore replaced.

Of the two genes, ATF3 is not a LINCS landmark and never enters the analysis.
CCL2 is landmark `6347` and a cohort target, so the defect reaches the registered
C1-RPE1 arm. The RPE1 axis matches 813 of the 978 landmarks, which is a property
of the deposit and not of this defect.

### The first implementation defect: a NaN routed into a threshold branch

The registered R0.4 rule fires when the median per-target cosine is below 0.99. For
C1-RPE1 the median was `NaN`, over 2,393 targets of which 1,943 returned `NaN` and
450 returned a finite cosine. A `NaN` median is not below 0.99, and
`float(np.median(cosines)) >= 0.99` is false for it, so the implementation labeled
C1-RPE1 as describing the single-cell construction rather than C1 on undefined
arithmetic rather than on measured disagreement. The registration specifies no
handling for non-finite input, and routing it into the failure branch is not what
it registers.

### The second implementation defect: the wrong gene space

`experiments/modal_03d_single_cell.py:stage_targets` computes the single-cell mean
and compares it with the released pseudobulk over the file's full gene axis, 8,749
genes for RPE1 and 8,563 for K562, and estimates split-half reliability on that
same full axis. The registration defines C1 on the matched landmarks: "The
direction is restricted to landmarks present in the file, matched by symbol through
the pinned gene-info file, and drug signatures are restricted to the same
landmarks." R0.4 verifies a cosine against C1 and R0.6 estimates the reliability of
the object R7c then qualifies, so both belong on each construction's analysis
basis. On the full axis, R0.4 verified one vector space while R0.6 assessed
another, and non-landmark ATF3 entered the verification at all only through this
defect.

### What had already been seen when the rule was written

The 450 finite RPE1 targets gave a median cosine of 1.0000000100052202 with a
minimum of 0.8097611585567198 and two targets below 0.99, and the 2,057 K562
targets gave 0.9999999986457878 with no `NaN`. The finite-subset numbers were
therefore in hand before any coordinate rule was frozen. They are recorded here and
are not reported as R0.4: the median over the targets that happened to survive on
the full gene axis is not an estimate of the median on a frozen landmark basis, and
R0.4 is recomputed once under Amendment 4 rather than read off this subset.

### What was found in the genome-wide release, as a separate defect

Applying the rule to every construction that reads an official pseudobulk file,
rather than to the gene the RPE1 diagnosis named, found a second affected source.
The genome-wide release `K562_gwps_normalized_bulk_01.h5ad` leaves six of its 721
matched landmarks undefined: ICAM1, MEST, PXN, SLC25A14, BAMBI and TCTN1, with 61
to 146 `+inf` values each over 11,258 rows, all `+inf` and no `NaN`. C1-GW and
C1-GW-phenotype-positive read that file, so both carried undefined coordinates into
R0.5 and R4, and no check had looked. The K562-essential release leaves none of its
728 matched landmarks undefined.

**No mechanism is attributed to these six.** The control-variance diagnostic that
settled the RPE1 case reads the raw single-cell file; the genome-wide single-cell
release is 66 GB, which is why the registration does not estimate genome-wide
reliability from it, and that diagnostic has not been run. They are recorded as
deposited `+inf` coordinates and nothing further is claimed about why.

| construction | source | matched | in basis | excluded |
|---|---|---|---|---|
| C1-K562 | `K562_essential_normalized_bulk_01.h5ad` | 728 | 728 | none |
| C1-RPE1 | `rpe1_normalized_bulk_01.h5ad` | 813 | 812 | CCL2 |
| C1-GW | `K562_gwps_normalized_bulk_01.h5ad` | 721 | 715 | ICAM1, MEST, PXN, SLC25A14, BAMBI, TCTN1 |
| C1-GW-phenotype-positive | the same file | 721 | 715 | the same six |

Each basis is above the registered 500-landmark floor. The counts, the positions,
the per-coordinate value counts and a hash of each ordered basis are in
`registry/frozen/analysis_bases.json`, measured by
`experiments/03i_freeze_analysis_bases.py` from the copies on the volume the
analyses read, and pinned from outside that file in
`registry/frozen/analysis_bases_pin.json`.

The registered comparisons run on the landmarks their two constructions share:
R0.5's C0 against C1 on 728, C1 against C1-RPE1 on 688, C1 against C1-GW on 712,
and R6a's C1 against C1-RPE1 on 688. CCL2 is not among the 728 landmarks the
K562-essential release matches, so excluding it changes no registered comparison;
the genome-wide exclusions cost six coordinates from C1 against C1-GW.

### What is not changed

No deposited value is altered, imputed, clipped or reconstructed. The affected
cells are not excluded. The two reconstruction attempts are recorded in
`results/03d_h3_reference_discordance/rpe1_reconstruction_not_justified.json` and
neither licenses replacing a published finite value.

The authors have been written to about the defect and asked for the exact
normalization procedure. A corrected release or the exact code, if supplied, is
analyzed as a separately identified sensitivity analysis and does not replace a
frozen run.

### Records

- `results/03d_h3_reference_discordance/rpe1_non_finite_diagnosis.json` — where the
  non-finite values are, and what R0.4 reported on them
- `results/03d_h3_reference_discordance/rpe1_normalization_mechanism.json` — the
  zero control standard deviation in the two affected groups
- `results/03d_h3_reference_discordance/rpe1_reconstruction_not_justified.json` —
  the two failed reconstructions
- `experiments/PREREG_H3_S1S3_AMENDMENT_4.md` — the prospective finite-coordinate
  rule written in response

## Deviation 14: one reproduction tolerance was carrying two different comparisons

**Date:** 2026-10-07
**Status:** the R0-R7 execution of 2026-10-07 refused at this gate and is void and
preserved. R0.4 and R0.6 have been computed under Amendment 4; no R0.5 or R1-R7
association, resampling result or criterion outcome exists.
**Applies to:** the reproduction gate of
`experiments/03d_h3_reference_discordance.py`, and through it every R0-R7
invocation.

### What refused

```
P_shrna reproduces to 1.14e-06, tolerance 1e-06
```

One drug of 795 crossed the flat tolerance: levofloxacin, deposited `proj_shrna`
8.6574612530, error 1.1385e-06, which is 2.21 float32 ULPs at that magnitude.
`P_shrna` spans 0.6424 to 16.95 across the cohort and one float32 ULP at 16.95 is
1.0106e-06. The maximum relative error is 1.96e-07 and the median error is 0.505
ULP. `E_shrna` and `D` are both well inside the flat tolerance, at 9.365e-08 and
7.496e-08.

### Two explanations, both refuted by their own measurements

The retired loader takes the mean of a float32 array and numpy accumulates that in
float32, so the first candidate was the accumulator. The float32 and float64
accumulations of the same signatures differ from each other by about 1e-9, and both
sit about 1e-6 from the deposited values. Not the cause.

The second was a float32 projection path, since the dot product and the mean over
context pairs would round at the magnitude of `P` rather than of a signature value.
No dtype route reproduces the deposited values below a 1e-7 relative floor. Not the
cause either.

What survives is a systematic relative offset of that order between the deposited
values and every route that can be constructed today. Context membership is
identical between routes for every drug examined, and among the nine drugs the
localization diagnostic examined the difference is spread rather than concentrated
in one context pair: 6 of levofloxacin's 21 pairs sit above half the maximum
per-pair difference of 5.658e-07.

So the diagnostics found no localized context-pair discrepancy among the examined
focus drugs and rejected the two tested dtype explanations; the remaining residual
has not been explained.

### A third candidate, tested and refuted, and what the test could establish

Deviation 9's pre-fix shRNA consensus indexing was the one remaining named candidate,
and the plan was frozen at `e684d16` before it ran so no variant of the indexing
could be searched. H2 holds. Over all 795 drugs the pre-fix reference differs from the
deposited records by up to 7.818 in `P_shrna` and 0.448 in `E_shrna`, against 9.755e-07
and 9.843e-08 for the corrected reference, and not one drug is closer under it.

The outcome was fixed by provenance rather than discovered. The records this gate
compares against are the corrected artifact, recomputed with the corrected function at
`f822fb1`, whose own commit title reads "the indexing defect was not the cause". A
defect whose repair produced the comparison target cannot account for a residual
measured against that target, and the 7.818 and 0.448 are the two figures this entry
already reported for the corrected reconstruction against the superseded artifact,
measured in the other direction. I proposed the candidate in the round-24 packet
without reading that far back in this log, and the review agreed to test it on the
strength of the description. The measurement is kept for what it does establish, that
the deposited reference directions are the corrected consensus. No further candidate
for the residual is named, and what remains unexplained is that a re-execution of the
corrected route does not return the corrected records bit for bit.

### Route to route, over the cohort rather than nine drugs

The earlier claim that the two modern routes agree with each other five to ten times
better than either agrees with the deposited values came from nine focus drugs and does
not survive the cohort. Over all 795:

| quantity | routes mutually closer | max ratio | median ratio |
|---|---|---|---|
| `P_shrna` | 568 of 795 | 2.5 | 2.5 |
| `E_shrna` | 611 of 795 | 5.4 | 4.5 |
| `D` | 778 of 795 | 31 | 38 |

The defensible statement is about the cohort's maxima and medians, where the routes are
mutually closer for every quantity. Per drug it is not universal, and the factor is 2.5
rather than five to ten for projected dispersion.

### The defect is that one rule was carrying two comparisons

The registration describes its flat absolute tolerance as holding between two
float64 computations from one extraction. Amendment 3's B1 then made the verified
GCTX rebuild the authoritative compound source, and the gate went on applying that
one rule to a comparison the registration does not describe: a reconstruction
against the deposited records.

Measured over all 795 drugs and all three invariant quantities:

| comparison | flat absolute 1e-06 | elementwise float32 | worst normalized residual |
|---|---|---|---|
| retired extraction against deposited | holds, 0 over | holds | 0.446 |
| authoritative rebuild against deposited | fails, 1 over | holds | 0.421 |

The frozen tolerance is correct for the comparison it was registered for and is
kept. Amendment 5 adds the elementwise float32 rule Amendment 2 already registers,
for the comparison Amendment 3 created.

### What this does not change

No registered hypothesis, statistic, criterion, interval, permutation, seed or
reading changes. No drug receives a tolerance of its own, no drug leaves the cohort,
and levofloxacin is not examined as a drug: it is neither the largest `P_shrna` in
the cohort nor the largest error in ULPs, and there is no reading under which its
identity matters.

### An earlier overstatement, corrected

The round-23 review packet said a flat absolute tolerance at the top of `P_shrna`'s
range is "unsatisfiable". That is too strong: float64 quantities derived from
identical float32 inputs can agree exactly. The defensible statement is narrower,
and is the one Amendment 5 makes: the flat absolute tolerance is incompatible with
the numerical equivalence expected between the deposited and the authoritative
routes, and it rejected the unperturbed verified computation for one of 795 drugs.

The same packet said no R0-R7 statistic exists, without qualification. R0.4 and R0.6
were recomputed under Amendment 4. What holds is that no R0.5 or R1-R7 association,
bootstrap, permutation or criterion outcome exists.

Amendment 5's draft described the legacy measurement as pinned while the gate
consuming it read whatever path it was handed, trusted three booleans, and computed
the file's hash after parsing and only recorded it. A file supplying those three
booleans would have satisfied the legacy layer, which is the defect Amendment 4's
own external-pinning rule exists to prevent. The audit is now hashed and compared
with a pin held outside it before it is parsed, and it carries the inputs it was
measured on rather than aggregate summaries alone.

### Records

- `results/03d_h3_reference_discordance/reproduction_tolerance_diagnosis.json`
- `results/03d_h3_reference_discordance/reproduction_tolerance_power.json`
- `results/03d_h3_reference_discordance/reproduction_error_localized.json`
- `results/03d_h3_reference_discordance/two_layer_gate_measurement.json`
- `experiments/PREREG_H3_S1S3_AMENDMENT_5.md`
- `experiments/PREREG_H3_S1S3_DIAGNOSTIC_D9_INDEXING.md`, frozen at `e684d16` before
  the diagnostic ran
- `results/03d_h3_reference_discordance/d9_indexing_diagnostic.json`
- `results/03d_h3_reference_discordance/legacy_reproduction_audit.json`, with
  `legacy_reproduction_per_drug.json` behind it and
  `registry/frozen/reproduction_legacy_pin.json` holding its expected identity
- `results/03d_h3_reference_discordance/route_to_route_all_drugs.json`
- `experiments/PREREG_H3_S1S3_DIAGNOSTIC_ENVIRONMENT.md`, the frozen plan for the
  library-version candidate. Its stage was added to `experiments/modal_h3_execute.py`
  after Amendment 5 pinned that file, so that pinned digest describes the file at the
  freeze and not afterwards; the reason is stated in the plan.

---

## Deviation 15: R4 did not compute the registered essential-minus-other difference

**Pre-registered:** "Split the genome-wide arm into essential-screen targets
(perturbed in the K562-essential file) and other targets. For each subset, compute
rho(P, E) and rho(D, E) under C1-GW with TCB intervals. **Estimate the difference
essential - other with a TCB that resamples targets within each subset
independently.**" (`PREREG_H3_REFERENCE_DISCORDANCE.md`, frozen `7f57136`, R4.)

**Actual:** the first complete R0-R7 execution, 2026-10-07, computed both subsets
with their intervals and no difference. There is no difference estimate anywhere in
`R4_essential_composition`.

**Why it matters:** both R4 rows of the interpretation grid turn on that comparison.
"R4: raw association depends on essential composition" reads as "the raw association
is specific to essential-screen targets", and "R4: H3 pattern outside the essential
set" reads as "H3's pattern replicates against a CRISPRi reference outside the
essential-gene set". Neither can be evaluated from two subset intervals read side by
side, which is the error an interval on the difference exists to prevent. Without it,
two registered readings go unreported.

**How it was identified:** reading the completed run's modules against the
registration, module by module, after the execution finished.

**What was done:** `independent_comparison` computes it as the registration specifies.
The subsets share no drug and no target, so `paired_comparison`, which resamples one
cohort and carries both estimates through the same replicate, does not apply. Each
subset is resampled over its own targets in its own stream and replicate `i` of one is
differenced against replicate `i` of the other. The difference is reported with the
registered sign, essential minus other. The reading follows the frozen comparison rule
with the essential subset as the reference, because the registered question is whether
the association survives outside the essential set.

**Why this is not an extension of the registration:** the quantity was registered at
`7f57136`, before any outcome existed. Computing it now executes the plan. Nothing was
added to the interpretation grid and no criterion was written after seeing a result;
Perplexity round 25 ruled against adding any, and this adds none.

**When:** after the first complete execution, before any R4 reading was reported.

## Deviation 16: the registration's "Reversal" is unsatisfied in the corrected data

**Pre-registered:** the follow-up is titled "Why does H3's pattern reverse between the
shRNA and CRISPRi references?" and defines its central term under Terms: "'The
qualitative pattern reverses' means that the projected association is positive only
under shRNA, and the raw association is positive only under CRISPRi."
(`PREREG_H3_REFERENCE_DISCORDANCE.md`, frozen `7f57136`.) Its Foreknowledge table
records the CRISPRi K562 arm under construction C0 at rho(P, E) = -0.1247 and
rho(D, E) = +0.2723, cited to `results/03b_h3_crispri/h3_crispri_results.json`,
sha256 `f00f8428`.

**Actual:** those two values are reproduced exactly from that file, which still hashes
to `f00f8428`. But the file predates the gene-axis mislabeling corrected under
Deviation 12 on 2026-10-05, two weeks after this registration was frozen on
2026-09-22. Recomputed on the corrected gene axis, the same C0 construction gives
rho(P, E) = +0.5410 and rho(D, E) = -0.1920 (R0.9 headline, run
`r0_to_r7_2026_10_07b`). Both signs are opposite to the frozen Foreknowledge.

The projected association is therefore positive under the shRNA reference (+0.3172),
under C0 (+0.5410), under C1 (+0.6074) and in RPE1 (+0.6090 to +0.7702). The raw
association is positive under none of them. By the registration's own definition, the
reversal is absent: the projected association is not positive *only* under shRNA, and
the raw association is not positive under CRISPRi at all.

**Why it matters:** three interpretation-grid rows name the reversal as their subject
("R0.9: raw reversal construction-dependent", "raw reversal replicates in RPE1", "raw
reversal is context-dependent"). Their premise term is unsatisfied, so the rows cannot
be read as written. The modules still return their mechanical readings, and those
readings are correct about what they measure: R0.9 returns "retained" for both halves
and R6a returns "retained", meaning the corrected pattern is stable across
constructions and across cell context. What is stable is a positive projected
association and a negative raw one, not a reversal. Reporting "the reversal is not
specific to K562" from R6a would assert the existence of something this data does not
contain.

**How it was identified:** recomputing rho(P, E) and rho(D, E) directly from the
pinned Foreknowledge file and comparing against the run's C0 headline, while placing
the completed readings against the grid row by row.

**What was done:** no registration is edited and no criterion is added or removed. The
three rows naming the reversal are recorded here as inapplicable, with the reason, and
the readings they would have produced are reported as what the modules measure:
construction-robustness and context-robustness of the corrected pattern. The
registered question the follow-up was built to answer does not arise, which is an
outcome of the analysis and is reported as one.

**When:** 2026-10-07, after the second complete execution, before any reading was
written into a manuscript.

## Deviation 17: the S1-S3 gate pinned two artifacts of a retired extraction

**Pre-registered:** `PREREG_H3_MAGNITUDE_AND_SHARED_AXIS.md` (frozen `f288507`) and the
amendment at `7f57136` pin the compound and shRNA extractions as
`lincs_subset.npz` sha256 `2ad0f5d3` and `lincs_shrna.npz` sha256 `4a990e50`.
`experiments/03c_h3_sensitivity.py` refuses unless the rebuild manifest records
both, and records the deposited artifact, at the top level and inside the
fingerprint.

**Actual:** the first attempt to run S1-T, S2-T and S3-T refused, correctly, with
`rebuild manifest is missing ['lincs_subset_sha256', 'lincs_shrna_sha256',
'deposited_sha256']`. The rebuild that produced the bundle on the `di-h3` volume
reads GSE92742 directly: `stage_fetch` downloads the GCTX from GEO and
`_fingerprint` pins `gctx`, `sig_info`, `shrna_sig_info`, `deposited`,
`landmark_order`, `loader` and `extract_code`. It never opens either `.npz`. Those
two files are products of the extraction retired when the gene-axis defect was
corrected under Deviation 12.

The locally tracked bundle is not an alternative. Its manifest pins
`deposited_sha256` `65e5d10e`, the artifact superseded by the correction, so that
bundle predates the fix.

**Why it matters:** the gate could be satisfied in three ways and two of them are
dishonest. Making the rebuild record hashes for files it does not read would
assert a dependency the run does not have. Hand-patching the manifest would assert
provenance nobody derived, which is the failure Amendment 4 exists to prevent.
Running on the local bundle would analyze pre-correction data.

**What was done:** the gate now pins what the rebuild actually reads. `deposited`
against the corrected `fd69e26f`, `sig_info` against `19da29c0` and
`shrna_sig_info` against `bd396fa0`, all three registered in advance, and the
presence of `landmark_order`, `loader` and `extract_code`. The GCTX is pinned to
`b293f3fb` and named in the source as trust on first use rather than a registered
pin, because `results/03c_h3_sensitivity/input_pin_check.json` records it as
"stamped on first retrieval, not registered in advance".

No inference criterion, seed, statistic or cohort rule is touched. The quantity the
gate protects is unchanged: the analysis still runs only against the corrected
artifact. What changed is which upstream objects the receipt names, because the
pipeline that writes the receipt was replaced.

**Why not regenerate the manifest instead:** `_fingerprint` includes
`extract_code`, the hash of the rebuild's own source, so any edit to the rebuild
invalidates every cached shard and forces a full re-read of GSE92742. That cost
would be acceptable on its own, but the rebuild would still have to name files it
does not read, so the expensive route does not buy an honest manifest.

**How it was identified:** the run refused on its own gate; the volume manifest was
retrieved and compared against the script's contract and against the tracked copy.

**When:** 2026-10-08, before any S1-T, S2-T or S3-T statistic existed.
