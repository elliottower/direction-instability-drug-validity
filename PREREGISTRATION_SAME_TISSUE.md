# Pre-registration: Same-Tissue Direction Instability (Experiment 9)

**Timestamp:** 2026-07-05T18:15:36Z
**Status:** FROZEN — do not edit after running Experiment 9

## Hypothesis

When direction instability is computed only across cell lines from the same
tissue, lineage-shared transcriptomic programs create false consistency:
drugs appear mechanism-conserving because their target cells respond
similarly due to shared biology, not because the drug has a conserved
mechanism. A "lineage correction" (removing genes that differentiate tissue
subtypes) should reorganize the ranking more than the all-tissue toxicity
correction did (rho = 0.83).

## Prediction (H_tissue)

For at least one tissue with >= 5 cell lines (large intestine, lung,
heme, breast, ovary, skin), Spearman rho between raw and
lineage-corrected direction instability will fall below 0.70 ("bites").

## Decision criteria (same as all prior experiments)

- rho < 0.70: "bites" — correction substantially reorganizes ranking
- 0.70 <= rho < 0.83: "marginal"
- rho >= 0.83: "doesn't bite"

## Method

1. For each tissue with >= 5 cell lines in LINCS L1000:
   - Restrict to drugs tested in ALL cell lines of that tissue
   - Compute direction instability using only within-tissue cell lines
   - Identify "lineage genes": genes with high variance across cell lines
     of that tissue in DMSO controls (top 100 by variance)
   - Compute lineage-corrected instability after removing lineage genes
   - Report Spearman rho (raw vs corrected)

2. Secondary analysis: compare within-tissue direction instability
   distributions to across-all-tissues direction instability for the same
   drugs. If within-tissue D is systematically lower, lineage similarity
   is creating apparent consistency.

## Tissues to test

| Tissue | N cell lines | Cell lines |
|--------|-------------|------------|
| Large intestine | 14 | CL34, HCT116, HT115, HT29, LOVO, MDST8, NCIH508, NCIH716, RKO, SNU1040, SNUC5, SW480, SW620, SW948 |
| Lung | 12 | A549, CORL23, DV90, H1299, HCC15, HCC515, NCIH1694, NCIH1836, NCIH2073, NCIH596, SKLU1, T3M10 |
| Heme | 8 | HL60, JURKAT, NOMO1, PL21, SKM1, THP1, U937, WSUDLCL2 |
| Breast | 6 | BT20, HS578T, MCF10A, MCF7, MDAMB231, SKBR3 |
| Ovary | 6 | COV644, EFO27, OV7, RMGI, RMUGS, TYKNU |
| Skin | 5 | A375, FIBRNPC, MCH58, SKMEL1, SKMEL28 |

## What would make this informative

- If rho < 0.70 in any tissue: the validity ladder bites when contexts
  share biology, proving the ladder is not vacuous.
- If rho >= 0.83 in all tissues: same-tissue restriction doesn't create
  sufficient confounding, further strengthening the "ranking is robust"
  thesis.
- Either outcome is publishable.
