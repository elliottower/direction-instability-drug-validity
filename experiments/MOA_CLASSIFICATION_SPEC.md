# MOA Classification Spec for Experiment 13f (Broad Definition)

**Status: FROZEN — committed before 13f_v2 execution. Do not modify after commit.**

**Purpose:** Define a principled classification of LINCS drugs into "machinery"
vs "receptor-mediated" categories for the post-hoc MOA classification AUROC
(Experiment 13f). This is broader than the narrow keyword list in Experiment 14
(pre-registered, JUMP-CP) and Experiment 13d (Cohen's d bootstrap). Those
experiments retain their original narrow definitions unchanged.

## Single classification principle

**Machinery** = drugs whose primary target is core genomic, epigenetic,
proteostatic, or cytoskeletal machinery — components of DNA replication/repair,
chromatin modification, cell-cycle regulation, protein degradation (proteasome),
or cytoskeletal dynamics. These are constitutively expressed structural and
enzymatic components present in essentially all dividing human cells.

**Receptor-mediated** = drugs whose primary target is a cell-surface receptor
(GPCR, receptor tyrosine kinase acting as receptor) or ion channel. These
targets are differentially expressed across cell types.

**Everything else is excluded.** No case-by-case discretionary exclusions. If a
target does not fall into either category above, it is automatically excluded.
This includes intracellular enzymes that are not part of core
genomic/proteostatic/cytoskeletal machinery (COX, PDE, MAO, AChE, carbonic
anhydrase, HMGCR, CYP450), nuclear hormone receptors (see below), bacterial
targets, transporters (reuptake inhibitors), and unclassifiable mechanisms.

### Nuclear hormone receptors: excluded

Nuclear hormone receptors (glucocorticoid, estrogen, androgen, progesterone,
PPAR, retinoid, vitamin D) are intracellular ligand-activated transcription
factors that act directly on chromatin. By strict intracellular/surface
classification, they belong on the machinery side. However, their expression is
tissue-specific (estrogen receptor is absent from many cell types), creating a
genuine edge case: mechanistically intracellular, phenotypically
context-dependent. Rather than making a discretionary call that could be tuned to
the result, we exclude them. A sensitivity analysis will report the AUROC with
nuclear receptors folded into machinery and separately into receptor, confirming
the result is robust to this boundary.

### Kinase inhibitors: separate third category

Kinase inhibitors are held out as a separate intermediate category, matching the
pre-registered Experiment 14 structure. Some kinases are constitutive (CDK),
some are receptor-like (EGFR). The binary AUROC is computed as machinery vs
receptor with kinases excluded; the ordinal Spearman uses all three groups
(machinery < kinase < receptor).

## Keyword lists

These lists follow mechanically from the principle above.

### Machinery (core genomic / epigenetic / proteostatic / cytoskeletal)

**DNA replication and repair:**
- DNA synthesis inhibitor
- DNA alkylating agent
- ribonucleotide reductase inhibitor
- topoisomerase inhibitor
- PARP inhibitor
- thymidylate synthase inhibitor
- dihydrofolate reductase inhibitor

**Chromatin and epigenetic modification:**
- HDAC inhibitor
- histone methyltransferase inhibitor
- histone demethylase inhibitor
- DNA methyltransferase inhibitor
- BET inhibitor
- sirtuin inhibitor

**Cell cycle and proliferation:**
- CDK inhibitor
- cell cycle inhibitor
- MDM inhibitor

**Proteostasis:**
- proteasome inhibitor
- ATPase inhibitor
- autophagy inhibitor
- HSP inhibitor

**Cytoskeleton:**
- tubulin inhibitor
- tubulin polymerization inhibitor

**Transcription and translation:**
- RNA synthesis inhibitor
- protein synthesis inhibitor
- mTOR inhibitor

**Signaling cascades acting on core machinery:**
- NFkB pathway inhibitor

### Receptor-mediated (cell-surface receptors and ion channels)

**GPCRs and related:**
- dopamine receptor (agonist / antagonist)
- serotonin receptor (agonist / antagonist)
- histamine receptor (agonist / antagonist)
- adrenergic receptor (agonist / antagonist)
- acetylcholine receptor (agonist / antagonist)
- glutamate receptor (agonist / antagonist)
- opioid receptor (agonist / antagonist)
- cannabinoid receptor (agonist / antagonist)
- adenosine receptor (agonist / antagonist)
- angiotensin receptor antagonist
- prostanoid receptor (agonist / antagonist)
- tachykinin antagonist
- leukotriene receptor antagonist
- endothelin receptor antagonist
- benzodiazepine receptor agonist

**Ion channels:**
- calcium channel blocker
- T-type calcium channel blocker
- sodium channel blocker
- potassium channel blocker
- potassium channel activator
- GABA receptor (agonist / antagonist)

### Kinase inhibitors (separate intermediate category)

- EGFR inhibitor
- MEK inhibitor
- MAPK inhibitor / p38 MAPK inhibitor
- SRC inhibitor
- PI3K inhibitor
- JAK inhibitor
- FLT3 inhibitor
- VEGFR inhibitor
- PDGFR inhibitor
- ABL inhibitor
- RAF inhibitor
- PKC inhibitor
- AKT inhibitor
- glycogen synthase kinase inhibitor
- IKK inhibitor
- rho associated kinase inhibitor
- tyrosine kinase inhibitor
- checkpoint kinase inhibitor
- ALK inhibitor
- BTK inhibitor
- KIT inhibitor

### Automatically excluded (by the principle, not by discretion)

The following do not meet either definition and are excluded mechanically:

- **Nuclear hormone receptors**: glucocorticoid, estrogen, androgen,
  progesterone, PPAR, retinoid, vitamin D, mineralocorticoid receptor
  agonists/antagonists (intracellular transcription factors — edge case
  documented above)
- **Intracellular enzymes outside core machinery**: cyclooxygenase inhibitor,
  phosphodiesterase inhibitor, monoamine oxidase inhibitor,
  acetylcholinesterase inhibitor, carbonic anhydrase inhibitor, HMGCR inhibitor,
  cytochrome P450 inhibitor, nitric oxide synthase inhibitor
- **Bacterial/fungal/viral targets**: bacterial cell wall synthesis inhibitor,
  bacterial DNA gyrase inhibitor, bacterial ribosomal subunit inhibitors, sterol
  demethylase inhibitor, HIV protease inhibitor, reverse transcriptase inhibitor
- **Transporter modulators**: SSRI, dopamine reuptake inhibitor, norepinephrine
  reuptake inhibitor
- **Unclassifiable**: local anesthetic, anthelmintic agent, progestogen hormone

## Priority order for compound MOA strings

LINCS MOA annotations are pipe-separated (e.g., "CDK inhibitor|cell cycle
inhibitor|MCL1 inhibitor"). Classification uses the first match in priority
order: **machinery > kinase > receptor > excluded**. This ensures specific
targets take precedence over vague annotations.

## Relationship to other definitions in the paper

| Analysis | Definition | Machinery n | Receptor n | Status |
|----------|-----------|-------------|------------|--------|
| Exp 14 (JUMP-CP) | Narrow, pre-registered | 19 (JUMP-CP) | 107 (JUMP-CP) | Frozen (commit 8d74bd0) |
| Exp 13d (Cohen's d) | Narrow, post-hoc | 114 (LINCS) | 481 (LINCS) | Completed, d=1.94 in paper |
| Exp 13f v2 (MOA AUROC) | Broad, post-hoc | TBD | TBD | **This spec** |

The narrow definition is a subset of the broad. The broad adds PARP, BET,
histone modifiers, DNA methyltransferase, sirtuin, HSP, autophagy, mTOR (if not
already matched), and additional cell-surface receptors (adenosine, angiotensin,
prostanoid, tachykinin, leukotriene, endothelin). The paper reports both and
notes their distinct purposes:
- Narrow (13d): conservative replication of the drug paper's per-class examples
- Broad (13f): principled test using the core-machinery-vs-surface-receptor rule

## Primary analysis and robustness commitment

**Primary analysis: nuclear receptors excluded.** This is the base case. The
AUROC reported as the headline 13f v2 result uses machinery (n=180) vs receptor
(n=586) with nuclear receptors, non-core intracellular enzymes, kinases,
bacterial targets, and transporters all removed from the evaluation set.

**Robustness requirement: DI > Frechet must hold in all three variants.** The
result is considered robust only if DI outperforms Frechet variance on the MOA
classification task under all three nuclear-receptor treatments:
1. Nuclear receptors excluded (primary)
2. Nuclear receptors folded into machinery
3. Nuclear receptors folded into receptor

If any variant shows Frechet >= DI, the robustness claim is withdrawn and the
result is reported as "directionally supportive but sensitive to the
nuclear-receptor boundary."

## Verified counts (pre-run)

Applying the broad classification to the LINCS 8,949-drug CSV:
- Machinery: 180
- Kinase: 153
- Receptor: 586
- Excluded (auto, by the rule): 863
- No MOA annotation: 7,167
- Total: 8,949

Aurora kinase inhibitors (n=11): all classified as kinase (Aurora kinase
inhibitor appears only in the kinase keyword list). No dual-match ambiguity.

Nuclear hormone receptors: all classified as excluded (verified for
glucocorticoid n=34, estrogen n=27, androgen n=10, progesterone n=15, PPAR
n=23+1 compound, retinoid n=12, vitamin D n=6, mineralocorticoid n=4). One PPAR
drug (oleoylethanolamide) classified as receptor due to compound MOA including
cannabinoid receptor — correct behavior under the priority rule.

## Implementation assertions

The 13f_v2 script must include the following runtime assertions:
1. Classes are disjoint: no drug appears in more than one of {machinery, kinase, receptor}
2. Excluded drugs are removed from the AUROC evaluation set (not silently relabeled)
3. machinery_n + kinase_n + receptor_n + excluded_n == total annotated drugs (1,782)
4. Binary AUROC computed on machinery-vs-receptor only; kinases held out
