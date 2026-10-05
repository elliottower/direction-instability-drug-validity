# Provenance: the gene-axis defect, and the audit of what it reaches

Deviation 11 and Deviation 12 in `DEVIATION_LOG.md` are the registration record.
This file carries the cross-repository consequence and the paper audit, which belong
to neither registration.

## Where the defect comes from

Both LINCS extractions this project pins, `lincs_subset.npz` and `lincs_shrna.npz`,
were written by `drug-perturbation-geometry`. Three sites there parse the pinned
GCTX with a requested row order and label the values with the request rather than
with the file, and nothing in that repository calls `reindex`. cmapPy ignores the
`rid` order it is given and returns the file's own, measured in
`results/03c_h3_sensitivity/reader_agreement.json`. Both files therefore carry the
same permutation of the gene axis: correct values, wrong labels.

`experiments/modal_h3_rebuild.py:218` is the one line that avoids it, which is why
the rebuild passes the gate and the extractions do not.
`drug-perturbation-geometry/docs/provenance.md` carries that side of the story.

## The audit of this project's papers

Every gene symbol in every paper version was classified by where its identity comes
from, because the two sources are affected differently:

- **a drug-target annotation**, from `frozen_drug_labels.json` — unaffected, the
  gene axis plays no part
- **an index into a signature column** — affected, the label is wrong

| paper | gene identities | verdict |
|---|---|---|
| `direction_instability_confound_audit_v1b.tex` (current) | PSMA3, PSMB1, PSMB6, PSMC4, TOP2A, SF3B1, SNRPG, DDX23 | **clean**: the rank-shift analysis ranks the 1,676 Perturb-seq knockdown genes (`:163`, `:429`), whose identities come from the h5ad's own labels, and the remainder are drug-target annotations |
| `direction_instability_zenodo.tex` (deposited) | the same set | **clean**, same reasoning |
| `direction_instability_bmcmrm_v1-v8`, `bioinfadv_v1`, `bmcbioinfo_v1`, `psb_v3`, `psb_v4` | target annotations only | **clean** |
| `direction_instability_v1-v8`, `psb.tex`, `psb_v2.tex` (superseded drafts) | SUV39H1, MYC, CDK6 named as "chromatin remodeling genes" of a conserved signature (`v8:170-173`) | **affected**, and the same claim as the drug-transport paper. Superseded, so no correction is owed; recorded so the sentence is not revived |

A scan by symbol alone is not sufficient and was not relied on. BIRC5, CDK6 and MYC
are both landmark genes and annotated drug targets, so membership in either list
says nothing about how a sentence uses them; each occurrence was read in context.

## Open

The rank-shift analysis names Perturb-seq genes, whose labels are sound. Whether the
*quantity* being ranked is computed in the LINCS landmark space, and so inherits the
permutation even though the names do not, has not been established. It does not
affect the gene identities and it may affect the shift magnitudes.

The CRISPRi arms are the live question and are tracked in Deviation 12:
`03d:604` builds the symbol list from the extraction's declared ids and uses it to
place five external references, and `03d:950` selects R7f's removed coordinate by
gene name.
