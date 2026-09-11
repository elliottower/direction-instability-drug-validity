# Experiment 7 / HP1 artifacts

Recovered from git history on 2026-09-11. Both files were committed in
`4a7f3fe` (2026-07-05) and later removed from the working tree; the generating
script, `experiments/07_perturbseq_correction.py`, was frozen in `249abaf` and
deleted in `8c64d56` (2026-07-07) and is recoverable the same way.

`corrected_distances.npz` holds `raw_dists`, `corr_dists` and `ess_mask` for
1,676 Perturb-seq knockdowns. Recomputed from it:

| quantity | value |
|---|---|
| genes | 1,676 |
| Spearman rho, raw vs. essential-corrected geodesic | 0.9066 |
| maximum absolute rank shift | 903 positions |

These back the HP1 row of the scorecard and the Perturb-seq bars of the
cross-domain figure, which reads this file directly rather than carrying
literals.
