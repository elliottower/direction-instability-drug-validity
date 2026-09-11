# Superseded: pre control-exclusion artifacts

**Do not use for scientific reporting.**

These were computed on a cohort that wrongly retained eight plate-level positive
control compounds (`JCP2022_012818`, `_025848`, `_035095`, `_037716`, `_046054`,
`_050797`, `_064022`, `_085227`). Each of those sits on 93.6-100% of eligible
COMPOUND plates and is an experimental control rather than a screened compound,
so their inclusion violates the intended unit population and inflates replicate
counts by three orders of magnitude against the cohort median.

Retained for provenance only. The corrected artifacts sit one directory up and
are the only ones the manuscript reports.

| | |
|---|---|
| cohort here | 25,000 compounds, controls included |
| corrected cohort | 24,992 compounds, controls excluded |
| exclusion rule | compound present on >= 50% of eligible COMPOUND plates |
| rule evidence | `results/qc_control_audit/control_audit.json` |

Hashes in `SHA256SUMS`. Produced before any commit of this work; the repository's
last commit at the time was `8b62908` (2026-07-09).
