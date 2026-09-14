"""Build the committed source table behind Figure 5.

Figure 5 previously hardcoded six direction-instability values and six
cell-line counts that reproduce no stored artifact. This script derives that
table from `results/01_toxicity_failure/toxicity_results.json`, the file this
paper's own toxicity-correction experiment writes, and records the selectivity
categories with the reasoning for the grouping.

Selectivity is recorded as three ordered categories rather than six untied
ranks. The three pan-HDAC inhibitors differ in potency rather than in isoform
breadth, and the relative selectivity of tubacin (HDAC6) and PCI-34051 (HDAC8)
is assay-dependent, so neither pair supports an internal ordering.

Usage:
    uv run --no-project --with numpy python experiments/build_fig5_source.py
"""
import json
from itertools import combinations
from math import sqrt
from pathlib import Path

REPO = Path("/Users/elliottower/Documents/GitHub/direction-instability-drug-validity")
SOURCE = REPO / "results/01_toxicity_failure/toxicity_results.json"
OUT_DIR = REPO / "results/fig5_hdac_source"

# drug key in the results file, display name, selectivity category, ordered rank
PANEL = [
    ("panobinostat", "Panobinostat", "Pan-HDAC", 1),
    ("vorinostat",   "Vorinostat",   "Pan-HDAC", 1),
    ("belinostat",   "Belinostat",   "Pan-HDAC", 1),
    ("entinostat",   "Entinostat",   "Class I-selective", 2),
    ("tubacin",      "Tubacin",      "Isoform-selective", 3),
    ("pci-34051",    "PCI-34051",    "Isoform-selective", 3),
]


def kendall_tau_b(categories, values):
    n = len(values)
    concordant = discordant = cross = 0
    for i, j in combinations(range(n), 2):
        if categories[i] == categories[j]:
            continue
        cross += 1
        s = (categories[j] - categories[i]) * (values[j] - values[i])
        concordant += s > 0
        discordant += s < 0
    n0 = n * (n - 1) // 2
    sizes = {}
    for c in categories:
        sizes[c] = sizes.get(c, 0) + 1
    ties_x = sum(k * (k - 1) // 2 for k in sizes.values())
    return {
        "cross_category_pairs": cross,
        "concordant": concordant,
        "discordant": discordant,
        "tau_b": round((concordant - discordant) / sqrt((n0 - ties_x) * n0), 4),
    }


def main():
    records = {str(r["drug"]).lower(): r for r in json.loads(SOURCE.read_text())}

    rows = []
    for key, display, category, rank in PANEL:
        r = records[key]
        rows.append({
            "drug": display,
            "results_key": key,
            "selectivity_category": category,
            "category_rank": rank,
            "corrected_instability": round(r["corrected_instability"], 4),
            "raw_instability": round(r["raw_instability"], 4),
            "n_celllines": r["n_celllines"],
        })

    stats = kendall_tau_b(
        [r["category_rank"] for r in rows],
        [r["corrected_instability"] for r in rows],
    )
    # exact two-sided p by enumerating category-label assignments
    from itertools import permutations
    cats = [r["category_rank"] for r in rows]
    vals = [r["corrected_instability"] for r in rows]
    observed = stats["tau_b"]
    seen, extreme, total = set(), 0, 0
    for perm in permutations(range(6)):
        key = tuple(cats[i] for i in perm)
        if key in seen:
            continue
        seen.add(key)
    for key in seen:
        total += 1
        if abs(kendall_tau_b(list(key), vals)["tau_b"]) >= abs(observed) - 1e-9:
            extreme += 1
    stats["exact_two_sided_p"] = round(extreme / total, 4)
    stats["n_distinct_category_assignments"] = total

    out = {
        "figure": "fig5_hdac",
        "plotted_quantity": "corrected_instability",
        "plotted_quantity_note": (
            "Toxicity-gene-corrected direction instability, computed in this "
            "repository by experiments/01_toxicity_failure. Lower means the "
            "response points the same direction across cell lines."
        ),
        "selectivity_note": (
            "Three ordered categories. The three pan-HDAC inhibitors differ in "
            "potency rather than isoform breadth; tubacin and PCI-34051 are both "
            "isoform-selective and their relative selectivity is assay-dependent. "
            "Neither pair supports an internal ordering, so both are tied."
        ),
        "source_file": str(SOURCE.relative_to(REPO)),
        "rows": rows,
        "statistics": stats,
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "fig5_hdac_source.json").write_text(json.dumps(out, indent=2))

    for r in rows:
        print(f"{r['drug']:14s} {r['selectivity_category']:18s} "
              f"D={r['corrected_instability']:.3f}  n={r['n_celllines']:2d}")
    print(f"\nKendall tau_b = {stats['tau_b']}  "
          f"({stats['concordant']}/{stats['cross_category_pairs']} concordant, "
          f"exact two-sided p = {stats['exact_two_sided_p']})")
    print(f"\nWritten to {OUT_DIR / 'fig5_hdac_source.json'}")


if __name__ == "__main__":
    main()
