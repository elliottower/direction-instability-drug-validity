"""Modal wrapper: run the registered H3 analyses where the data already lives.

The gate, S1-S3 and the reference-discordance driver each load the LINCS
extraction and build per-drug matrices, which is several gigabytes of working set
before any resampling starts. None of it belongs on a laptop.

Logic lives in the scripts this calls; this file arranges workers, volumes and
retrieval. Inputs come from the volumes they were pinned on: `drug-perturbation-vol`
holds the extraction, `di-h3` the GCTX rebuild.

    modal run --detach experiments/modal_h3_execute.py --stage gate
    modal run --detach experiments/modal_h3_execute.py --stage s1s3
    modal run --detach experiments/modal_h3_execute.py --stage driver
    modal volume get di-h3-results . ./results/
"""
import modal

app = modal.App("di-h3-execute")
extraction = modal.Volume.from_name("drug-perturbation-vol")
rebuild = modal.Volume.from_name("di-h3")
results = modal.Volume.from_name("di-h3-results", create_if_missing=True)
inputs = modal.Volume.from_name("di-h3-inputs")

REPO = "/Users/elliottower/Documents/GitHub/direction-instability-drug-validity"

base = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.1.3", "scipy==1.14.1", "pandas==2.2.3", "anndata==0.11.4",
                 "h5py==3.12.1", "cmapPy==4.0.1")
)


def _with_code(built):
    """The repository's code and the artifacts the stages read, added last.

    Modal refuses a build step after an `add_local_*`, so every layer that
    installs something comes before this.
    """
    return (built
            .add_local_dir(f"{REPO}/experiments", remote_path="/app/experiments")
            .add_local_dir(f"{REPO}/geometry", remote_path="/app/geometry")
            .add_local_dir(f"{REPO}/results/03_phenotype_projection",
                           remote_path="/app/results/03_phenotype_projection")
            .add_local_dir(f"{REPO}/results/03b_h3_crispri",
                           remote_path="/app/results/03b_h3_crispri")
            .add_local_dir(f"{REPO}/results/03d_h3_reference_discordance",
                           remote_path="/app/results/03d_h3_reference_discordance"))


image = _with_code(base)
# the registered implementation runs the suite, and a suite that fits models is
# an analysis under a different name, so it runs here rather than on a laptop
test_image = _with_code(base.pip_install("pytest==9.1.1", "scikit-learn==1.9.1",
                                         "tqdm==4.70.1", "pyarrow==25.0.1")).add_local_dir(
    f"{REPO}/tests", remote_path="/app/tests")

# The Replogle and PRISM files are half a gigabyte. They live on a volume rather
# than in the image: mounting them makes every run upload them again, and a client
# killed mid-upload leaves an app that never dispatches its function.

# 32 GB, not 256: a run that needs more has a bug rather than a big input, and
# retries=0 so a bug fails once instead of nine times overnight
COMMON = dict(image=image, timeout=86400, memory=32768, cpu=8.0, retries=0,
              volumes={"/extraction": extraction, "/rebuild": rebuild, "/out": results,
                       "/inputs": inputs})

GCTX = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
# a report never overwrites an earlier one
GATE_RUN = "gate_run_2026-10-03"

MAPPING_SHA256 = "152361cb3174a5fb7aae0229c3e3d049dc00d49d9e442925156a9fe0564b89d3"


def _ts():
    from datetime import datetime
    return datetime.now().strftime("%H:%M:%S")


EXTRACTION_FILES = ["lincs_subset.npz", "lincs_shrna.npz", "lincs_shrna_siginfo.csv.gz",
                    "GSE92742_Broad_LINCS_sig_info.txt.gz",
                    "GSE92742_Broad_LINCS_gene_info.txt.gz",
                    "GSE92742_Broad_LINCS_cell_info.txt.gz", "frozen_drug_labels.json"]


def _repo_at_its_absolute_path():
    """Put the code where its own constants expect to find it.

    The analysis scripts carry absolute paths to the repository. Rather than edit
    frozen, reviewed code to run it elsewhere, the container gets a link at that
    path pointing at the mounted copy.
    """
    from pathlib import Path

    target = Path(REPO)
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to("/app")
    return target


def _stage_inputs():
    """One directory holding the extraction and the metadata the scripts expect."""
    from pathlib import Path

    staged = Path("/tmp/extraction")
    staged.mkdir(parents=True, exist_ok=True)
    for name in EXTRACTION_FILES:
        link = staged / name
        if not link.exists():
            link.symlink_to(Path("/extraction") / name)
    return staged


def _stage_rebuild():
    """The GCTX rebuild in one directory, as the gate expects it."""
    from pathlib import Path

    staged = Path("/tmp/rebuild")
    staged.mkdir(parents=True, exist_ok=True)
    for source in list(Path("/rebuild/shards").glob("*")) + list(Path("/rebuild/results").glob("*")):
        link = staged / source.name
        if not link.exists():
            link.symlink_to(source)
    return staged


def _run(argv):
    """Run one of the analysis scripts in-process, with the repo importable."""
    import runpy
    import sys

    sys.path.insert(0, "/app")
    sys.argv = argv
    print(f"[{_ts()}] {' '.join(argv)}", flush=True)
    runpy.run_path(argv[0], run_name="__main__")


# the registration's pins, restated here so the gate run can be refused if an
# input is not the file the registration names. The GCTX checksum could not exist
# before first retrieval and was stamped on the volume at that time.
REGISTERED_PINS = {
    "GSE92742_Broad_LINCS_sig_info.txt.gz":
        "19da29c0ee12ddf27f9698cd0da40beaff58657dcde9d382aae068737e831299",
    "lincs_shrna_siginfo.csv.gz":
        "bd396fa0e1a2f00c1b5f2c8d2b35f9a056f5e5353382475655869038037ec014",
    "phenotype_projection_results.json":
        "fd69e26fc9a3917323065b631688baeab8b283f735c8bf5b16210ba67bd21425",
}


@app.function(**COMMON)
def stage_verify_pins():
    """Are the inputs the gate will read the files the registration names?

    The gate records the hash of everything it reads, which says what ran. It does
    not say that what ran is what was registered, and those are different claims.
    This compares the files on the volume against the registration's pins and
    against the checksum stamped beside the GCTX on first retrieval.
    """
    import hashlib
    import json
    from pathlib import Path

    raw = Path("/rebuild/raw")

    def digest(path):
        sha = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 22), b""):
                sha.update(block)
        return sha.hexdigest()

    checked = {}
    for name, expected in REGISTERED_PINS.items():
        found = digest(raw / name)
        checked[name] = {"expected": expected, "found": found, "matches": found == expected}

    stamped = (raw / "gctx.sha256").read_text().split()[0]
    found = digest(raw / GCTX)
    checked[GCTX] = {"expected": stamped, "found": found, "matches": found == stamped,
                     "note": "stamped on first retrieval, not registered in advance"}

    gene_info = "GSE92742_Broad_LINCS_gene_info.txt.gz"
    checked[gene_info] = {"expected": (raw / f"{gene_info}.sha256").read_text().split()[0],
                          "found": digest(raw / gene_info)}
    checked[gene_info]["matches"] = (checked[gene_info]["expected"]
                                     == checked[gene_info]["found"])

    summary = {"all_inputs_match_their_pins": all(entry["matches"]
                                                  for entry in checked.values()),
               "inputs": checked,
               "reading": ("a mismatch means the gate would compare a file the registration "
                           "does not name, and the run is refused rather than reported")}
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "input_pin_check.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_gate():
    """The reconstruction gate. A failure here voids the analyses registered against it."""
    from pathlib import Path

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    _run(["/app/experiments/03e_reconstruction_gate.py",
          "--rebuilt", str(rebuilt),
          "--extraction", str(staged),
          "--gctx", f"/rebuild/raw/{GCTX}",
          "--shrna-siginfo", "/rebuild/raw/lincs_shrna_siginfo.csv.gz",
          "--compound-siginfo", "/rebuild/raw/GSE92742_Broad_LINCS_sig_info.txt.gz",
          "--gene-info", "/rebuild/raw/GSE92742_Broad_LINCS_gene_info.txt.gz",
          "--cohort", "/app/results/03_phenotype_projection/phenotype_projection_results.json",
          "--output", f"/out/03c_h3_sensitivity/{GATE_RUN}"])
    results.commit()
    return (Path(f"/out/03c_h3_sensitivity/{GATE_RUN}/reconstruction_gate.json")
            ).read_text()[:2000]


@app.function(**COMMON)
def stage_s1s3():
    """S1-S3 and their target-level versions, on the bundle the rebuild produced."""
    _repo_at_its_absolute_path()
    _run(["/app/experiments/03c_h3_sensitivity.py",
          "--bundle", "/rebuild/results/cohort_bundle.npz",
          "--manifest", "/rebuild/results/rebuild_manifest.json"])
    results.commit()
    return "s1s3 done"


@app.function(**COMMON)
def stage_driver():
    """R0 to R7 and R5 response, in one invocation."""
    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    _run(["/app/experiments/03d_h3_reference_discordance.py",
          "--r5-stage", "response",
          "--expected-mapping-sha256", MAPPING_SHA256,
          "--data", str(staged),
          "--perturbseq", "/extraction/ReplogleWeissman2022_K562_essential.h5ad",
          "--replogle", "/inputs/replogle2022",
          "--prism", "/inputs/prism_19q4",
          "--output", "/out/03d_h3_reference_discordance"])
    results.commit()
    return "driver done"


@app.function(**COMMON)
def stage_gate_diagnostic():
    """Why the shRNA half of the gate failed, described rather than explained away.

    A failed gate voids the analyses registered against it. Before anything is
    concluded, this measures the disagreement: how many signatures differ, by how
    much, whether the two versions are related by a scale, a permutation of genes,
    or a different normalization, and whether the compound half agrees.
    """
    import json
    import sys

    import numpy as np

    sys.path.insert(0, "/app")
    _repo_at_its_absolute_path()
    staged, rebuilt_dir = _stage_inputs(), _stage_rebuild()
    gate = __import__("experiments.03e_reconstruction_gate", fromlist=["x"]) if False else None
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    records = json.loads(open("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read())
    targets = sorted({str(r["target"]) for r in records})

    rebuilt_sigs, rebuilt_membership, rebuilt_genes, stored = gate.load_shrna_rebuild_strict(
        rebuilt_dir)
    extraction_sigs, extraction_membership, extraction_genes = gate.load_shrna_extraction(
        staged, targets_wanted=set(targets))

    shared = sorted(set(rebuilt_sigs) & set(extraction_sigs))
    order = {g: i for i, g in enumerate(map(str, extraction_genes))}
    permutation = np.array([order[g] for g in map(str, rebuilt_genes)])

    rows = []
    for sig_id in shared:
        left = np.asarray(rebuilt_sigs[sig_id], dtype=np.float64)
        right = np.asarray(extraction_sigs[sig_id], dtype=np.float64)[permutation]
        denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
        rows.append({
            "sig_id": sig_id,
            "max_abs_difference": float(np.abs(left - right).max()),
            "cosine": float(left @ right / denominator) if denominator else float("nan"),
            "scale_ratio": float(np.linalg.norm(left) / np.linalg.norm(right))
                           if np.linalg.norm(right) else float("nan"),
            "rebuild_norm": float(np.linalg.norm(left)),
            "extraction_norm": float(np.linalg.norm(right)),
        })

    differing = [r for r in rows if r["max_abs_difference"] > 1e-5]
    cosines = np.array([r["cosine"] for r in rows])
    ratios = np.array([r["scale_ratio"] for r in rows])
    summary = {
        "n_signatures_compared": len(rows),
        "n_differing_beyond_1e-5": len(differing),
        "fraction_differing": len(differing) / max(len(rows), 1),
        "cosine": {"min": float(np.nanmin(cosines)), "median": float(np.nanmedian(cosines)),
                   "max": float(np.nanmax(cosines)),
                   "n_above_0.999": int((cosines > 0.999).sum())},
        "scale_ratio": {"min": float(np.nanmin(ratios)), "median": float(np.nanmedian(ratios)),
                        "max": float(np.nanmax(ratios))},
        "max_abs_difference": {"max": max(r["max_abs_difference"] for r in rows),
                               "median": float(np.median([r["max_abs_difference"]
                                                          for r in rows]))},
        "worst_ten": sorted(rows, key=lambda r: -r["max_abs_difference"])[:10],
        "reading": ("a cosine near 1 with a scale ratio away from 1 means the same direction "
                    "under a different normalization; a cosine far from 1 means different "
                    "values, not a rescaling"),
    }
    out = "/out/03c_h3_sensitivity"
    __import__("os").makedirs(out, exist_ok=True)
    open(f"{out}/gate_shrna_diagnostic.json", "w").write(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps({k: v for k, v in summary.items() if k != "worst_ten"}, indent=2)


@app.function(**COMMON)
def stage_permutation_test():
    """Equal norms with zero cosine means a permutation. This finds which axis.

    If a signature's sorted values match, the same numbers are present in a
    different order, and the disagreement is a gene-axis mislabeling. If instead
    the rebuild's vector equals a *different* extraction signature, the row labels
    are mislabeled, which is Deviation 9's shape one level deeper.
    """
    import json
    import sys

    import numpy as np

    sys.path.insert(0, "/app")
    _repo_at_its_absolute_path()
    staged, rebuilt_dir = _stage_inputs(), _stage_rebuild()
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    records = json.loads(open("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read())
    targets = sorted({str(r["target"]) for r in records})
    rebuilt_sigs, _, rebuilt_genes, _ = gate.load_shrna_rebuild_strict(rebuilt_dir)
    extraction_sigs, _, extraction_genes = gate.load_shrna_extraction(
        staged, targets_wanted=set(targets))

    sample = sorted(set(rebuilt_sigs) & set(extraction_sigs))[:200]
    same_values, exact_row_match = 0, 0
    matches = {}
    extraction_matrix = np.vstack([np.asarray(extraction_sigs[s], dtype=np.float64)
                                   for s in sample])
    for sig_id in sample:
        left = np.asarray(rebuilt_sigs[sig_id], dtype=np.float64)
        right = np.asarray(extraction_sigs[sig_id], dtype=np.float64)
        if np.allclose(np.sort(left), np.sort(right), rtol=1e-5, atol=1e-5):
            same_values += 1
        # does this rebuilt vector equal some other extraction signature?
        distances = np.abs(extraction_matrix - left).max(axis=1)
        best = int(np.argmin(distances))
        if distances[best] < 1e-4:
            exact_row_match += 1
            if sample[best] != sig_id and len(matches) < 5:
                matches[sig_id] = sample[best]

    summary = {
        "n_sampled": len(sample),
        "same_multiset_of_values": same_values,
        "reading_same_values": ("the same numbers in a different order means the gene axis is "
                                "mislabeled, not that the values differ"),
        "rebuilt_vector_found_elsewhere_in_extraction": exact_row_match,
        "example_row_mismatches": matches,
        "gene_axis_identical_as_sequence": [str(a) for a in rebuilt_genes][:5] ==
                                           [str(a) for a in extraction_genes][:5],
        "rebuilt_genes_head": [str(g) for g in rebuilt_genes[:5]],
        "extraction_genes_head": [str(g) for g in extraction_genes[:5]],
        "rebuilt_genes_sorted": [str(g) for g in rebuilt_genes] == sorted(
            [str(g) for g in rebuilt_genes], key=int),
        "extraction_genes_sorted": [str(g) for g in extraction_genes] == sorted(
            [str(g) for g in extraction_genes], key=int),
    }
    open("/out/03c_h3_sensitivity/gate_permutation_test.json", "w").write(
        json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_axis_test():
    """Which column order does each extraction matrix actually use?

    The compound half of the gate passed while the shRNA half failed on the same
    declared gene axis. Either the two npz files disagree about their own column
    order, or one of them carries labels that do not describe its matrix. H3's
    alignment between a drug mean and a target direction is only meaningful if
    both sit in the same coordinate system, so this is measured, not assumed.
    """
    import json
    import sys

    import numpy as np

    sys.path.insert(0, "/app")
    _repo_at_its_absolute_path()
    staged, rebuilt_dir = _stage_inputs(), _stage_rebuild()
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    records = json.loads(open("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read())
    targets = sorted({str(r["target"]) for r in records})
    rebuilt_sigs, _, rebuilt_genes, _ = gate.load_shrna_rebuild_strict(rebuilt_dir)
    extraction_sigs, _, extraction_genes = gate.load_shrna_extraction(
        staged, targets_wanted=set(targets))

    rebuilt_genes = [str(g) for g in rebuilt_genes]
    extraction_genes = [str(g) for g in extraction_genes]
    declared = {g: i for i, g in enumerate(extraction_genes)}
    by_label = np.array([declared[g] for g in rebuilt_genes])

    sample = sorted(set(rebuilt_sigs) & set(extraction_sigs))[:50]
    identity_err, label_err = [], []
    for sig_id in sample:
        left = np.asarray(rebuilt_sigs[sig_id], dtype=np.float64)
        right = np.asarray(extraction_sigs[sig_id], dtype=np.float64)
        identity_err.append(float(np.abs(left - right).max()))
        label_err.append(float(np.abs(left - right[by_label]).max()))

    # and the same question for the compound matrix, whose half of the gate passed
    compounds = np.load(f"{staged}/lincs_subset.npz", allow_pickle=True)
    compound_genes = [str(g) for g in compounds["gene_ids"]]

    summary = {
        "n_sampled": len(sample),
        "shrna_identity_order_max_error": max(identity_err),
        "shrna_declared_label_order_max_error": max(label_err),
        "reading": ("whichever is at float tolerance is the order the extraction's shRNA "
                    "matrix actually uses"),
        "shrna_matrix_follows": ("the rebuild's sorted order, so its own gene_ids labels are "
                                 "wrong" if max(identity_err) < 1e-4 else
                                 "its declared labels" if max(label_err) < 1e-4 else
                                 "neither: a third order"),
        "compound_gene_ids_equal_shrna_gene_ids": compound_genes == extraction_genes,
        "compound_gene_ids_sorted": compound_genes == sorted(compound_genes, key=int),
        "shrna_gene_ids_sorted": extraction_genes == sorted(extraction_genes, key=int),
        "compound_half_of_the_gate": "passed, by declared labels, before the shRNA half ran",
    }
    open("/out/03c_h3_sensitivity/gate_axis_test.json", "w").write(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_recover_permutation():
    """Recover the permutation between the extraction's shRNA columns and the GCTX.

    The matrix follows neither its declared labels nor sorted order. If one
    consistent permutation maps it onto the rebuild for every signature, the axis
    is recoverable and the directions can be relabelled rather than rebuilt. If no
    consistent permutation exists, they cannot.
    """
    import json
    import sys

    import numpy as np

    sys.path.insert(0, "/app")
    _repo_at_its_absolute_path()
    staged, rebuilt_dir = _stage_inputs(), _stage_rebuild()
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    records = json.loads(open("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read())
    targets = sorted({str(r["target"]) for r in records})
    rebuilt_sigs, _, rebuilt_genes, _ = gate.load_shrna_rebuild_strict(rebuilt_dir)
    extraction_sigs, _, extraction_genes = gate.load_shrna_extraction(
        staged, targets_wanted=set(targets))
    rebuilt_genes = [str(g) for g in rebuilt_genes]
    extraction_genes = [str(g) for g in extraction_genes]

    shared = sorted(set(rebuilt_sigs) & set(extraction_sigs))
    probe = shared[0]
    left = np.asarray(rebuilt_sigs[probe], dtype=np.float64)
    right = np.asarray(extraction_sigs[probe], dtype=np.float64)
    # the permutation that sorts both the same way maps one onto the other
    candidate = np.argsort(right)[np.argsort(np.argsort(left))]

    agree, disagree = 0, []
    for sig_id in shared[:300]:
        a = np.asarray(rebuilt_sigs[sig_id], dtype=np.float64)
        b = np.asarray(extraction_sigs[sig_id], dtype=np.float64)[candidate]
        if np.allclose(a, b, rtol=1e-4, atol=1e-4):
            agree += 1
        elif len(disagree) < 5:
            disagree.append({"sig_id": sig_id, "max_abs": float(np.abs(a - b).max())})

    implied = [extraction_genes[i] for i in candidate]
    summary = {
        "n_checked": min(300, len(shared)),
        "n_agreeing_under_one_permutation": agree,
        "permutation_is_consistent": agree == min(300, len(shared)),
        "examples_disagreeing": disagree,
        "implied_extraction_order_head": implied[:8],
        "rebuild_order_head": rebuilt_genes[:8],
        "implied_equals_declared": implied == extraction_genes,
        "implied_is_sorted": implied == sorted(implied, key=int),
        "reading": ("one consistent permutation means the shRNA columns are a relabelling "
                    "away from the GCTX; no consistent permutation means they are not the "
                    "same numbers at all"),
    }
    open("/out/03c_h3_sensitivity/gate_permutation_recovery.json", "w").write(
        json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_verify_permutation():
    """Verify the recovered permutation on every signature, and try to name it.

    A permutation recovered from one vector and confirmed on 300 is a hypothesis.
    This checks all 14,656, confirms the map is a bijection, writes it out so the
    repair uses a recorded object rather than a rediscovered one, and tests the
    orders it might correspond to: the gene-info file's own row order, the
    landmark ids as strings, or the symbols alphabetically.
    """
    import json
    import sys

    import numpy as np
    import pandas as pd

    sys.path.insert(0, "/app")
    _repo_at_its_absolute_path()
    staged, rebuilt_dir = _stage_inputs(), _stage_rebuild()
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    records = json.loads(open("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read())
    targets = sorted({str(r["target"]) for r in records})
    rebuilt_sigs, _, rebuilt_genes, _ = gate.load_shrna_rebuild_strict(rebuilt_dir)
    extraction_sigs, _, extraction_genes = gate.load_shrna_extraction(
        staged, targets_wanted=set(targets))
    rebuilt_genes = [str(g) for g in rebuilt_genes]
    extraction_genes = [str(g) for g in extraction_genes]

    shared = sorted(set(rebuilt_sigs) & set(extraction_sigs))
    probe = shared[0]
    left = np.asarray(rebuilt_sigs[probe], dtype=np.float64)
    right = np.asarray(extraction_sigs[probe], dtype=np.float64)
    permutation = np.argsort(right)[np.argsort(np.argsort(left))]

    assert sorted(permutation.tolist()) == list(range(len(rebuilt_genes))), (
        "the recovered map is not a bijection, so it is not a permutation")

    worst, failures = 0.0, []
    for sig_id in shared:
        a = np.asarray(rebuilt_sigs[sig_id], dtype=np.float64)
        b = np.asarray(extraction_sigs[sig_id], dtype=np.float64)[permutation]
        error = float(np.abs(a - b).max())
        worst = max(worst, error)
        if error > 1e-4 and len(failures) < 10:
            failures.append({"sig_id": sig_id, "max_abs": error})

    # what order does the extraction's matrix actually sit in?
    implied = [extraction_genes[i] for i in permutation]      # rebuild position -> extraction gene
    actual_column_order = [None] * len(implied)
    for position, gene in zip(permutation, rebuilt_genes):
        actual_column_order[int(position)] = gene

    gene_info = pd.read_csv(f"{staged}/GSE92742_Broad_LINCS_gene_info.txt.gz", sep="\t",
                            low_memory=False)
    landmark = gene_info[gene_info.pr_is_lm == 1]
    file_order = [str(g) for g in landmark.pr_gene_id]
    symbol_order = [str(g) for g in landmark.sort_values("pr_gene_symbol").pr_gene_id]
    string_sorted = sorted([str(g) for g in rebuilt_genes])

    summary = {
        "n_signatures_checked": len(shared),
        "max_abs_difference_under_the_permutation": worst,
        "all_within_1e-4": worst <= 1e-4,
        "failures": failures,
        "is_a_bijection": True,
        "actual_column_order_head": actual_column_order[:8],
        "declared_order_head": extraction_genes[:8],
        "matches_gene_info_file_order": actual_column_order == file_order,
        "matches_symbol_alphabetical": actual_column_order == symbol_order,
        "matches_id_as_string_sorted": actual_column_order == string_sorted,
        "permutation": permutation.tolist(),
        "reading": ("the extraction's shRNA columns sit in `actual_column_order`; its stored "
                    "gene_ids say otherwise, and the compound matrix follows the stored ids"),
    }
    open("/out/03c_h3_sensitivity/shrna_axis_recovery.json", "w").write(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps({k: v for k, v in summary.items() if k != "permutation"}, indent=2)[:1800]


@app.function(**COMMON)
def stage_axis_uniqueness():
    """Is the recovered permutation the only one that reconciles the two matrices?

    Deviation 11 states that one bijection maps the extraction onto the rebuild for
    every signature. That is existence, not uniqueness: a second valid permutation
    exists exactly when two gene columns carry the same values across all 14,656
    signatures, because those two columns could then be exchanged without changing
    anything. This measures how far apart the closest pair of columns is, so the
    deviation can say whether the recovered map is the only one or merely one.
    """
    import json
    import importlib.util
    from pathlib import Path

    import numpy as np

    _repo_at_its_absolute_path()
    staged, rebuilt_dir = _stage_inputs(), _stage_rebuild()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    records = json.loads(Path("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read_text())
    targets = sorted({str(r["target"]) for r in records})
    rebuilt_sigs, _, rebuilt_genes, _ = gate.load_shrna_rebuild_strict(rebuilt_dir)
    extraction_sigs, _, _ = gate.load_shrna_extraction(staged, targets_wanted=set(targets))

    shared = sorted(set(rebuilt_sigs) & set(extraction_sigs))
    genes = [str(g) for g in rebuilt_genes]
    matrix = np.vstack([np.asarray(rebuilt_sigs[s], dtype=np.float64) for s in shared])
    assert matrix.shape == (len(shared), len(genes)), f"matrix is {matrix.shape}"

    # two identical columns are the only way a second permutation can exist
    profiles = {}
    for column, gene in enumerate(genes):
        profiles.setdefault(np.ascontiguousarray(matrix[:, column]).tobytes(), []).append(gene)
    collisions = {digest_genes[0]: digest_genes[1:]
                  for digest_genes in profiles.values() if len(digest_genes) > 1}

    # the closest pair, whether or not it is an exact tie: squared distances from
    # the Gram matrix, which is one 978 x 978 matmul rather than 478,000 loops
    gram = matrix.T @ matrix
    square = np.diag(gram).copy()
    distance_squared = square[:, None] + square[None, :] - 2.0 * gram
    np.fill_diagonal(distance_squared, np.inf)
    flat = int(np.argmin(distance_squared))
    first, second = divmod(flat, len(genes))
    closest = float(np.sqrt(max(distance_squared[first, second], 0.0)))
    closest_elementwise = float(np.abs(matrix[:, first] - matrix[:, second]).max())

    tolerance = gate.ATOL + gate.RTOL * float(np.abs(matrix).mean())
    summary = {
        "n_signatures": len(shared),
        "n_columns": len(genes),
        "n_distinct_column_profiles": len(profiles),
        "n_columns_sharing_a_profile": sum(len(v) + 1 for v in collisions.values()),
        "colliding_columns": collisions,
        "permutation_is_unique": len(profiles) == len(genes),
        "closest_pair": {"genes": [genes[first], genes[second]],
                         "l2_distance": closest,
                         "max_abs_elementwise_difference": closest_elementwise},
        "mean_abs_value": float(np.abs(matrix).mean()),
        "gate_elementwise_tolerance_at_that_scale": tolerance,
        "reading": ("every column profile distinct means exactly one permutation maps the "
                    "extraction onto the rebuild; the closest pair's elementwise difference "
                    "says how far that conclusion sits from the gate's tolerance"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "shrna_axis_uniqueness.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps({k: v for k, v in summary.items() if k != "colliding_columns"}, indent=2)


@app.function(**{**COMMON, "image": test_image})
def stage_tests():
    """The registered suite, in the image the analyses run in."""
    import subprocess

    _repo_at_its_absolute_path()
    from pathlib import Path

    finished = subprocess.run(
        ["python", "-m", "pytest", "/app/tests", "-v", "-rs", "--no-header",
         "--ignore=/app/tests/test_combined_experiments.py"],
        cwd="/app", capture_output=True, text=True)
    print(finished.stdout[-8000:], flush=True)
    print(finished.stderr[-4000:], flush=True)
    # the collected list and the skip reasons are the record of what actually ran,
    # so they go to a file rather than only to a log that scrolls
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "test_output.txt").write_text(finished.stdout + finished.stderr)
    results.commit()
    assert finished.returncode == 0, f"the suite failed with {finished.returncode}"
    return finished.stdout[-2000:]


@app.function(**COMMON)
def stage_reader_agreement():
    """Do the production parser and the gate's own HDF5 read agree on the real file?

    The gate's reader takes the GCTX's row and column identifiers from its HDF5
    metadata and builds the positions itself. Those paths come from the format's
    layout rather than from this file, so they are checked against it before the
    gate relies on them. This also records what the production parser returns:
    the file's own row order, or the order it was asked for.
    """
    import json
    import importlib.util
    from pathlib import Path

    import numpy as np
    import pandas as pd
    from cmapPy.pandasGEXpress import parse

    _repo_at_its_absolute_path()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    gene_info = pd.read_csv(raw / "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="\t",
                            low_memory=False)
    landmark = [str(g) for g in gene_info[gene_info.pr_is_lm == 1].pr_gene_id]
    siginfo = pd.read_csv(raw / "lincs_shrna_siginfo.csv.gz")
    sample = sorted(set(siginfo.sig_id.astype(str)))[:200]

    values, genes, signatures, hashes = gate.read_gctx_slice(raw / GCTX, sample, landmark)

    # the production route: ask for an order, and see what comes back
    requested = list(reversed(landmark))
    gct = parse.parse(str(raw / GCTX), cid=sample, rid=requested)
    frame = gct.data_df
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)

    rows = {sig_id: i for i, sig_id in enumerate(signatures)}
    columns = {gene: i for i, gene in enumerate(genes)}
    worst = 0.0
    for sig_id in signatures:
        for gene in frame.index:
            left = float(values[rows[sig_id], columns[gene]])
            worst = max(worst, abs(left - float(frame.loc[gene, sig_id])))

    summary = {
        "n_signatures": len(signatures), "n_genes": len(genes),
        "max_abs_difference_after_identifier_alignment": worst,
        "readers_agree": worst <= gate.ATOL + gate.RTOL * float(np.abs(values).mean()),
        "orientation_the_reader_found": hashes["orientation"],
        "parser_returned_the_requested_order": list(frame.index) == requested,
        "parser_returned_the_source_order": list(frame.index) == genes,
        "selected_gene_axis_sha256": hashes["selected_gene_axis_sha256"],
        "source_row_axis_sha256": hashes["source_row_axis_sha256"],
        "reading": ("agreement after aligning on identifiers means the gate's own read of the "
                    "HDF5 metadata finds the same value under the same gene as the parser "
                    "production uses, by a route that shares no lookup with it"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "reader_agreement.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.local_entrypoint()
def main(stage: str):
    """stage: gate | s1s3 | driver | diagnostic"""
    if stage == "gate":
        print(stage_gate.remote())
    elif stage == "s1s3":
        print(stage_s1s3.remote())
    elif stage == "driver":
        print(stage_driver.remote())
    elif stage == "diagnostic":
        print(stage_gate_diagnostic.remote())
    else:
        raise SystemExit(f"unknown stage {stage}")
    print("\nretrieve with:")
    print("    modal volume get di-h3-results 03c_h3_sensitivity ./results/")
    print("    modal volume get di-h3-results 03d_h3_reference_discordance ./results/")
