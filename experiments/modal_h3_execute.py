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
                           remote_path="/app/results/03d_h3_reference_discordance")
            # the pinned cohort manifests and the frozen identities gate v2 reads
            .add_local_dir(f"{REPO}/registry", remote_path="/app/registry")
            .add_local_dir(f"{REPO}/results/03c_h3_sensitivity",
                           remote_path="/app/results/03c_h3_sensitivity",
                           ignore=["gctx_rebuild/**", "gate_*/**", "*.npz"]))


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
# a report never overwrites an earlier one, so the directory carries the time the
# run started rather than a date typed by hand
GATE_RUN = None          # set per invocation in stage_gate

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


def _sha256_of(path):
    """The digest of one file, for the stages that check a frozen artifact."""
    import hashlib
    from pathlib import Path

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _run(argv):
    """Run one of the analysis scripts, with the repo importable.

    A separate process rather than `runpy` in this one, so the script receives a
    real argv instead of one assigned onto the interpreter, and so a refusal
    arrives as an exit code this function raises on rather than as an exception
    crossing a boundary it was not raised across.
    """
    import os
    import subprocess

    print(f"[{_ts()}] {' '.join(argv)}", flush=True)
    environment = {**os.environ, "PYTHONPATH": "/app"}
    finished = subprocess.run(["python", *argv], cwd="/app", env=environment)
    if finished.returncode != 0:
        raise RuntimeError(
            f"{argv[0]} exited {finished.returncode}; its output is above and the stage "
            "does not continue past a script that refused")


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
def stage_gate(cohort: str = "/app/results/03_phenotype_projection/phenotype_projection_results.json"):
    """The reconstruction gate. A failure here voids the analyses registered against it."""
    from pathlib import Path

    from datetime import datetime

    run = f"gate_run_{datetime.now():%Y-%m-%dT%H%M%S}"
    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    _run(["/app/experiments/03e_reconstruction_gate.py",
          "--rebuilt", str(rebuilt),
          "--extraction", str(staged),
          "--gctx", f"/rebuild/raw/{GCTX}",
          "--shrna-siginfo", "/rebuild/raw/lincs_shrna_siginfo.csv.gz",
          "--compound-siginfo", "/rebuild/raw/GSE92742_Broad_LINCS_sig_info.txt.gz",
          "--gene-info", "/rebuild/raw/GSE92742_Broad_LINCS_gene_info.txt.gz",
          "--cohort", cohort,
          "--output", f"/out/03c_h3_sensitivity/{run}"])
    results.commit()
    return (Path(f"/out/03c_h3_sensitivity/{run}/reconstruction_gate.json")
            ).read_text()[:2000]


@app.function(**COMMON)
def stage_s1s3():
    """S1-S3 and their target-level versions, on the bundle the rebuild produced."""
    import shutil
    from pathlib import Path

    _repo_at_its_absolute_path()
    _run(["/app/experiments/03c_h3_sensitivity.py",
          "--bundle", "/rebuild/results/cohort_bundle.npz",
          "--manifest", "/rebuild/results/rebuild_manifest.json"])
    # The script writes beside its own source. The symlink puts that in the image
    # layer, not on the results volume, so the outputs are copied across before the
    # commit; without this the run completes and is lost when the container exits.
    written = Path("/app/results/03c_h3_sensitivity")
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    carried = []
    for name in ("h3_sensitivity_results.json", "h3_sensitivity_results.json.sha256",
                 "h3_sensitivity_draws.npz"):
        source = written / name
        if not source.exists():
            raise RuntimeError(
                f"{name} is absent after the script returned; the stage does not "
                "commit a partial run")
        shutil.copy2(source, out / name)
        carried.append(name)
    results.commit()
    return f"s1s3 done; carried {', '.join(carried)}"


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
    import site

    import numpy as np

    site.addsitedir("/app")
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
    import site

    import numpy as np

    site.addsitedir("/app")
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
    import site

    import numpy as np

    site.addsitedir("/app")
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
    import site

    import numpy as np

    site.addsitedir("/app")
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
    import site

    import numpy as np
    import pandas as pd

    site.addsitedir("/app")
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

    if not (sorted(permutation.tolist()) == list(range(len(rebuilt_genes)))):
        raise AssertionError("the recovered map is not a bijection, so it is not a permutation")

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
    if not (matrix.shape == (len(shared), len(genes))):
        raise AssertionError(f"matrix is {matrix.shape}")

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
def stage_tests(selector: str = "/app/tests", optimized: bool = False):
    """The registered suite, in the image the analyses run in.

    `selector` narrows it to one or more files while a change is being worked on;
    the release check runs it with the default, which is the whole suite.
    `optimized` runs it under `python -O`, which strips `assert` statements: a
    contract that only holds without that flag is not a contract.
    """
    import subprocess

    _repo_at_its_absolute_path()
    from pathlib import Path

    interpreter = ["python", "-O"] if optimized else ["python"]
    finished = subprocess.run(
        interpreter + ["-m", "pytest", *selector.split(), "-v", "-rs", "--no-header",
                       "--ignore=/app/tests/test_combined_experiments.py"],
        cwd="/app", capture_output=True, text=True)
    print(finished.stdout[-8000:], flush=True)
    print(finished.stderr[-4000:], flush=True)
    # the collected list and the skip reasons are the record of what actually ran,
    # so they go to a file rather than only to a log that scrolls
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    suffix = "_optimized" if optimized else ""
    name = f"test_output{suffix}.txt" if selector == "/app/tests" else \
        f"test_output_{Path(selector.split()[0]).stem}{suffix}.txt"
    (out / name).write_text(finished.stdout + finished.stderr)
    results.commit()
    if finished.returncode != 0:
        raise AssertionError(
            f"the suite failed with {finished.returncode} under "
            f"{' '.join(interpreter)}; see {name}")
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


@app.function(**COMMON)
def stage_diagnose_compound():
    """Why does the retained extraction disagree with the source by exactly 20.0?

    The gate compares drug-by-cell aggregates. A constant maximum of 20.0 across
    every drug is not what a grouping difference looks like, and the rebuild agrees
    with the source to 1e-7, so this compares the retained extraction against the
    source at the level of single signatures, joined on signature id and gene id,
    where no aggregation can be responsible. It also hashes both copies of the
    compound metadata, because two exist on two volumes.
    """
    import hashlib
    import importlib.util
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")

    def digest(path):
        sha = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 22), b""):
                sha.update(block)
        return sha.hexdigest()

    metadata = {
        "beside the extraction": digest(Path("/extraction")
                                        / "GSE92742_Broad_LINCS_sig_info.txt.gz"),
        "pinned on the rebuild volume": digest(raw / "GSE92742_Broad_LINCS_sig_info.txt.gz"),
    }

    records = json.loads(Path("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read_text())
    drugs = sorted({r["drug"] for r in records})
    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    siginfo = siginfo[siginfo.pert_iname.isin(set(drugs))]

    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    extraction_genes = [str(g) for g in extraction["gene_ids"]]
    matrix = extraction["signatures"]

    shared = sorted(set(siginfo.sig_id.astype(str)) & set(held))
    sample = shared[:50]
    values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, sample, extraction_genes)
    rows = {s: i for i, s in enumerate(signatures)}
    column = {g: i for i, g in enumerate(genes)}
    order = np.array([column[g] for g in extraction_genes])

    rows_out = []
    for sig_id in sample:
        left = np.asarray(matrix[held[sig_id]], dtype=np.float64)
        right = values[rows[sig_id]][order]           # source, on the extraction's axis
        denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
        rows_out.append({
            "sig_id": sig_id,
            "max_abs_difference": float(np.abs(left - right).max()),
            "cosine": float(left @ right / denominator) if denominator else None,
            "max_abs_difference_if_negated": float(np.abs(left + right).max()),
            "extraction_range": [float(left.min()), float(left.max())],
            "source_range": [float(right.min()), float(right.max())],
            "n_at_plus_ten_extraction": int((left >= 9.999).sum()),
            "n_at_minus_ten_extraction": int((left <= -9.999).sum()),
            "n_at_plus_ten_source": int((right >= 9.999).sum()),
            "n_at_minus_ten_source": int((right <= -9.999).sum()),
        })

    worst = max(r["max_abs_difference"] for r in rows_out)
    cosines = [r["cosine"] for r in rows_out if r["cosine"] is not None]
    summary = {
        "metadata_copies": metadata,
        "metadata_copies_identical": len(set(metadata.values())) == 1,
        "n_cohort_signatures_in_metadata": len(set(siginfo.sig_id.astype(str))),
        "n_of_those_held_by_the_extraction": len(shared),
        "n_sampled": len(sample),
        "max_abs_difference_at_signature_level": worst,
        "median_cosine": float(np.median(cosines)) if cosines else None,
        "n_explained_by_negation": sum(1 for r in rows_out
                                       if r["max_abs_difference_if_negated"] < 1e-3),
        "worst_five": sorted(rows_out, key=lambda r: -r["max_abs_difference"])[:5],
        "reading": ("a disagreement at the level of single signatures rules out grouping and "
                    "aggregation; a cosine near -1 or agreement under negation would mean a "
                    "sign convention, and saturation counts say whether clipping is involved"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "compound_discrepancy_diagnostic.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps({k: v for k, v in summary.items() if k != "worst_five"}, indent=2)


@app.function(**COMMON)
def stage_compound_axis_recovery():
    """Is the compound extraction's axis permuted too, and is it the same permutation?

    A median cosine of 0.003 at single-signature level, with matching value ranges,
    is what Deviation 11 looked like. This asks the same questions of the compound
    matrix: does one bijection reconcile it with the source for every sampled
    signature, does that order match anything the pipeline uses, and is it the same
    permutation the shRNA matrix carries.
    """
    import importlib.util
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    records = json.loads(Path("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read_text())
    drugs = sorted({r["drug"] for r in records})
    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    siginfo = siginfo[siginfo.pert_iname.isin(set(drugs))]

    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    declared = [str(g) for g in extraction["gene_ids"]]
    matrix = extraction["signatures"]

    sample = sorted(set(siginfo.sig_id.astype(str)) & set(held))[:300]
    values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, sample, declared)
    rows = {s: i for i, s in enumerate(signatures)}
    column = {g: i for i, g in enumerate(genes)}
    on_declared = np.array([column[g] for g in declared])

    probe = sample[0]
    left = np.asarray(matrix[held[probe]], dtype=np.float64)       # extraction
    right = values[rows[probe]][on_declared]                       # source, declared order
    permutation = np.argsort(left)[np.argsort(np.argsort(right))]
    is_bijection = sorted(permutation.tolist()) == list(range(len(declared)))

    worst, agree = 0.0, 0
    for sig_id in sample:
        a = values[rows[sig_id]][on_declared]
        b = np.asarray(matrix[held[sig_id]], dtype=np.float64)[permutation]
        error = float(np.abs(a - b).max())
        worst = max(worst, error)
        agree += error <= 1e-4

    # what order does the extraction's matrix actually sit in?
    actual = [None] * len(declared)
    for position, gene in zip(permutation, declared):
        actual[int(position)] = gene
    gene_info = pd.read_csv(raw / "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="\t",
                            low_memory=False)
    landmark = gene_info[gene_info.pr_is_lm == 1]
    frozen = [str(g) for g in landmark.sort_values("pr_gene_id").pr_gene_id]

    shrna_recovery = Path("/app/results/03c_h3_sensitivity/shrna_axis_recovery.json")
    shrna_order = (json.loads(shrna_recovery.read_text()).get("actual_column_order_head")
                   if shrna_recovery.exists() else None)

    summary = {
        "n_sampled": len(sample),
        "recovered_map_is_a_bijection": is_bijection,
        "n_agreeing_under_one_permutation": agree,
        "max_abs_difference_under_the_permutation": worst,
        "one_permutation_reconciles_every_sampled_signature": agree == len(sample),
        "declared_order_head": declared[:8],
        "actual_column_order_head": actual[:8],
        "declared_equals_frozen_landmark_order": declared == frozen,
        "actual_equals_frozen_landmark_order": actual == frozen,
        "actual_equals_gene_info_file_order": actual == [str(g) for g in landmark.pr_gene_id],
        "shrna_actual_order_head_for_comparison": shrna_order,
        "reading": ("one bijection reconciling every signature means the compound matrix "
                    "carries correct values under a wrong axis, as the shRNA matrix does; "
                    "comparing the two recovered orders says whether it is the same wrong axis"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "compound_axis_recovery.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_compound_row_test():
    """Does each extraction row hold the signature its label claims?

    No single gene permutation reconciles the compound extraction with the source,
    which rules out Deviation 11's shape. Matching value ranges with a near-zero
    cosine and no consistent column map is what row mislabeling looks like instead:
    the row is a real signature, just not the one the label names. This searches the
    source for the signature each row actually holds.
    """
    import importlib.util
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    records = json.loads(Path("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read_text())
    drugs = sorted({r["drug"] for r in records})
    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    cohort_ids = sorted(set(siginfo[siginfo.pert_iname.isin(set(drugs))].sig_id.astype(str)))

    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    stored_ids = [str(s) for s in extraction["sig_ids"]]
    held = {sig_id: i for i, sig_id in enumerate(stored_ids)}
    declared = [str(g) for g in extraction["gene_ids"]]
    matrix = extraction["signatures"]

    # read a block of the source large enough that a mislabeled row's true owner is
    # likely inside it, on the extraction's declared gene order
    block = cohort_ids[:3000]
    values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, block, declared)
    column = {g: i for i, g in enumerate(genes)}
    on_declared = np.array([column[g] for g in declared])
    source = values[:, on_declared]
    source_unit = source / np.linalg.norm(source, axis=1, keepdims=True)
    position_in_source = {sig_id: i for i, sig_id in enumerate(signatures)}

    probes = [sig_id for sig_id in block[:120] if sig_id in held]
    found_elsewhere, exact_self, offsets, examples = 0, 0, [], []
    for sig_id in probes:
        row = np.asarray(matrix[held[sig_id]], dtype=np.float64)
        norm = np.linalg.norm(row)
        if norm == 0:
            continue
        similarity = source_unit @ (row / norm)
        best = int(np.argmax(similarity))
        if similarity[best] > 0.999:
            owner = signatures[best]
            if owner == sig_id:
                exact_self += 1
            else:
                found_elsewhere += 1
                offsets.append(position_in_source[owner] - position_in_source[sig_id])
                if len(examples) < 8:
                    examples.append({"label": sig_id, "actually_holds": owner,
                                     "cosine": float(similarity[best])})

    summary = {
        "n_probed": len(probes),
        "rows_holding_the_signature_their_label_names": exact_self,
        "rows_holding_a_different_signature": found_elsewhere,
        "rows_matching_nothing_in_the_block": len(probes) - exact_self - found_elsewhere,
        "examples": examples,
        "offset_is_constant": len(set(offsets)) == 1 if offsets else None,
        "distinct_offsets": sorted(set(offsets))[:10],
        "stored_ids_are_sorted": stored_ids == sorted(stored_ids),
        "n_signatures_in_the_extraction": len(stored_ids),
        "reading": ("a row matching a different source signature at cosine ~1 means the "
                    "values are real and the labels are wrong, which is Deviation 9's shape "
                    "rather than Deviation 11's; a constant offset would name the mechanism"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "compound_row_test.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_compound_column_match():
    """Recover the compound axis by column profile rather than by sorting values.

    The earlier recovery derived a permutation from one signature's `argsort`, which
    cannot work here: tens of values per signature sit pinned at the clipping bounds,
    so the order among ties is arbitrary and a map derived from one signature does not
    transfer. Identical minima, maxima and saturation counts on both sides say the
    values are the same multiset, so this matches each declared column against the
    source columns by its profile across many signatures, where ties do not arise.
    """
    import importlib.util
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    records = json.loads(Path("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read_text())
    drugs = sorted({r["drug"] for r in records})
    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    cohort = sorted(set(siginfo[siginfo.pert_iname.isin(set(drugs))].sig_id.astype(str)))

    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    declared = [str(g) for g in extraction["gene_ids"]]
    sample = [s for s in cohort if s in held][:400]

    # the array is decompressed once. Indexing an NpzFile inside a comprehension
    # re-inflates the whole 167,266 by 978 matrix on every access, which is the
    # mistake this file already carries a comment about in load_shrna_extraction.
    stored = extraction["signatures"]
    left = np.vstack([np.asarray(stored[held[s]], dtype=np.float64) for s in sample])
    del stored
    values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, sample, declared)
    row_of = {s: i for i, s in enumerate(signatures)}
    right = np.vstack([values[row_of[s]] for s in sample])        # source, its own order

    profiles = {}
    for column in range(right.shape[1]):
        profiles.setdefault(np.ascontiguousarray(right[:, column]).tobytes(), []).append(column)
    collisions = sum(1 for cols in profiles.values() if len(cols) > 1)

    matched, unmatched = {}, []
    for column, gene in enumerate(declared):
        key = np.ascontiguousarray(left[:, column]).tobytes()
        where = profiles.get(key)
        if where and len(where) == 1:
            matched[gene] = genes[where[0]]
        else:
            unmatched.append(gene)

    images = list(matched.values())
    summary = {
        "n_signatures_used": len(sample),
        "n_declared_columns": len(declared),
        "n_source_columns_with_a_unique_profile": len(profiles),
        "n_source_profile_collisions": collisions,
        "n_declared_columns_matched_exactly": len(matched),
        "n_unmatched": len(unmatched),
        "unmatched_examples": unmatched[:5],
        "match_is_a_bijection": len(set(images)) == len(images) == len(declared),
        "n_columns_already_in_the_right_place": sum(1 for gene, image in matched.items()
                                                    if gene == image),
        "example_relabelings": [{"declared": gene, "actually": image}
                                for gene, image in list(matched.items())[:8]],
        "reading": ("every declared column matching exactly one source column by profile, as a "
                    "bijection, means the compound extraction carries correct values under a "
                    "wrong gene axis: the same defect class as Deviation 11, on the drug data"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "compound_column_match.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_scope_of_the_compound_defect():
    """Does raw direction instability itself differ between the extraction and the source?

    Deviation 11 scoped the axis defect to quantities that pair the two matrices, on
    the ground that a permutation applied to every signature alike leaves the cosines
    between signatures unchanged. The compound disagreement is not a shared
    permutation, so that reasoning does not carry, and D is computed from this file.
    This recomputes D for a sample of cohort drugs from both sources and compares them
    against the deposited values. It also tests the multiset directly, by sorting each
    signature's values, which the earlier summary statistics only suggested.
    """
    import importlib.util
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    records = json.loads(Path("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read_text())
    deposited = {r["drug"]: r for r in records}
    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    siginfo = siginfo[siginfo.pert_iname.isin(set(deposited))]

    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    declared = [str(g) for g in extraction["gene_ids"]]
    stored = extraction["signatures"]            # decompressed once

    drugs = sorted(deposited)[:20]
    wanted = sorted(set(siginfo[siginfo.pert_iname.isin(drugs)].sig_id.astype(str))
                    & set(held))
    values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, wanted, declared)
    column = {g: i for i, g in enumerate(genes)}
    on_declared = np.array([column[g] for g in declared])
    source_row = {s: i for i, s in enumerate(signatures)}

    def instability(matrix):
        unit = matrix / np.linalg.norm(matrix, axis=1, keepdims=True)
        upper = np.triu_indices(matrix.shape[0], 1)
        return float(1.0 - (unit @ unit.T)[upper].mean())

    # the multiset test the summary statistics only hinted at
    multiset_matches, probes = 0, wanted[:40]
    for sig_id in probes:
        left = np.sort(np.asarray(stored[held[sig_id]], dtype=np.float64))
        right = np.sort(values[source_row[sig_id]][on_declared])
        multiset_matches += bool(np.allclose(left, right, rtol=1e-5, atol=1e-5))
    monotone = sum(1 for sig_id in probes
                   if np.all(np.diff(np.asarray(stored[held[sig_id]])) >= 0)
                   or np.all(np.diff(np.asarray(stored[held[sig_id]])) <= 0))

    rows = []
    for drug, group in siginfo[siginfo.pert_iname.isin(drugs)].groupby("pert_iname"):
        per_cell_left, per_cell_right = [], []
        for cell, members in group.groupby("cell_id"):
            ids = sorted({s for s in members.sig_id.astype(str) if s in held})
            if not ids:
                continue
            per_cell_left.append(np.vstack([np.asarray(stored[held[s]], dtype=np.float64)
                                            for s in ids]).mean(axis=0))
            per_cell_right.append(np.vstack([values[source_row[s]][on_declared]
                                             for s in ids]).mean(axis=0))
        if len(per_cell_left) < 2:
            continue
        left, right = np.vstack(per_cell_left), np.vstack(per_cell_right)
        rows.append({"drug": drug, "n_cell_lines": left.shape[0],
                     "D_extraction": instability(left), "D_source": instability(right),
                     "D_deposited": deposited[drug]["raw_instability"]})

    for entry in rows:
        entry["extraction_minus_deposited"] = entry["D_extraction"] - entry["D_deposited"]
        entry["source_minus_deposited"] = entry["D_source"] - entry["D_deposited"]

    summary = {
        "n_drugs": len(rows),
        "multiset_test": {"n_probed": len(probes), "n_matching_sorted_values": multiset_matches,
                          "n_rows_monotonic": monotone},
        "D_extraction_reproduces_deposited": bool(rows) and all(
            abs(e["extraction_minus_deposited"]) < 1e-6 for e in rows),
        "D_source_reproduces_deposited": bool(rows) and all(
            abs(e["source_minus_deposited"]) < 1e-6 for e in rows),
        "max_abs_extraction_minus_deposited": max((abs(e["extraction_minus_deposited"])
                                                   for e in rows), default=None),
        "max_abs_source_minus_deposited": max((abs(e["source_minus_deposited"])
                                               for e in rows), default=None),
        "per_drug": rows,
        "reading": ("if D from the extraction reproduces the deposited values and D from the "
                    "source does not, every quantity the paper reports was computed from this "
                    "file and the disagreement reaches the primary construct, not only H3"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "compound_defect_scope.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps({k: v for k, v in summary.items() if k != "per_drug"}, indent=2)


@app.function(**COMMON)
def stage_compare_the_two_permutations():
    """Is the compound axis permuted, and is it the same permutation the shRNA axis carries?

    The earlier column match compared a float64 array from the npz against a float32
    array from HDF5 by their bytes, which cannot match whatever the values are, so its
    zero of 978 said nothing. Both sides are cast here before comparison.

    The question this settles matters more than the mechanism. H3 pairs a drug
    signature with a target direction. If both files carry the same wrong axis, the
    pairing was consistent and the cosine between them is unaffected; if they carry
    different wrong axes, it is not.
    """
    import importlib.util
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    records = json.loads(Path("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read_text())
    drugs = sorted({r["drug"] for r in records})
    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    cohort = sorted(set(siginfo[siginfo.pert_iname.isin(set(drugs))].sig_id.astype(str)))

    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    declared = [str(g) for g in extraction["gene_ids"]]
    stored = extraction["signatures"]
    sample = [s for s in cohort if s in held][:400]

    left = np.vstack([np.asarray(stored[held[s]]) for s in sample]).astype(np.float32)
    values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, sample, declared)
    row_of = {s: i for i, s in enumerate(signatures)}
    right = np.vstack([values[row_of[s]] for s in sample]).astype(np.float32)

    profiles = {}
    for index in range(right.shape[1]):
        profiles.setdefault(np.ascontiguousarray(right[:, index]).tobytes(), []).append(index)

    matched, unmatched = {}, []
    for index, gene in enumerate(declared):
        where = profiles.get(np.ascontiguousarray(left[:, index]).tobytes())
        if where and len(where) == 1:
            matched[gene] = genes[where[0]]
        else:
            unmatched.append(gene)

    # the shRNA permutation, as Deviation 11 recovered it
    shrna = json.loads(Path("/app/results/03c_h3_sensitivity/shrna_axis_recovery.json")
                       .read_text()) if Path(
        "/app/results/03c_h3_sensitivity/shrna_axis_recovery.json").exists() else {}
    shrna_actual = shrna.get("actual_column_order_head")

    # the compound matrix's actual order: declared position -> the gene it really holds
    actual = [matched.get(gene) for gene in declared]
    agreement = None
    if shrna_actual and all(a is not None for a in actual[:len(shrna_actual)]):
        agreement = actual[:len(shrna_actual)] == shrna_actual

    summary = {
        "n_signatures_used": len(sample),
        "n_declared_columns": len(declared),
        "n_matched_exactly": len(matched),
        "n_unmatched": len(unmatched),
        "match_is_a_bijection": len(set(matched.values())) == len(matched) == len(declared),
        "n_columns_already_correct": sum(1 for gene, image in matched.items()
                                         if gene == image),
        "declared_head": declared[:8],
        "compound_actual_order_head": actual[:8],
        "compound_actual_order": actual,
        "declared_order": declared,
        "shrna_actual_order_head": shrna_actual,
        "the_two_files_share_the_same_wrong_axis": agreement,
        "reading": ("a bijection with few columns in place means the compound matrix is "
                    "permuted like the shRNA one; whether the two orders agree decides "
                    "whether H3 paired two matrices in one coordinate system or two"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "two_permutations_compared.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_all_cohort_permutation_check():
    """Does the recovered permutation hold for every cohort signature, not a sample?

    The permutation was recovered from 400 signatures and the invariance claim now
    sits in a deviation record, so it is checked against all 41,643 cohort
    signatures: apply the inverse to the extraction and compare against the source,
    signature by signature, on the identifiers.
    """
    import importlib.util
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    records = json.loads(Path("/app/results/03_phenotype_projection/"
                              "phenotype_projection_results.json").read_text())
    drugs = sorted({r["drug"] for r in records})
    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    cohort = sorted(set(siginfo[siginfo.pert_iname.isin(set(drugs))].sig_id.astype(str)))

    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    declared = [str(g) for g in extraction["gene_ids"]]
    stored = extraction["signatures"]
    wanted = [s for s in cohort if s in held]

    # the permutation, recovered from a 400-signature block and then applied here
    block = wanted[:400]
    left = np.vstack([np.asarray(stored[held[s]]) for s in block]).astype(np.float32)
    values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, block, declared)
    row_of = {s: i for i, s in enumerate(signatures)}
    right = np.vstack([values[row_of[s]] for s in block]).astype(np.float32)
    profiles = {np.ascontiguousarray(right[:, i]).tobytes(): i for i in range(right.shape[1])}
    onto = [profiles.get(np.ascontiguousarray(left[:, i]).tobytes()) for i in range(len(declared))]
    if not (all(index is not None for index in onto)):
        raise AssertionError("the permutation did not recover")
    onto = np.array(onto)
    genes_in_source_order = list(genes)

    worst, failures, compared = 0.0, [], 0
    for start in range(0, len(wanted), 2000):
        chunk = wanted[start:start + 2000]
        values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, chunk,
                                                            genes_in_source_order)
        if not (genes == genes_in_source_order):
            raise AssertionError("the source returned a different gene order")
        row_of = {s: i for i, s in enumerate(signatures)}
        for sig_id in chunk:
            mine = np.asarray(stored[held[sig_id]], dtype=np.float64)
            theirs = values[row_of[sig_id]].astype(np.float64)
            error = float(np.abs(theirs[onto] - mine).max())
            worst = max(worst, error)
            compared += 1
            if error > gate.ATOL + gate.RTOL * 10.0 and len(failures) < 10:
                failures.append({"sig_id": sig_id, "max_abs_difference": error})

    summary = {
        "n_cohort_signatures_compared": compared,
        "n_recovered_from": len(block),
        "max_abs_difference_under_the_permutation": worst,
        "n_failures": len(failures),
        "failures": failures,
        "one_permutation_holds_for_every_cohort_signature": not failures,
        "reading": ("the permutation was recovered from a 400-signature block and applied to "
                    "every cohort signature; holding throughout makes the invariance claim "
                    "general rather than sampled"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "all_cohort_permutation_check.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps({k: v for k, v in summary.items() if k != "failures"}, indent=2)


@app.function(**COMMON)
def stage_core_gene_identities():
    """Which genes are the HDAC core, under the declared labels and under the true ones?

    `drug-perturbation-geometry/experiments/03_core_defenses.py:321-325` selects a
    core by column-wise statistics and then names it with `gene_ids[i]` from the
    extraction. The statistics travel with the column, so the core set is right and
    the names are not. This recomputes the core for the pan-HDAC drugs the paper
    describes and reports both namings, so the corrected identities exist rather
    than only the knowledge that the printed ones are wrong.
    """
    import importlib.util
    import json
    from collections import Counter
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    labels = json.loads((Path(staged) / "frozen_drug_labels.json").read_text())
    hdac = sorted({entry["pert_iname"] for entry in labels["drugs"]
                   if "HDAC inhibitor" in str(entry.get("moa", ""))})

    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    siginfo = siginfo[siginfo.pert_iname.isin(set(hdac))]
    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    declared = [str(g) for g in extraction["gene_ids"]]
    stored = extraction["signatures"]

    # the permutation, recovered from the source exactly as before
    wanted = sorted(set(siginfo.sig_id.astype(str)) & set(held))
    block = wanted[:400]
    left = np.vstack([np.asarray(stored[held[s]]) for s in block]).astype(np.float32)
    values, source_genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, block, declared)
    row_of = {s: i for i, s in enumerate(signatures)}
    right = np.vstack([values[row_of[s]] for s in block]).astype(np.float32)
    profiles = {np.ascontiguousarray(right[:, i]).tobytes(): i for i in range(right.shape[1])}
    truth = [source_genes[profiles[np.ascontiguousarray(left[:, i]).tobytes()]]
             for i in range(len(declared))]

    symbols = pd.read_csv(raw / "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="\t",
                          low_memory=False)
    symbol = dict(zip(symbols.pr_gene_id.astype(str), symbols.pr_gene_symbol))

    # the paper's core rule, per drug, over its per-cell mean signatures
    counts_declared, counts_true, per_drug = Counter(), Counter(), {}
    for drug, group in siginfo.groupby("pert_iname"):
        per_cell = []
        for cell, members in group.groupby("cell_id"):
            ids = sorted({s for s in members.sig_id.astype(str) if s in held})
            if ids:
                per_cell.append(np.vstack([np.asarray(stored[held[s]], dtype=np.float64)
                                           for s in ids]).mean(axis=0))
        if len(per_cell) < 5:
            continue
        sigs = np.vstack(per_cell)
        consistency = np.mean(np.sign(sigs) == np.sign(sigs.mean(axis=0, keepdims=True)), axis=0)
        cv = np.std(np.abs(sigs), axis=0) / np.maximum(np.mean(np.abs(sigs), axis=0), 1e-10)
        effect = np.abs(sigs.mean(axis=0))
        mask = (consistency > 0.8) & (cv < 1.0) & (effect > 0.5)
        positions = np.where(mask)[0]
        per_drug[drug] = int(mask.sum())
        for i in positions:
            counts_declared[declared[i]] += 1
            counts_true[truth[i]] += 1

    n_drugs = len(per_drug)
    shared_declared = sorted(g for g, c in counts_declared.items() if c == n_drugs)
    shared_true = sorted(g for g, c in counts_true.items() if c == n_drugs)
    named_in_the_paper = ["SUV39H1", "MYC", "CDK6", "BIRC5", "ORC1"]
    declared_symbols = sorted(symbol.get(g, g) for g in shared_declared)
    true_symbols = sorted(symbol.get(g, g) for g in shared_true)

    summary = {
        "n_hdac_drugs_used": n_drugs,
        "core_size_per_drug": per_drug,
        "n_shared_core_declared": len(shared_declared),
        "n_shared_core_true": len(shared_true),
        "core_size_is_invariant": len(shared_declared) == len(shared_true),
        "shared_core_under_declared_labels": declared_symbols,
        "shared_core_under_the_true_axis": true_symbols,
        "genes_named_in_the_paper": named_in_the_paper,
        "named_genes_in_the_declared_core": [g for g in named_in_the_paper
                                             if g in declared_symbols],
        "named_genes_in_the_true_core": [g for g in named_in_the_paper if g in true_symbols],
        "n_names_unchanged": len(set(declared_symbols) & set(true_symbols)),
        "reading": ("the core is selected by column statistics, so its size is invariant and "
                    "its membership as a set of columns is too; only the names differ. The "
                    "paper's named genes are read off the declared labels."),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "core_gene_identities.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps({k: v for k, v in summary.items() if k != "core_size_per_drug"}, indent=2)


EXTENSION_COHORT = [
    "ABT-751", "CYT-997", "D-64131", "PJ-34", "SB-334867", "SB-408124",
    "cycloheximide", "emetine", "fenbendazole", "flubendazole", "homoharringtonine",
    "ketoconazole", "oxibendazole", "parbendazole", "salubrinal", "tipifarnib",
    "vindesine",
]


@app.function(**COMMON)        # cmapPy is already in the analysis image
def stage_extend_rebuild():
    """Add the 17 drugs the CRISPRi arm needs and the rebuild lacks.

    The rebuild covers the H3 and H4 cohort completely, 795 drugs over 258 targets,
    and is short of the 812-drug CRISPRi arm by these seventeen. They are built the
    way `modal_h3_rebuild.py` builds every other shard: the pinned GCTX, the frozen
    landmark order, reindexed after the parse because `rid=` selects rows without
    ordering them, and the per-drug per-cell mean the deposited artifact takes.
    """
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd
    from cmapPy.pandasGEXpress import parse

    raw = Path("/rebuild/raw")
    shards = Path("/rebuild/shards")
    gene_info = pd.read_csv(raw / "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="\t",
                            low_memory=False)
    landmark = gene_info[gene_info.pr_is_lm == 1].sort_values("pr_gene_id")
    ids = [str(g) for g in landmark.pr_gene_id]
    if not (len(ids) == 978):
        raise AssertionError(f"{len(ids)} landmark genes")

    existing = set()
    for shard in sorted(shards.glob("shard_*.npz")):
        with np.load(shard, allow_pickle=True) as data:
            existing |= {k for k in data.files if k not in ("fingerprint", "__cells__")}
    wanted = [d for d in EXTENSION_COHORT if d not in existing]
    print(f"{len(existing)} drugs present; extending by {len(wanted)}", flush=True)
    if not wanted:
        return "nothing to add"

    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    sub = siginfo[siginfo.pert_iname.isin(set(wanted)) & siginfo.pert_iname.notna()]
    gct = parse.parse(str(raw / GCTX), cid=sorted(set(sub.sig_id.astype(str))), rid=ids)
    if not (gct.data_df.shape[0] == 978):
        raise AssertionError(f"parsed {gct.data_df.shape[0]} rows")
    if not (gct.data_df.columns.is_unique):
        raise AssertionError("parsed matrix has duplicate signature ids")
    frame = gct.data_df
    frame.index = frame.index.astype(str)
    frame = frame.reindex(index=ids)               # rid= selects, it does not order
    if not (not frame.isna().any().any()):
        raise AssertionError("a landmark gene is missing from the parse")
    if not (list(frame.index) == ids):
        raise AssertionError("landmark order is not the frozen order")
    mat = frame.T

    payload, cells_by_drug = {}, {}
    for drug, group in sub.groupby("pert_iname"):
        per_cell, cell_ids = [], []
        for cell_id, rows in group.groupby("cell_id"):
            sids = [s for s in sorted(set(rows.sig_id.astype(str))) if s in mat.index]
            if sids:
                per_cell.append(np.vstack([mat.loc[s].to_numpy(np.float64)
                                           for s in sids]).mean(axis=0))
                cell_ids.append(str(cell_id))
        if per_cell:
            payload[drug] = np.vstack(per_cell)
            cells_by_drug[drug] = cell_ids

    out = shards / "shard_ext_000.npz"
    np.savez_compressed(out, fingerprint=np.array("crispri-extension-2026-10-05"),
                        __cells__=np.array(json.dumps(cells_by_drug)), **payload)
    rebuild.commit()

    summary = {"n_requested": len(wanted), "n_written": len(payload),
               "not_found_in_the_metadata": sorted(set(wanted) - set(payload)),
               "cells_per_drug": {d: len(c) for d, c in cells_by_drug.items()},
               "shard": str(out), "gene_axis": "the frozen landmark order, reindexed",
               "reading": ("these drugs are built by the same route as every other shard, "
                           "so the gate can verify them against the source alongside the rest")}
    result_dir = Path("/out/03c_h3_sensitivity")
    result_dir.mkdir(parents=True, exist_ok=True)
    (result_dir / "rebuild_extension.json").write_text(json.dumps(summary, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_gate_on_cohort(cohort_path: str):
    """The gate against a cohort other than H3's, so the extension can be verified.

    `stage_gate` runs the registered H3 cohort. The CRISPRi arm is a wider set, 812
    drugs, and the seventeen added to the rebuild sit outside H3's 795, so they are
    not covered by a run of the registered gate. This runs the same gate against the
    cohort given.
    """
    from datetime import datetime
    from pathlib import Path

    run = f"gate_cohort_{datetime.now():%Y-%m-%dT%H%M%S}"
    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    _run(["/app/experiments/03e_reconstruction_gate.py",
          "--rebuilt", str(rebuilt),
          "--extraction", str(staged),
          "--gctx", f"/rebuild/raw/{GCTX}",
          "--shrna-siginfo", "/rebuild/raw/lincs_shrna_siginfo.csv.gz",
          "--compound-siginfo", "/rebuild/raw/GSE92742_Broad_LINCS_sig_info.txt.gz",
          "--gene-info", "/rebuild/raw/GSE92742_Broad_LINCS_gene_info.txt.gz",
          "--cohort", cohort_path,
          "--output", f"/out/03c_h3_sensitivity/{run}"])
    results.commit()
    return (Path(f"/out/03c_h3_sensitivity/{run}/reconstruction_gate.json")).read_text()[:2500]


@app.function(**COMMON)
def stage_crispri_routes_checked():
    """Is the shRNA-versus-CRISPRi reversal an artifact of the mislabeled gene axis?

    The CRISPRi reference is placed by symbol: position i receives the Perturb-seq
    value for whichever gene the extraction's declared labels name at i. The drug
    signature at position i holds the gene that is actually there. Where those
    differ, the two operands are in different coordinate systems.

    Four routes, on one fixed cohort, with the permutation recovered from the source
    before any CRISPRi quantity is computed:

    1. deposited      declared symbols, extraction signatures. Must reproduce.
    2. corrected      true symbols, extraction signatures. Both operands aligned.
    3. from the source  true symbols, rebuild signatures. Must agree with 2.
    4. control        declared symbols, signatures permuted to the declared axis.
                      Applying the permutation to both operands must return route 1.
    """
    import importlib.util
    import json
    from pathlib import Path

    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import stats

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)

    raw = Path("/rebuild/raw")
    extraction = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    declared = [str(g) for g in extraction["gene_ids"]]
    stored = extraction["signatures"]

    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          low_memory=False)
    labels = json.loads((Path(staged) / "frozen_drug_labels.json").read_text())
    targets = {e["pert_iname"]: str(e["target"]).split("|")[0].strip()
               for e in labels["drugs"] if e.get("target")}

    # ---- the true axis, recovered from the source, before any CRISPRi value -----
    cohort_sigs = sorted(set(siginfo[siginfo.pert_iname.isin(set(targets))]
                             .sig_id.astype(str)) & set(held))
    block = cohort_sigs[:400]
    left = np.vstack([np.asarray(stored[held[s]]) for s in block]).astype(np.float32)
    values, source_genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, block, declared)
    row_of = {s: i for i, s in enumerate(signatures)}
    right = np.vstack([values[row_of[s]] for s in block]).astype(np.float32)
    profile = {np.ascontiguousarray(right[:, i]).tobytes(): i for i in range(right.shape[1])}
    truth = [source_genes[profile[np.ascontiguousarray(left[:, i]).tobytes()]]
             for i in range(len(declared))]
    if not (len(set(truth)) == len(truth)):
        raise AssertionError("the recovered axis repeats a gene")
    # position in the declared axis of each true gene, for the control route
    declared_position = {gene: i for i, gene in enumerate(declared)}
    to_declared = np.array([declared_position[g] for g in truth])

    info = pd.read_csv(raw / "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="\t",
                       low_memory=False)
    symbol = dict(zip(info.pr_gene_id.astype(str), info.pr_gene_symbol.astype(str)))
    declared_symbols = [symbol.get(g, g) for g in declared]
    true_symbols = [symbol.get(g, g) for g in truth]

    # ---- the CRISPRi reference, placed on a given symbol list -------------------
    adata = ad.read_h5ad("/extraction/ReplogleWeissman2022_K562_essential.h5ad")
    perturbseq_genes = list(adata.var_names)
    control = adata.obs["gene"] == "non-targeting"
    if not (control.sum() > 0):
        raise AssertionError("no non-targeting cells")
    control_mean = np.asarray(adata[control].X.mean(axis=0)).ravel()

    def crispri_reference(symbols):
        order = {s: i for i, s in enumerate(symbols)}
        shared = [g for g in perturbseq_genes if g in order]
        columns = [perturbseq_genes.index(g) for g in shared]
        rows = [order[g] for g in shared]
        out = {}
        for gene in sorted(set(adata.obs["gene"]) - {"non-targeting"}):
            cells = adata.obs["gene"] == gene
            if cells.sum() < 10:
                continue
            delta = np.asarray(adata[cells].X.mean(axis=0)).ravel() - control_mean
            vector = np.zeros(len(symbols))
            vector[rows] = delta[columns]
            if np.linalg.norm(vector) > 0:
                out[gene] = vector
        return out, len(shared)

    # ---- the per-drug matrices, four ways ---------------------------------------
    def per_drug(reorder=None, from_rebuild=False):
        matrices = {}
        if from_rebuild:
            shards = Path("/rebuild/shards")
            cells = {}
            for shard in sorted(shards.glob("shard_*.npz")):
                with np.load(shard, allow_pickle=True) as data:
                    blob = json.loads(str(data["__cells__"]))
                    for key in data.files:
                        if key not in ("fingerprint", "__cells__"):
                            matrices[key] = np.asarray(data[key], dtype=np.float64)
                    cells.update(blob)
            return matrices
        sub = siginfo[siginfo.pert_iname.isin(set(targets))]
        for drug, group in sub.groupby("pert_iname"):
            per_cell = []
            for cell, rows_ in group.groupby("cell_id"):
                ids = [s for s in sorted(set(rows_.sig_id.astype(str))) if s in held]
                if ids:
                    per_cell.append(np.vstack([np.asarray(stored[held[s]], dtype=np.float64)
                                               for s in ids]).mean(axis=0))
            if per_cell:
                matrix = np.vstack(per_cell)
                matrices[drug] = matrix if reorder is None else matrix[:, reorder]
        return matrices

    def instability(matrix):
        unit = matrix / np.linalg.norm(matrix, axis=1, keepdims=True)
        upper = np.triu_indices(matrix.shape[0], 1)
        return float(1.0 - (unit @ unit.T)[upper].mean())

    def projected(matrix, direction):
        upper = np.triu_indices(matrix.shape[0], 1)
        unit = direction / np.linalg.norm(direction)
        return float(np.abs((matrix[upper[0]] - matrix[upper[1]]) @ unit).mean())

    def cosine_squared(mean_sig, direction):
        a, b = np.linalg.norm(mean_sig), np.linalg.norm(direction)
        if a < 1e-10 or b < 1e-10:
            return 0.0
        return float(mean_sig @ direction / (a * b)) ** 2

    def run(symbols, matrices, name, both=None):
        """`both` applies one permutation to the drug matrices and the reference
        alike, which cannot change a cosine and is therefore the control."""
        reference, n_shared = crispri_reference(symbols)
        if both is not None:
            reference = {t: v[both] for t, v in reference.items()}
            matrices = {d: m[:, both] for d, m in matrices.items()}
        records = []
        for drug, matrix in matrices.items():
            target = targets.get(drug)
            if target is None or matrix.shape[0] < 5 or target not in reference:
                continue
            direction = reference[target]
            records.append({"drug": drug, "target": target,
                            "raw": instability(matrix),
                            "proj": projected(matrix, direction),
                            "enrich": cosine_squared(matrix.mean(axis=0), direction)})
        proj = stats.spearmanr([r["proj"] for r in records], [r["enrich"] for r in records])
        raw_r = stats.spearmanr([r["raw"] for r in records], [r["enrich"] for r in records])
        return {"route": name, "n": len(records), "n_shared_genes": n_shared,
                "n_targets": len({r["target"] for r in records}),
                "projected_vs_enrichment": [float(proj.statistic), float(proj.pvalue)],
                "raw_vs_enrichment": [float(raw_r.statistic), float(raw_r.pvalue)],
                "records": records}

    # the rebuild sits in the frozen landmark order: the landmark ids ascending by
    # Entrez id, which is what modal_h3_rebuild.py parses against
    frozen = [str(g) for g in info[info.pr_is_lm == 1].sort_values("pr_gene_id").pr_gene_id]
    if not (len(frozen) == 978):
        raise AssertionError(f"{len(frozen)} landmark genes")
    frozen_symbols = [symbol.get(g, g) for g in frozen]

    rng = np.random.default_rng()
    shuffle = rng.permutation(len(declared))
    routes = [
        run(declared_symbols, per_drug(), "1 deposited: declared symbols, extraction"),
        run(true_symbols, per_drug(), "2 corrected: true symbols, extraction"),
        run(frozen_symbols, per_drug(from_rebuild=True), "3 from the source: rebuild"),
        run(true_symbols, per_drug(), "4 control: route 2 with one permutation on both",
            both=shuffle),
    ]

    deposited = {"projected": -0.12470506592643997, "raw": 0.2722681898254417, "n": 131}
    summary = {
        "deposited_values": deposited,
        "routes": [{k: v for k, v in r.items() if k != "records"} for r in routes],
        "route_1_reproduces_the_deposited_result": bool(
            abs(routes[0]["projected_vs_enrichment"][0] - deposited["projected"]) < 0.01),
        "control_returns_route_2": bool(
            abs(routes[3]["projected_vs_enrichment"][0]
                - routes[1]["projected_vs_enrichment"][0]) < 1e-9),
        "routes_2_and_3_agree": bool(
            abs(routes[1]["projected_vs_enrichment"][0]
                - routes[2]["projected_vs_enrichment"][0]) < 0.05),
        "reading": ("route 1 must reproduce the deposited number before anything is "
                    "concluded; route 4 applies one random permutation to both operands "
                    "of route 2 and must return route 2 exactly, since no cosine can "
                    "change under that; routes 2 and 3 are the corrected result reached "
                    "by two independent paths and must agree"),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "crispri_three_routes.json").write_text(json.dumps(
        {**summary, "per_route_records": {r["route"]: r["records"] for r in routes}}, indent=2))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_preflight_checked():
    """Load and validate every input R0-R7 will read, and compute no registered statistic.

    Amendment 3 pins the rebuild production reads, and the pin table is empty until
    the values exist. This reports them from the volume, which is what production
    reads rather than a retrieved copy, and exercises every loader and contract so
    a failure surfaces here rather than part way through the registered run.
    """
    import hashlib
    import importlib.util
    import json
    import site
    from pathlib import Path

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    site.addsitedir("/app")          # 03d imports geometry, as `_run` arranges
    spec = importlib.util.spec_from_file_location(
        "discordance", "/app/experiments/03d_h3_reference_discordance.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    def digest(path):
        sha = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 22), b""):
                sha.update(block)
        return sha.hexdigest()

    rebuild_hashes = {}
    for name in ("landmark_gene_ids.json", "shrna_consensus.npz", "shrna_signatures.npz",
                 "rebuild_manifest.json", "cohort_bundle.npz"):
        candidate = Path(rebuilt) / name
        if candidate.exists():
            rebuild_hashes[name] = digest(candidate)
    shards = sorted(Path(rebuilt).glob("shard_*.npz"))
    rebuild_hashes["shards"] = {s.name: digest(s) for s in shards}

    # every loader and contract, with nothing computed from them
    per_drug, gene_axis = mod.build_drug_signatures_from_rebuild(rebuilt)
    frozen = mod.frozen_landmark_order(
        Path(staged) / "GSE92742_Broad_LINCS_gene_info.txt.gz")
    mod.require(gene_axis == frozen, "the rebuild's axis is not the frozen landmark order")
    shrna = mod.build_shrna_reference_from_rebuild(rebuilt, gene_axis)

    # Amendment 4: the bases artifact, through the shared validator, and every
    # construction actually built on it, so a basis that does not describe its
    # source fails here rather than part way through the registered run
    symbols = mod.landmark_symbols(
        Path(staged) / "GSE92742_Broad_LINCS_gene_info.txt.gz", gene_axis)
    bases = mod.validated_analysis_bases(
        "/app/registry/frozen/analysis_bases.json", gene_axis, mod.analysis_bases_pin(),
        generator="/app/experiments/03i_freeze_analysis_bases.py")
    built = {}
    for name, filename in (("C1-K562", "K562_essential_normalized_bulk_01.h5ad"),
                           ("C1-RPE1", "rpe1_normalized_bulk_01.h5ad"),
                           ("C1-GW", "K562_gwps_normalized_bulk_01.h5ad")):
        reference = mod.load_replogle_bulk(
            Path("/inputs/replogle2022") / filename, symbols, name, basis=bases[name])
        built[name] = {"n_landmarks": int(len(reference.positions)),
                       "n_directions": len(reference.directions),
                       "basis_sha256": bases[name]["basis_sha256"],
                       "source_sha256": bases[name]["file_sha256"]}
    records = json.loads(Path("/app/results/03b_h3_crispri/h3_crispri_results.json").read_text())
    cohort = mod.check_cohort_identity(records, per_drug)

    cells = sorted({cell for drug in per_drug.values() for cell in drug})
    summary = {
        "rebuild_sha256": rebuild_hashes,
        "analysis_bases": {
            "file_sha256": digest("/app/registry/frozen/analysis_bases.json"),
            "constructions": built,
            "registered_comparisons": json.loads(
                Path("/app/registry/frozen/analysis_bases.json").read_text()
            )["registered_comparisons"]},
        "cohort_identity": cohort,
        "n_drugs_loaded": len(per_drug),
        "n_cell_lines_seen": len(cells),
        "n_shrna_targets": len(shrna.directions),
        "gene_axis_is_the_frozen_order": gene_axis == frozen,
        "n_shards": len(shards),
        "computed_no_registered_statistic": True,
        "reading": ("every loader and contract R0-R7 depends on ran without computing a "
                    "registered quantity. The hashes here are the values Amendment 3's pin "
                    "table takes, read from the volume production reads."),
    }
    out = Path("/out/03c_h3_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "preflight.json").write_text(json.dumps(summary, indent=2))
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


@app.function(**COMMON)
def stage_pin_legacy_axis_map_r4():
    """Build the legacy coordinate map Amendment 3's B6 pins, direction verified here.

    The map says, for each column of the retained extraction, which index of its
    declared axis holds the gene that column actually carries. Its direction is the
    part I have already got wrong once, by zipping a reconstruction against the
    extraction's declared genes instead of the rebuild's, so nothing here rests on a
    field name: the recovery is the one `stage_all_cohort_permutation_check` already
    ran to max 0.0 over every cohort signature, and the map is then applied and
    required to reproduce the source exactly. A reversed map fails that step rather
    than passing quietly.

    Float32 on both sides before the byte comparison. A float64 view compared against
    float32 source bytes returns zero matches whatever the data.
    """
    import hashlib
    import importlib.util
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = Path(_stage_inputs())     # the extraction and its metadata
    raw = Path("/rebuild/raw")         # the pinned GCTX, which is not in that directory
    spec = importlib.util.spec_from_file_location(
        "gate", "/app/experiments/03e_reconstruction_gate.py")
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    out = Path("/out/03c_h3_sensitivity/legacy_axis_map")
    out.mkdir(parents=True, exist_ok=True)

    extraction_path = staged / "lincs_subset.npz"
    extraction = np.load(extraction_path, allow_pickle=True)
    held = {str(s): i for i, s in enumerate(extraction["sig_ids"])}
    declared = [str(g) for g in extraction["gene_ids"]]
    stored = extraction["signatures"]
    if len(set(declared)) != len(declared):
        raise AssertionError("the extraction's declared axis repeats an identifier")

    siginfo = pd.read_csv(staged / "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="\t",
                          usecols=["sig_id", "pert_iname"], low_memory=False)
    drugs = {r["drug"] for r in json.loads(
        Path("/app/registry/cohorts/cohort_812_compound.json").read_text())["records"]}
    cohort = sorted(set(siginfo[siginfo.pert_iname.isin(drugs)].sig_id.astype(str)))
    wanted = [s for s in cohort if s in held]
    if not wanted:
        raise AssertionError("no cohort signature is present in the extraction")

    block = wanted[:400]
    left = np.vstack([np.asarray(stored[held[s]]) for s in block]).astype(np.float32)
    values, genes, signatures, _ = gate.read_gctx_slice(raw / GCTX, block, declared)
    row_of = {s: i for i, s in enumerate(signatures)}
    right = np.vstack([values[row_of[s]] for s in block]).astype(np.float32)
    profiles = {np.ascontiguousarray(right[:, i]).tobytes(): i
                for i in range(right.shape[1])}
    onto = [profiles.get(np.ascontiguousarray(left[:, i]).tobytes())
            for i in range(len(declared))]
    if any(index is None for index in onto):
        raise AssertionError(
            f"{sum(1 for i in onto if i is None)} extraction columns matched no source "
            "column by exact equality, so the map is not recoverable")
    if sorted(onto) != list(range(len(declared))):
        raise AssertionError("the recovered map is not a total bijection")

    # onto[i] indexes the source's own gene order, so the gene truly in extraction
    # column i is genes[onto[i]]; the map wants that gene's index in the declared axis
    declared_index = {gene: i for i, gene in enumerate(declared)}
    order = [declared_index[str(genes[onto[i]])] for i in range(len(declared))]
    if sorted(order) != list(range(len(declared))):
        raise AssertionError("the map over the declared axis is not a total bijection")
    inverse = [0] * len(order)
    for actual, j in enumerate(order):
        inverse[j] = actual
    if [inverse[order[i]] for i in range(len(order))] != list(range(len(order))):
        raise AssertionError("the map does not round-trip to the identity")

    genes_in_source_order = list(genes)
    worst, compared = 0.0, 0
    for start in range(0, len(wanted), 2000):
        chunk = wanted[start:start + 2000]
        values, chunk_genes, signatures, _ = gate.read_gctx_slice(
            raw / GCTX, chunk, genes_in_source_order)
        if chunk_genes != genes_in_source_order:
            raise AssertionError("the source returned a different gene order")
        row_of = {s: i for i, s in enumerate(signatures)}
        for sig_id in chunk:
            mine = np.asarray(stored[held[sig_id]], dtype=np.float64)
            theirs = values[row_of[sig_id]].astype(np.float64)
            worst = max(worst, float(np.abs(theirs[onto] - mine).max()))
            compared += 1
        json.dump({"compared": compared, "max_abs_difference": worst},
                  open(out / "progress.json", "w"))      # RULE ONE: inside the loop
        results.commit()
    if worst != 0.0:
        raise AssertionError(f"the map does not reproduce the source exactly: {worst}")

    canonical = json.dumps(order, sort_keys=True, separators=(",", ":"))
    payload = {
        "name": "legacy_axis_map", "schema_version": 1,
        "declared_index_of_each_actual_column": order,
        "map_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "legacy_sha256_observed": hashlib.sha256(extraction_path.read_bytes()).hexdigest(),
        "recovered_from": {"n_signatures_in_the_recovery_block": len(block),
                           "n_cohort_signatures_confirmed": compared,
                           "max_abs_difference_under_the_map": worst,
                           "source": "the pinned GSE92742 GCTX, read by this stage"},
        "direction": ("declared_index_of_each_actual_column[i] == j means column i of "
                      "the extraction's matrix holds the gene its declared axis lists "
                      "at index j"),
    }
    # no expected hash is written here: an artifact that carries its own expected hash
    # attests to itself. The expected values go into registry/frozen/, which Amendment 3
    # pins, and the gate reads them from there.
    (out / "legacy_axis_map.json").write_text(json.dumps(payload, indent=2) + "\n")
    results.commit()
    print(json.dumps({k: v for k, v in payload.items()
                      if k != "declared_index_of_each_actual_column"}, indent=2))
    return payload["map_sha256"]


@app.function(**COMMON)
def stage_gate_v2_on_scope(scope: str):
    """Gate v2 on one cohort scope, into its own immutable directory.

    Amendment 3's B5 requires both scopes. The report filename is fixed, so each run
    gets its own directory and neither overwrites the other; the gate names the
    directory itself from the scope and the time it started.

    Every path here is copied from `stage_gate_on_cohort`, which runs, rather than
    written from what seemed reasonable.
    """
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
          "--cohort-manifest", "/app/registry/cohorts/cohort_812_compound.json",
          "--paired-manifest", "/app/registry/cohorts/cohort_795_shrna_paired.json",
          "--eligible-targets", "/app/registry/cohorts/shrna_eligible_targets.json",
          "--excluded-records", "/app/registry/cohorts/shrna_excluded_records.json",
          "--axis-map", "/app/registry/frozen/legacy_axis_map.json",
          "--frozen-identities", "/app/registry/frozen/expected_identities.json",
          "--coverage-report",
          "/app/results/03c_h3_sensitivity/shrna_coverage_identity_v2.json",
          "--deposited-records", "/app/results/03b_h3_crispri/h3_crispri_results.json",
          "--cohort-scope", scope,
          "--output", "/out/03c_h3_sensitivity"])
    results.commit()
    reports = sorted(Path("/out/03c_h3_sensitivity").glob(
        f"gate_v2_{scope}_*/reconstruction_gate.json"))
    return reports[-1].read_text()[:3000] if reports else "no report was written"


@app.function(**COMMON)
def stage_r0_to_r7():
    """R0 to R7 and the R5 response, in one invocation of the registered driver.

    The driver runs once, as its registration requires, so no module is computed
    twice under different invocation states. Amendment 3 is frozen at 03175d0 and
    gate v2 holds all five predicates on both scopes, which is what B5 requires
    before this runs.

    Paths are copied from the stages that already run. `_repo_at_its_absolute_path`
    is what makes the driver's absolute CRISPRI_RECORDS path resolve.
    """
    from pathlib import Path

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    bases = Path("/app/registry/frozen/analysis_bases.json")
    if not bases.exists():
        raise AssertionError(
            "registry/frozen/analysis_bases.json is absent from the image, and Amendment 4 "
            "evaluates every C1 construction on a frozen basis. Run --stage bases, fetch the "
            "result into registry/frozen/, commit it, and redeploy.")
    # Amendment 4's recomputation, under its own name: the full-axis audit it
    # supersedes is pinned by sha256 in the amendment's foreknowledge section and
    # stays where that section names it
    audit = Path("/app/results/03d_h3_reference_discordance/"
                 "single_cell_audit_amendment_4.json")
    if not audit.exists():
        raise AssertionError(
            "single_cell_audit_amendment_4.json is absent from the image, and R0.3, R0.4 "
            "and R0.6 are registered analyses; R0.4 also decides how R0.6 and R7c may be "
            "read. Run experiments/modal_03d_single_cell.py --stage preflight, then "
            "--stage targets and --stage merge, fetch the result into "
            "results/03d_h3_reference_discordance/ under that name, and redeploy.")

    # Amendment 1's A4 splits R5 in two: stage one freezes the drug-to-treatment
    # mapping before any response value is read, and the response stage consumes that
    # table rather than building one. The reviewed table is committed in the
    # repository, so it is placed where the response stage looks and is never rebuilt
    # here: rebuilding it inside the run is the thing A4 exists to prevent.
    frozen_mapping = Path("/app/results/03d_h3_reference_discordance/r5_prism_mapping.json")
    if not frozen_mapping.exists():
        raise AssertionError(
            "r5_prism_mapping.json is absent from the image. Run --r5-stage mapping, "
            "review the table, commit it to results/03d_h3_reference_discordance/, and "
            "redeploy. The response stage does not build a mapping.")
    digest = _sha256_of(frozen_mapping)
    if digest != MAPPING_SHA256:
        raise AssertionError(
            f"the committed R5 mapping hashes to {digest} and the registered table is "
            f"{MAPPING_SHA256}; the response stage reads the frozen table or none")
    out_dir = Path("/out/03d_h3_reference_discordance")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "r5_prism_mapping.json").write_bytes(frozen_mapping.read_bytes())
    print(f"[{_ts()}] R5 mapping placed, sha256 {digest}", flush=True)

    _run(["/app/experiments/03d_h3_reference_discordance.py",
          "--r5-stage", "response",
          "--data", str(staged),
          "--rebuilt", str(rebuilt),
          "--perturbseq", "/inputs/scperturb/ReplogleWeissman2022_K562_essential.h5ad",
          "--replogle", "/inputs/replogle2022",
          "--prism", "/inputs/prism_19q4",
          "--single-cell", str(audit),
          "--analysis-bases", "/app/registry/frozen/analysis_bases.json",
          # the sealed audit, authenticated against registry/frozen/reproduction_legacy_pin.json
          # before it is parsed. The two-layer measurement that preceded it carries no
          # provenance and no per-drug table, so the gate refuses it by design.
          "--legacy-reproduction",
          "/app/results/03d_h3_reference_discordance/legacy_reproduction_audit.json",
          "--expected-mapping-sha256", MAPPING_SHA256,
          "--output", "/out/03d_h3_reference_discordance"])
    results.commit()
    written = sorted(Path("/out/03d_h3_reference_discordance").glob("*.json"))
    return [p.name for p in written]


@app.function(**COMMON)
def stage_freeze_analysis_bases():
    """Measure each C1 construction's finite landmark basis, on the volume's own files.

    Amendment 4's rule is measured rather than asserted, and it is measured on the
    copies production reads: a coverage report computed against a different copy of
    the same artifact is how Deviation 12's neighbourhood was found.
    """
    from pathlib import Path

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    out = Path("/out/registry")
    out.mkdir(parents=True, exist_ok=True)
    _run(["/app/experiments/03i_freeze_analysis_bases.py",
          "--replogle", "/inputs/replogle2022",
          "--gene-info", str(staged / "GSE92742_Broad_LINCS_gene_info.txt.gz"),
          "--perturbseq", "/inputs/scperturb/ReplogleWeissman2022_K562_essential.h5ad",
          "--output", str(out / "analysis_bases.json")])
    results.commit()
    return (out / "analysis_bases.json").read_text()


@app.function(**COMMON)
def stage_diagnose_reproduction_tolerance():
    """Why P_shrna reproduces to 1.14e-06 against a flat 1e-06, per drug.

    The driver's own gate refused and this asks what the refusal is made of. It
    computes no registered statistic: the per-drug P, E and D it compares are the
    ones the driver builds before its gate fires, and what is reported here is the
    reproduction error against the deposited records, disaggregated from the
    maximum the gate already took. No correlation, interval or reading is computed.

    The question is whether the error is proportional to the magnitude of the
    quantity, which is what a float32 round trip looks like, or scattered, which
    is what a wrong value looks like.
    """
    import importlib.util
    import json
    import site
    from pathlib import Path

    import numpy as np

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    site.addsitedir("/app")
    spec = importlib.util.spec_from_file_location(
        "discordance", "/app/experiments/03d_h3_reference_discordance.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    per_drug, gene_axis = mod.build_drug_signatures_from_rebuild(rebuilt)
    shrna = mod.build_shrna_reference_from_rebuild(rebuilt, gene_axis)
    targets, _, _ = mod.drug_targets(staged)
    records = json.loads(Path(
        "/app/results/03b_h3_crispri/h3_crispri_results.json").read_text())
    arm = mod.assemble_arm("shRNA arm", per_drug, targets, shrna,
                           drugs={r["drug"] for r in records if "proj_shrna" in r})
    values = mod.quantities(arm, shrna)

    by_drug = {r["drug"]: r for r in records}
    rows = []
    for i, drug in enumerate(arm.drugs):
        record = by_drug[drug]
        for quantity, computed, deposited in (
                ("P_shrna", values["P"][i], record["proj_shrna"]),
                ("E_shrna", values["E"][i], record["enrich_shrna"]),
                ("D", values["D"][i], record["raw_instability"])):
            error = abs(float(computed) - float(deposited))
            rows.append({"drug": drug, "quantity": quantity,
                         "computed": float(computed), "deposited": float(deposited),
                         "absolute_error": error,
                         "relative_error": error / abs(float(deposited))
                         if float(deposited) != 0 else None,
                         "float32_ulp_here": abs(float(deposited)) * 2.0 ** -24,
                         "error_in_ulps": error / (abs(float(deposited)) * 2.0 ** -24)
                         if float(deposited) != 0 else None})

    summary = {}
    for quantity in ("P_shrna", "E_shrna", "D"):
        group = [r for r in rows if r["quantity"] == quantity]
        errors = np.array([r["absolute_error"] for r in group])
        magnitudes = np.array([abs(r["deposited"]) for r in group])
        ulps = np.array([r["error_in_ulps"] for r in group if r["error_in_ulps"] is not None])
        over = [r for r in group if r["absolute_error"] >= 1e-6]
        finite = np.isfinite(errors) & np.isfinite(magnitudes) & (magnitudes > 0)
        summary[quantity] = {
            "n_drugs": len(group),
            "max_absolute_error": float(errors.max()),
            "max_relative_error": float((errors[finite] / magnitudes[finite]).max()),
            "max_error_in_float32_ulps": float(ulps.max()) if ulps.size else None,
            "median_error_in_float32_ulps": float(np.median(ulps)) if ulps.size else None,
            "n_at_or_over_flat_1e-6": len(over),
            "drugs_at_or_over_flat_1e-6": [
                {"drug": r["drug"], "deposited": r["deposited"],
                 "absolute_error": r["absolute_error"], "error_in_ulps": r["error_in_ulps"]}
                for r in sorted(over, key=lambda r: -r["absolute_error"])[:20]],
            "magnitude_range": [float(magnitudes.min()), float(magnitudes.max())],
            # a float32 round trip gives an error that tracks magnitude; a wrong
            # value does not, so the correlation is the discriminating measurement
            "pearson_error_against_magnitude": float(
                np.corrcoef(errors[finite], magnitudes[finite])[0, 1]),
            "elementwise_float32_rule_holds": bool(np.all(
                errors[finite] <= 1e-8 + 1e-6 * magnitudes[finite])),
        }

    out = {
        "what_this_is": (
            "why the driver's reproduction gate refused. RECON_TOL is a flat absolute "
            "1e-06 and the registration states the deposited records reproduce 'to within "
            "1e-6 for D, P and E'. P_shrna spans 0.7 to 17, and one float32 ULP at 17 is "
            "1.01e-06, so a flat 1e-06 is below float32 resolution at the top of that "
            "range. This measures whether the error tracks magnitude, as a float32 round "
            "trip does, or is scattered, as a wrong value is."),
        "computed_no_registered_statistic": True,
        "flat_tolerance_in_code": 1e-6,
        "summary": summary,
        "per_drug": rows,
    }
    directory = Path("/out/03d_h3_reference_discordance")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "reproduction_tolerance_diagnosis.json").write_text(
        json.dumps(out, indent=2, allow_nan=False))
    results.commit()
    return json.dumps(summary, indent=2)


@app.function(**COMMON)
def stage_reproduction_tolerance_power():
    """What size of real error each candidate tolerance would catch.

    A tolerance that admits the observed float32 discrepancy is only defensible if
    it still refuses a wrong value. This injects a systematic relative scaling
    error into the drug signatures, the shape a unit or normalization mistake
    takes, recomputes the reproduction errors, and reports the smallest injected
    error each rule rejects.

    No registered statistic is computed: the perturbed quantities are thrown away
    and only the pass or fail of each rule is reported.
    """
    import importlib.util
    import json
    import site
    from pathlib import Path

    import numpy as np

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    site.addsitedir("/app")
    spec = importlib.util.spec_from_file_location(
        "discordance", "/app/experiments/03d_h3_reference_discordance.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    per_drug, gene_axis = mod.build_drug_signatures_from_rebuild(rebuilt)
    shrna = mod.build_shrna_reference_from_rebuild(rebuilt, gene_axis)
    targets, _, _ = mod.drug_targets(staged)
    records = json.loads(Path(
        "/app/results/03b_h3_crispri/h3_crispri_results.json").read_text())
    eligible = {r["drug"] for r in records if "proj_shrna" in r}
    by_drug = {r["drug"]: r for r in records}

    def worst_errors(scale):
        scaled = {drug: {cell: value * scale for cell, value in cells.items()}
                  for drug, cells in per_drug.items()}
        arm = mod.assemble_arm("shRNA arm", scaled, targets, shrna, drugs=eligible)
        values = mod.quantities(arm, shrna)
        rows = []
        for i, drug in enumerate(arm.drugs):
            record = by_drug[drug]
            for quantity, computed, deposited in (
                    ("P_shrna", values["P"][i], record["proj_shrna"]),
                    ("E_shrna", values["E"][i], record["enrich_shrna"]),
                    ("D", values["D"][i], record["raw_instability"])):
                rows.append((quantity, abs(float(computed) - float(deposited)),
                             abs(float(deposited))))
        out = {}
        for quantity in ("P_shrna", "E_shrna", "D"):
            errors = np.array([e for q, e, _ in rows if q == quantity])
            magnitudes = np.array([m for q, _, m in rows if q == quantity])
            out[quantity] = {
                "max_absolute_error": float(errors.max()),
                "flat_1e-6_holds": bool(errors.max() < 1e-6),
                "elementwise_float32_holds": bool(
                    np.all(errors <= 1e-8 + 1e-6 * magnitudes)),
            }
        return out

    def worst_errors_one_context(epsilon):
        """One context per drug scaled, which is not scale-invariant.

        D and E are built from cosines and are unchanged by scaling every
        signature of a drug by the same factor, so a uniform ladder cannot
        exercise them. Perturbing a single context changes the pairwise cosines
        and reaches all three quantities.
        """
        scaled = {}
        for drug, cells in per_drug.items():
            order = sorted(cells)
            scaled[drug] = {cell: (value * (1.0 + epsilon) if cell == order[0] else value)
                            for cell, value in cells.items()}
        arm = mod.assemble_arm("shRNA arm", scaled, targets, shrna, drugs=eligible)
        values = mod.quantities(arm, shrna)
        rows = []
        for i, drug in enumerate(arm.drugs):
            record = by_drug[drug]
            for quantity, computed, deposited in (
                    ("P_shrna", values["P"][i], record["proj_shrna"]),
                    ("E_shrna", values["E"][i], record["enrich_shrna"]),
                    ("D", values["D"][i], record["raw_instability"])):
                rows.append((quantity, abs(float(computed) - float(deposited)),
                             abs(float(deposited))))
        out = {}
        for quantity in ("P_shrna", "E_shrna", "D"):
            errors = np.array([e for q, e, _ in rows if q == quantity])
            magnitudes = np.array([m for q, _, m in rows if q == quantity])
            out[quantity] = {
                "max_absolute_error": float(errors.max()),
                "flat_1e-6_holds": bool(errors.max() < 1e-6),
                "elementwise_float32_holds": bool(
                    np.all(errors <= 1e-8 + 1e-6 * magnitudes)),
            }
        return out

    def worst_errors_one_coordinate(epsilon):
        """One gene coordinate of one context moved, which changes its direction.

        A cosine is invariant to the scale of either argument, so D, the mean
        pairwise cosine, is unchanged by any rescaling of a signature. Only a
        perturbation that turns a signature reaches it.
        """
        scaled = {}
        for drug, cells in per_drug.items():
            order = sorted(cells)
            scaled[drug] = {}
            for cell, value in cells.items():
                if cell == order[0]:
                    moved = np.array(value, dtype=np.float64, copy=True)
                    moved[0] = moved[0] * (1.0 + epsilon) + epsilon
                    scaled[drug][cell] = moved
                else:
                    scaled[drug][cell] = value
        arm = mod.assemble_arm("shRNA arm", scaled, targets, shrna, drugs=eligible)
        values = mod.quantities(arm, shrna)
        rows = []
        for i, drug in enumerate(arm.drugs):
            record = by_drug[drug]
            for quantity, computed, deposited in (
                    ("P_shrna", values["P"][i], record["proj_shrna"]),
                    ("E_shrna", values["E"][i], record["enrich_shrna"]),
                    ("D", values["D"][i], record["raw_instability"])):
                rows.append((quantity, abs(float(computed) - float(deposited)),
                             abs(float(deposited))))
        out = {}
        for quantity in ("P_shrna", "E_shrna", "D"):
            errors = np.array([e for q, e, _ in rows if q == quantity])
            magnitudes = np.array([m for q, _, m in rows if q == quantity])
            out[quantity] = {
                "max_absolute_error": float(errors.max()),
                "flat_1e-6_holds": bool(errors.max() < 1e-6),
                "elementwise_float32_holds": bool(
                    np.all(errors <= 1e-8 + 1e-6 * magnitudes)),
            }
        return out

    ladders = {}
    for name, perturb in (("every_context_scaled", lambda e: worst_errors(1.0 + e)),
                          ("one_context_scaled", worst_errors_one_context),
                          ("one_coordinate_moved", worst_errors_one_coordinate)):
        ladder = {}
        for exponent in (0, 8, 7, 6, 5, 4, 3):
            epsilon = 0.0 if exponent == 0 else 10.0 ** -exponent
            label = "unperturbed" if epsilon == 0 else f"relative {epsilon:.0e}"
            ladder[label] = {"injected_relative_error": epsilon,
                             **_as_plain(perturb(epsilon))}
            print(f"[{_ts()}] {name} {label}: "
                  + ", ".join(f"{q} flat={v['flat_1e-6_holds']} elementwise="
                              f"{v['elementwise_float32_holds']}"
                              for q, v in ladder[label].items()
                              if isinstance(v, dict)), flush=True)
        ladders[name] = ladder

    out = {
        "what_this_is": (
            "the smallest injected relative error in the drug signatures that each "
            "candidate reproduction tolerance rejects. A rule that admits the observed "
            "float32 discrepancy has to refuse a wrong value, or it is not a check."),
        "computed_no_registered_statistic": True,
        "rules": {"flat": "max |a-b| < 1e-6",
                  "elementwise_float32": "|a-b| <= 1e-8 + 1e-6 |b| for every drug"},
        "note_on_the_ladders": (
            "a cosine does not change when either argument is rescaled, so D, the mean "
            "pairwise cosine, is invariant to every rescaling and neither scaling ladder "
            "can move it. E is a squared cosine against the mean signature, so scaling "
            "one context of a drug does move it but scaling all of them does not. Only "
            "the third ladder, which moves one gene coordinate and so turns a signature, "
            "reaches all three quantities. A ladder that cannot move the quantity it "
            "claims to test establishes nothing about the tolerance on it."),
        "ladders": ladders,
    }
    directory = Path("/out/03d_h3_reference_discordance")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "reproduction_tolerance_power.json").write_text(
        json.dumps(out, indent=2, allow_nan=False))
    results.commit()
    return json.dumps(ladder, indent=2)


def _as_plain(mapping):
    return {key: dict(value) for key, value in mapping.items()}


@app.function(**COMMON)
def stage_localize_reproduction_error():
    """Where the reproduction difference comes from, before any tolerance is amended.

    Perplexity's round 23 ruling blocks an amendment behind this: a refusal that
    hides a localized implementation discrepancy is not a tolerance problem. The
    question is whether the difference between the deposited values and the
    authoritative rebuild is a property of the aggregation arithmetic or of the
    data.

    Three routes are computed for the same drugs:

      legacy float32   the retired loader as it stands, `signatures[rows].mean(0)`
                       on a float32 array, which numpy accumulates in float32
      legacy float64   the same signatures and the same membership, accumulated in
                       float64
      rebuild          the authoritative per-drug per-cell means from the shards

    If the deposited values match the legacy float32 route and the rebuild matches
    the legacy float64 route, the difference is the float32 accumulator in the
    retired aggregation and the deposited value is the less exact of the two. If
    the membership or the per-pair structure differs instead, that is an
    implementation discrepancy and no tolerance should be amended.

    No registered statistic is computed: the per-pair projections reported here are
    intermediates of a reproduction check, and no correlation, interval or reading
    is produced.
    """
    import importlib.util
    import json
    import site
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    site.addsitedir("/app")
    spec = importlib.util.spec_from_file_location(
        "discordance", "/app/experiments/03d_h3_reference_discordance.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    FOCUS = ["BMS-387032", "JNJ-38877605", "KU-60019", "amonafide",
             "hexamethylenebisacetamide", "latanoprost", "levofloxacin",
             "pefloxacin", "targinine"]

    records = json.loads(Path(
        "/app/results/03b_h3_crispri/h3_crispri_results.json").read_text())
    by_drug = {r["drug"]: r for r in records}

    # the retired route, and the same route with a float64 accumulator
    compounds = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    raw = compounds["signatures"]
    stored_dtype = str(raw.dtype)
    position = {str(sig_id): i for i, sig_id in enumerate(compounds["sig_ids"])}
    siginfo = pd.read_csv(Path(staged) / "GSE92742_Broad_LINCS_sig_info.txt.gz",
                          sep="\t", low_memory=False)
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(position)
                      & siginfo.pert_iname.notna()].copy()
    siginfo["_row"] = siginfo.sig_id.astype(str).map(position)

    legacy32, legacy64, membership = {}, {}, {}
    for (drug, cell), group in siginfo.groupby(["pert_iname", "cell_id"]):
        rows = group._row.values
        legacy32.setdefault(drug, {})[cell] = raw[rows].mean(axis=0)
        legacy64.setdefault(drug, {})[cell] = raw[rows].astype(np.float64).mean(axis=0)
        membership.setdefault(drug, {})[cell] = sorted(
            str(s) for s in group.sig_id.astype(str))

    shrna_legacy = mod.build_shrna_reference(staged)
    rebuild, gene_axis = mod.build_drug_signatures_from_rebuild(rebuilt)
    shrna_rebuild = mod.build_shrna_reference_from_rebuild(rebuilt, gene_axis)
    targets, _, _ = mod.drug_targets(staged)

    def projection_pairs(cells, reference, target, dtype=np.float64):
        """|(s_i - s_j) . u| per ordered context pair, keyed by the two cell names.

        `dtype` carries the whole projection: the difference, the matrix product
        and the mean. A float32 path is what the deposited pipeline would have
        taken on a float32 extraction, and its rounding lands at the magnitude of
        P rather than of a signature value.
        """
        names = sorted(cells)
        matrix = np.vstack([np.asarray(cells[name], dtype=dtype) for name in names])
        direction = reference.directions[target][reference.positions] \
            if len(reference.positions) != matrix.shape[1] else reference.directions[target]
        direction = np.asarray(direction, dtype=dtype)
        upper = np.triu_indices(len(names), k=1)
        differences = matrix[upper[0]] - matrix[upper[1]]
        values = np.abs(differences @ direction)
        return {f"{names[a]}|{names[b]}": float(v)
                for a, b, v in zip(upper[0], upper[1], values)}, float(values.mean(dtype=dtype))

    findings = {}
    for drug in FOCUS:
        target = targets.get(drug)
        if target is None or drug not in rebuild or drug not in legacy32:
            findings[drug] = {"skipped": "absent from a route or carries no target"}
            continue
        entry = {"target": target, "deposited_P": float(by_drug[drug]["proj_shrna"])}
        legacy_cells = set(legacy32[drug])
        rebuild_cells = set(rebuild[drug])
        entry["membership"] = {
            "n_cells_legacy": len(legacy_cells), "n_cells_rebuild": len(rebuild_cells),
            "cell_sets_identical": legacy_cells == rebuild_cells,
            "only_in_legacy": sorted(legacy_cells - rebuild_cells),
            "only_in_rebuild": sorted(rebuild_cells - legacy_cells),
            "signatures_per_cell": {c: len(v) for c, v in membership[drug].items()},
        }
        if legacy_cells != rebuild_cells:
            findings[drug] = entry
            continue

        routes = (("legacy_float32", legacy32[drug], shrna_legacy, np.float64),
                  ("legacy_float64", legacy64[drug], shrna_legacy, np.float64),
                  ("rebuild", rebuild[drug], shrna_rebuild, np.float64),
                  # the projection itself in float32, which is what a float32
                  # extraction carried end to end
                  ("legacy_float32_projection", legacy32[drug], shrna_legacy, np.float32),
                  ("rebuild_float32_projection", rebuild[drug], shrna_rebuild, np.float32))
        for label, cells, reference, dtype in routes:
            pairs, mean = projection_pairs(cells, reference, target, dtype=dtype)
            entry[label] = {"P": mean, "n_pairs": len(pairs), "projection_dtype": str(dtype.__name__)}
            entry.setdefault("_pairs", {})[label] = pairs

        pairs = entry.pop("_pairs")
        entry["differences_from_deposited"] = {
            label: abs(entry[label]["P"] - entry["deposited_P"])
            for label, _, _, _ in routes}
        entry["rebuild_against_legacy_float64"] = abs(
            entry["rebuild"]["P"] - entry["legacy_float64"]["P"])
        entry["rebuild_against_legacy_float32"] = abs(
            entry["rebuild"]["P"] - entry["legacy_float32"]["P"])
        shared = sorted(set(pairs["rebuild"]) & set(pairs["legacy_float32"]))
        per_pair = np.array([abs(pairs["rebuild"][k] - pairs["legacy_float32"][k])
                             for k in shared])
        entry["per_pair_against_legacy_float32"] = {
            "n_pairs": len(shared),
            "max": float(per_pair.max()), "mean": float(per_pair.mean()),
            "worst_pair": shared[int(per_pair.argmax())],
            # one pair carrying the whole difference is a localized discrepancy;
            # a difference spread over every pair is accumulation
            "n_pairs_above_half_the_max": int((per_pair > per_pair.max() / 2).sum()),
        }
        findings[drug] = entry
        differences = entry["differences_from_deposited"]
        print(f"[{_ts()}] {drug}: deposited {entry['deposited_P']:.10f} | "
              f"f64 routes: legacy {differences['legacy_float32']:.3g} "
              f"rebuild {differences['rebuild']:.3g} | f32 projection: legacy "
              f"{differences['legacy_float32_projection']:.3g} rebuild "
              f"{differences['rebuild_float32_projection']:.3g}", flush=True)

    out = {
        "what_this_is": (
            "where the reproduction difference comes from, measured before any tolerance "
            "is amended. The retired loader takes the mean of a float32 array, which "
            "numpy accumulates in float32; the rebuild's shards carry a mean computed "
            "elsewhere. This asks whether the deposited values track the float32 "
            "accumulation and the rebuild tracks the float64 one, or whether the "
            "membership or per-pair structure differs, which would be an implementation "
            "discrepancy rather than arithmetic."),
        "computed_no_registered_statistic": True,
        "lincs_subset_signatures_dtype": stored_dtype,
        "focus_drugs": FOCUS,
        "why_these_drugs": ("the drug that crossed the flat tolerance, and the five "
                            "largest by absolute error, by relative error and by error "
                            "in float32 ULPs"),
        "findings": findings,
    }
    directory = Path("/out/03d_h3_reference_discordance")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "reproduction_error_localized.json").write_text(
        json.dumps(out, indent=2, allow_nan=False))
    results.commit()
    return json.dumps({d: v.get("differences_from_deposited") for d, v in findings.items()},
                      indent=2)


@app.function(**COMMON)
def stage_two_layer_gate_measurement():
    """Does the two-layer gate round 23 proposes actually hold, on all 795 drugs?

    Layer one, legacy reproduction: the retired extraction and the retired loader
    against the deposited records, under the frozen flat absolute 1e-06. If the
    legacy route cannot reproduce the deposited analysis to its own registered
    tolerance, that is a finding and no amendment should paper over it.

    Layer two, production equivalence: the authoritative rebuild against the
    deposited records under the elementwise float32 rule Amendment 2 registers.

    Both layers are measured here and neither is enforced. No registered statistic
    is computed; these are the reproduction errors a gate would compare.
    """
    import importlib.util
    import json
    import site
    from pathlib import Path

    import numpy as np

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    site.addsitedir("/app")
    spec = importlib.util.spec_from_file_location(
        "discordance", "/app/experiments/03d_h3_reference_discordance.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    records = json.loads(Path(
        "/app/results/03b_h3_crispri/h3_crispri_results.json").read_text())
    eligible = {r["drug"] for r in records if "proj_shrna" in r}
    by_drug = {r["drug"]: r for r in records}
    targets, _, _ = mod.drug_targets(staged)

    def errors_for(per_drug, reference, label):
        arm = mod.assemble_arm(label, per_drug, targets, reference, drugs=eligible)
        values = mod.quantities(arm, reference)
        rows = []
        for i, drug in enumerate(arm.drugs):
            record = by_drug[drug]
            for quantity, computed, deposited in (
                    ("P_shrna", values["P"][i], record["proj_shrna"]),
                    ("E_shrna", values["E"][i], record["enrich_shrna"]),
                    ("D", values["D"][i], record["raw_instability"])):
                error = abs(float(computed) - float(deposited))
                allowed = 1e-8 + 1e-6 * abs(float(deposited))
                rows.append({"drug": drug, "quantity": quantity,
                             "computed": float(computed), "deposited": float(deposited),
                             "absolute_error": error, "allowed_elementwise": allowed,
                             "normalized_residual": error / allowed})
        return arm, rows

    # The legacy layer reads the retired extraction on its own declared terms.
    # `build_drug_signatures` refuses it, because Amendment 3's B2 put a
    # frozen-axis check there so a known-bad declaration cannot reach production.
    # That check is exactly right for production and exactly wrong here: the
    # deposited analysis consumed this file with that declaration, and reproducing
    # the deposited analysis means consuming it the same way. The aggregation
    # below is the retired loader's, line for line, without that check.
    import pandas as pd

    compounds = np.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    signatures = compounds["signatures"]
    position = {str(sig_id): i for i, sig_id in enumerate(compounds["sig_ids"])}
    siginfo = pd.read_csv(Path(staged) / "GSE92742_Broad_LINCS_sig_info.txt.gz",
                          sep="\t", low_memory=False)
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(position)
                      & siginfo.pert_iname.notna()].copy()
    siginfo["_row"] = siginfo.sig_id.astype(str).map(position)
    legacy_per_drug = {}
    for (drug, cell), group in siginfo.groupby(["pert_iname", "cell_id"]):
        legacy_per_drug.setdefault(drug, {})[cell] = \
            signatures[group._row.values].mean(axis=0)
    legacy_axis = [str(g) for g in compounds["gene_ids"]]
    legacy_reference = mod.build_shrna_reference(staged)
    rebuild_per_drug, rebuild_axis = mod.build_drug_signatures_from_rebuild(rebuilt)
    rebuild_reference = mod.build_shrna_reference_from_rebuild(rebuilt, rebuild_axis)

    layers = {}
    for name, per_drug, reference in (("legacy_reproduction", legacy_per_drug, legacy_reference),
                                      ("production_equivalence", rebuild_per_drug,
                                       rebuild_reference)):
        arm, rows = errors_for(per_drug, reference, name)
        summary = {"n_drugs": len(arm.drugs)}
        for quantity in ("P_shrna", "E_shrna", "D"):
            group = [r for r in rows if r["quantity"] == quantity]
            errors = np.array([r["absolute_error"] for r in group])
            residuals = np.array([r["normalized_residual"] for r in group])
            over_flat = [r for r in group if r["absolute_error"] >= 1e-6]
            over_elementwise = [r for r in group if r["normalized_residual"] > 1.0]
            summary[quantity] = {
                "max_absolute_error": float(errors.max()),
                "max_normalized_residual": float(residuals.max()),
                "worst_drug_by_residual": group[int(residuals.argmax())]["drug"],
                "worst_drug_by_absolute": group[int(errors.argmax())]["drug"],
                "flat_1e-6_holds": bool(errors.max() < 1e-6),
                "n_over_flat_1e-6": len(over_flat),
                "drugs_over_flat_1e-6": sorted(r["drug"] for r in over_flat),
                "elementwise_float32_holds": len(over_elementwise) == 0,
                "n_over_elementwise": len(over_elementwise),
                "drugs_over_elementwise": sorted(r["drug"] for r in over_elementwise),
            }
        layers[name] = summary
        print(f"[{_ts()}] {name}: " + ", ".join(
            f"{q} flat={summary[q]['flat_1e-6_holds']}({summary[q]['n_over_flat_1e-6']}) "
            f"elementwise={summary[q]['elementwise_float32_holds']} "
            f"maxres={summary[q]['max_normalized_residual']:.3f}"
            for q in ("P_shrna", "E_shrna", "D")), flush=True)

    out = {
        "what_this_is": (
            "both layers of the gate round 23 proposes, measured and not enforced. Layer "
            "one asks whether the retired route reproduces the deposited analysis to the "
            "frozen flat absolute 1e-06. Layer two asks whether the authoritative rebuild "
            "satisfies the elementwise float32 rule Amendment 2 registers. A legacy layer "
            "that fails its own frozen tolerance is a finding, not something to amend."),
        "computed_no_registered_statistic": True,
        "rules": {"flat": "max |a-b| < 1e-6",
                  "elementwise_float32": "|a-b| <= 1e-8 + 1e-6|b| for every drug"},
        "layers": layers,
    }
    directory = Path("/out/03d_h3_reference_discordance")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "two_layer_gate_measurement.json").write_text(
        json.dumps(out, indent=2, allow_nan=False))
    results.commit()
    return json.dumps(layers, indent=2)


def _discordance_module():
    """03d, importable in a container, without putting `sys` in reach."""
    import importlib.util
    import site

    site.addsitedir("/app")
    spec = importlib.util.spec_from_file_location(
        "discordance", "/app/experiments/03d_h3_reference_discordance.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _retired_drug_signatures(staged, numpy, pandas):
    """The retired loader's aggregation, line for line, without the axis check.

    `build_drug_signatures` refuses this extraction, because Amendment 3's B2 put a
    frozen-axis check there so a known-bad declaration cannot reach production. That
    check is right for production and wrong here: the deposited analysis consumed
    this file with that declaration, and reproducing the deposited analysis means
    consuming it the same way. This runs once, outside production, and production
    consumes the sealed audit rather than this path.
    """
    from pathlib import Path

    compounds = numpy.load(Path(staged) / "lincs_subset.npz", allow_pickle=True)
    signatures = compounds["signatures"]
    position = {str(sig_id): i for i, sig_id in enumerate(compounds["sig_ids"])}
    siginfo = pandas.read_csv(Path(staged) / "GSE92742_Broad_LINCS_sig_info.txt.gz",
                              sep="\t", low_memory=False)
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(position)
                      & siginfo.pert_iname.notna()].copy()
    siginfo["_row"] = siginfo.sig_id.astype(str).map(position)
    per_drug = {}
    for (drug, cell), group in siginfo.groupby(["pert_iname", "cell_id"]):
        per_drug.setdefault(drug, {})[cell] = signatures[group._row.values].mean(axis=0)
    return per_drug


@app.function(**COMMON)
def stage_legacy_reproduction_audit():
    """The sealed legacy audit production consumes, with what it was measured on.

    Round 24 found the gate trusting three booleans from whatever file it was
    handed: the record was described as pinned and the hash was computed after
    parsing and only reported. A record that attests to itself establishes nothing,
    so this writes an audit that carries its schema, its rule, the inputs it was
    measured on, and the hash of a complete per-drug table, and `verified_legacy_audit`
    authenticates it against an external pin before parsing it.

    The route-to-route comparison is computed here over all 795 drugs, because the
    claim that the two reconstructable routes agree with each other more closely
    than either agrees with the deposited values was drawn from nine focus drugs.
    """
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    mod = _discordance_module()

    records = json.loads(Path(
        "/app/results/03b_h3_crispri/h3_crispri_results.json").read_text())
    eligible = {r["drug"] for r in records if "proj_shrna" in r}
    by_drug = {r["drug"]: r for r in records}
    targets, _, _ = mod.drug_targets(staged)

    legacy_per_drug = _retired_drug_signatures(staged, np, pd)
    legacy_reference = mod.build_shrna_reference(staged)
    rebuild_per_drug, rebuild_axis = mod.build_drug_signatures_from_rebuild(rebuilt)
    rebuild_reference = mod.build_shrna_reference_from_rebuild(rebuilt, rebuild_axis)

    computed = {}
    for route, per_drug, reference in (
            ("legacy", legacy_per_drug, legacy_reference),
            ("authoritative", rebuild_per_drug, rebuild_reference)):
        arm = mod.assemble_arm(route, per_drug, targets, reference, drugs=eligible)
        values = mod.quantities(arm, reference)
        computed[route] = {"drugs": [str(d) for d in arm.drugs], "values": values}

    drugs = computed["legacy"]["drugs"]
    if drugs != computed["authoritative"]["drugs"]:
        raise RuntimeError(
            "the two routes do not hold the same drugs, so a paired comparison "
            f"between them is not defined: {len(drugs)} against "
            f"{len(computed['authoritative']['drugs'])}")

    quantities = (("P_shrna", "P", "proj_shrna"), ("E_shrna", "E", "enrich_shrna"),
                  ("D", "D", "raw_instability"))
    rows, paired = [], []
    for i, drug in enumerate(drugs):
        record = by_drug[drug]
        for name, key, field in quantities:
            deposited = float(record[field])
            allowed = 1e-8 + 1e-6 * abs(deposited)
            for route in ("legacy", "authoritative"):
                value = float(computed[route]["values"][key][i])
                error = abs(value - deposited)
                rows.append({"drug": drug, "quantity": name, "route": route,
                             "computed": value, "deposited": deposited,
                             "absolute_error": error,
                             "allowed_elementwise": allowed,
                             "normalized_residual": error / allowed,
                             "within_flat_1e-6": error < 1e-6})
            a = float(computed["legacy"]["values"][key][i])
            b = float(computed["authoritative"]["values"][key][i])
            paired.append({"drug": drug, "quantity": name,
                           "legacy": a, "authoritative": b, "deposited": deposited,
                           "route_to_route": abs(a - b),
                           "legacy_to_deposited": abs(a - deposited),
                           "authoritative_to_deposited": abs(b - deposited)})

    directory = Path("/out/03d_h3_reference_discordance")
    directory.mkdir(parents=True, exist_ok=True)

    # the complete table is written and hashed first, because the audit carries its
    # digest and aggregate summaries are not accepted without it
    legacy_rows = [r for r in rows if r["route"] == "legacy"]
    table_path = directory / "legacy_reproduction_per_drug.json"
    table_path.write_text(json.dumps(
        {"what_this_is": "every legacy-route comparison behind the sealed audit's summaries",
         "rule": "max |a-b| < 1e-6", "n_rows": len(legacy_rows), "rows": legacy_rows},
        indent=2, allow_nan=False))
    table_sha = mod.sha256_file(table_path)

    layer = {}
    for name, _, _ in quantities:
        group = [r for r in legacy_rows if r["quantity"] == name]
        errors = np.array([r["absolute_error"] for r in group])
        residuals = np.array([r["normalized_residual"] for r in group])
        over = [r for r in group if r["absolute_error"] >= 1e-6]
        layer[name] = {
            "max_absolute_error": float(errors.max()),
            "max_normalized_residual": float(residuals.max()),
            "worst_drug_by_residual": group[int(residuals.argmax())]["drug"],
            "worst_drug_by_absolute": group[int(errors.argmax())]["drug"],
            "flat_1e-6_holds": bool(errors.max() < 1e-6),
            "n_over_flat_1e-6": len(over),
            "drugs_over_flat_1e-6": sorted(r["drug"] for r in over)}

    provenance = {
        "retired_extraction_sha256": mod.sha256_file(Path(staged) / "lincs_subset.npz"),
        "retired_shrna_sha256": mod.sha256_file(Path(staged) / "lincs_shrna.npz"),
        "retired_shrna_siginfo_sha256": mod.sha256_file(
            Path(staged) / "lincs_shrna_siginfo.csv.gz"),
        "compound_siginfo_sha256": mod.sha256_file(
            Path(staged) / "GSE92742_Broad_LINCS_sig_info.txt.gz"),
        "deposited_records_sha256": mod.sha256_file(
            Path("/app/results/03b_h3_crispri/h3_crispri_results.json")),
        "analysis_code_sha256": mod.sha256_file(
            Path("/app/experiments/03d_h3_reference_discordance.py")),
        "measurement_code_sha256": mod.sha256_file(Path("/root/modal_h3_execute.py")),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
    }

    audit = {
        "schema": "legacy_reproduction_audit/1",
        "what_this_is": (
            "the quarantined legacy route, executed once outside production and sealed. "
            "Production never reads the retired extraction; it verifies this audit "
            "against its external pin and consumes it."),
        "computed_no_registered_statistic": True,
        "flat_rule": "max |a-b| < 1e-6",
        "flat_tolerance": 1e-6,
        "quantities": [name for name, _, _ in quantities],
        "n_drugs": len(drugs),
        "layers": {"legacy_reproduction": layer},
        "per_drug_table": {"file": table_path.name, "sha256": table_sha,
                           "n_rows": len(legacy_rows)},
        "provenance": provenance,
    }
    audit_path = directory / "legacy_reproduction_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, allow_nan=False))

    # the route-to-route comparison, over the whole cohort rather than nine drugs
    route_summary = {}
    for name, _, _ in quantities:
        group = [p for p in paired if p["quantity"] == name]
        between = np.array([p["route_to_route"] for p in group])
        to_legacy = np.array([p["legacy_to_deposited"] for p in group])
        to_auth = np.array([p["authoritative_to_deposited"] for p in group])
        closer = int(sum(p["route_to_route"] < min(p["legacy_to_deposited"],
                                                   p["authoritative_to_deposited"])
                         for p in group))
        route_summary[name] = {
            "n_drugs": len(group),
            "max_route_to_route": float(between.max()),
            "median_route_to_route": float(np.median(between)),
            "max_legacy_to_deposited": float(to_legacy.max()),
            "max_authoritative_to_deposited": float(to_auth.max()),
            "median_legacy_to_deposited": float(np.median(to_legacy)),
            "median_authoritative_to_deposited": float(np.median(to_auth)),
            "n_drugs_where_routes_agree_more_closely_than_either_agrees_with_deposited":
                closer,
            "holds_for_every_drug": closer == len(group)}

    (directory / "route_to_route_all_drugs.json").write_text(json.dumps(
        {"what_this_is": (
            "the paired legacy-against-authoritative comparison over every drug of the "
            "registered cohort. The nine-drug localization diagnostic could not "
            "establish a statement about the cohort."),
         "n_drugs": len(drugs), "summary": route_summary, "rows": paired},
        indent=2, allow_nan=False))
    results.commit()

    print(f"[{_ts()}] audit sha256       {mod.sha256_file(audit_path)}", flush=True)
    print(f"[{_ts()}] per-drug sha256    {table_sha}", flush=True)
    for name, _, _ in quantities:
        print(f"[{_ts()}] {name}: legacy flat holds={layer[name]['flat_1e-6_holds']} "
              f"max={layer[name]['max_absolute_error']:.4g}; routes agree more closely "
              f"for {route_summary[name]['n_drugs_where_routes_agree_more_closely_than_either_agrees_with_deposited']}"
              f"/{route_summary[name]['n_drugs']}", flush=True)
    return json.dumps({"audit_sha256": mod.sha256_file(audit_path),
                       "per_drug_sha256": table_sha,
                       "provenance": provenance,
                       "n_drugs": len(drugs),
                       "route_to_route": route_summary}, indent=2)


@app.function(**COMMON)
def stage_d9_indexing_diagnostic():
    """Does the pre-fix shRNA consensus indexing explain the residual?

    The plan is frozen in `experiments/PREREG_H3_S1S3_DIAGNOSTIC_D9_INDEXING.md`.
    One reconstruction of the defect Deviation 9 names, no variants, no tuning, and
    `D` excluded because it does not depend on the reference direction at all.

    This is quarantined forensic work. The defective indexing is not a production
    route under any outcome, and Amendment 5's gate rule does not change on it.
    """
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    mod = _discordance_module()

    records = json.loads(Path(
        "/app/results/03b_h3_crispri/h3_crispri_results.json").read_text())
    eligible = {r["drug"] for r in records if "proj_shrna" in r}
    by_drug = {r["drug"]: r for r in records}
    targets, _, _ = mod.drug_targets(staged)

    # the defect, reconstructed exactly: each signature id mapped to its row
    # position in the metadata file, and that number used to index the matrix
    shrna = np.load(Path(staged) / "lincs_shrna.npz", allow_pickle=True)
    signatures = shrna["signatures"]
    matrix_position = {str(s): i for i, s in enumerate(shrna["sig_ids"])}
    siginfo = pd.read_csv(Path(staged) / "lincs_shrna_siginfo.csv.gz")
    metadata_position = {str(s): i for i, s in enumerate(siginfo.sig_id.astype(str))}
    n_before = len(siginfo)
    siginfo = siginfo[siginfo.sig_id.astype(str).isin(matrix_position)].copy()
    if len(siginfo) != n_before:
        raise RuntimeError(
            f"the metadata filter drops {n_before - len(siginfo)} rows, so a row "
            "position in the file and a position after filtering are not the same "
            "number and the reconstruction would be ambiguous")

    directions, counts = {}, {}
    for gene, group in siginfo.groupby("pert_iname"):
        ids = sorted(set(group.sig_id.astype(str)))
        if len(ids) < mod.MIN_SIGNATURES:
            continue
        consensus = signatures[[metadata_position[i] for i in ids]].mean(axis=0)
        if np.linalg.norm(consensus) == 0:
            continue
        directions[gene] = mod.unit(consensus)
        counts[gene] = len(ids)
    prefix_reference = mod.Reference("shRNA pre-fix indexing", directions,
                                     np.arange(mod.N_LANDMARK),
                                     {"n_targets": len(directions),
                                      "indexing": "metadata row position, the Deviation 9 defect"})

    legacy_per_drug = _retired_drug_signatures(staged, np, pd)
    corrected_reference = mod.build_shrna_reference(staged)
    rebuild_per_drug, rebuild_axis = mod.build_drug_signatures_from_rebuild(rebuilt)
    rebuild_reference = mod.build_shrna_reference_from_rebuild(rebuilt, rebuild_axis)

    arms = {}
    for label, per_drug, reference in (
            ("prefix_indexing", legacy_per_drug, prefix_reference),
            ("corrected_legacy", legacy_per_drug, corrected_reference),
            ("authoritative", rebuild_per_drug, rebuild_reference)):
        arm = mod.assemble_arm(label, per_drug, targets, reference, drugs=eligible)
        arms[label] = (arm, mod.quantities(arm, reference))
        print(f"[{_ts()}] {label}: {len(arm.drugs)} drugs, "
              f"{len(reference.directions)} directions", flush=True)

    shared = sorted(set.intersection(*[{str(d) for d in arm.drugs}
                                       for arm, _ in arms.values()]))
    index = {label: {str(d): i for i, d in enumerate(arm.drugs)}
             for label, (arm, _) in arms.items()}

    rows = []
    for drug in shared:
        record = by_drug[drug]
        for name, key, field in (("P_shrna", "P", "proj_shrna"),
                                 ("E_shrna", "E", "enrich_shrna")):
            deposited = float(record[field])
            value = {label: float(arms[label][1][key][index[label][drug]])
                     for label in arms}
            rows.append({
                "drug": drug, "quantity": name, "deposited": deposited,
                **{f"{label}": value[label] for label in arms},
                "prefix_to_deposited": abs(value["prefix_indexing"] - deposited),
                "corrected_to_deposited": abs(value["corrected_legacy"] - deposited),
                "authoritative_to_deposited": abs(value["authoritative"] - deposited)})

    summary, verdicts = {}, {}
    for name in ("P_shrna", "E_shrna"):
        group = [r for r in rows if r["quantity"] == name]
        prefix = np.array([r["prefix_to_deposited"] for r in group])
        corrected = np.array([r["corrected_to_deposited"] for r in group])
        summary[name] = {
            "n_drugs": len(group),
            "max_prefix_to_deposited": float(prefix.max()),
            "max_corrected_to_deposited": float(corrected.max()),
            "median_prefix_to_deposited": float(np.median(prefix)),
            "median_corrected_to_deposited": float(np.median(corrected)),
            "max_improvement_factor": float(corrected.max() / prefix.max())
            if prefix.max() > 0 else None,
            "median_improvement_factor": float(np.median(corrected) / np.median(prefix))
            if np.median(prefix) > 0 else None,
            "n_drugs_prefix_closer": int(sum(prefix < corrected)),
            "reaches_1e-8": bool(prefix.max() < 1e-8)}
        # the criterion the frozen plan states, applied as written
        verdicts[name] = {
            "H1": bool(summary[name]["reaches_1e-8"]
                       and summary[name]["median_improvement_factor"] is not None
                       and summary[name]["median_improvement_factor"] >= 10),
            "improved_tenfold_at_the_maximum": bool(
                summary[name]["max_improvement_factor"] is not None
                and summary[name]["max_improvement_factor"] >= 10),
            "H2": bool(summary[name]["max_improvement_factor"] is not None
                       and summary[name]["max_improvement_factor"] < 2)}

    holds = all(verdicts[name]["H1"] for name in ("P_shrna", "E_shrna"))
    rejected = all(verdicts[name]["H2"] for name in ("P_shrna", "E_shrna"))
    reading = ("H1 holds: the pre-fix indexing accounts for the residual on both "
               "reference-dependent quantities across the cohort" if holds else
               "H2 holds: the pre-fix indexing does not account for the residual"
               if rejected else
               "neither H1 nor H2 as the plan states them; the partial outcome is "
               "reported and the residual is not described as explained")

    out = {
        "registered_by": "experiments/PREREG_H3_S1S3_DIAGNOSTIC_D9_INDEXING.md",
        "question": ("does the exact historical indexing defect explain the residual "
                     "across the complete cohort and both reference-dependent quantities?"),
        "what_this_is": (
            "quarantined forensic work. The defective indexing is reconstructed once to "
            "ask how the deposited numbers were produced, and is not a production route."),
        "computed_no_registered_statistic": True,
        "D_excluded_because": (
            "D is one minus the mean pairwise cosine of a drug's own signatures and does "
            "not depend on the reference direction, so it cannot move under any change of "
            "reference"),
        "n_drugs": len(shared),
        "summary": summary,
        "verdicts": verdicts,
        "H1_holds": holds,
        "H2_holds": rejected,
        "reading": reading,
        "rows": rows,
    }
    directory = Path("/out/03d_h3_reference_discordance")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "d9_indexing_diagnostic.json").write_text(
        json.dumps(out, indent=2, allow_nan=False))
    results.commit()
    print(f"[{_ts()}] {reading}", flush=True)
    return json.dumps({"summary": summary, "verdicts": verdicts, "reading": reading},
                      indent=2)


# The library versions `uv.lock` resolves at f822fb1, the commit that produced the
# deposited records by running 03b locally. cmapPy is absent on purpose: it is
# imported only by the GCTX rebuild stages in this file, never by the geometry
# chain the legacy route loads, and pinning it here would constrain numpy for a
# reason unrelated to the hypothesis.
matched_env = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.4.6", "scipy==1.18.0", "pandas==2.3.3", "anndata==0.12.19")
)
MATCHED = dict(COMMON, image=_with_code(matched_env))


@app.function(**MATCHED)
def stage_environment_diagnostic():
    """Does the library-version difference account for the residual?

    The plan is frozen in `experiments/PREREG_H3_S1S3_DIAGNOSTIC_ENVIRONMENT.md`.
    One image, pinned to the versions `uv.lock` resolves at `f822fb1`, which is the
    commit that produced the deposited records by running 03b locally. The legacy
    route and its inputs are the sealed audit's; only the libraries move.

    No second version combination is tried. If the residual survives, that is the
    answer: fitting an environment to an observed difference is what the frozen plan
    forbids.
    """
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd

    _repo_at_its_absolute_path()
    staged = _stage_inputs()
    mod = _discordance_module()

    print(f"[{_ts()}] numpy {np.__version__}, pandas {pd.__version__}", flush=True)
    audit = json.loads(Path(
        "/app/results/03d_h3_reference_discordance/legacy_reproduction_audit.json").read_text())
    # the inputs must be the ones the audit was measured on, or this is a comparison
    # between two different things rather than between two environments
    for field, path in (("retired_extraction_sha256", Path(staged) / "lincs_subset.npz"),
                        ("retired_shrna_sha256", Path(staged) / "lincs_shrna.npz"),
                        ("deposited_records_sha256",
                         Path("/app/results/03b_h3_crispri/h3_crispri_results.json"))):
        digest = mod.sha256_file(path)
        if digest != audit["provenance"][field]:
            raise RuntimeError(
                f"{field}: this run reads {digest} and the audit was measured on "
                f"{audit['provenance'][field]}, so a difference between the two would "
                "not be a difference between environments")

    records = json.loads(Path(
        "/app/results/03b_h3_crispri/h3_crispri_results.json").read_text())
    eligible = {r["drug"] for r in records if "proj_shrna" in r}
    by_drug = {r["drug"]: r for r in records}
    targets, _, _ = mod.drug_targets(staged)

    per_drug = _retired_drug_signatures(staged, np, pd)
    reference = mod.build_shrna_reference(staged)
    arm = mod.assemble_arm("legacy under matched libraries", per_drug, targets,
                           reference, drugs=eligible)
    values = mod.quantities(arm, reference)

    rows = []
    for i, drug in enumerate(arm.drugs):
        record = by_drug[drug]
        for name, key, field in (("P_shrna", "P", "proj_shrna"),
                                 ("E_shrna", "E", "enrich_shrna"),
                                 ("D", "D", "raw_instability")):
            computed, deposited = float(values[key][i]), float(record[field])
            rows.append({"drug": str(drug), "quantity": name, "computed": computed,
                         "deposited": deposited,
                         "absolute_error": abs(computed - deposited)})

    # the criteria the frozen plan states, applied by the code rather than by eye
    summary, verdicts = {}, {}
    for name in ("P_shrna", "E_shrna", "D"):
        group = [r for r in rows if r["quantity"] == name]
        errors = np.array([r["absolute_error"] for r in group])
        was = audit["layers"]["legacy_reproduction"][name]["max_absolute_error"]
        now = float(errors.max())
        summary[name] = {
            "n_drugs": len(group),
            "max_absolute_error_under_matched_libraries": now,
            "max_absolute_error_in_the_sealed_audit": was,
            "improvement_factor": (was / now) if now > 0 else None,
            "median_absolute_error": float(np.median(errors)),
            "reaches_1e-8": bool(now < 1e-8),
            "n_drugs_exact": int((errors == 0.0).sum()),
        }
        verdicts[name] = {
            "reaches_1e-8": summary[name]["reaches_1e-8"],
            "fell_tenfold": bool(now > 0 and was / now >= 10),
            "fell_twofold": bool(now > 0 and was / now >= 2),
        }

    holds = all(verdicts[q]["reaches_1e-8"] for q in verdicts)
    partly = (not holds) and any(verdicts[q]["fell_tenfold"] for q in verdicts)
    rejected = not any(verdicts[q]["fell_twofold"] for q in verdicts)
    reading = ("H1 holds: under the recorded library versions the legacy route "
               "reproduces the deposited records below 1e-08 on every quantity"
               if holds else
               "H1 holds in part: at least one quantity fell tenfold without reaching "
               "1e-08, and the residual is not described as a library-version effect"
               if partly else
               "H2 holds: no quantity's maximum absolute error fell by as much as a "
               "factor of two, so the residual is not a library-version effect"
               if rejected else
               "neither H1 nor H2 as the plan states them; the partial outcome is "
               "reported and the residual is not described as explained")

    out = {
        "registered_by": "experiments/PREREG_H3_S1S3_DIAGNOSTIC_ENVIRONMENT.md",
        "question": ("does running the legacy route under the library versions that "
                     "produced the deposited records remove the residual?"),
        "what_this_is": (
            "one image, pinned to the versions uv.lock resolves at f822fb1. No second "
            "version combination is tried: fitting an environment to an observed "
            "difference is what the frozen plan forbids."),
        "computed_no_registered_statistic": True,
        "environment": {
            "numpy": np.__version__, "pandas": pd.__version__,
            "pinned_from": "uv.lock at f822fb1",
            "interpreter_not_recoverable": (
                "uv.lock records requires-python >= 3.10 and pins no interpreter, and "
                "nothing records which Python ran on 2026-09-15. The interpreter is held "
                "at the audit's, so only the libraries move."),
            "platform_differs": (
                "the deposited records were produced on an Intel Mac and this runs in a "
                "Debian container, so the BLAS a wheel carries is not reconstructed"),
        },
        "inputs_match_the_sealed_audit": True,
        "audit_compared_against":
            "results/03d_h3_reference_discordance/legacy_reproduction_audit.json",
        "n_drugs": len(arm.drugs),
        "summary": summary,
        "verdicts": verdicts,
        "H1_holds": holds,
        "H1_holds_in_part": partly,
        "H2_holds": rejected,
        "reading": reading,
        "rows": rows,
    }
    directory = Path("/out/03d_h3_reference_discordance")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "environment_diagnostic.json").write_text(
        json.dumps(out, indent=2, allow_nan=False))
    results.commit()
    for name in ("P_shrna", "E_shrna", "D"):
        s = summary[name]
        print(f"[{_ts()}] {name}: audit {s['max_absolute_error_in_the_sealed_audit']:.4g} "
              f"-> matched {s['max_absolute_error_under_matched_libraries']:.4g} "
              f"(factor {s['improvement_factor']:.3g}), exact for {s['n_drugs_exact']}"
              f"/{s['n_drugs']}", flush=True)
    print(f"[{_ts()}] {reading}", flush=True)
    return json.dumps({"summary": summary, "verdicts": verdicts, "reading": reading},
                      indent=2)
