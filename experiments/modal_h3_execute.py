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

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.1.3", "scipy==1.14.1", "pandas==2.2.3", "anndata==0.11.4",
                 "h5py==3.12.1")
    .add_local_dir(f"{REPO}/experiments", remote_path="/app/experiments")
    .add_local_dir(f"{REPO}/geometry", remote_path="/app/geometry")
    .add_local_dir(f"{REPO}/results/03_phenotype_projection",
                   remote_path="/app/results/03_phenotype_projection")
    .add_local_dir(f"{REPO}/results/03b_h3_crispri", remote_path="/app/results/03b_h3_crispri")
    .add_local_dir(f"{REPO}/results/03d_h3_reference_discordance",
                   remote_path="/app/results/03d_h3_reference_discordance")
)
# The Replogle and PRISM files are half a gigabyte. They live on a volume rather
# than in the image: mounting them makes every run upload them again, and a client
# killed mid-upload leaves an app that never dispatches its function.

# 32 GB, not 256: a run that needs more has a bug rather than a big input, and
# retries=0 so a bug fails once instead of nine times overnight
COMMON = dict(image=image, timeout=86400, memory=32768, cpu=8.0, retries=0,
              volumes={"/extraction": extraction, "/rebuild": rebuild, "/out": results,
                       "/inputs": inputs})

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


@app.function(**COMMON)
def stage_gate():
    """The reconstruction gate. A failure here voids the analyses registered against it."""
    from pathlib import Path

    _repo_at_its_absolute_path()
    staged, rebuilt = _stage_inputs(), _stage_rebuild()
    _run(["/app/experiments/03e_reconstruction_gate.py",
          "--rebuilt", str(rebuilt),
          "--extraction", str(staged),
          "--cohort", "/app/results/03_phenotype_projection/phenotype_projection_results.json",
          "--output", "/out/03c_h3_sensitivity"])
    results.commit()
    return (Path("/out/03c_h3_sensitivity/reconstruction_gate.json")).read_text()[:2000]


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


@app.local_entrypoint()
def main(stage: str):
    """stage: gate | s1s3 | driver"""
    if stage == "gate":
        print(stage_gate.remote())
    elif stage == "s1s3":
        print(stage_s1s3.remote())
    elif stage == "driver":
        print(stage_driver.remote())
    else:
        raise SystemExit(f"unknown stage {stage}")
    print("\nretrieve with:")
    print("    modal volume get di-h3-results 03c_h3_sensitivity ./results/")
    print("    modal volume get di-h3-results 03d_h3_reference_discordance ./results/")
