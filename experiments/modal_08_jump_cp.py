"""Modal wrapper for Experiment 8 on JUMP-CP. Logic stays in the 08* scripts.

    modal run --detach experiments/modal_08_jump_cp.py

Runs the three stages in order on one CPU worker, writing everything to a Modal
volume and committing after each stage so a killed app resumes rather than
restarts. Nothing here computes; it calls the same functions that run locally.

Afterwards:
    modal volume get di-exp8 results /path/to/results
"""
import modal

app = modal.App("di-exp8-jump-cp")
vol = modal.Volume.from_name("di-exp8", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "numpy==2.1.3",
        "pandas==2.2.3",
        "pyarrow==18.1.0",
        "scipy==1.14.1",
    )
    .add_local_dir("experiments", remote_path="/app/experiments")
)

COMMON = dict(
    image=image,
    volumes={"/vol": vol},
    timeout=86400,          # 24h
    memory=65536,           # the consensus groupby is the memory peak
    cpu=8.0,
)


@app.function(**COMMON)
def stage_eligible():
    import importlib, sys
    from pathlib import Path
    sys.path.insert(0, "/app/experiments")
    m = importlib.import_module("08_jump_cp_feature_ablation")
    m.REPO = Path("/vol"); m.CACHE = Path("/vol/cache"); m.OUT_DIR = Path("/vol/results")
    m.eligible()
    vol.commit()
    return "eligible done"


@app.function(**COMMON)
def stage_consensus():
    import importlib, sys
    from pathlib import Path
    sys.path.insert(0, "/app/experiments")
    m = importlib.import_module("08_jump_cp_feature_ablation")
    m.REPO = Path("/vol"); m.CACHE = Path("/vol/cache"); m.OUT_DIR = Path("/vol/results")
    m.consensus()
    vol.commit()
    blocks = list(Path("/vol/cache/consensus_blocks").glob("block_*.parquet"))
    return f"consensus done: {len(blocks)} blocks"


@app.function(**COMMON)
def stage_analyze():
    import importlib, sys
    from pathlib import Path
    sys.path.insert(0, "/app/experiments")
    m = importlib.import_module("08_jump_cp_feature_ablation")
    m.REPO = Path("/vol"); m.CACHE = Path("/vol/cache"); m.OUT_DIR = Path("/vol/results")
    m.analyze()
    vol.commit()
    return "ablation done"


@app.function(**COMMON)
def stage_matched():
    import importlib, sys
    from pathlib import Path
    sys.path.insert(0, "/app/experiments")
    m = importlib.import_module("08c_matched_well_contrast")
    m.REPO = Path("/vol"); m.OUT_DIR = Path("/vol/results")
    m.CACHE = Path("/vol/cache")                 # where --eligible wrote well_metadata
    m.CKPT = Path("/vol/results/08c_checkpoints")
    m.main()
    vol.commit()
    return "matched contrast done"


@app.local_entrypoint()
def main(stage: str = "all"):
    """stage: all | eligible | consensus | analyze | matched"""
    if stage in ("all", "eligible"):
        print(stage_eligible.remote())
    if stage in ("all", "consensus"):
        print(stage_consensus.remote())
    if stage in ("all", "analyze"):
        print(stage_analyze.remote())
    if stage in ("all", "matched"):
        print(stage_matched.remote())
    print("\nall stages complete. retrieve with:")
    print("    modal volume get di-exp8 results ./results_from_modal")
