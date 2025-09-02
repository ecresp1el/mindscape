DIV90 Cell Ranger – Run Guide
=============================

How the pieces fit
- Driver: `scripts/cellranger/div90/drive_div90.sh` orchestrates runs (dry|local|slurm).
- Config: `scripts/cellranger/scripts/config_div90.sh` provides defaults (paths, IDs, cores).
- Wrapper: `scripts/cellranger/scripts/create_project_cellranger_div90.sh` prepares CSV and launches Snakemake.
- Snakefile: `scripts/cellranger/cellranger.smk` runs `cellranger multi` for the chosen `--id`.
- Slurm: `scripts/cellranger/slurm/submit_cellranger_div90.sh` submits, and `job_cellranger_div90.sh` executes the wrapper.

Flow (end‑to‑end)
1) You run `drive_div90.sh` with a mode.
2) It sources config and exports env used downstream (including `WRAPPER_PATH`). See:
   - `scripts/cellranger/div90/drive_div90.sh:63`
   - `scripts/cellranger/div90/drive_div90.sh:65`
   - `scripts/cellranger/div90/drive_div90.sh:67`
3) Local/dry: it calls the wrapper directly. Slurm: it calls the submit script, which writes logs under `$TEST_DIR/logs`. See:
   - `scripts/cellranger/slurm/submit_cellranger_div90.sh:46`
4) The job script resolves the wrapper reliably from Slurm spool using `WRAPPER_PATH`/`SLURM_SUBMIT_DIR`. See:
   - `scripts/cellranger/slurm/job_cellranger_div90.sh:10`
5) The wrapper sources config, copies your source CSV to `$TEST_DIR/multi_config.csv`, injects `create-bam,true`, and patches `reference` and `probe-set`. See:
   - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:73`
   - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:77`
   - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:103`
   - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:106`
6) Snakemake runs the single rule in `cellranger.smk`, which executes `cellranger multi` to produce `$TEST_DIR/<OUTPUT_ID>`.

Defaults (safe to start with)
- Work dir: `TEST_DIR` → `/nfs/turbo/umms-parent/$USER/mindscape_div90/<timestamp>`
  - `scripts/cellranger/scripts/config_div90.sh:20`
- Source CSV: `TURBO_CONFIG_SOURCE` → prior DIV90 path (contains spaces; already quoted)
  - `scripts/cellranger/scripts/config_div90.sh:25`
- Probe set: `PROBE_PATH` → 10X Human Refs 2020‑A
  - `scripts/cellranger/scripts/config_div90.sh:30`
- Reference: `REF_GENOME` → 10X Human Refs 2020‑A
  - `scripts/cellranger/scripts/config_div90.sh:35`
- Output ID: `OUTPUT_ID=div90-reanalysis`
  - `scripts/cellranger/scripts/config_div90.sh:46`

Requirements
- Tools on PATH or via modules (wrapper tries to load if `module` exists):
  - cellranger, snakemake
- Read access to NFS turbo paths in config.

Run examples
- Dry-run (no execution):
  - `scripts/cellranger/div90/drive_div90.sh dry`
- Local run:
  - `scripts/cellranger/div90/drive_div90.sh local`
- Slurm run (example resources):
  - `TIME=48:00:00 CPUS=32 MEM=128G ACCOUNT=parent0 PARTITION=standard JOB_NAME=div90_full scripts/cellranger/div90/drive_div90.sh slurm`

Monitor progress
- Queue state: `squeue -j <JOBID>`
- Logs (stdout/stderr): `$TEST_DIR/logs/<JOB_NAME>_<JOBID>.out|.err`
  - The driver prints a ready-to-tail command after submission. See:
  - `scripts/cellranger/div90/drive_div90.sh:95`

What gets created
- `$TEST_DIR/multi_config.csv` (copied and patched)
- `$TEST_DIR/<OUTPUT_ID>` (Cell Ranger output directory)
- `$TEST_DIR/logs/*` (Slurm logs when running via Slurm)

Advanced / env overrides
- Override any default before running, e.g. switch to alternative probe/ref:
  - `export PROBE_PATH=/nfs/turbo/.../Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv`
  - `export REF_GENOME=/nfs/turbo/.../refdata-gex-GRCh38-2020-A`
- Dry-run via Slurm (planning only):
  - `DRY_RUN=1 TIME=00:10:00 CPUS=1 MEM=8G JOB_NAME=test_div90_dry scripts/cellranger/div90/drive_div90.sh slurm`
  - Wrapper honors `DRY_RUN` to pass `-n` to Snakemake. See:
    - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:128`

Troubleshooting
- PENDING (Priority): lower `CPUS/MEM`, use a different `PARTITION`/`QOS`, or wait.
- snakemake not found: ensure your module environment provides it or add to PATH.
- Reference folder not found: update `REF_GENOME` or ensure the NFS path is mounted.
- CSV path contains spaces: keep quotes when overriding `TURBO_CONFIG_SOURCE` (defaults already quoted).
- Wrapper not found from Slurm spool: the driver exports `WRAPPER_PATH`; job resolves it. See:
  - `scripts/cellranger/div90/drive_div90.sh:67`
  - `scripts/cellranger/slurm/job_cellranger_div90.sh:10`

Reproduce our run
- Submitted:
  - `TIME=48:00:00 CPUS=32 MEM=128G ACCOUNT=parent0 PARTITION=standard JOB_NAME=div90_full scripts/cellranger/div90/drive_div90.sh slurm`
- Job ID printed on submit; monitor with:
  - `squeue -j <JOBID>`
  - `tail -f $TEST_DIR/logs/div90_full_<JOBID>.out`

