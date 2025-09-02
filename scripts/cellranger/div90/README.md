DIV90 Cell Ranger – How To Run
================================

Overview
- Entry script: `scripts/cellranger/div90/drive_div90.sh` (dry, local, slurm)
- Config: `scripts/cellranger/scripts/config_div90.sh`
- Wrapper: `scripts/cellranger/scripts/create_project_cellranger_div90.sh`
- Snakefile: `scripts/cellranger/cellranger.smk`
- Slurm: `scripts/cellranger/slurm/{submit_cellranger_div90.sh,job_cellranger_div90.sh}`

Quick Start
- Dry-run (plan only):
  - `scripts/cellranger/div90/drive_div90.sh dry`
- Local run:
  - `scripts/cellranger/div90/drive_div90.sh local`
- Slurm submission:
  - `TIME=24:00:00 CPUS=32 MEM=128G ACCOUNT=parent0 scripts/cellranger/div90/drive_div90.sh slurm`

What Gets Used
- `config_div90.sh`: sets defaults for `TEST_DIR`, `TURBO_CONFIG_SOURCE`, `PROBE_PATH`, `REF_GENOME`.
  - TEST_DIR defaults under `/nfs/turbo/umms-parent/$USER/mindscape_div90/<timestamp>`
  - TURBO_CONFIG_SOURCE defaults to the prior DIV90 config.csv path.
  - Probe/reference default to 10X Human Refs 2020-A paths on turbo.
- `create_project_cellranger_div90.sh`: copies and patches `multi_config.csv`, unlocks dir, runs Snakemake.
  - Honors `DRY_RUN=1` to pass `-n` and skip reference existence check.
- `cellranger.smk`: one rule calling `cellranger multi --id {id} --csv multi_config.csv`.
- Slurm wrappers: submit and job scripts that call the same wrapper.

Environment Overrides (optional)
- `TEST_DIR`, `TURBO_CONFIG_SOURCE`, `PROBE_PATH`, `REF_GENOME`, `OUTPUT_ID`, `CORES`.
- For Slurm: `ACCOUNT`, `PARTITION`, `QOS`, `TIME`, `CPUS`, `MEM`, `MAIL_USER`, `JOB_NAME`, `LOG_DIR`.

Sanity Checks
- Inspect patched CSV after dry-run:
  - `grep -E '^(reference|probe-set|create-bam),' "$TEST_DIR/multi_config.csv"`
- Confirm planned output dir:
  - `echo "$TEST_DIR/$OUTPUT_ID"`

Notes
- No files were moved; these wrappers only organize the entry points to avoid conflicts with other pipelines.
- To switch probe/reference defaults, edit `config_div90.sh` or set env vars before running.

