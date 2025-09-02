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
   - `scripts/cellranger/div90/drive_div90.sh:77`
   - `scripts/cellranger/div90/drive_div90.sh:79`
   - `scripts/cellranger/div90/drive_div90.sh:81`
3) Local/dry: it calls the wrapper directly. Slurm: it calls the submit script, which writes logs under `$TEST_DIR/logs`. See:
   - `scripts/cellranger/slurm/submit_cellranger_div90.sh:46`
4) The job script resolves the wrapper reliably from Slurm spool using `WRAPPER_PATH`/`SLURM_SUBMIT_DIR`. See:
   - `scripts/cellranger/slurm/job_cellranger_div90.sh:13`
5) The wrapper sources config, copies your source CSV to `$TEST_DIR/multi_config.csv`, injects `create-bam,true`, and patches `reference` and `probe-set`. See:
   - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:83`
   - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:89`
   - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:114`
   - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:117`
6) Snakemake runs the single rule in `cellranger.smk`, which executes `cellranger multi` to produce `$TEST_DIR/<OUTPUT_ID>`.

Defaults (safe to start with)
- Work dir: `TEST_DIR` → `/nfs/turbo/umms-parent/$USER/mindscape_div90/<timestamp>`
  - `scripts/cellranger/scripts/config_div90.sh:29`
- Source CSV: `TURBO_CONFIG_SOURCE` → prior DIV90 path (contains spaces; already quoted)
  - `scripts/cellranger/scripts/config_div90.sh:34`
- Probe set: `PROBE_PATH` → 10X Human Refs 2020‑A
  - `scripts/cellranger/scripts/config_div90.sh:39`
- Reference: `REF_GENOME` → 10X Human Refs 2020‑A
  - `scripts/cellranger/scripts/config_div90.sh:44`
- Output ID: `OUTPUT_ID=div90-reanalysis`
  - `scripts/cellranger/scripts/config_div90.sh:55`

FASTQ paths (important)
- Some historical config files point FASTQs to restricted prefixes (e.g., `/nfs/turbo/agc-data/...`). If those paths
  are not accessible from your compute nodes, normalize the `fastqs` column in `[libraries]` before running.
- Options (set one):
  - `FASTQS_DIR`: force the `fastqs` column for all rows to a single directory.
  - `FASTQS_REPLACE_FROM` + `FASTQS_REPLACE_TO`: replace a path prefix inside the `fastqs` column.
- Where it’s applied: `scripts/cellranger/scripts/create_project_cellranger_div90.sh` (after reference/probe patching).
- Example:
  - `export FASTQS_DIR=/nfs/turbo/umms-parent/Accessible_fastqs/10496-MW/fastqs_10496-MW`
  - Or
  - `export FASTQS_REPLACE_FROM=/nfs/turbo/agc-data/processing`
  - `export FASTQS_REPLACE_TO=/nfs/turbo/umms-parent/Accessible_fastqs`
  - `scripts/cellranger/div90/drive_div90.sh slurm`

Enforcement
- By default `ENFORCE_ACCESSIBLE=1` ensures:
  - `TURBO_CONFIG_SOURCE` comes from `Accessible_multi-config_csvs`.
  - All `[libraries].fastqs` paths start with `ALLOWED_FASTQS_PREFIX` (default `/nfs/turbo/umms-parent`).
- Override (not recommended): set `ENFORCE_ACCESSIBLE=0` to bypass checks, or change `ALLOWED_FASTQS_PREFIX`.

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
  - Email: the driver auto-sets `MAIL_USER` to `$USER@umich.edu` (override domain via `DEFAULT_EMAIL_DOMAIN`).
    - Override domain: `DEFAULT_EMAIL_DOMAIN=your.edu ... drive_div90.sh slurm`
    - Override email explicitly: `MAIL_USER=you@your.edu ... drive_div90.sh slurm`

Monitor progress
- Queue state: `squeue -j <JOBID>`
- Logs (stdout/stderr): `$TEST_DIR/logs/<JOB_NAME>_<JOBID>.out|.err`
  - The driver prints a ready-to-tail command after submission. See:
  - `scripts/cellranger/div90/drive_div90.sh:111`

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
    - `scripts/cellranger/scripts/create_project_cellranger_div90.sh:141`

Troubleshooting
- PENDING (Priority): lower `CPUS/MEM`, use a different `PARTITION`/`QOS`, or wait.
- snakemake not found: ensure your module environment provides it or add to PATH.
- Reference folder not found: update `REF_GENOME` or ensure the NFS path is mounted.
- CSV path contains spaces: keep quotes when overriding `TURBO_CONFIG_SOURCE` (defaults already quoted).
- Wrapper not found from Slurm spool: the driver exports `WRAPPER_PATH`; job resolves it. See:
  - `scripts/cellranger/div90/drive_div90.sh:81`
  - `scripts/cellranger/slurm/job_cellranger_div90.sh:13`

Reproduce Example (elcrespo)
- Submit (same settings we used):
  - `TIME=48:00:00 CPUS=32 MEM=128G ACCOUNT=parent0 PARTITION=standard JOB_NAME=div90_accessible_full2 scripts/cellranger/div90/drive_div90.sh slurm`
- The driver prints the Job ID and ready-to-tail path. Example from our run:
  - Job: `31730391`
  - Workdir: `/nfs/turbo/umms-parent/elcrespo/mindscape_div90/20250902_144025`
- Monitor:
  - `squeue -j 31730391`
  - `tail -f /nfs/turbo/umms-parent/elcrespo/mindscape_div90/20250902_144025/logs/div90_accessible_full2_31730391.out`
  - `tail -f /nfs/turbo/umms-parent/elcrespo/mindscape_div90/20250902_144025/logs/div90_accessible_full2_31730391.err`
- View UI (tunnel from your laptop):
  - Find the node:port from the log line: `grep -m1 'Serving UI at' $TEST_DIR/logs/div90_accessible_full2_31730391.out`
  - Tunnel: `ssh -N -L 127.0.0.1:8000:<node>:<port> elcrespo@greatlakes.arc-ts.umich.edu`
  - Open: `http://localhost:8000/?auth=<token>` (token is in the same log line)
- Completion check:
  - `sacct -j 31730391 -o JobID,State,ExitCode,Elapsed`
  - `test -d /nfs/turbo/umms-parent/elcrespo/mindscape_div90/20250902_144025/div90-reanalysis/outs && echo outs ready`
