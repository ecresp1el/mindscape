#!/usr/bin/env bash
set -euo pipefail
# Ensures proper execution of script without errors

# ─────────────────────────────────────────────────────────────────────────────
# scripts/cellranger/scripts/create_project_cellranger_div90.sh
#
# - Sources parameterized config_div90.sh
# - Validates required vars (no hardcoded defaults)
# - Copies and patches multi_config.csv (create-bam, reference, probe-set)
# - Unlocks Snakemake dir if needed
# - Launches Snakemake for the specified OUTPUT_ID target
# ─────────────────────────────────────────────────────────────────────────────

# Resolve this script's directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # Dynamically creates the name of the script directory 

# Load config (env vars can override values inside)
# Inputs expected (via config/env):
#   - TEST_DIR: working directory for the run (created if missing)
#   - TURBO_CONFIG_SOURCE: source CSV (10x multi config)
#   - PROBE_PATH: Flex probe-set CSV path
#   - REF_GENOME or TURBO_REF_BASE + REF_SUBPATH: reference folder
#   - OUTPUT_ID: Cell Ranger --id / Snakemake target (default in config)
#   - SNAKEFILE: relative or absolute path to the Snakefile to use
CONFIG_FILE="$SCRIPT_DIR/config_div90.sh" # Creates config file path
# Ensures config file is present and loads it in
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "❌ Missing config file: $CONFIG_FILE" >&2
  exit 1
fi
source "$CONFIG_FILE"

# Determine cellranger folder root (one up from this scripts dir)
CELLRANGER_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Validate required inputs (ensure all inputs exist)
_fail_if_empty() { local n="$1" v="$2"; [[ -n "$v" ]] || { echo "❌ Required var $n is empty" >&2; exit 1; }; }
_fail_if_empty TEST_DIR "${TEST_DIR:-}"
_fail_if_empty TURBO_CONFIG_SOURCE "${TURBO_CONFIG_SOURCE:-}"
_fail_if_empty PROBE_PATH "${PROBE_PATH:-}"

# Resolve reference path
# - Allow REF_GENOME directly, or derive from TURBO_REF_BASE/REF_SUBPATH for clarity, depending on whether the reference genome is initially provided.
if [[ -z "${REF_GENOME:-}" ]]; then
  if [[ -n "${TURBO_REF_BASE:-}" && -n "${REF_SUBPATH:-}" ]]; then
    REF_GENOME="$TURBO_REF_BASE/$REF_SUBPATH"
  else
    echo "❌ Provide REF_GENOME or TURBO_REF_BASE+REF_SUBPATH" >&2
    exit 1
  fi
fi

# In normal mode, require the reference dir to exist. In DRY_RUN, allow missing.
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  if [[ ! -d "$REF_GENOME" ]]; then
    echo "❌ Reference folder not found at $REF_GENOME" >&2
    exit 1
  fi
else
  if [[ ! -d "$REF_GENOME" ]]; then
    echo "ℹ️ [DRY_RUN] Skipping reference existence check: $REF_GENOME" >&2
  fi
fi

# Snakefile absolute path (prefer inside cellranger folder)
# - SNAKEFILE can be relative (using /*) to scripts/cellranger/ or absolute.
SNAKEFILE_REL_OR_ABS="${SNAKEFILE:-cellranger.smk}"
if [[ "$SNAKEFILE_REL_OR_ABS" = /* ]]; then
  SNAKEFILE_ABS="$SNAKEFILE_REL_OR_ABS"
else
  SNAKEFILE_ABS="$CELLRANGER_ROOT/$SNAKEFILE_REL_OR_ABS"
fi
if [[ ! -f "$SNAKEFILE_ABS" ]]; then
  echo "❌ Snakefile not found at $SNAKEFILE_ABS" >&2
  exit 1
fi

# Prepare working directory and stage CSV into it (can run even if directory exists)
mkdir -p "$TEST_DIR"
CONFIG_DEST="$TEST_DIR/multi_config.csv"

# Copy original config
echo "📄 Copying multi-config → $CONFIG_DEST"
cp -f "$TURBO_CONFIG_SOURCE" "$CONFIG_DEST"

# Inject 'create-bam,true' after [gene-expression] if not already present
# - 10x Flex requires BAMs; this ensures the setting is present in the CSV.
if ! grep -q '^create-bam,' "$CONFIG_DEST"; then
  echo "🛠  Injecting create-bam after [gene-expression]…"
  awk '
    BEGIN { ins=0 }
    {
      print
      if (!ins && $0=="[gene-expression]") {
        getline; print "create-bam,true"; print $0
        ins=1
      }
    }
  ' "$CONFIG_DEST" > "${CONFIG_DEST}.tmp" && mv "${CONFIG_DEST}.tmp" "$CONFIG_DEST"
fi
# Inserted into the multi_config that sets up cellranger

# Cross-platform in-place sed helper (macOS/Linux)
# - Uses GNU sed -i if available; falls back to BSD sed behavior on macOS.
_sed_inplace() {
  local expr="$1" file="$2"
  if sed --version >/dev/null 2>&1; then
    sed -i "$expr" "$file"
  else
    # macOS/BSD sed
    sed -i '' "$expr" "$file" 2>/dev/null || { sed -i.bak "$expr" "$file" && rm -f "${file}.bak"; }
  fi
}

echo "🧬 Patching reference → $REF_GENOME"
_sed_inplace "s|^reference,.*|reference,$REF_GENOME|" "$CONFIG_DEST"

echo "🧬 Patching probe-set → $PROBE_PATH"
_sed_inplace "s|^probe-set,.*|probe-set,$PROBE_PATH|" "$CONFIG_DEST"

# Enforce curated accessible CSV source if requested
if [[ "${ENFORCE_ACCESSIBLE:-1}" == "1" ]]; then
  if [[ "$TURBO_CONFIG_SOURCE" != *"/Accessible_multi-config_csvs/"* ]]; then
    echo "❌ TURBO_CONFIG_SOURCE must come from an Accessible_multi-config_csvs path for DIV90." >&2
    echo "   Current: $TURBO_CONFIG_SOURCE" >&2
    echo "   Tip: set TURBO_CONFIG_SOURCE to the curated accessible CSV under /nfs/turbo/umms-parent/Accessible_multi-config_csvs" >&2
    exit 1
  fi
fi

# Optionally normalize FASTQ paths in [libraries], if not a dry run and FASTQS_DIR does not exist
if [[ -n "${FASTQS_DIR:-}" ]]; then
  if [[ "${DRY_RUN:-0}" != "1" && ! -d "$FASTQS_DIR" ]]; then
    echo "❌ FASTQS_DIR not found: $FASTQS_DIR" >&2; exit 1
  fi
  echo "📦 Patching [libraries] fastqs → $FASTQS_DIR"
  awk -v dir="$FASTQS_DIR" 'BEGIN{FS=OFS=","} /^\[/ {sec=$0; print; next} sec=="[libraries]" { if ($1=="fastq_id") {print; next} $2=dir; print; next } {print}' "$CONFIG_DEST" > "${CONFIG_DEST}.tmp" && mv "${CONFIG_DEST}.tmp" "$CONFIG_DEST"
elif [[ -n "${FASTQS_REPLACE_FROM:-}" && -n "${FASTQS_REPLACE_TO:-}" ]]; then
  echo "📦 Rewriting [libraries] fastqs prefix: $FASTQS_REPLACE_FROM → $FASTQS_REPLACE_TO"
  awk -v from="$FASTQS_REPLACE_FROM" -v to="$FASTQS_REPLACE_TO" 'BEGIN{FS=OFS=","} /^\[/ {sec=$0; print; next} sec=="[libraries]" { if ($1=="fastq_id") {print; next} sub(from,to,$2); print; next } {print}' "$CONFIG_DEST" > "${CONFIG_DEST}.tmp" && mv "${CONFIG_DEST}.tmp" "$CONFIG_DEST"
fi

# Validate FASTQ accessibility if enforcement is enabled
if [[ "${ENFORCE_ACCESSIBLE:-1}" == "1" ]]; then
  allowed_prefix="${ALLOWED_FASTQS_PREFIX:-/nfs/turbo/umms-parent}"
  # Only evaluate non-empty data rows inside [libraries]
  bad_count=$(awk -v pfx="$allowed_prefix" '
    BEGIN{FS=","; sec=""}
    /^\[/ {sec=$0; next}
    sec=="[libraries]" && NF>0 {
      if ($1=="fastq_id") next;
      # Skip pure blank lines
      if ($0 ~ /^\s*$/) next;
      # Require column 2 (fastqs) to start with allowed prefix
      if (index($2, pfx) != 1) bad++
    }
    END{print bad+0}
  ' "$CONFIG_DEST")
  if [[ "$bad_count" -gt 0 ]]; then
    echo "❌ Detected $bad_count FASTQ path(s) outside allowed prefix: $allowed_prefix" >&2
    echo "   Use FASTQS_DIR or FASTQS_REPLACE_FROM/FASTQS_REPLACE_TO to normalize the [libraries] fastqs paths." >&2
    echo "   CSV: $CONFIG_DEST" >&2
    exit 1
  fi
fi

# Load modules if available (optional)
# - On module-based HPC, try to load cellranger + snakemake for convenience.
if command -v module >/dev/null 2>&1; then
  echo "⏳ Loading modules…"
  set +u
  module purge || true
  module load Bioinformatics cellranger || true
  module load snakemake || true
  set -u
fi

# Unlock if needed
# - Clears any stale Snakemake locks in TEST_DIR.
echo "🔓 Ensuring Snakemake working directory is unlocked…"
snakemake --unlock \
  --snakefile "$SNAKEFILE_ABS" \
  --directory "$TEST_DIR" || true

# Run Snakemake target
start_time=$(date +%s)

# Respect DRY_RUN=1 to pass -n to Snakemake (plan only)
SMK_DRY=""
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  SMK_DRY="-n"
  echo "🧪 DRY_RUN enabled: Snakemake will not execute Cell Ranger."
fi

echo "🚀 Running Snakemake (Cell Ranger multi)… target=$OUTPUT_ID cores=$CORES ${SMK_DRY:+[dry-run]}"
snakemake --rerun-incomplete -j "$CORES" $SMK_DRY \
  --snakefile "$SNAKEFILE_ABS" \
  --directory "$TEST_DIR" \
  "$OUTPUT_ID"

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "✅ Workflow complete for '$OUTPUT_ID'."
echo "⏱ Total runtime: $(( elapsed / 60 ))m $(( elapsed % 60 ))s"
