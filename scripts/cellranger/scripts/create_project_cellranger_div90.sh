#!/usr/bin/env bash
set -euo pipefail

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
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load config (env vars can override values inside)
CONFIG_FILE="$SCRIPT_DIR/config_div90.sh"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "❌ Missing config file: $CONFIG_FILE" >&2
  exit 1
fi
source "$CONFIG_FILE"

# Determine cellranger folder root (one up from this scripts dir)
CELLRANGER_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Validate required inputs
_fail_if_empty() { local n="$1" v="$2"; [[ -n "$v" ]] || { echo "❌ Required var $n is empty" >&2; exit 1; }; }

_fail_if_empty TEST_DIR "${TEST_DIR:-}"
_fail_if_empty TURBO_CONFIG_SOURCE "${TURBO_CONFIG_SOURCE:-}"
_fail_if_empty PROBE_PATH "${PROBE_PATH:-}"

# Resolve reference path
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

# Prepare working directory
mkdir -p "$TEST_DIR"
CONFIG_DEST="$TEST_DIR/multi_config.csv"

# Copy original config
echo "📄 Copying multi-config → $CONFIG_DEST"
cp -f "$TURBO_CONFIG_SOURCE" "$CONFIG_DEST"

# Inject 'create-bam,true' after [gene-expression] if not already present
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

# Cross-platform in-place sed helper (macOS/Linux)
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

# Load modules if available (optional)
if command -v module >/dev/null 2>&1; then
  echo "⏳ Loading modules…"
  set +u
  module purge || true
  module load Bioinformatics cellranger || true
  module load snakemake || true
  set -u
fi

# Unlock if needed
echo "🔓 Ensuring Snakemake working directory is unlocked…"
snakemake --unlock \
  --snakefile "$SNAKEFILE_ABS" \
  --directory "$TEST_DIR" || true

# Run Snakemake target
start_time=$(date +%s)

# Respect DRY_RUN=1 to pass -n to Snakemake
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
