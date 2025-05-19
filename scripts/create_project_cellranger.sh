#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# scripts/create_project_cellranger.sh
#
# - Source config.sh
# - Auto-detect refdata-gex* under TURBO_REF_BASE
# - Copy & patch multi_config.csv
# - Unlock Snakemake dir if locked
# - Launch Snakemake with --rerun-incomplete
# ─────────────────────────────────────────────────────────────────────────────

# 1) Load config
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$SCRIPT_DIR/config.sh"

# 2) Prep paths
REPO_ROOT="$( dirname "$SCRIPT_DIR" )"
mkdir -p "$TEST_DIR"
CONFIG_DEST="$TEST_DIR/multi_config.csv"
PROBE_FILE="$PROBE_PATH"
SNAKEFILE_ABS="$REPO_ROOT/$SNAKEFILE"

# 3) Use manually set reference folder
REF_GENOME="$TURBO_REF_BASE/$REF_SUBPATH"
if [[ ! -d "$REF_GENOME" ]]; then
  echo "❌ ERROR: Reference folder not found at $REF_GENOME"
  exit 1
fi
echo "🧬 Using reference: $REF_GENOME"

# 4) Load modules (disable nounset to avoid CellRanger's unbound vars)
echo "⏳ Loading modules…"
set +u
module purge
module load Bioinformatics cellranger
module load snakemake
set -u

# 5) Copy & patch the multi_config.csv
echo "📄 Copying accessible multi_config file → $CONFIG_DEST"
cp -f "$TURBO_CONFIG_SOURCE" "$CONFIG_DEST"

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

echo "🧬 Patching reference → $REF_GENOME"
sed -i "s|^reference,.*|reference,$REF_GENOME|" "$CONFIG_DEST"

echo "🧬 Patching probe-set → $PROBE_FILE"
sed -i "s|^probe-set,.*|probe-set,$PROBE_FILE|" "$CONFIG_DEST"

# 5.5) Unlock Snakemake directory if needed
echo "🔓 Ensuring Snakemake working directory is unlocked…"
snakemake --unlock \
  --snakefile "$SNAKEFILE_ABS" \
  --directory "$TEST_DIR" || true

# 6) Run Snakemake → Cell Ranger multi with --rerun-incomplete and timing
start_time=$(date +%s)

echo "🚀 Running Snakemake (Cell Ranger multi)…"
snakemake --rerun-incomplete -j "$CORES" \
  --snakefile "$SNAKEFILE_ABS" \
  --directory "$TEST_DIR" \
  "$OUTPUT_ID"

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))

echo "✅ Workflow complete for '$OUTPUT_ID'."
echo "⏱ Finished at: $(date)"
echo "⏱ Total runtime: $(( elapsed / 60 )) minutes and $(( elapsed % 60 )) seconds"
