#!/usr/bin/env bash
set -euo pipefail

# ————————————————————————
# Load user config
# ————————————————————————
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$SCRIPT_DIR/config.sh"

# ————————————————————————
# Resolve absolute paths
# ————————————————————————
REPO_ROOT="$( dirname "$SCRIPT_DIR" )"

TEST_DIR_ABS="$TEST_DIR"
mkdir -p "$TEST_DIR_ABS"

CONFIG_DEST="$TEST_DIR_ABS/multi_config.csv"
PROBE_FILE="$PROBE_PATH"
REF_GENOME="$TURBO_REF_BASE/$REF_SUBPATH"
SNAKEFILE_ABS="$REPO_ROOT/$SNAKEFILE"

# ————————————————————————
# Load modules (workaround for unbound‐var in Cell Ranger module)
# ————————————————————————
echo "⏳ Loading modules…"
set +u
module purge
module load Bioinformatics cellranger
module load snakemake
set -u

# ————————————————————————
# Copy & patch multi_config.csv
# ————————————————————————
echo "📄 Copying original config → $CONFIG_DEST"
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

# ————————————————————————
# Run Snakemake → Cell Ranger multi
# ————————————————————————
echo "🚀 Running Snakemake (Cell Ranger multi)…"
snakemake -j "$CORES" \
  --snakefile "$SNAKEFILE_ABS" \
  --directory "$TEST_DIR_ABS" \
  "$OUTPUT_ID"

echo "✅ Workflow complete for '$OUTPUT_ID'."