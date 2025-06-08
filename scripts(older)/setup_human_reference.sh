#!/bin/bash

###############################################################################
# Build a 10x Cell Ranger Reference (GRCh38-2020-A) for Flex
# Based on 10x Genomics instructions: https://support.10xgenomics.com
#
# Builds reference in /nfs/turbo/umms-parent/Manny_human_ref
###############################################################################

set -e  # Exit on any error

GENOME="GRCh38"
VERSION="2020-A"
BUILD="GRCh38-GENCODEv32_build"
BASE_DIR="/nfs/turbo/umms-parent/Manny_human_ref"
BUILD_DIR="${BASE_DIR}/${BUILD}"
SOURCE_DIR="${BASE_DIR}/reference_sources"
FASTA_URL="https://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
FASTA_IN="${SOURCE_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_IN="${SOURCE_DIR}/gencode.v32.primary_assembly.annotation.gtf"

echo "Creating necessary folders..."
mkdir -p "$BUILD_DIR" "$SOURCE_DIR"

echo "Downloading FASTA and GTF if needed..."
[ ! -f "$FASTA_IN" ] && curl -sS "$FASTA_URL" | zcat > "$FASTA_IN"
[ ! -f "$GTF_IN" ] && curl -sS "$GTF_URL" | zcat > "$GTF_IN"

echo "Modifying FASTA headers..."
FASTA_MODIFIED="${BUILD_DIR}/$(basename "$FASTA_IN").modified"
cat "$FASTA_IN" \
  | sed -E 's/^>(\S+).*/>\1 \1/' \
  | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
  | sed -E 's/^>MT />chrM /' \
  > "$FASTA_MODIFIED"

echo "Stripping version suffixes from GTF IDs..."
GTF_MODIFIED="${BUILD_DIR}/$(basename "$GTF_IN").modified"
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$GTF_IN" \
  | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
  | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
  | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
  > "$GTF_MODIFIED"

echo "Filtering GTF by gene biotypes..."
BIOTYPE_PATTERN="(protein_coding|protein_coding_LoF|lncRNA|IG_.*|TR_.*|TR_.*_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
ALLOWLIST="${BUILD_DIR}/gene_allowlist"

cat "$GTF_MODIFIED" \
  | awk '$3 == "transcript"' \
  | grep -E "$GENE_PATTERN" \
  | grep -E "$TX_PATTERN" \
  | grep -Ev "$READTHROUGH_PATTERN" \
  | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
  | sort \
  | uniq \
  > "$ALLOWLIST"

echo "Filtering final GTF file..."
GTF_FILTERED="${BUILD_DIR}/$(basename "$GTF_IN").filtered"
grep -E "^#" "$GTF_MODIFIED" > "$GTF_FILTERED"
grep -Ff "$ALLOWLIST" "$GTF_MODIFIED" \
  | awk -F "\t" '$1 != "chrY" || $4 >= 2752083 && $4 < 56887903 && !/ENSG00000290840/' \
  >> "$GTF_FILTERED"

echo "Loading Bioinformatics + cellranger module..."
module purge
module load Bioinformatics cellranger

echo "Running cellranger mkref..."

cellranger mkref \
  --genome="$GENOME" \
  --fasta="$FASTA_MODIFIED" \
  --genes="$GTF_FILTERED" \
  --nthreads=16 \
  --output-dir="$BUILD_DIR/refdata-gex-${GENOME}-${VERSION}"

echo "Reference build complete at:"
echo "$BUILD_DIR/refdata-gex-${GENOME}-${VERSION}"