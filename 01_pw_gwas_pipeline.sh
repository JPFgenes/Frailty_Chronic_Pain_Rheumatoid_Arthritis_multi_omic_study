#!/usr/bin/env bash
# Run GWAS-PW, summarise blocks, extract per-SNP posteriors, build CS90,
# and optionally annotate credible-set variants with BEDtools.
# Usage: bash scripts/01_pw_gwas_pipeline.sh INPUT.tsv BLOCKS.bed TRAIT1 TRAIT2 OUTPUT_PREFIX [GENES.bed] [EXCLUDE_APOE]
set -euo pipefail
if (( $# < 5 || $# > 7 )); then
 echo 'Usage: 01_pw_gwas_pipeline.sh INPUT BLOCKS TRAIT1 TRAIT2 PREFIX [GENES_BED] [EXCLUDE_APOE=FALSE]' >&2; exit 2
fi
input=$1; blocks=$2; trait1=$3; trait2=$4; prefix=$5
genes=${6:-NONE}; exclude_apoe=${7:-FALSE}
for f in "$input" "$blocks"; do [[ -s "$f" ]] || { echo "Missing $f" >&2; exit 1; }; done
command -v gwas-pw >/dev/null || { echo 'Missing gwas-pw (tested source version 0.21)' >&2; exit 1; }
command -v Rscript >/dev/null || { echo 'Missing Rscript' >&2; exit 1; }
scriptdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$(dirname "$prefix")"
# Study-specific input preparation (Z-scores, variance and sample overlap)
# must be performed and documented before running this script.
gwas-pw -i "$input" -phenos "$trait1" "$trait2" -bed "$blocks" -o "${prefix}_post"
seg="${prefix}_post.segbfs.gz"; bfs="${prefix}_post.bfs.gz"
[[ -s "$seg" && -s "$bfs" ]] || { echo 'Expected GWAS-PW outputs not found' >&2; exit 1; }
gwas-pw --extract-snp-posteriors "$bfs" -i "$input" -phenos "$trait1" "$trait2" -bed "$blocks" -o "${prefix}_extracted"
extracted="${prefix}_extracted.bfs.gz"
[[ -s "$extracted" ]] || { echo 'SNP posterior extraction failed' >&2; exit 1; }
Rscript "$scriptdir/02_pw_gwas_postprocess.R" "$seg" "$extracted" "$prefix" "$exclude_apoe"
if [[ "$genes" != NONE ]]; then
 [[ -s "$genes" ]] || { echo "Gene BED not found: $genes" >&2; exit 1; }
 command -v bedtools >/dev/null || { echo 'Gene annotation requires bedtools' >&2; exit 1; }
 variants="${prefix}_CS90_variants.tsv"
 tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
 awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++)h[$i]=i;next} {p=$(h["pos"]); if(p>0) print $(h["chr"])"\t"p-1"\t"p"\t"$(h["block"])"\t"$(h["PPA_3"])}' "$variants" | LC_ALL=C sort -k1,1 -k2,2n > "$tmp/variants.bed"
 LC_ALL=C sort -k1,1 -k2,2n "$genes" > "$tmp/genes.bed"
 if [[ -s "$tmp/variants.bed" ]]; then
   bedtools closest -d -a "$tmp/variants.bed" -b "$tmp/genes.bed" > "${prefix}_CS90_nearest_genes.tsv"
   # For 4-column gene BED (chr, start, end, gene name): prioritise by max SNP PPA_3.
   awk 'BEGIN{OFS="\t"}{key=$4 OFS $9; if(!(key in best) || $5>best[key])best[key]=$5}END{print "block","gene","max_variant_PPA3";for(k in best)print k,best[k]}' "${prefix}_CS90_nearest_genes.tsv" > "${prefix}_gene_priorities.tsv"
 fi
fi
echo "PW-GWAS pipeline complete: ${prefix}_*"
