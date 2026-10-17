#!/bin/bash
# usage: lookup_chr.sh CHROM KIND KEYS_DIR OUT_DIR      e.g. lookup_chr.sh chr21 snv keys scores
# Streams one SpliceAI2 precomputed file (KIND snv or indel; there is no chrY indel file) from Hugging Face and
# keeps the rows matching KEYS_DIR/CHROM.KIND.keys.tsv (make_keys.py): OUT_DIR/CHROM.KIND.tsv.gz, and the keys the
# release lacks in OUT_DIR/CHROM.KIND.missing.tsv (gpu_set.py picks those worth scoring). Nothing large is stored:
# the SNV file of chr1 alone is ~30 GB compressed. A finished output is skipped, so rerun after any failure.
# Each chromosome is independent: run them as an array job (47 tasks: 24 SNV, 23 indel).
set -euo pipefail
c=$1; kind=$2; keys=$3/$c.$kind.keys.tsv; out=$4/$c.$kind.tsv.gz
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$4"
[[ -s $out ]] && { echo "$c $kind: done already"; exit 0; }
[[ -s $keys ]] || { echo "$c $kind: no keys"; exit 0; }
B=https://huggingface.co/datasets/illumina-ai/SpliceAI2-data/resolve/main/precomputed_scores_v2.0/GRCh38
t0=$(date +%s)
# column 84 of the release is spliceai2_summary_score (the max of the four site deltas, v1's DS_MAX analogue).
# Written under a .part name and renamed only when every stage of the pipe succeeded: a failed download never
# leaves a file that a rerun would skip.
bash "$here/hf_ranged_cat.sh" "$B/$kind/$c.tsv.gz" | gzip -dc | cut -f1-5,84 \
  | python3 "$here/lookup.py" "$keys" "${out%.tsv.gz}.part.tsv.gz"
mv "${out%.tsv.gz}.part.tsv.gz" "$out"; mv "${out%.tsv.gz}.part.missing.tsv" "${out%.tsv.gz}.missing.tsv"
echo "$c $kind: $(( $(date +%s) - t0 )) s"
