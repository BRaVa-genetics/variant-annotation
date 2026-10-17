#!/bin/bash
#SBATCH -J brava_annot
#SBATCH -o brava_annot-%A_%a.out
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH --array=1-23
# Example: steps 4 and 5 for every chromosome as a SLURM array (task 23 = chrX). Edit the paths below.
set -euo pipefail

repo=/path/to/variant-annotation
us=/path/to/universal-saige                                    # resources/ from download_resources.sh --alphamissense
vep_dir=/path/to/vep105_loftee/out                             # step 1 tables
sai2=/path/to/spliceai2                                        # step 2: scores/ and gpu/
out_dir=/path/to/group_files

c=${SLURM_ARRAY_TASK_ID}; [[ $c -eq 23 ]] && c=X
mkdir -p "$out_dir"
g=$out_dir/brava.chr$c.saige_group.txt.gz
python3 "$repo/SAIGE_annotations/scripts/brava_create_annot.py" \
    -v "$vep_dir/sites_only_output_chr${c}_vep.gnomad_popmax_0.01_processed.txt" \
    --spliceai2 $(ls "$sai2"/scores/chr$c.{snv,indel}.tsv.gz "$sai2"/gpu/chunk_*.scores.tsv.gz 2>/dev/null) \
    --gtf "$us/resources/gencode.v39.annotation.gtf.gz" \
    -w "$g"
bash "$repo/SAIGE_annotations/scripts/alphamissense_files.sh" "$us" "$g" "$c"
