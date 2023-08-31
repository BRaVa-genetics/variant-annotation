#!/usr/bin/env bash
#SBATCH --account=lindgren.prj
#SBATCH --job-name=annotate
#SBATCH --chdir=/well/lindgren/barney/brava-ground-truth-annotations/SAIGE-annotations-for-BRaVa/
#SBATCH --output=logs/annotate.log
#SBATCH --error=logs/annotate.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-23

#set -o errexit
#set -o nounset

get_current_cluster() {
  if [ ! -z "${SGE_ACCOUNT:-}" ]; then
    echo "sge"
  elif [ ! -z "${SLURM_JOB_ID:-}" ]; then
    echo "slurm"
  else
    raise_error "Could not find SGE/SLURM variables!"
  fi
}

get_array_task_id() {
  local _cluster=$( get_current_cluster )
  if [ "${_cluster}" == "sge" ]; then
    echo "${SGE_TASK_ID}"
  elif [ "${_cluster}" == "slurm" ]; then
    echo "${SLURM_ARRAY_TASK_ID}"
  else
    raise_error "${_cluster} is not a valid cluster!"
  fi
}

get_chr() {
  if [ $1 -eq 23 ]; then
    echo "X"
  elif [ $1 -eq 24 ]; then
    echo "Y"
  elif [ $1 -ge 1 ] && [ $1 -le 22 ]; then
    echo $1
  else
    raise_error "chromosome number must be between 1 and 24"
  fi
}

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly vep="/well/lindgren/barney/variant-annotation/run-vep/out/ukb_wes_450k.qced.chr${chr}.vep_processed.txt"

readonly spliceai_dir="/well/lindgren/barney/spliceai/out"
readonly spliceai_path="${spliceai_dir}/ukb_wes_450k.qced.v6.sites_only.${chr}.all.vcf"

readonly out_dir="data/out"
readonly out_file="${out_dir}/ukb_wes_450k.july.qced.brava_common_rare.v7.chr${chr}.saige.txt"

readonly vep_snp_id_col="SNP_ID"
readonly vep_gene_col="GENE"
readonly vep_lof_col="LOF"
readonly vep_max_pop_col="MAX_AF"
readonly vep_mane_select_col="MANE_SELECT"
readonly vep_revel_col="REVEL_SCORE"
readonly vep_cadd_phred_col="CADD_PHRED"
readonly vep_consequence_col="CSQ"
readonly vep_canonical_col="CANONICAL"
readonly vep_biotype_col="BIOTYPE"

mkdir -p ${out_dir}

. ~/job-scripts/conda.sh base

python3 scripts/brava_create_annot.py \
    --chr=$chr \
    --vep=$vep \
    --spliceai=$spliceai_path \
    --out_file=$out_file \
    --maf_cutoff=1.0 \
    --vep_snp_id_col=$vep_snp_id_col \
    --vep_gene_col=$vep_gene_col \
    --vep_lof_col=$vep_lof_col \
    --vep_max_pop_col=$vep_max_pop_col \
    --vep_mane_select_col=$vep_mane_select_col \
    --vep_revel_col=$vep_revel_col \
    --vep_cadd_phred_col=$vep_cadd_phred_col \
    --vep_consequence_col=$vep_consequence_col \
    --vep_canonical_col=$vep_canonical_col \
    --vep_biotype_col=$vep_biotype_col

gzip $out_file
