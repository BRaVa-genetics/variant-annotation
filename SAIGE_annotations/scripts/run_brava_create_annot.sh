#!/bin/bash

#SBATCH -J annot_long
#SBATCH -o annot_long-%A_%a.out
#SBATCH -e annot_long-%A_%a.err
#SBATCH -p short
#SBATCH --array=1-23

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

# Load R module
module purge
module load BCFtools

spliceai_dir="/well/lindgren/barney/spliceai/out"
vep_annotations_dir="/well/lindgren/barney/variant-annotation/run-vep/out"
out_dir="/well/lindgren/UKBIOBANK/dpalmer/annotation_munging"
annotation_repo_dir="/well/lindgren/UKBIOBANK/dpalmer/variant-annotation"
brava_create_annot="${annotation_repo_dir}/SAIGE_annotations/scripts/brava_create_annot.py"

SLURM_ARRAY_TASK_ID=21

mkdir -p ${out_dir}

# Map SLURM_ARRAY_TASK_ID to chromosome
if [ "$SLURM_ARRAY_TASK_ID" -eq 23 ]
then
    chr="X"
else
    chr="$SLURM_ARRAY_TASK_ID"
fi

vep="${vep_annotations_dir}/ukb_wes_450k.qced.chr${chr}.vep_processed.txt"
spliceai="${spliceai_dir}/ukb_wes_450k.qced.v6.sites_only.${chr}.all.vcf"
spliceai_no_chr="${out_dir}/ukb_wes_450k.qced.v6.sites_only.no_chr.${chr}.all.vcf.gz"
out="${out_dir}/ukb_wes_450k.july.qced.brava_common_rare.v7.chr${chr}.saige.txt"

# First, remove the chr from the CHROM column
bcftools annotate --rename-chrs <(echo -e "chr${chr}\t${chr}") ${spliceai} -O z -o ${spliceai_no_chr}
python3 ${brava_create_annot} --vep ${vep} --spliceai ${spliceai_no_chr} --out_file $out